# GapFill - Gap Filling Pipeline

## Quick Reference

```bash
# 单倍体
python -m gapfill -a assembly.fa --hifi hifi.fq --ont ont.fq -o output

# 多倍体 (自动检测)
python -m gapfill -a hap1.fa hap2.fa --hifi hifi.fq --ont ont.fq -o output

# 多倍体 + 优化模式 (减少75%比对)
python -m gapfill -a hap1.fa hap2.fa --hifi hifi.fq --optimized -o output

# 断点续跑 (中断后恢复)
python -m gapfill -a assembly.fa --hifi hifi.fq -o output --resume

# 清除断点，重新开始
python -m gapfill -a assembly.fa --hifi hifi.fq -o output --clear-checkpoint
```

## 核心原则

1. **N 长度无意义** - N 序列只是占位符，长度不代表真实 gap 大小
2. **Gap 标准化** - 所有 N 占位符标准化为 500bp，确保 spanning reads 检测正确
3. **预处理阶段** - 在迭代前完成 normalize、pre-assess、polish
4. **本地验证 (validate-before-apply)** - 填充后立即本地验证，只有通过验证的填充才应用到 assembly
5. **HiFi 优先** - 7层策略优先使用高准确性 HiFi，ONT 提供长度优势
6. **多倍体先标准化** - 在 SNP 检测之前先标准化所有 haplotype 的 gap
7. **Alignment-based SNP 检测** - 使用 minimap2 比对检测 SNP，正确处理 haplotype 独有的 gap

## 性能优化 (v2.0)

### 优化策略一览

| 优化 | 效果 | 原理 |
|------|------|------|
| **Reads 过滤** | 减少 80% 比对时间 | 过滤掉锚定在非gap区域的reads，后续只比对有用的reads |
| **并行填充** | 3-5x 加速 | 多进程同时处理多个gap |
| **共识优先** | 减少 50% wtdbg2 调用 | 高质量 spanning reads 直接取共识，跳过组装 |
| **本地验证** | 消除回退迭代 | 填充后立即本地验证，失败的填充不会被应用 |

**预期整体加速：4-5x** (例如 74h → 16h)

### Reads 过滤策略 (utils/reads_cache.py)

```
原始流程 (慢):
  每轮迭代: 全量比对所有 reads → assembly
  10 轮 × 50GB HiFi = 数十小时

优化流程 (快):
  预处理: 全量比对一次 → 过滤掉锚定在非gap区域的reads
  后续迭代: 只比对保留的 reads (可能只有 10-20%)
```

**ReadsCache 类方法：**
- `set_gap_regions()` - 设置 gap 区域，建立索引
- `filter_bam()` - 过滤 BAM，提取有用 reads 到 FASTA
- `get_filtered_reads_path()` - 获取过滤后的 reads 文件路径

**过滤规则：**

| Read 类型 | 判定 | 说明 |
|----------|------|------|
| Spanning | 保留 | 跨越 gap 的 reads |
| Overlap | 保留 | 与 gap 区域重叠的 reads |
| Near-gap | 保留 | 靠近 gap 边界的 reads (1000bp 内) |
| Soft-clipped | 保留 | 有 soft-clip 的 reads (可能跨越 gap) |
| Anchored | 过滤 | 完全锚定在非 gap 区域的 reads |

**关键实现逻辑：**
```python
def _is_useful_read(self, read) -> Tuple[bool, str]:
    # Case 1: Read spans the gap
    if read_start <= gap_start and read_end >= gap_end:
        return True, 'spanning'
    # Case 2: Read overlaps with gap region
    if read_start < gap_end and read_end > gap_start:
        return True, 'overlap'
    # Case 3: Soft-clips near gap
    # Case 4: Near gap boundary
    # Case 5: Anchored in non-gap region → filter out
    return False, 'anchored'
```

### 共识优先策略 (core/consensus.py)

**触发条件：**
- HiFi spanning reads >= 3, or ONT spanning reads >= 5
- K-mer 一致性 > 85% (POA) or > 95% (direct)

**ConsensusBuilder 类：**
- `build_consensus(reads, gap_info, read_type)` - 主入口，选择最佳策略
- `_estimate_identity()` - K-mer based 一致性估计
- `_direct_consensus()` - 直接 majority vote (>95% 一致性, HiFi only)
- `_poa_consensus()` - POA 共识 (>85% 一致性，使用 spoa/abpoa)

**共识策略选择：**

| Read 类型 | 一致性 | Reads 数 | 策略 |
|-----------|--------|---------|------|
| HiFi | >95% | >=3 | Direct consensus (majority vote) |
| HiFi | >85% | >=3 | POA consensus (spoa/abpoa) |
| ONT | >85% | >=5 | POA consensus only |
| any | <85% | any | Fall back to wtdbg2 |

### 本地验证 (validate-before-apply)

**机制：**
- 每个 tier 组装成功后，立即构建本地参考序列 (flank + fill + flank)
- 从当前 BAM 中提取附近的 reads，比对到本地参考
- 使用与延迟验证相同的标准 (spanning reads, coverage) 检查填充质量
- 验证失败的填充不会被应用，而是尝试下一个 tier

**效果：**
- 消除了整个延迟验证/回退机制 (~350 行代码)
- 不再产生 assembly_reverted.fasta
- 不再需要 validation.ok 标记
- 每轮迭代都有效推进，不会浪费在回退上

### 并行 Gap 填充 (core/parallel.py)

**实现方式：**
- 使用 `ProcessPoolExecutor` 多进程并行
- `spawn` context 避免 pysam fork 问题
- 每个 worker 独立创建 GapFiller 实例

**关键函数：**
- `fill_gaps_parallel()` - 并行填充多个 gap
- `fill_gaps_sequential()` - 顺序填充 (fallback)
- `GapBatcher.create_batches()` - 按优先级批量处理

**Worker 数量：**
```python
max_workers = min(threads, len(gaps), 16)  # 最多 16 个进程
```

**异常处理：**
- 单个 gap 超时 (5分钟) → 标记失败
- 并行执行异常 → 自动 fallback 到顺序处理

### 命令行参数

```bash
# 默认启用所有优化
python -m gapfill -a assembly.fa --hifi hifi.fq -o output

# 禁用 reads 过滤 (内存不足时)
python -m gapfill -a assembly.fa --hifi hifi.fq -o output --no-filter-reads

# 禁用并行填充 (调试时)
python -m gapfill -a assembly.fa --hifi hifi.fq -o output --no-parallel
```

## 包结构

```
gapfill/
├── __init__.py
├── __main__.py
├── cli.py                    # 命令行解析
├── core/
│   ├── filler.py             # GapFiller - 7层策略 (含共识优先)
│   ├── validator.py          # GapValidator + GapStatusTracker
│   ├── polisher.py           # FlankPolisher - 侧翼序列抛光
│   ├── consensus.py          # ConsensusBuilder - 共识序列构建
│   └── parallel.py           # 并行gap填充
├── engines/
│   ├── haploid.py            # HaploidEngine - 单倍体引擎 + 优化
│   ├── polyploid.py          # PolyploidEngine - 多倍体引擎
│   └── optimized_polyploid.py # OptimizedPolyploidEngine - 批量比对优化
└── utils/
    ├── indexer.py            # AssemblyIndexer
    ├── scanner.py            # GapScanner
    ├── tempfiles.py          # TempFileManager
    ├── checkpoint.py         # CheckpointManager - 断点续跑
    └── reads_cache.py        # ReadsCache - reads过滤与缓存
```

## 工作流程

### 单倍体工作流程 (优化版)

```
┌─────────────────────────────────────────────────────────────────┐
│  PREPROCESSING PHASE (一次性，在迭代前)                          │
├─────────────────────────────────────────────────────────────────┤
│  0a. Find initial gaps                                          │
│  0b. Normalize gaps → 500N                                      │
│  0c. Initial alignment → BAM                                    │
│  0d. [优化] Filter reads → 过滤锚定在非gap区域的reads           │
│  0e. Pre-assess gap flanks                                      │
│      ├─ PENDING: 侧翼正常                                       │
│      └─ NEEDS_POLISH: 侧翼有问题                                │
│  0f. Polish problematic flanks → polished assembly              │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  ITERATION LOOP (使用过滤后的reads)                              │
├─────────────────────────────────────────────────────────────────┤
│  1. [优化] Align FILTERED reads → BAM (而非全量reads)           │
│  2. Find remaining gaps                                         │
│  3. [优化] Fill gaps in PARALLEL (每个 tier 含本地验证)         │
│      ├─ TIER 0: 共识优先 (跳过wtdbg2)                           │
│      └─ TIER 1-6: 常规策略                                      │
│      每个 tier: 组装 → 本地验证 → 通过则返回，失败则下一 tier    │
│  4. Process results:                                            │
│      ├─ 已验证完全填充 → FILLED_COMPLETE                        │
│      ├─ 已验证部分填充 → FILLED_PARTIAL (应用前截断失败侧)      │
│      └─ 组装但验证全失败 → FAILED                               │
│  5. Apply validated fills to assembly                           │
└─────────────────────────────────────────────────────────────────┘
```

### 多倍体工作流程 (优化版)

```
┌─────────────────────────────────────────────────────────────────┐
│  STEP 0: Normalize ALL haplotypes → 500N                        │
│  STEP 0b: [优化] Filter reads                                   │
│      ├─ 全量比对一次到 ref_haplotype                            │
│      └─ 过滤掉锚定在非gap区域的reads                            │
│  STEP 1: SNP detection (alignment-based)                        │
│  STEP 2: Phase reads (使用过滤后的reads)                        │
│      ├─ HiFi reads → phased_hap1_hifi.fa, phased_hap2_hifi.fa   │
│      └─ ONT reads → phased_hap1_ont.fa, phased_hap2_ont.fa      │
│  STEP 3: Gap filling per haplotype                              │
└─────────────────────────────────────────────────────────────────┘
```

### Gap 状态流转

```
初始: PENDING
       │
预处理 ├──→ NEEDS_POLISH ──→ [Polish] ──→ PENDING
       │
迭代 N ├──→ 组装 + 本地验证通过 (完全) ──→ FILLED_COMPLETE
       │
       ├──→ 组装 + 本地验证通过 (部分) ──→ FILLED_PARTIAL
       │     (应用前截断验证失败的一侧)
       │
       ├──→ 组装成功但所有 tier 验证失败 ──→ FAILED (下轮重试)
       │
       └──→ 填充失败 ──→ UNFILLABLE (永久跳过)
                    └──→ FAILED (下轮重试)
```

**部分填充截断逻辑 (应用前)：**
```
原始: [左填充][500N][右填充]

左失败+右通过: [500N][右填充]
左通过+右失败: [左填充][500N]
两侧都失败:    [500N]
```

### GapStatus 枚举

| 状态 | 说明 | 下轮操作 |
|------|------|---------|
| PENDING | 待处理 | 尝试填充 |
| NEEDS_POLISH | 侧翼需要 Polish | 预处理阶段 Polish |
| FILLED_COMPLETE | 完全填充，已验证 | 跳过 |
| FILLED_PARTIAL | 部分填充，已验证 | 跳过 |
| UNFILLABLE | 确认无法填充 | 永久跳过 |
| FAILED | 填充/验证失败 | 下轮重试 |

## 核心组件

### 1. 7层 HiFi/ONT 分层策略 (core/filler.py)

```
TIER 0: Direct consensus       → [优化] 高质量reads直接取共识，跳过wtdbg2
TIER 1: HiFi-only spanning     → 最高准确性
TIER 2: ONT-only spanning      → 长度优势 + 可选 HiFi 抛光
TIER 3: Hybrid spanning        → 混合跨越reads
TIER 4: HiFi flanking + merge  → HiFi 侧翼组装 + 智能合并
TIER 5: ONT flanking + merge   → ONT 侧翼组装 + 可选抛光
TIER 6: Hybrid flanking + 500N → 最后备选

**Flanking 填充说明：**
- 左侧和右侧 reads 分别组装
- 尝试合并左右序列（寻找重叠区域）
- 若无法合并：`left_seq + 500N + right_seq`（标记为 `*_flanking_500N`）
- `_500N` 后缀表示两侧无法合并，不是只填充了一侧
```

**TIER 0 共识优先策略：**
- HiFi: 当 spanning reads >= 3 且一致性 > 85% 时触发
- ONT: 当 spanning reads >= 5 且一致性 > 85% 时触发 (POA only)
- 直接计算共识序列（POA 或 majority vote），跳过 wtdbg2
- 减少约 50% 的 wtdbg2 调用

**关键方法：**
- `_collect_reads_by_type()` - 分别收集 HiFi/ONT reads
- `_validate_fill_locally()` - 本地验证填充 (构建 local ref + align)
- `try_consensus_fill()` - 共识优先填充 (core/consensus.py)
- `_assemble_spanning_reads()` - 类型特定的 wtdbg2 参数
- `_cluster_reads()` - 组装前 reads 聚类，减少内部 N
- `_polish_with_hifi()` - 用 racon 抛光 ONT 组装
- `_run_wtdbg2_assembly()` - HiFi 用 ccs preset，ONT 用 ont preset

**Reads 聚类策略 (v2.1)：**
- 组装前使用 `minimap2 -x ava-pb` 进行 all-vs-all 比对
- 基于序列相似性 (>90% identity) 聚类 reads
- 每个 cluster 独立组装，选择最佳结果
- 避免混合不同来源的 reads 导致内部 N

**内部 N 处理 (保守策略)：**
- wtdbg2 组装结果若包含内部 N，直接拒绝
- 不尝试拆分或部分使用含 N 的序列
- 确保所有填充序列的完整性和真实性

### 2. 本地验证 (core/validator.py)

**validate-before-apply 架构：**
- 每个 tier 组装成功后，立即在本地验证
- 构建本地参考序列: `left_flank (5kb) + fill_sequence + right_flank (5kb)`
- 从当前 BAM 提取区域内 reads，比对到本地参考
- 使用与原有验证相同的标准检查质量

**GapValidator 类：**
- `validate_fill_locally()` - 本地验证填充 (构建 local ref + align + validate)
- `validate_complete_fill()` - 验证完全填充
- `validate_partial_fill()` - 验证部分填充（独立左右判断）
- `pre_assess_gap()` - 预评估单个 gap 侧翼
- `pre_assess_gaps()` - 批量预评估
- `analyze_failed_gap()` - 分析失败原因
- `_extract_reads_from_bams()` - 从 BAM 提取区域内 reads

**验证标准（完全填充 vs 部分填充）：**

| 类型 | 验证标准 | 说明 |
|------|----------|------|
| **完全填充** | Spanning reads ≥1-2 | 必须有 reads 跨越整个填充区域 |
| **部分填充** | Junction coverage ≥5, Insert coverage ≥5, 无断点 | 只检查连接点质量 |

**部分填充独立左右验证：**
- 左侧和右侧独立判断通过/失败
- 通过的一侧保留，失败的一侧在应用前截断
- 使用 `PartialFillValidationResult` 返回独立结果

```python
@dataclass
class PartialFillValidationResult:
    valid: bool              # 至少一侧通过
    status: GapStatus
    left_valid: bool         # 左侧是否通过
    left_fill_length: int
    left_junction_coverage: float
    right_valid: bool        # 右侧是否通过
    right_fill_length: int
    right_junction_coverage: float
```

**GapStatusTracker 类：**
- `should_attempt()` - 判断是否应该尝试填充 (跳过 FILLED_COMPLETE, FILLED_PARTIAL, UNFILLABLE)
- `set_status()` - 设置 gap 状态
- `from_dict()` - 反序列化时自动将旧 FILLED_PENDING 状态转为 FAILED

### 3. 侧翼 Polish (core/polisher.py)

**FlankPolisher 类：**
- `polish_gap_flanks()` - Polish 单个 gap 的侧翼
- `polish_assembly_flanks()` - 批量 Polish 并更新 assembly

**Polish 方法：**
1. **Racon polish**（首选）- 使用 racon 进行 consensus 校正
2. **Pileup consensus**（备选）- 基于 pileup 的多数投票

**触发条件：**
- Clip accumulation: ≥5 reads 在同一位置 clip
- High mismatch density: > 5% 错配率

### 4. 多倍体 SNP 检测 (engines/polyploid.py)

**流程：**
```
STEP 0: Gap 标准化 ALL haplotypes (统一为 500N)
STEP 1: Alignment-based SNP 检测 (minimap2 asm5)
STEP 2: Phasing reads
STEP 3: Gap filling (skip_normalization=True)
```

**关键类和方法：**
- `ReadPhaser.normalize_all_assemblies()` - 标准化所有 haplotype
- `ReadPhaser._detect_snps_alignment()` - alignment-based SNP 检测
- `ReadPhaser._assign_read_to_haplotype_detailed()` - SNP-based 分配

### 5. 断点续跑 (.ok 文件机制)

**使用方法：**
```bash
# 首次运行
python -m gapfill -a assembly.fa --hifi hifi.fq -o output

# 中断后恢复
python -m gapfill -a assembly.fa --hifi hifi.fq -o output --resume

# 清除断点，重新开始
python -m gapfill -a assembly.fa --hifi hifi.fq -o output --clear-checkpoint
```

**进度标记 (.ok 文件)：**

```
output/
├── normalized.ok              # gap 标准化完成
├── preprocessing.ok           # 预处理完成
├── gap_tracker.json           # GapStatusTracker 状态
├── preprocessing/
│   ├── alignment.ok           # 初始比对完成
│   ├── filter.ok              # reads 过滤完成
│   └── polish.ok              # flank polish 完成
├── iteration_1/
│   ├── alignment.ok           # 比对完成
│   └── iteration.ok           # 本轮完成 (含填充+应用)
└── complete.ok                # 全部完成
```

**关键改进：**
- 使用 `.ok` 文件而非 checkpoint.json 判断进度
- `gap_tracker.json` 单独保存 GapStatusTracker 状态
- 每步完成后立即创建 `.ok` 文件
- resume 时检查 `.ok` 文件决定跳过哪些步骤

**状态跟踪：**

| 阶段 | 单倍体 | 多倍体 | .ok 文件 | 可复用文件 |
|------|--------|--------|----------|-----------|
| normalization | ✓ | ✓ | normalized.ok | assembly_normalized.fasta |
| polish | ✓ | ✓ | preprocessing/polish.ok | assembly_polished.fasta |
| snp_detection | - | ✓ | snp_detection.ok | snp_database.json |
| phasing | - | ✓ | phasing.ok | phased_*_hifi.fasta |
| iteration N | ✓ | ✓ | iteration_N/iteration.ok | assembly_filled.fasta |

## 命令行参数

```
-a, --assembly FILE(s)    # 1个=单倍体，2+个=多倍体
--hifi FILE               # HiFi reads
--ont FILE                # ONT reads
-o, --output DIR          # 输出目录
-t, --threads N           # 线程数 (default: 8)
--max-iterations N        # 最大迭代 (default: 10)
--min-gap-size N          # 最小 gap 大小 (default: 100)
--phasing METHOD          # builtin | whatshap
--no-ambiguous-reads      # 多倍体不使用 ambiguous reads
--optimized               # 使用批量比对优化 (多倍体)
--resume                  # 从断点恢复运行
--clear-checkpoint        # 清除断点，从头开始
--no-filter-reads         # 禁用 reads 过滤优化 (单倍体)
--no-parallel             # 禁用并行填充 (单倍体)
-v, --verbose             # 详细日志
```

## 输出结构

**单倍体：**
```
output/
├── checkpoint.json              # 断点状态文件
├── preprocessing/               # 预处理阶段文件
│   ├── hifi.bam, ont.bam
│   └── polish_work/
├── assembly_normalized.fasta
├── assembly_polished.fasta      # (如果有 polish)
├── final_assembly.fasta
├── final_stats.json
└── iteration_N/
    ├── assembly_filled.fasta
    ├── hifi.bam, ont.bam
    └── work/gap_*/
```

**多倍体：**
```
output/
├── hap1_normalized.fasta
├── hap2_normalized.fasta
├── snp_database.json
├── phased_*_hifi.fasta
├── phased_*_ont.fasta
├── hap1/
│   └── final_assembly.fasta
├── hap2/
│   └── final_assembly.fasta
└── polyploid_summary.json
```

## 外部工具参数

### wtdbg2 (gap 组装)

| 参数 | HiFi | ONT | 说明 |
|------|------|-----|------|
| `-x` | ccs | ont | 读取类型 preset |
| `-L` | 1000 | 2000 | 最小读长 |
| `-e` | 2 | 2 | 边缘覆盖度阈值 |
| `-S` | 1 | 1 | 救回低覆盖度边缘 |
| `-g` | auto | auto | 基于 reads 统计估计 |

**组装前 reads 聚类：**
```bash
# all-vs-all alignment
minimap2 -x ava-pb -t {threads} reads.fa reads.fa > overlaps.paf
# 根据 overlap 信息聚类，每个 cluster 独立组装
```

**内部 N 检测与处理：**
- 组装完成后检查 `'N' in sequence`
- 若包含内部 N：直接拒绝该组装结果，返回 None
- 不拆分含 N 序列，确保数据完整性

### minimap2 (比对)

| 用途 | Preset | 额外参数 | 说明 |
|------|--------|----------|------|
| HiFi → Assembly | map-hifi | - | HiFi 读取比对 |
| ONT → Assembly | map-ont | - | ONT 读取比对 |
| Hap → Hap (SNP) | asm5 | - | 同种 haplotype 比对 |
| Phasing | map-hifi/ont | --secondary=no | 唯一分配 |

### racon (Polish)

```
racon -t {threads} reads.fa align.paf reference.fa
```

### samtools sort

```
samtools sort -@ {threads} -m 2G -o {output} -
```

- `-m 2G`: 每线程内存限制，防止 OOM

## 依赖

**Python:**
- biopython, pysam, numpy

**External:**
- minimap2, samtools, wtdbg2, wtpoa-cns
- racon (Polish 功能)
- bcftools, whatshap (多倍体 whatshap 模式)

## 版本历史

### v3.0 (Validate-Before-Apply Architecture)

**架构重构：**
- 从"延迟验证"切换到"本地验证 (validate-before-apply)"架构
- 每个 tier 组装成功后立即构建本地参考 (flank + fill + flank) 并验证
- 验证失败的 tier 自动 fallthrough 到下一个 tier
- 只有通过验证的填充才会被应用到 assembly

**代码精简 (~250 行净减少)：**
- 删除 `_validate_previous_fills()` (~165 行)
- 删除 `_revert_fills()` (~46 行)
- 删除 `_partial_revert_fills()` (~103 行)
- 删除 `PendingFill` 跟踪和坐标更新逻辑
- 简化 `_finalize_result()` (从 ~130 行到 ~15 行)
- 移除 `FILLED_PENDING` 状态流 (保留枚举值用于反序列化兼容)

**新增本地验证：**
- `GapValidator.validate_fill_locally()` - 本地验证主入口
- `GapValidator._extract_reads_from_bams()` - BAM reads 提取
- `GapFiller._validate_fill_locally()` - 填充流程中的验证包装器
- `HaploidEngine._truncate_partial_fill()` - 部分填充截断 (替代复杂回退逻辑)

**共识策略扩展：**
- HiFi 共识阈值降低: 5 reads → 3 reads
- POA 一致性阈值降低: 0.90 → 0.85
- 新增 ONT 共识支持 (POA only, ≥5 reads)

**并行填充改进：**
- Worker 进程启用本地验证 (`enable_validation=True`)
- 每个 worker 独立验证，无需全局协调

### v2.1 (Bug Fixes)

**参数传递修复：**
- `min_mapq` 参数现正确传递到 `OptimizedPolyploidEngine`
- 修复 `filter_bam()` 中硬编码的 `min_mapq=20`，改用 `self.min_mapq`
- CLI → Engine 调用现包含所有必要参数 (`min_gap_size`, `min_mapq`)

**Resume 机制修复：**
- 修复 `haploid.py:490` 的占位符 bug
- Resume 时检测不一致状态（apply.ok 不存在但 filling.ok 存在）
- 自动重新运行填充以获取真实序列

**Reads 聚类与内部 N 处理：**
- 新增 `_cluster_reads()` 函数，组装前聚类 reads
- 使用 `minimap2 -x ava-pb` all-vs-all 比对
- 保守策略：含内部 N 的组装结果直接拒绝
