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
4. **延迟验证** - 填充后下一轮迭代时用新 BAM 验证
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
| **高置信度同轮验证** | 减少 30% 迭代 | 高置信度填充直接确认，不等下一轮 |

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
- HiFi spanning reads >= 5
- K-mer 一致性 > 95%

**ConsensusBuilder 类：**
- `build_consensus()` - 主入口，选择最佳策略
- `_estimate_identity()` - K-mer based 一致性估计
- `_direct_consensus()` - 直接 majority vote (>95% 一致性)
- `_poa_consensus()` - POA 共识 (>90% 一致性，使用 spoa/abpoa)

**共识策略选择：**

| 一致性 | Reads 数 | 策略 |
|--------|---------|------|
| >95% | >=3 | Direct consensus (majority vote) |
| >90% | >=5 | POA consensus (spoa/abpoa) |
| <90% | any | Fall back to wtdbg2 |

### 高置信度同轮验证

**触发条件：**
- Tier 0 或 Tier 1 填充成功
- Spanning reads >= 5
- 覆盖度 >= 5x

**效果：**
- 满足条件的填充直接标记为 FILLED_COMPLETE
- 跳过延迟验证，减少约 30% 迭代次数

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
│  2. Validate previous fills (iteration 2+)                      │
│      ├─ 完全填充验证：检查 spanning reads                       │
│      ├─ 部分填充验证：独立判断左右两侧                          │
│      │   ├─ 左侧通过 + 右侧失败 → 保留左侧，右侧回退 500N       │
│      │   ├─ 左侧失败 + 右侧通过 → 左侧回退 500N，保留右侧       │
│      │   └─ 两侧都失败 → 完全回退为 500N gap                    │
│      └─ 2b. 如有回退 → 重新比对 reads 到 reverted assembly      │
│  3. Find remaining gaps (基于 reverted assembly 坐标)           │
│  4. [优化] Fill gaps in PARALLEL                                │
│      ├─ TIER 0: 共识优先 (跳过wtdbg2)                           │
│      └─ TIER 1-6: 常规策略                                      │
│  5. Process results:                                            │
│      ├─ 高置信度完全填充 → FILLED_COMPLETE (立即确认)           │
│      └─ 其他 → FILLED_PENDING (延迟验证)                        │
│  6. Apply fills to assembly (更新后续填充坐标)                  │
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
迭代 N ├──→ 完全填充成功 ──→ FILLED_PENDING ──→ 迭代 N+1 验证
       │                                            │
       │                                    ├─ spanning reads 足够 → FILLED_COMPLETE
       │                                    └─ 验证失败 → FAILED (完全回退 500N)
       │
       ├──→ 部分填充成功 ──→ FILLED_PENDING ──→ 迭代 N+1 独立验证
       │                                            │
       │                            ├─ 两侧都通过 → FILLED_PARTIAL
       │                            ├─ 一侧通过 → FILLED_PARTIAL (部分回退)
       │                            └─ 两侧都失败 → FAILED (完全回退 500N)
       │
       └──→ 填充失败 ──→ UNFILLABLE (永久跳过)
                    └──→ FAILED (下轮重试)
```

**部分回退逻辑：**
```
原始: [左侧翼][左填充][500N][右填充][右侧翼]

左失败+右通过: [左侧翼][500N][右填充][右侧翼]
左通过+右失败: [左侧翼][左填充][500N][右侧翼]
两侧都失败:    [左侧翼][500N][右侧翼]
```

### GapStatus 枚举

| 状态 | 说明 | 下轮操作 |
|------|------|---------|
| PENDING | 待处理 | 尝试填充 |
| NEEDS_POLISH | 侧翼需要 Polish | 预处理阶段 Polish |
| FILLED_PENDING | 已填充，待验证 | 下轮验证 |
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
```

**TIER 0 共识优先策略：**
- 当 HiFi spanning reads >= 5 且一致性 > 95% 时触发
- 直接计算共识序列（POA），跳过 wtdbg2
- 减少约 50% 的 wtdbg2 调用

**关键方法：**
- `_collect_reads_by_type()` - 分别收集 HiFi/ONT reads
- `try_consensus_fill()` - [新] 共识优先填充 (core/consensus.py)
- `_assemble_spanning_reads()` - 类型特定的 wtdbg2 参数
- `_polish_with_hifi()` - 用 racon 抛光 ONT 组装
- `_run_wtdbg2_assembly()` - HiFi 用 ccs preset，ONT 用 ont preset

### 2. 延迟验证 (core/validator.py)

**为什么需要延迟验证：**
- 填充时的 BAM 是基于旧 assembly（gap 区域是 N）
- 验证需要 reads 比对到填充后的序列
- 下一轮迭代的 BAM 基于 filled assembly，可以准确验证

**GapValidator 类：**
- `validate_complete_fill()` - 验证完全填充
- `validate_partial_fill()` - 验证部分填充（独立左右判断）
- `pre_assess_gap()` - 预评估单个 gap 侧翼
- `pre_assess_gaps()` - 批量预评估
- `analyze_failed_gap()` - 分析失败原因

**验证标准（完全填充 vs 部分填充）：**

| 类型 | 验证标准 | 说明 |
|------|----------|------|
| **完全填充** | Spanning reads ≥3-5 | 必须有 reads 跨越整个填充区域 |
| **部分填充** | Junction coverage ≥5, Insert coverage ≥5, 无断点 | 只检查连接点质量 |

**部分填充独立左右验证：**
- 左侧和右侧独立判断通过/失败
- 通过的一侧保留，失败的一侧回退为 500N
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
- `add_pending_fill()` - 添加待验证的填充
- `get_pending_fills()` - 获取待验证列表
- `should_attempt()` - 判断是否应该尝试填充
- `set_status()` - 设置 gap 状态

**PendingFill 数据结构：**
```python
@dataclass
class PendingFill:
    gap_id: str
    chrom: str
    original_start: int
    original_end: int
    filled_start: int
    filled_end: int
    sequence: str
    is_complete: bool
    source: str  # e.g., "hifi_spanning"
    tier: int
```

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
├── gap_tracker.json           # GapStatusTracker 状态 (含 pending_fills)
├── preprocessing/
│   ├── alignment.ok           # 初始比对完成
│   ├── filter.ok              # reads 过滤完成
│   └── polish.ok              # flank polish 完成
├── iteration_1/
│   ├── alignment.ok           # 比对完成
│   ├── validation.ok          # 验证完成
│   └── iteration.ok           # 本轮完成
└── complete.ok                # 全部完成
```

**关键改进：**
- 使用 `.ok` 文件而非 checkpoint.json 判断进度
- `gap_tracker.json` 单独保存，包含 pending_fills
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
    ├── assembly_reverted.fasta  # (如果有验证失败)
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
