# GapFill - Gap Filling Pipeline

## Quick Reference

```bash
# 单倍体
python -m gapfill -a assembly.fa --hifi hifi.fq --ont ont.fq -o output

# 多倍体 (自动检测)
python -m gapfill -a hap1.fa hap2.fa --hifi hifi.fq --ont ont.fq -o output

# 多倍体 + 优化模式 (减少75%比对)
python -m gapfill -a hap1.fa hap2.fa --hifi hifi.fq --optimized -o output

# 使用 Hi-C 数据
python -m gapfill -a assembly.fa --hifi hifi.fq --hic hic_R1.fq hic_R2.fq -o output

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
5. **HiFi 优先** - 6层策略优先使用高准确性 HiFi，ONT 提供长度优势
6. **多倍体先标准化** - 在 SNP 检测之前先标准化所有 haplotype 的 gap
7. **Alignment-based SNP 检测** - 使用 minimap2 比对检测 SNP，正确处理 haplotype 独有的 gap

## 包结构

```
gapfill/
├── __init__.py
├── __main__.py
├── cli.py                    # 命令行解析
├── core/
│   ├── filler.py             # GapFiller - 6层 HiFi/ONT 策略
│   ├── validator.py          # GapValidator + GapStatusTracker
│   └── polisher.py           # FlankPolisher - 侧翼序列抛光
├── engines/
│   ├── haploid.py            # HaploidEngine - 单倍体引擎
│   ├── polyploid.py          # PolyploidEngine - 多倍体引擎
│   └── optimized_polyploid.py # OptimizedPolyploidEngine - 批量比对优化
└── utils/
    ├── indexer.py            # AssemblyIndexer
    ├── scanner.py            # GapScanner
    ├── tempfiles.py          # TempFileManager
    ├── hic.py                # HiCAnalyzer - Hi-C 数据分析
    └── checkpoint.py         # CheckpointManager - 断点续跑
```

## 工作流程

### 单倍体工作流程

```
┌─────────────────────────────────────────────────────────────────┐
│  PREPROCESSING PHASE (一次性，在迭代前)                          │
├─────────────────────────────────────────────────────────────────┤
│  0a. Find initial gaps                                          │
│  0b. Prepare Hi-C data (optional)                               │
│  0c. Normalize gaps → 500N                                      │
│  0d. Initial alignment → BAM                                    │
│  0e. Pre-assess gap flanks                                      │
│      ├─ PENDING: 侧翼正常                                       │
│      └─ NEEDS_POLISH: 侧翼有问题                                │
│  0f. Polish problematic flanks → polished assembly              │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  ITERATION 1                                                     │
├─────────────────────────────────────────────────────────────────┤
│  1. Align reads → BAM                                           │
│  2. (skip validation - no previous fills)                       │
│  3. Find gaps                                                   │
│  4. Filter gaps (skip UNFILLABLE, FILLED_*)                     │
│  5. Fill gaps → mark as FILLED_PENDING                          │
│  6. Apply fills to assembly                                     │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  ITERATION 2+                                                    │
├─────────────────────────────────────────────────────────────────┤
│  1. Align reads → NEW BAM (基于上轮 filled assembly)             │
│  2. Validate previous fills (用新 BAM)                          │
│      ├─ 验证通过 → FILLED_COMPLETE/PARTIAL                      │
│      └─ 验证失败 → 回退填充，恢复 500N gap                       │
│  3. Find gaps                                                   │
│  4. Filter gaps                                                 │
│  5. Fill gaps → FILLED_PENDING                                  │
│  6. Apply fills                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Gap 状态流转

```
初始: PENDING
       │
预处理 ├──→ NEEDS_POLISH ──→ [Polish] ──→ PENDING
       │
迭代 N ├──→ 填充成功 ──→ FILLED_PENDING
       │                      │
       │              迭代 N+1 验证
       │                      │
       │              ├─ 通过 → FILLED_COMPLETE
       │              └─ 失败 → FAILED (回退 gap)
       │
       └──→ 填充失败 ──→ UNFILLABLE (永久跳过)
                    └──→ FAILED (下轮重试)
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

### 1. 6层 HiFi/ONT 分层策略 (core/filler.py)

```
TIER 1: HiFi-only spanning     → 最高准确性
TIER 2: ONT-only spanning      → 长度优势 + 可选 HiFi 抛光
TIER 3: Hybrid spanning        → 混合跨越reads
TIER 4: HiFi flanking + merge  → HiFi 侧翼组装 + 智能合并
TIER 5: ONT flanking + merge   → ONT 侧翼组装 + 可选抛光
TIER 6: Hybrid flanking + 500N → 最后备选
```

**关键方法：**
- `_collect_reads_by_type()` - 分别收集 HiFi/ONT reads
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
- `validate_partial_fill()` - 验证部分填充
- `pre_assess_gap()` - 预评估单个 gap 侧翼
- `pre_assess_gaps()` - 批量预评估
- `analyze_failed_gap()` - 分析失败原因

**验证标准：**
- Spanning reads 数量（短 gap ≥3，长 gap ≥5）
- 覆盖度连续性（零覆盖比例 < 10%）
- 平均覆盖度（≥ 3x）
- Junction 质量（clip 比例 < 30%）

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

### 4. Hi-C 数据整合 (utils/hic.py)

| 阶段 | Hi-C 作用 | 方法 |
|------|----------|------|
| **填充前** | 估计 gap 真实大小 | `estimate_gap_sizes()` |
| **填充中** | 指导 wtdbg2 -g 参数 | `gap_size_estimates` |
| **填充中** | 多候选序列选择 | `_select_best_candidate_with_hic()` |
| **填充后** | 验证填充正确性 | `validate_fills()` |
| **Phasing** | 增强 reads 分配 | `enhance_phasing()` |

**多倍体 Hi-C 合并参考比对：**
```
1. 创建合并参考 (hap1__Chr1, hap2__Chr1, hap3__Chr1, ...)
2. Hi-C 比对一次到合并参考
3. 按前缀分流 BAM 到各 haplotype
4. 每个 haplotype 获得独立的 HiCAnalyzer
```

### 5. 多倍体 SNP 检测 (engines/polyploid.py)

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

### 6. 断点续跑 (utils/checkpoint.py)

**使用方法：**
```bash
# 首次运行
python -m gapfill -a assembly.fa --hifi hifi.fq -o output

# 中断后恢复
python -m gapfill -a assembly.fa --hifi hifi.fq -o output --resume

# 清除断点，重新开始
python -m gapfill -a assembly.fa --hifi hifi.fq -o output --clear-checkpoint
```

**状态跟踪：**

| 阶段 | 单倍体 | 多倍体 | 可复用文件 |
|------|--------|--------|-----------|
| normalization | ✓ | ✓ | assembly_normalized.fasta |
| polish | ✓ | ✓ | assembly_polished.fasta |
| snp_detection | - | ✓ | snp_database.json |
| phasing | - | ✓ | phased_*_hifi.fasta |
| filling | ✓ | ✓ | iteration_N/*.bam |

## 命令行参数

```
-a, --assembly FILE(s)    # 1个=单倍体，2+个=多倍体
--hifi FILE               # HiFi reads
--ont FILE                # ONT reads
--hic R1 R2               # Hi-C paired-end reads
--hic-bam FILE            # 预比对的 Hi-C BAM
-o, --output DIR          # 输出目录
-t, --threads N           # 线程数 (default: 8)
--max-iterations N        # 最大迭代 (default: 10)
--min-gap-size N          # 最小 gap 大小 (default: 100)
--phasing METHOD          # builtin | whatshap
--no-ambiguous-reads      # 多倍体不使用 ambiguous reads
--optimized               # 使用批量比对优化 (多倍体)
--resume                  # 从断点恢复运行
--clear-checkpoint        # 清除断点，从头开始
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
| `-g` | auto | auto | 基于 reads 或 Hi-C 估计 |

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
- bwa-mem2 (Hi-C 比对)
- bcftools, whatshap (多倍体 whatshap 模式)
