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
```

## 核心原则

1. **N 长度无意义** - N 序列只是占位符，长度不代表真实 gap 大小
2. **Gap 标准化** - 所有 N 占位符标准化为 500bp，确保 spanning reads 检测正确
3. **多轮迭代** - 逐步填充，每轮重新比对 reads
4. **Supplementary 检测** - 利用 supplementary alignments 发现真正跨 gap 的 reads
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
│   ├── filler.py             # GapFiller - 6层 HiFi/ONT 策略 + Hi-C 辅助
│   └── validator.py          # GapValidator - 填充验证
├── engines/
│   ├── haploid.py            # HaploidEngine - 单倍体引擎 + Hi-C 整合
│   ├── polyploid.py          # PolyploidEngine - 多倍体引擎
│   └── optimized_polyploid.py # OptimizedPolyploidEngine - 批量比对优化
└── utils/
    ├── indexer.py            # AssemblyIndexer
    ├── scanner.py            # GapScanner
    ├── tempfiles.py          # TempFileManager
    └── hic.py                # HiCAnalyzer - Hi-C 数据分析
```

## 核心优化

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

### 2. Hi-C 数据整合 (utils/hic.py)

| 阶段 | Hi-C 作用 | 方法 |
|------|----------|------|
| **填充前** | 估计 gap 真实大小 | `estimate_gap_sizes()` |
| **填充中** | 指导 wtdbg2 -g 参数 | `gap_size_estimates` |
| **填充中** | 多候选序列选择 | `_select_best_candidate_with_hic()` |
| **填充后** | 验证填充正确性 | `validate_fills()` |
| **Phasing** | 增强 reads 分配 | `enhance_phasing()` |

**HiCAnalyzer 类：**
- `estimate_gap_sizes()` - Hi-C pairs 跨越 gap 推断真实大小
- `validate_fills()` - 检查接触图连续性
- `enhance_phasing()` - 长程信息救回 ambiguous reads
- `get_contact_matrix()` - 获取区域接触矩阵

### 3. 多倍体 SNP 检测优化 (engines/polyploid.py)

**新流程（解决 gap N 长度不一致问题）：**
```
原始流程（有bug）：
  1. SNP 检测 (逐位置比较，gap长度不同导致坐标偏移)
  2. Gap 标准化
  3. Phasing
  4. Gap filling

新流程（已修复）：
  STEP 0: Gap 标准化 ALL haplotypes (统一为 500N)
  STEP 1: Alignment-based SNP 检测 (minimap2 asm5)
  STEP 2: Phasing reads
  STEP 3: Gap filling (skip_normalization=True)
```

**Alignment-based SNP 检测：**
- 使用 `minimap2 -ax asm5` 比对 hap2→hap1, hap3→hap1, ...
- 从 CIGAR 字符串中提取 SNP（跳过 gap 区域）
- 正确处理 haplotype 独有的 gap

**关键类和方法：**
- `ReadPhaser.normalize_all_assemblies()` - 标准化所有 haplotype
- `ReadPhaser._detect_snps_alignment()` - alignment-based SNP 检测
- `ReadPhaser._identify_haplotype_specific_gaps()` - 识别单倍型独有的 gap
- `GapRegion` dataclass - gap 元数据

### 4. 批量比对优化 (engines/optimized_polyploid.py)

**原理：**
```
原始：每个 haplotype 每次迭代独立比对
  4倍体 × 10迭代 × 2类型 = 80 次比对

优化：合并所有 haplotypes，只比对一次
  10迭代 × 2类型 = 20 次比对 (减少 75%)
```

**为什么不稀释覆盖度：**
- Phased reads 带有 SNP 特征
- 比对到合并参考时，自然选择正确的 haplotype
- Primary alignment 落在正确目标上

**关键方法：**
- `_create_merged_reference()` - 合并参考 (hap1__Chr1, hap2__Chr1, ...)
- `_merge_phased_reads()` - 合并所有 phased reads
- `_split_bam_by_prefix()` - 按前缀分流 BAM

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
-v, --verbose             # 详细日志
```

## HaploidEngine 参数

```python
HaploidEngine(
    assembly_file: str,
    hifi_reads: str = None,
    ont_reads: str = None,
    hic_reads: List[str] = None,
    hic_bam: str = None,
    output_dir: str = "output",
    threads: int = 8,
    max_iterations: int = 10,
    min_gap_size: int = 100,
    min_mapq: int = 20,
    skip_normalization: bool = False  # 跳过 gap 标准化（多倍体已预处理）
)
```

## 填充策略结果

| Strategy | Source | is_complete | has_placeholder | 说明 |
|----------|--------|-------------|-----------------|------|
| spanning | hifi_spanning | True | False | HiFi 直接跨越 |
| spanning | ont_spanning | True | False | ONT 跨越 |
| spanning | ont_spanning_hifi_polished | True | False | ONT + HiFi 抛光 |
| spanning | hybrid_ccs/hybrid_ont | True | False | 混合跨越 |
| flanking_merged | *_flanking_merged | True | False | 侧翼合并成功 |
| flanking | *_flanking_500N | False | True | 侧翼 + 500N |

## 输出结构

**单倍体：**
```
output/
├── assembly_normalized.fasta
├── final_assembly.fasta
├── final_stats.json
├── hic_aligned.bam              # (如果使用 Hi-C)
└── iteration_N/
    ├── assembly_filled.fasta
    ├── hifi.bam, ont.bam
    └── work/gap_*/
```

**多倍体：**
```
output/
├── hap1_normalized.fasta        # 标准化后的 haplotype assemblies
├── hap2_normalized.fasta
├── snp_database.json            # Alignment-based SNP 数据库
├── phased_*_hifi.fasta
├── phased_*_ont.fasta
├── hap1/
│   ├── assembly_normalized.fasta  # 符号链接
│   └── final_assembly.fasta
├── hap2/
│   └── final_assembly.fasta
└── polyploid_summary.json
```

**多倍体 (优化模式)：**
```
output/
├── phase_hifi.bam, phase_ont.bam   # 一次性 phasing
├── iteration_N/
│   ├── merged_reference.fasta      # 合并参考
│   ├── merged_hifi.bam             # 单次比对！
│   ├── merged_ont.bam
│   └── hap*/                       # 分流后的各 haplotype
├── hap1_filled.fasta
├── hap2_filled.fasta
└── summary.json
```

## 依赖

**Python:**
- biopython, pysam, numpy
- pyfaidx (可选，强烈推荐)

**External:**
- minimap2, samtools, wtdbg2, wtpoa-cns
- racon (可选，用于 HiFi 抛光)
- bwa-mem2 (可选，用于 Hi-C 比对)
- bcftools, whatshap (多倍体 whatshap 模式)

## 外部工具参数

### wtdbg2 (gap 组装)

| 参数 | HiFi | ONT | 说明 |
|------|------|-----|------|
| `-x` | ccs | ont | 读取类型 preset |
| `-L` | 1000 | 2000 | 最小读长 (减少噪音) |
| `-e` | 2 | 2 | 边缘覆盖度阈值 |
| `-S` | 1 | 1 | **关键**: 救回低覆盖度边缘 |
| `-g` | auto | auto | 基于 reads 或 Hi-C 估计 |

**Genome size 估计逻辑：**
```
1. 优先使用 Hi-C 估计值
2. 否则：avg_read_len * 0.8 (如果 reads >= 5)
3. 否则：total_bases / 10 (保守估计)
```

### minimap2 (比对)

| 用途 | Preset | 额外参数 | 说明 |
|------|--------|----------|------|
| HiFi → Assembly | map-hifi | - | HiFi 读取比对 |
| ONT → Assembly | map-ont | - | ONT 读取比对 |
| Hap → Hap (SNP) | asm5 | - | 同种 haplotype 比对 (~0.1% div) |
| Phasing | map-hifi/ont | --secondary=no | 唯一分配用于 phasing |

### bwa-mem2 (Hi-C)

```
bwa-mem2 mem -5SPM -t {threads}
```

| 参数 | 说明 |
|------|------|
| `-5` | split alignments 标记较短为 secondary |
| `-S` | 跳过 mate rescue |
| `-P` | 跳过 pairing (Hi-C 特性) |
| `-M` | Picard 兼容性 |

### samtools sort

```
samtools sort -@ {threads} -m 8G -o {output} -
```

- `-m 8G`: 每线程内存限制，防止大文件 OOM

## 比对次数对比 (多倍体)

| 倍性 | 迭代 | 原始 | 优化 | 减少 |
|------|------|------|------|------|
| 2n | 10 | 40 | 20 | 50% |
| 4n | 10 | 80 | 20 | 75% |
| 6n | 10 | 120 | 20 | 83% |
| 8n | 10 | 160 | 20 | 87.5% |
