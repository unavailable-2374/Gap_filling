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

### 3. 批量比对优化 (engines/optimized_polyploid.py)

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
├── snp_database.json
├── phased_*_hifi.fasta
├── phased_*_ont.fasta
├── hap1/final_assembly.fasta
├── hap2/final_assembly.fasta
└── summary.json
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

## 比对次数对比 (多倍体)

| 倍性 | 迭代 | 原始 | 优化 | 减少 |
|------|------|------|------|------|
| 2n | 10 | 40 | 20 | 50% |
| 4n | 10 | 80 | 20 | 75% |
| 6n | 10 | 120 | 20 | 83% |
| 8n | 10 | 160 | 20 | 87.5% |
