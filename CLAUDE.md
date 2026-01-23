# GapFill - Gap Filling Pipeline

## Quick Reference

```bash
# 单倍体
python -m gapfill -a assembly.fa --hifi reads.fq.gz -o output

# 多倍体 (自动检测)
python -m gapfill -a hap1.fa hap2.fa --hifi reads.fq.gz -o output
```

## 核心原则

1. **N 长度无意义** - N 序列只是占位符，长度不代表真实 gap 大小
2. **Gap 标准化** - 所有 N 占位符标准化为 500bp，确保 spanning reads 检测正确
3. **多轮迭代** - 逐步填充，每轮重新比对 reads
4. **Supplementary 检测** - 利用 supplementary alignments 发现真正跨 gap 的 reads

## 包结构

```
gapfill/
├── __init__.py           # 导出: HaploidEngine, PolyploidEngine, GapFiller, GapValidator
├── __main__.py           # python -m gapfill 入口
├── cli.py                # 命令行解析，自动判断单/多倍体模式
├── core/
│   ├── filler.py         # GapFiller - 核心填充逻辑 (三步策略)
│   └── validator.py      # GapValidator - 填充验证
├── engines/
│   ├── haploid.py        # HaploidEngine - 单倍体迭代引擎
│   └── polyploid.py      # PolyploidEngine + ReadPhaser - 多倍体引擎
└── utils/
    ├── indexer.py        # AssemblyIndexer - pyfaidx 快速序列访问
    ├── scanner.py        # GapScanner - 扫描 N-runs
    └── tempfiles.py      # TempFileManager - 临时文件管理
```

## 核心类详解

### GapFiller (core/filler.py)

**三步填充策略：**

```python
def fill_gap(self, gap: Dict) -> Dict:
    # Step 1: 尝试 spanning reads (直接 + supplementary-linked)
    result = self._try_spanning_reads_assembly(...)
    if result['success']: return result  # is_complete=True

    # Step 2: 尝试 flanking reads + 智能合并
    result = self._try_flanking_reads_with_merge(...)
    if result['success']: return result  # is_complete=True/False

    # Step 3: 失败
    return {'success': False, ...}
```

**关键方法：**
- `_get_direct_spanning_reads()` - 直接跨越 gap 的 reads
- `_get_supplementary_spanning_reads()` - 利用 supplementary 发现跨 gap reads
- `_try_merge_sequences()` - minimap2 asm5 检测 left/right 重叠并合并
- `_get_flanking_reads()` - 收集两侧 flanking reads (包含 supplementary)
- `_run_wtdbg2_assembly()` - wtdbg2 局部组装

**Supplementary 检测原理：**
```
当 read 真正跨越 gap 时，minimap2 产生:
- Primary alignment: 结束于 gap 左边界
- Supplementary alignment: 开始于 gap 右边界
→ 同一个 read_name 出现在 gap 两侧 = 真正的 spanning read

代码逻辑:
left_reads[read.query_name] = {...}   # 结束于 gap_start
right_reads[read.query_name] = {...}  # 开始于 gap_end
common = set(left_reads) & set(right_reads)  # 真正跨越的 reads
```

### HaploidEngine (engines/haploid.py)

**工作流程：**
```
STEP 0: _normalize_gaps() → 所有 gap 标准化为 500N
   ↓
ITERATION LOOP:
  Step 1: _align_reads() → minimap2 + samtools
  Step 2: _find_gaps() → regex 扫描 N+
  Step 3: GapFiller.fill_gap() → 三步策略
  Step 4: _apply_fills() → 更新 assembly
  Step 5: 检查进度，决定是否继续
   ↓
OUTPUT: final_assembly.fasta, final_stats.json
```

**关键方法：**
- `_normalize_gaps()` - 将所有 N 占位符替换为 500N
- `_align_reads()` - minimap2 比对 (map-hifi/map-ont)
- `_find_gaps()` - regex 扫描 `N+` (min_gap_size 过滤)
- `_apply_fills()` - 从后向前应用填充 (避免坐标偏移)

### PolyploidEngine (engines/polyploid.py)

**工作流程：**
```
STEP 1: 比对 reads 到参考单倍型
STEP 2: ReadPhaser.detect_haplotype_snps() → 检测单倍型特异 SNPs
STEP 3: ReadPhaser.phase_reads() → 根据 SNP 分配 reads 到单倍型
STEP 4: 对每个单倍型独立运行 HaploidEngine
STEP 5: 输出多单倍型填充结果
```

**ReadPhaser 类：**
- `_detect_snps_builtin()` - 比较单倍型序列找 SNPs
- `_detect_snps_whatshap()` - 使用 WhatsHap (需要 bcftools)
- `phase_reads()` - 根据 SNP profile 分配 reads
- `_assign_read_to_haplotype()` - 计算 hap_scores，选最高分

### GapValidator (core/validator.py)

**验证逻辑：**
- 检查填充区域的 read 覆盖度
- 检查 junction 处覆盖度是否连续
- 对于 partial fill，只验证 left junction (right 连接 N's)

**关键类：**
- `BAMPool` - BAM 句柄池，避免重复打开

### AssemblyIndexer (utils/indexer.py)

**优化：**
- pyfaidx 快速序列访问 (10-1000x faster)
- 全局缓存，跨模块共享 index
- 自动 fallback 到 BioPython SeqIO

## 填充策略结果

| Strategy | is_complete | has_placeholder | 说明 |
|----------|-------------|-----------------|------|
| `spanning` | True | False | 直接或 supplementary 跨越 |
| `flanking_merged` | True | False | left+right 成功合并 |
| `flanking` | False | True | left + 500N + right |

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
--no-ambiguous-reads      # 多倍体模式不使用 ambiguous reads
-v, --verbose             # 详细日志
```

## 输出结构

**单倍体：**
```
output/
├── assembly_normalized.fasta   # 标准化后的 assembly
├── final_assembly.fasta        # 最终填充结果
├── final_stats.json            # 统计信息
└── iteration_N/
    ├── assembly_filled.fasta
    ├── hifi.bam, ont.bam
    └── work/gap_*/             # 每个 gap 的工作文件
```

**多倍体：**
```
output/
├── snp_database.json           # SNP 数据库
├── phased_hap1_reads.fasta     # 分型后的 reads
├── hap1/final_assembly.fasta   # 各单倍型结果
├── hap2/final_assembly.fasta
└── polyploid_summary.json      # 汇总报告
```

## 依赖

**Python:**
- biopython, pysam, numpy
- pyfaidx (可选，强烈推荐)

**External:**
- minimap2, samtools, wtdbg2
- bcftools, whatshap (多倍体 whatshap 模式)
