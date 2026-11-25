# GFFmerge

GFFmerge 是一个独立的基因结构注释合并工具。它将多套 GFF3/GTF 按基因组座位和 CDS 剪接链聚类，结合同源蛋白和二代 RNA-seq 选择主代表；有 Iso-seq 时，还可进行局部 UTR 共识修正并输出可变剪切。

## 默认输入

常规运行只需要：

- 一套或多套 GFF3/GTF；
- 基因组 FASTA；
- 二代 RNA-seq BAM；
- 同源蛋白 FASTA。

二代 BAM 会自动用于 junction 提取、CDS 覆盖度和平均深度计算。蛋白 FASTA 会自动完成基因组比对并参与模型裁决，不需要准备中间文件。

**Iso-seq 是可选输入，不是运行前提。** 没有 Iso-seq 时，GFFmerge 正常完成合并、蛋白/RNA-seq 裁决和主代表输出。

## 安装

```bash
conda env create -f environment.yml
conda activate gffmerge
```

环境包含 Python、`pysam` 和 `miniprot`，不需要额外安装 Python 包。

## 快速开始

```bash
python3 GFFmerge.py \
  --set setA=A.gff3 \
  --set setB=B.gff3 \
  --set setC=C.gff3 \
  --genome genome.fa \
  --bam RNA.sorted.bam \
  --proteins related_proteins.faa \
  --threads 16 \
  --outdir GFFmerge_out
```

`--set` 和 `--bam` 都可以重复使用。多套 BAM 的证据会合并计算。

有 Iso-seq 时，在默认命令上增加：

```bash
  --isoseq isoseq.bam \
  --utr-consensus
```

没有 Iso-seq 时不需要这两个参数，也不影响主流程。

## 依赖

- Python 3.8 或更高版本；
- 读取 BAM 时需要 `pysam`；
- 蛋白 FASTA 自动比对时需要 `miniprot`；
- 其他输入处理使用 Python 标准库。

## 推荐比对参数

以下 STAR 和 minimap2 是可选的上游步骤；如果已有符合要求的 BAM，可直接运行 GFFmerge。

### 二代 RNA-seq：STAR

```bash
STAR \
  --genomeDir star_index \
  --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
  --readFilesCommand zcat \
  --runThreadN 16 \
  --twopassMode Basic \
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --outFilterMultimapNmax 200 \
  --outFilterMultimapScoreRange 50 \
  --outFilterMismatchNmax 50 \
  --outSAMattributes NH HI AS nM NM MD jM jI XS MC \
  --outSAMprimaryFlag AllBestScore \
  --outSAMtype BAM SortedByCoordinate \
  --genomeLoad NoSharedMemory

samtools index -c Aligned.sortedByCoord.out.bam
```

然后使用：

```bash
--bam Aligned.sortedByCoord.out.bam
```

关键参数说明：

- `NH` 用于区分唯一比对和多重比对 reads；唯一比对用于 CDS 深度，junction 则分别累计 unique/multi 支持。
- `XS` 提供转录链方向；缺少或不可靠时，GFFmerge 会根据基因组剪接 motif 补充判定。
- `--outSAMtype BAM SortedByCoordinate` 是必需的，因为程序按基因区间查询 BAM。
- `samtools index -c` 生成 CSI 索引，适合长 scaffold/contig；普通染色体规模也可使用默认 BAI 索引。
- `--twopassMode Basic` 有助于发现未注释 junction，适合给基因结构纠错提供证据。

### 可选 Iso-seq：minimap2

HiFi/高质量 Iso-seq 推荐：

```bash
minimap2 \
  -ax splice:hq \
  -ub \
  --secondary=no \
  -t 16 \
  genome.fa isoseq.hifi.fastq.gz \
  | samtools sort -@ 8 -o isoseq.bam -

samtools index -c isoseq.bam
```

较低质量的三代 reads 可将 `-x splice:hq` 改为 `-x splice`。

其中 `-ub` 很重要：不要替换为 `-uf`，否则可能使负链转录本的链向判定产生偏差。GFFmerge 会进一步用参考基因组上的 GT-AG、GC-AG、AT-AC 及其反向互补 motif 校验链向。

有 Iso-seq 时才增加：

```bash
--isoseq isoseq.bam --utr-consensus
```

### 蛋白比对

不需要手动准备蛋白比对结果。提供 `--proteins related_proteins.faa` 后，程序自动以每个蛋白保留一个最优基因组命中的方式运行，并将结果用于 CDS 覆盖度计算和模型裁决；`--threads` 控制比对线程数。

## 保留与选优规则

### 多套注释

同一 CDS 剪接链得到至少两个输入集合支持时，该链可以保留。每个链组先选最高评分模型，再解决同座位、同链的 CDS 冲突。

### 重复基因与同座位冲突

GFFmerge 使用以下规则识别重复基因和同座位冲突：

1. 多外显子模型要求完整 CDS 内含子链一致，才归为同一结构组。
2. 单外显子模型没有内含子可比，CDS 重叠达到较短模型的 90% 时归为同组。
3. 不同结构组若位于同一座位、同一转录链，且 CDS 双向重叠任一比例超过 30%，判为重复/冲突模型；保留证据评分最高者，其余写入 `dropped.gff3`，原因为 `locus_conflict_lower_score`。
4. 反义重叠基因不会因这一条规则互相删除，低重叠的相邻基因也可分别保留。

对应的核心 Python 片段如下，`blocks_coverage` 计算两组 CDS 区块的碱基重叠数：

```python
def reciprocal_cds_overlap(a, b):
    value = blocks_coverage(a.cds, b.cds)
    return value / max(a.cds_len, 1), value / max(b.cds_len, 1)

# 同链且任一方向 CDS 重叠超过阈值时，低分模型进入 dropped
if best.strand == other.strand and \
        max(reciprocal_cds_overlap(best, other)) > args.conflict_overlap:
    dropped[best.key] = "locus_conflict_lower_score"
```

双向比例是必要的：一个短 CDS 可能几乎完全落入长 CDS，但只占长 CDS 的一小部分。任一方向超过阈值即冲突，可避免同一基因的长短版本同时进入最终注释。阈值可用 `--conflict-overlap` 调整，默认 `0.3`。

### 串联重复与二联体

串联重复是基因组中的真实组织结构，不能与同一座位的冗余注释混为一谈。GFFmerge 在完成主代表选择后，对最终基因按同染色体、同链和相邻位置进行聚类，**只标记，不删除串联成员**。

相邻基因默认需同时满足：

- 基因组间距不超过 20 kb；
- CDS 长度比不超过 1.8；
- CDS 区块数差不超过 1。

相似关系可以沿相邻基因传递，形成一个重复簇。默认至少 3 个成员标为 `tandem`，恰好 2 个成员标为 `duplicate`。核心判据为：

```python
gap = max(0, model.cds[0][0] - previous.cds[-1][1] - 1)
length_ratio = max(model.cds_len, previous.cds_len) / min(model.cds_len, previous.cds_len)
similar = (
    gap <= args.tandem_max_gap
    and length_ratio <= args.tandem_cds_ratio
    and abs(len(model.cds) - len(previous.cds)) <= 1
)
```

最终 GFF 的基因行会增加：

```text
duplicate_type=tandem
duplicate_cluster=TANDEM000001
duplicate_size=3
```

`decision_report.tsv` 同样记录这三个字段，`summary.json` 提供 `tandem_clusters`、`tandem_genes`、`duplicate_clusters` 和 `duplicate_genes` 统计。可用以下参数调整：

```text
--tandem-max-gap 20000
--tandem-cds-ratio 1.8
--tandem-min-members 3
```

### 单套注释

非可信集合中的单套模型满足以下任一条件即可保留：

1. 有达到阈值的二代 junction、CDS 深度或 Iso-seq 完整链支持；
2. 有合格的同源蛋白支持；
3. 输入集合通过 `--trusted-set NAME` 指定为可信集合。

蛋白支持的默认门槛：

- 命中与模型同链，`Identity ≥ 0.3`；
- **同一条蛋白命中**对模型 CDS 碱基覆盖比例 ≥0.5；
- 命中区间对查询蛋白长度覆盖比例 ≥0.5；
- 模型没有内部终止子。

可通过以下参数调整：

```text
--min-protein-identity 0.3
--min-protein-coverage 0.5
--min-protein-query-coverage 0.5
```

蛋白证据还会按模型 CDS 覆盖比例加入代表评分。多条不同蛋白的局部命中不会拼接成“完整覆盖”。蛋白比对用于评估已有注释，不自动创建新基因。

### 二代深度

二代 BAM 的整链 CDS 平均深度只计一次，避免多内含子模型重复获得深度奖励。代表选择同时考虑 junction、覆盖度、蛋白和结构完整性，不等同于选择最长 CDS。

### 可选 Iso-seq 增强

提供 `--isoseq` 后，完整长读链会独立参与评分。

`--utr-consensus` 会对 CDS 兼容的长读进行局部链聚类。主导链达到分子数和主导度阈值后，可调整 UTR 外显子和末端；只有共识外显子完整覆盖 CDS 时才应用，CDS 不变。

有 Iso-seq 时，每基因默认最多输出 3 条转录本：1 条主代表和最多 2 条替代链。替代链须有不少于 2 条 Iso-seq 分子精确支持完整 CDS 内含子链，且 CDS 覆盖 ≥90%。蛋白或二代证据不能替代这一要求。`--max-isoforms 1` 可关闭替代链输出。

## 常用参数

```text
--set NAME=GFF              输入集合，可重复
--genome GENOME.fa          基因组 FASTA，必需
--bam RNA.bam               二代 BAM，可重复
--proteins PROTEINS.faa     同源蛋白 FASTA
--isoseq ISO.bam            可选 Iso-seq BAM 或 BED12，可重复
--utr-consensus             可选 Iso-seq UTR 共识修正
--threads 16                蛋白比对线程数
--set-weight NAME=FLOAT     覆盖集合默认权重
--trusted-set NAME          允许可信集合的单套模型直接保留
--max-isoforms 3            每基因总转录本上限，范围 1–3
--id-prefix GFFM            输出基因 ID 前缀
--outdir DIR                输出目录，必需
```

蛋白比对结果会在输出目录缓存并自动复用。输入蛋白更新后，建议使用新的输出目录重新运行。

## 输出

| 文件 | 内容 |
|---|---|
| `final.gff3` | 最终基因、mRNA、exon、CDS 和 UTR，可载入 IGV |
| `final.proteins.faa` | 按最终 CDS 翻译的全部转录本蛋白 |
| `decision_report.tsv` | 来源、评分、RNA 深度、Iso-seq 分子数及蛋白覆盖指标 |
| `dropped.gff3` | 未保留模型及原因 |
| `summary.json` | 基因、转录本、蛋白保留和 UTR 修正统计 |

`decision_report.tsv` 中 `reason=unique_with_protein` 表示该主代表在没有足够 RNA 证据的情况下由蛋白支持保留。`protein_coverage` 是命中对模型 CDS 的覆盖，`protein_query_coverage` 是命中对查询蛋白的覆盖。
