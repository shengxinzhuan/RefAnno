# RefAnno

**RefAnno: multi-evidence gene structure refinement and representative annotation.**

RefAnno 是一个独立的基因结构注释修订工具。它整合多套 GFF3/GTF、二代 RNA-seq、同源蛋白和可选 Iso-seq，分别完成结构评估、基因保留和代表转录本选择。目标不是拼接最长模型，而是得到结构可信、证据来源清楚、表达具有代表性的基因注释。

RefAnno 评估和修订输入模型，不从零创建新基因。局部共识目前用于可选的 UTR 修订，不改变 CDS；蛋白证据用于评估和保留模型，也不自动创建新基因。

## 默认输入

常规运行只需要：

- 一套或多套 GFF3/GTF；
- 基因组 FASTA；
- 二代 RNA-seq BAM；
- 同源蛋白 FASTA。

二代 BAM 会自动用于 junction 提取、CDS 覆盖度和平均深度计算。蛋白 FASTA 会自动完成基因组比对并参与模型裁决，不需要准备中间文件。

**Iso-seq 是可选输入，不是运行前提。** 没有 Iso-seq 时，RefAnno 正常完成注释整合、蛋白/RNA-seq 裁决和主代表输出。

## 安装

```bash
conda env create -f environment.yml
conda activate refanno
```

环境包含 Python、`pysam` 和 `miniprot`，不需要额外安装 Python 包。

## 快速开始

```bash
python3 RefAnno.py \
  --set setA=A.gff3 \
  --set setB=B.gff3 \
  --set setC=C.gff3 \
  --genome genome.fa \
  --bam RNA.sorted.bam \
  --proteins related_proteins.faa \
  --threads 16 \
  --outdir RefAnno_out
```

`--set` 和 `--bam` 都可以重复使用。多套 BAM 的证据会合并计算。

有 Iso-seq 时，在默认命令上增加：

```bash
  --isoseq isoseq.bam \
  --utr-consensus
```

没有 Iso-seq 时不需要这两个参数，也不影响主流程。

## 三层裁决

RefAnno 不再把蛋白、RNA、外显子数和 UTR 混入一个无法解释的总分，而是分别回答三个问题。

### 1. 结构是否可信

程序逐个评估 CDS 内含子和转录本完整链：

- 二代 unique junction 达到 `--min-sj`，或 Iso-seq junction 达到 `--min-isoform-molecules`，该 junction 记为已支持。
- `junctions_supported / junctions_total` 给出已支持 CDS junction 的比例。
- `structure_coverage` 多外显子模型为已支持 junction 比例，单外显子模型为 CDS 覆盖度。
- `canonical_junctions` 记录符合基因组剪接 motif 的 junction 数；它与 reads 支持分开报告。
- `complete_cds` 检查完整密码子框架、起始和终止密码子以及内部终止子。

`structure_status` 的含义：

| 状态 | 含义 |
|---|---|
| `exact_isoseq` | 完整 CDS 剪接链达到 Iso-seq 分子数要求 |
| `supported` | 全部 CDS junction 已支持，或单外显子模型达到覆盖度和深度要求 |
| `partial` | 有部分 RNA 证据，但完整结构未全部得到支持 |
| `unsupported` | 有 RNA 数据，但该模型未得到结构支持 |
| `untested` | 没有 RNA 数据，无法评价转录结构 |

结构状态用于主代表排序，不会单独删除基因。一个模型可以结构仍为 `partial`，但因多集合一致或同源蛋白支持而保留，并在结果中明确标记。

### 2. 基因是否保留

多套注释中，同一 CDS 剪接链得到至少两个输入集合支持时可以保留。非可信集合中的单套模型满足以下任一条件也可保留：

1. 有达到阈值的二代 junction、CDS 深度或 Iso-seq 完整链支持；
2. 有合格的同源蛋白支持；
3. 输入集合通过 `--trusted-set NAME` 指定为可信集合；
4. 符合重复家族的间接证据保留规则。

`egapx`、`eviann` 和 `e` 默认可信，其他集合可用 `--trusted-set` 加入。可信集合只影响保留，不会直接压过具有实验证据的冲突模型。

蛋白支持的默认门槛：

- 命中与模型同链，`Identity ≥ 0.3`；
- 同一条蛋白命中对模型 CDS 碱基覆盖比例 ≥0.5；
- 命中区间对查询蛋白长度覆盖比例 ≥0.5；
- 模型没有内部终止子。

可通过以下参数调整：

```text
--min-protein-identity 0.3
--min-protein-coverage 0.5
--min-protein-query-coverage 0.5
```

多条不同蛋白的局部命中不会拼接成“完整覆盖”。蛋白证据可以保留一个缺乏 RNA 支持的基因，但不能替错误 junction 提供结构支持。

### 3. 哪条是代表转录本

对同座位的冲突结构，RefAnno 优先使用直接实验证据，再考虑间接证据：

1. Iso-seq 完整链；
2. 二代 RNA 支持；
3. 直接蛋白支持；
4. 重复家族支持；
5. 多集合结构一致；
6. 可信输入来源。

同一证据层级内，依次比较 CDS 完整性、结构覆盖、Iso-seq 分子数、CDS 深度和覆盖度。蛋白保留通道使用蛋白 CDS 覆盖、查询蛋白覆盖和一致性；输入集合数量和权重只用于后续排序。

**外显子数、UTR 长度和 CDS 长度不再额外加分，CDS 长度也不作为平分时的生物学偏好。** 最终仍用模型键排序只是为了让结果可重复。

输出中的两个字段分别记录不同决定：

- `reason`：为什么保留这个基因，如 `unique_with_evidence`、`unique_with_protein` 或 `2set_chain_consensus`。
- `selection`：为什么选它作主代表，如 `rna_structure_depth`、`isoseq_chain_abundance` 或 `protein_coverage`。

## 重复基因与同座位冲突

RefAnno 使用以下规则识别重复基因和同座位冲突：

1. 多外显子模型要求完整 CDS 内含子链一致，才归为同一结构组。
2. 单外显子模型没有内含子可比，CDS 重叠达到较短模型的 90% 时归为同组。
3. 不同结构组若位于同一座位、同一转录链，且 CDS 双向重叠任一比例超过 30%，判为重复/冲突模型；按上述代表规则保留优先者，其余写入 `dropped.gff3`，原因为 `locus_conflict_lower_rank`。
4. 反义重叠基因不会因这一条规则互相删除，低重叠的相邻基因也可分别保留。

双向 CDS 重叠计算如下：

```python
def reciprocal_cds_overlap(a, b):
    value = blocks_coverage(a.cds, b.cds)
    return value / max(a.cds_len, 1), value / max(b.cds_len, 1)

if best.strand == other.strand and \
        max(reciprocal_cds_overlap(best, other)) > args.conflict_overlap:
    dropped[best.key] = "locus_conflict_lower_rank"
```

双向比例是必要的：一个短 CDS 可能几乎完全落入长 CDS，但只占长 CDS 的一小部分。任一方向超过阈值即冲突，可避免同一基因的长短版本同时进入最终注释。阈值可用 `--conflict-overlap` 调整，默认 `0.3`。

## 串联重复与二联体

串联重复是基因组中的真实组织结构，不能与同一座位的冗余注释混为一谈。RefAnno 在完成主代表选择后，对最终基因按同染色体、同链和相邻位置进行聚类，**只标记，不删除串联成员**。

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
--tandem-min-multi-depth 5
```

这些条件识别的是位置和结构相似的**候选重复簇**，不等同于经过序列系统发育验证的同源家族。

### 重复家族参与的保留裁决

重复信息不仅用于输出标签，也会在主代表选择前参与候选模型的保留判断。处理原则是：**结构相似的独立拷贝可以共享家族级证据，但不能共享精确到某个拷贝的证据。**

- unique read 深度、unique junction 和 Iso-seq 完整链仍只归属于实际匹配的模型，不在家族成员之间传递。
- 任一成员有合格蛋白支持时，可形成家族级蛋白证据。这对每个蛋白只保留一个最优命中的情况尤其有用，其他拷贝不会被误判为完全无同源支持。
- 多重比对 CDS 覆盖 ≥0.5 且平均深度 ≥5 时，形成 `multi_depth` 家族证据。默认阈值可用 `--tandem-min-multi-depth` 调整。
- 多重比对 junction 支持达到普通 junction 阈值的 3 倍时，形成 `multi_junction` 家族证据。
- 单套注释中的候选重复成员可凭上述家族证据保留，原因是 `tandem_family_support`，并同时进入复核标签，避免与精确 unique 支持混淆。
- 家族证据不增加直接证据等级，不改变同座位直接证据模型优先的规则，也不救援含内部终止子的模型。

**二联体 `duplicate` 与多成员串联簇 `tandem` 享有相同的家族救援规则。** 一个成员有合格蛋白或家族级多重比对证据时，另一个成员也可保留。两类关系统一使用 `tandem_family_support` 作为保留原因。

输出中的 `family_support` 记录具体来源，如 `protein`、`multi_depth`、`multi_junction` 或其组合；`cds_multi_depth` 和 `cds_multi_coverage` 记录模型自身的多重比对覆盖。`summary.json` 的 `tandem_family_retained_genes` 同时统计二联体和多成员串联簇。

## 二代深度

二代 BAM 的整链 CDS 平均深度只计一次，避免多内含子模型重复获得深度奖励。深度用于代表选择，不与蛋白支持混合成一个保留总分。unique/multi 深度和 junction 分开统计，多重比对证据不会伪装成某个拷贝的精确支持。

## 可选 Iso-seq 增强

提供 `--isoseq` 后，完整长读链独立参与结构和代表裁决。

`--utr-consensus` 会对 CDS 兼容的长读进行局部链聚类。主导链达到分子数和主导度阈值后，可调整 UTR 外显子和末端；只有共识外显子完整覆盖 CDS 时才应用，CDS 不变。

有 Iso-seq 时，每基因默认最多输出 3 条转录本：1 条主代表和最多 2 条替代链。替代链须有不少于 2 条 Iso-seq 分子精确支持完整 CDS 内含子链，且 CDS 覆盖 ≥90%。蛋白或二代证据不能替代这一要求。`--max-isoforms 1` 可关闭替代链输出。

## 依赖

- Python 3.8 或更高版本；
- 读取 BAM 时需要 `pysam`；
- 蛋白 FASTA 自动比对时需要 `miniprot`；
- 其他输入处理使用 Python 标准库。

## 推荐比对参数

以下 STAR 和 minimap2 是可选的上游步骤；如果已有符合要求的 BAM，可直接运行 RefAnno。

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
- `XS` 提供转录链方向；缺少或不可靠时，RefAnno 会根据基因组剪接 motif 补充判定。
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

其中 `-ub` 很重要：不要替换为 `-uf`，否则可能使负链转录本的链向判定产生偏差。RefAnno 会进一步用参考基因组上的 GT-AG、GC-AG、AT-AC 及其反向互补 motif 校验链向。

有 Iso-seq 时才增加：

```bash
--isoseq isoseq.bam --utr-consensus
```

### 蛋白比对

不需要手动准备蛋白比对结果。提供 `--proteins related_proteins.faa` 后，程序自动以每个蛋白保留一个最优基因组命中的方式运行，并将结果用于 CDS 覆盖度计算和模型裁决；`--threads` 控制比对线程数。

## 常用参数

```text
--set NAME=GFF              输入集合，可重复
--genome GENOME.fa          基因组 FASTA，必需
--bam RNA.bam               二代 BAM，可重复
--proteins PROTEINS.faa     同源蛋白 FASTA
--isoseq ISO.bam            可选 Iso-seq BAM 或 BED12，可重复
--utr-consensus             可选 Iso-seq UTR 共识修正
--threads 16                蛋白比对线程数
--set-weight NAME=FLOAT     覆盖集合默认权重，仅用于后续代表排序
--trusted-set NAME          允许可信集合的单套模型直接保留
--max-isoforms 3            每基因总转录本上限，范围 1–3
--id-prefix REFA            输出基因 ID 前缀
--outdir DIR                输出目录，必需
```

蛋白比对结果会在输出目录缓存并自动复用。输入蛋白更新后，建议使用新的输出目录重新运行。

## 输出

| 文件 | 内容 |
|---|---|
| `final.gff3` | 最终基因、mRNA、exon、CDS 和 UTR，可载入 IGV |
| `final.proteins.faa` | 按最终 CDS 翻译的全部转录本蛋白 |
| `decision_report.tsv` | 结构状态、保留原因、代表依据及各类证据指标 |
| `dropped.gff3` | 未保留模型及原因 |
| `summary.json` | 基因、转录本、结构状态、蛋白保留和 UTR 修正统计 |

`decision_report.tsv` 不再提供混合 `score`。主要字段包括：

- `structure_status`、`structure_coverage`、`junctions_supported`、`junctions_total`：结构可信程度；
- `reason`：基因保留原因；
- `selection`：主代表选择依据；
- `cds_depth`、`cds_coverage`、`isoform_molecules`：转录证据；
- `protein_target`、`protein_identity`、`protein_coverage`、`protein_query_coverage`：同源蛋白证据；
- `family_support`、`duplicate_type/cluster/size`：重复家族证据与标记；
- `review`：内部终止子、家族间接保留或 UTR 共识不确定等复核标签。
