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
python -m pip install --no-deps .
refanno --version
```

环境包含 Python、`pysam` 和 `miniprot`。最后一步安装 RefAnno 包和 `refanno` 命令，`--no-deps` 避免 pip 重复安装 conda 已管理的依赖。

## 快速开始

```bash
refanno \
  --set setA=A.gff3 \
  --set setB=B.gff3 \
  --set setC=C.gff3 \
  --genome genome.fa \
  --bam RNA.sorted.bam \
  --proteins related_proteins.faa \
  --threads 16 \
  --outdir RefAnno_out
```

`--set` 和 `--bam` 都可以重复使用。多套 BAM 的证据会合并计算。也可以不安装包，直接在源码目录使用 `python3 RefAnno.py` 或 `python3 -m refanno`，三者调用同一套实现。

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
| `sequence_supported` | 没有 RNA 数据，但全部 CDS 内含子的供体/受体位点通过 minisplice 序列校验 |
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

**`--proteins` 应使用独立来源的蛋白**（近缘物种、Swiss-Prot 等），不要用同一基因组自身或其另一单倍型预测得到的蛋白。后者与输入模型完全相同，蛋白支持会退化为"预测器自己确认自己"，无法区分真实基因和假基因/转座子 ORF。若蛋白命中中一致性为 1.0 的比例超过 30%，RefAnno 会在运行时给出警告。

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
--tandem-min-protein-similarity 0.7
```

这些条件识别的是位置和结构相似的**候选重复簇**，不等同于经过序列系统发育验证的同源家族。

### 重复家族参与的保留裁决

重复信息不仅用于输出标签，也会在主代表选择前参与候选模型的保留判断。处理原则是：**真正相关的独立拷贝可以共享家族级证据，但不能共享精确到某个拷贝的证据。**

位置和结构相似只用于形成候选簇。跨拷贝共享证据前，两个成员还须具有相同的蛋白靶标，或蛋白序列相似度达到 `--tandem-min-protein-similarity`（默认 0.7）。序列相似度使用蛋白全局编辑距离按较长序列归一化；没有相关成员的模型不能形成家族救援。

- unique read 深度、unique junction 和 Iso-seq 完整链仍只归属于实际匹配的模型，不在家族成员之间传递。
- 相关成员有合格蛋白支持时，可形成家族级蛋白证据。这对每个蛋白只保留一个最优命中的情况尤其有用，其他拷贝不会被误判为完全无同源支持。
- 多重比对 CDS 覆盖 ≥0.5 且平均深度 ≥5 时，形成 `multi_depth` 家族证据。默认阈值可用 `--tandem-min-multi-depth` 调整。
- 多重比对 junction 支持达到普通 junction 阈值的 3 倍时，形成 `multi_junction` 家族证据。
- 单套注释中的候选重复成员可凭上述家族证据保留，原因是 `tandem_family_support`，并同时进入复核标签，避免与精确 unique 支持混淆。
- 家族证据不增加直接证据等级，不改变同座位直接证据模型优先的规则，也不救援含内部终止子的模型。

**二联体 `duplicate` 与多成员串联簇 `tandem` 享有相同的家族救援规则。** 一个成员有合格蛋白或家族级多重比对证据时，另一个成员也可保留。两类关系统一使用 `tandem_family_support` 作为保留原因。

输出中的 `family_support` 记录具体来源，如 `protein`、`multi_depth`、`multi_junction` 或其组合；`family_donors` 和 `family_rescues.tsv` 进一步记录证据来自哪个模型、哪个最终基因和哪个重复簇。`cds_multi_depth` 和 `cds_multi_coverage` 记录模型自身的多重比对覆盖。`summary.json` 的 `tandem_family_retained_genes` 同时统计二联体和多成员串联簇。

## 二代深度

RefAnno 对每个二代 BAM 顺序扫描一遍，在同一趟里提取 junction 并按 CDS 块累计深度，不按基因做随机查询，因此在大基因组和深度较高的 BAM 上耗时与 BAM 总大小成正比，而不随模型数量倍增。整链 CDS 平均深度只计一次，避免多内含子模型重复获得深度奖励。深度用于代表选择，不与蛋白支持混合成一个保留总分。unique/multi 深度和 junction 分开统计，多重比对证据不会伪装成某个拷贝的精确支持。

**`--bam` 只应提供 RNA-seq 比对。** Hi-C、ATAC、基因组重测序等数据的覆盖度不是表达信号，却会被当作 CDS 深度参与裁决，使单外显子模型被误判为有 RNA 支持。RefAnno 在扫描时统计每个 BAM 的比对数和含剪接 reads 数并写入 `run_manifest.json`；若一个 BAM 的剪接 reads 比例低于 1%，会在运行时给出警告。

## 无转录组模式：minisplice 剪接位点序列校验

有些类群没有转录组数据，只能提供 GFF + 近缘蛋白。此时不带 `--bam` 运行即可（所有模型记为 `untested`），并建议接入 [minisplice](https://github.com/lh3/minisplice) 做剪接位点序列校验，它有两方面作用：

1. **校验模型结构**：minisplice 用 1D-CNN 估计基因组中每个 GT/AG 位点是真实剪接位点的几率比分。RefAnno 检查每个 CDS 内含子的供体和受体分数，全部达到 `--min-splice-score`（默认 2）的多外显子模型记为 `sequence_supported`；未全部通过且缺乏直接证据的模型加 `weak_splice_sites` 复查标记。
2. **改善远缘蛋白比对**：miniprot ≥0.16 时自动把分数文件通过 `--spsc` 传入，提高远缘同源比对的剪接位点准确性，从而改善无转录组模式下最关键的蛋白证据质量。

minisplice 不在 conda 仓库，需要从源码编译后用 `--minisplice` 指定路径：

```bash
git clone https://github.com/lh3/minisplice && cd minisplice && make
```

```bash
refanno --set A=a.gff3 --set B=b.gff3 --genome genome.fa --proteins related.faa \
  --minisplice /path/to/minisplice --threads 16 --outdir RefAnno_out
```

模型来源有三种方式：

- 不提供 `--splice-model` 时，RefAnno 用**多套输入共同支持的转录本**（同一座位内 ≥2 个集合 CDS 剪接链完全一致的链，每条链取 CDS 最长者）自动训练和校准；一致链不足 500 条时退回可信集合（`egapx`/`eviann`/`e` 加 `--trusted-set`）的每基因最长转录本，再不足则用全部输入。训练数据、模型和分数都写入输出目录并按指纹缓存。
- 脊椎动物/昆虫基因组可以下载预训练模型（如 Zenodo 上的 `vi2-7k.kan` 及 `.cali`）后用 `--splice-model vi2-7k.kan` 指定；旁边没有 `.cali` 时会用输入注释自动校准。
- 已有全基因组分数文件（如 Zenodo 预计算的 human/fly 文件）时用 `--splice-scores score.tsv.gz` 直接提供，跳过训练与预测。

**序列校验是注释自洽性检查，不是独立证据**：ab initio 预测器本就会产生规范剪接位点，因此 `sequence_supported` 只提升代表排序、参与结构评价，不能单独保留基因；它真正的价值是反向使用——多外显子模型的位点未通过校验（`weak_splice_sites`）时值得人工查看。

## 可选 Iso-seq 增强
提供 `--isoseq` 后，完整长读链独立参与结构和代表裁决。

`--utr-consensus` 会对 CDS 兼容的长读进行局部链聚类。主导链达到分子数和主导度阈值后，可调整 UTR 外显子和末端；只有共识外显子完整覆盖 CDS 时才应用，CDS 不变。

有 Iso-seq 时，每基因默认最多输出 3 条转录本：1 条主代表和最多 2 条替代链。替代链须有不少于 2 条 Iso-seq 分子精确支持完整 CDS 内含子链，且 CDS 覆盖 ≥90%。蛋白或二代证据不能替代这一要求。`--max-isoforms 1` 可关闭替代链输出。

## 依赖

- Python 3.8 或更高版本；
- 读取 BAM 时需要 `pysam`；
- 蛋白 FASTA 自动比对时需要 `miniprot`；
- 可选：剪接位点序列校验需要自行编译的 `minisplice`（conda 无此包），用 `--minisplice` 指定路径；
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
- RefAnno 对每个 BAM **顺序扫描一遍**，同时提取 junction 和 CDS 深度，不按基因区间随机查询，因此不要求 BAM 已排序或建索引；排序和索引仍建议保留，供 IGV 等工具使用。
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
--minisplice PATH           minisplice 可执行文件（源码编译，conda 无此包）
--splice-model FILE.kan     预训练 minisplice 模型；缺省自动训练
--splice-scores FILE.tsv    已有基因组剪接位点分数，跳过训练与预测
--min-splice-score 2        供体/受体位点通过校验的最低分数
--utr-consensus             可选 Iso-seq UTR 共识修正
--threads 16                蛋白比对线程数
--set-weight NAME=FLOAT     覆盖集合默认权重，仅用于后续代表排序
--trusted-set NAME          允许可信集合的单套模型直接保留
--max-isoforms 3            每基因总转录本上限，范围 1–3
--id-prefix REFA            输出基因 ID 前缀
--outdir DIR                输出目录，必需
```

蛋白比对结果会在输出目录缓存。缓存按基因组、蛋白、miniprot 版本、线程数和选项的 SHA256 指纹校验，任一项改变即自动失效并重新比对，不需要手动换输出目录。

## 工程结构与测试

核心代码按职责拆分：

```text
refanno/
├── cli.py             命令行与主流程
├── models.py          转录本、长读链和蛋白命中数据模型
├── evidence.py        GFF/BAM/BED12、蛋白与深度证据
├── splice.py          minisplice 训练、预测、缓存与位点索引
├── adjudication.py    结构评估、保留和代表裁决
├── consensus.py       UTR 共识与可变剪切
├── output.py          GFF3、蛋白、报告和摘要输出
└── audit.py           运行清单与输入输出校验和
```

`RefAnno.py` 是源码兼容入口；安装后的 `refanno` 和 `python -m refanno` 使用同一个 `refanno.cli` 主流程。输出摘要有独立的 `output_schema_version`，报告列由端到端测试固定。

开发模式安装和测试：

```bash
python -m pip install --no-deps -e .
python -m unittest discover -v
```

测试共 7 个端到端用例，覆盖真实合成 BAM 的深度与 junction（含深度累加器按唯一/多重分开、多 CDS 块求和的直接校验）、蛋白保留、重复家族救援及其 donor 记录、非同源相邻模型不得互相救援的反例、同座位冲突关系、代表优先级、Iso-seq 可变剪切、真实 UTR 修订的 `changes.tsv` 记录，以及 miniprot 缓存复用与输入改变后的缓存失效。环境中有 miniprot 时还会运行蛋白 FASTA 自动比对集成测试；也可用 `REFANNO_MINIPROT=/path/to/miniprot` 指定可执行文件。GitHub Actions 对推送和 Pull Request 运行同一套测试。

## 输出

| 文件 | 内容 |
|---|---|
| `final.gff3` | 高置信最终基因：有直接证据（Iso-seq/二代/同源蛋白/序列校验）或多集合一致共识；含 mRNA、exon、CDS 和 UTR，可载入 IGV |
| `final.proteins.faa` | 高置信层全部转录本按最终 CDS 翻译的蛋白 |
| `weak.gff3` | 弱支持基因：`family_indirect`、`annotation_only`、`untested` 且非共识；属性与 final 相同但明确隔离，不建议直接进入下游分析 |
| `weak.proteins.faa` | 弱支持层转录本的蛋白 |
| `decision_report.tsv` | 结构状态、证据等级、保留原因、代表依据及各类证据指标 |
| `decision_events.tsv` | 每个输入模型的去向：主代表、替代转录本或淘汰，以及原因和相关模型 |
| `junction_evidence.tsv` | 逐 junction 的 unique/multi reads、Iso-seq 分子、motif 和支持状态 |
| `supporting_junctions.bed` | 最终主代表和替代转录本的支持 junction，可直接载入 IGV |
| `family_rescues.tsv` | 被家族证据救援的基因、donor 模型、重复簇、证据类型和数值 |
| `changes.tsv` | UTR 共识修订前后 exon、CDS 是否改变、支持分子和主导度 |
| `review.gff3` | 需要人工复查的模型：部分结构支持、间接救援或无结构支持 |
| `gene_confidence.bed` | 全部最终基因的 BED9 置信度轨道：itemRgb 按证据等级着色，score 为 0–1000 置信度，名称带证据等级 |
| `curation_sites.bed` | 需要人工矫正的位点轨道，按风险分从高到低排序，名称标注主要问题（internal_stop、UTR 歧义、证据等级等） |
| `conflicting_models.gff3` | 同座位冲突中被淘汰的模型及其获胜基因 |
| `dropped.gff3` | 其余未保留模型及原因 |
| `input_fates.tsv` | 每套输入的模型数、完整 CDS 数、结构支持数和最终去向 |
| `qc_report.html` | 离线 QC 报告：运行概览、证据等级与保留原因分布、输入集合命运、最终基因集结构、逐 BAM 与 junction 统计、Iso-seq/蛋白证据分布、重复家族与救援、复查优先级 Top 50 |
| `summary.json` | 软件与输出格式版本，以及基因、转录本、结构状态、修订和证据等级统计 |
| `run_manifest.json` | 命令、参数、运行环境版本、起止时间、输入输出文件大小与 SHA256；每个 BAM 另记录比对数、MAPQ 通过数和含剪接 reads 数 |
| `miniprot.cache.json` | 蛋白比对缓存指纹，输入或参数改变时自动失效 |

`decision_report.tsv` 不再提供混合 `score`。主要字段包括：

- `structure_status`、`structure_coverage`、`junctions_supported`、`junctions_total`：结构可信程度；
- `evidence_class`：模型最终依据的最高证据等级，含义见下节；
- `reason`：基因保留原因；
- `selection`：主代表选择依据；
- `cds_depth`、`cds_coverage`、`isoform_molecules`：转录证据；
- `splice_sites_supported`、`splice_min_score`：minisplice 序列校验的内含子通过数与最弱位点分数；
- `protein_target`、`protein_identity`、`protein_coverage`、`protein_query_coverage`：同源蛋白证据；
- `family_support`、`family_donors`、`duplicate_type/cluster/size`：重复家族证据、donor 与标记；
- `review`：内部终止子、家族间接保留或 UTR 共识不确定等复核标签；
- `output_tier`：`final`（高置信）或 `weak`（弱支持）。

**final/weak 分层**：`final.gff3` 只放高置信基因——证据等级为 `isoseq_exact`、`short_read_complete`、`short_read_partial`、`homology_only`、`sequence_only`，或保留原因是多集合一致链共识（`2set/3set_chain_consensus`）。其余基因（`family_indirect`、`annotation_only`、`untested` 且非共识）写入 `weak.gff3`。这样在没有转录组数据的物种上，下游不会把"凭可信来源或家族间接证据留下来的基因"和"有独立证据筛选过的基因"混成一个合集使用。两层的基因 ID 连续编号，`decision_report.tsv` 的 `output_tier` 列记录每个基因的去向。

`summary.json` 中的 `review_genes` 和 `structural_changes` 分别给出需要复查的基因数和实际发生 UTR 修订的基因数，`evidence_classes` 按等级统计最终基因。

## 证据链与复查

每个最终基因都带有 `evidence_class`，表示它依据的最高证据等级，从强到弱为：

| 等级 | 含义 |
|---|---|
| `isoseq_exact` | 完整 CDS 剪接链有 Iso-seq 分子精确支持 |
| `short_read_complete` | 全部 CDS junction 有二代支持，或单外显子达到覆盖和深度要求 |
| `short_read_partial` | 有部分二代 RNA 证据，但结构未全部支持 |
| `homology_only` | 没有可达阈值的 RNA 结构支持，凭合格同源蛋白保留 |
| `family_indirect` | 凭重复家族的间接证据保留，自身无直接证据 |
| `sequence_only` | 无 RNA 与蛋白证据，但全部 CDS 内含子通过 minisplice 序列校验 |
| `untested` | 没有 RNA 数据，无法评价转录结构 |
| `annotation_only` | 有 RNA 数据但未获结构支持，也未进入其他通道 |

推荐的 IGV 复查方式是同时载入 `final.gff3`、`review.gff3` 和 `supporting_junctions.bed`：`final.gff3` 查看保留结果，`review.gff3` 高亮需要复核的模型，`supporting_junctions.bed` 直接对照 BAM 中的剪接位点。发现可疑基因后，再到 `decision_report.tsv` 查看该基因的证据指标和 `review` 标签，在 `decision_events.tsv` 追溯同座位被淘汰的模型，在 `junction_evidence.tsv` 查看逐 junction 的 reads 支持。`family_indirect` 的基因可沿 `family_rescues.tsv` 找到提供证据的 donor 模型。

`qc_report.html` 不依赖外部资源，可直接用浏览器打开。报告分九节：运行概览；证据链（证据等级、保留原因、主代表依据、结构状态分布，证据等级 × 单外显子交叉表，final/weak 分层 × 证据等级，可选的 minisplice 序列校验小节）；输入集合命运（保留率、final/weak 贡献、淘汰原因分布、集合间一致链矩阵）；最终基因集的外显子数/CDS 长度/转录本数分布；RNA 证据（逐 BAM 统计并对剪接率异常的 BAM 给出警告、主代表 junction 支持与弱支持 junction 列表、CDS 深度与覆盖度分箱及二维交叉）；Iso-seq 证据（含 UTR 修订的长度变化分布）；蛋白同源证据（一致性分布，并在命中高度一致时提示自我确认风险）；重复与基因家族（含家族救援明细）；复查优先级（review 标记计数与按风险分排序的 Top 50 基因）。基因数量和 BUSCO 不能单独证明结构正确，QC 报告因此以证据等级分布、复查清单和结构变更为主线，而不是只给一个通过比例。

## 对接 GSAman 人工矫正

RefAnno 的输出可以直接进入 [GSAman](https://github.com/CJ-Chen/GSAman) 的"人工介入"矫正流程。推荐按以下方式组织轨道：

1. `cp final.gff3 final.gsaman.gff3` 后导入，即进入可编辑模式（GSAman 通过 `.gsaman` 后缀开启编辑），矫正完成后从 GSAman 导出 GFF3；
2. 载入基因组 FASTA 和二代 BAM（同一次运行使用的比对文件）；
3. 载入 `supporting_junctions.bed`，直接对照每个保留 junction 的 reads 支持；
4. 载入 `curation_sites.bed` 作为矫正任务清单：文件按风险分从高到低排序，从上到下逐个跳转处理即可；名称中的后缀说明该位点的主要问题（`internal_stop`、`utr_consensus_ambiguous`、`homology_only` 等）；
5. 可选载入 `gene_confidence.bed` 作为全局底色轨道，一眼区分可信区与风险区。

颜色与置信度约定：

| 颜色 | 含义 |
|---|---|
| 红 | 含内部终止子，优先级最高 |
| 绿 | Iso-seq 精确支持或二代完整支持，且无复查标记 |
| 黄 | 二代部分支持 |
| 橙 | 仅凭同源蛋白或重复家族间接证据保留 |
| 灰 | 无有效证据 |

BED 第 5 列 score 是 0–1000 的置信度：以证据等级为基准分（isoseq_exact 900、short_read_complete 800、short_read_partial 550、homology_only 400、family_indirect 300、annotation_only 200），每有一个不支持的 junction 扣 30、每有一个复查标记扣 50，内部终止子直接压到 50。`curation_sites.bed` 的排序规则与 `qc_report.html` 第 9 节的风险分一致。

## 许可证

RefAnno 使用 MIT License，详见 `LICENSE`。
