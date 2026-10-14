# E1 · 外部数据准入核查报告（修订版 v2）

**核查日期**：2026-09-18
**修订日期**：2026-09-18（作者核查后更正 SRA 一项，并补充 GSE237846）
**核查对象**：GSE217845（scRNA-seq）、GSE310352（CosMx）、GSE277782（CosMx）、GSE237846（Molecular Cartography）
**核查方式**：GEO E-utilities（esearch / esummary / elink）+ FTP 目录 + 实样下载（基因面板、barcode、元数据）
**核查性质**：只读元数据与文件内容检查，未运行任何表达分析

---

## 〇、对 v1 的更正

### 更正一：GSE217845 **存在 SRA 原始测序数据**

**v1 的结论"三个系列在 GEO 均未关联 SRA/BioProject"，对 GSE217845 是错的。**

我的判断来自 ESummary 的 `relation` 字段为空，**没有使用 `elink`**。用 `elink` 从 GDS 连到 SRA 后：

```
dbfrom=gds  db=sra  id=200217845  ->  16 个链接
```

逐样本核对也确认，例如 `GSM6727546`：

```
!Sample_type      = SRA
!Sample_relation  = BioSample: https://www.ncbi.nlm.nih.gov/biobsample/SAMN31700947
!Sample_relation  = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX18243182
```

**因此撤回以下判断**：

> ~~GSE217845 无 SRA，因此无法从原始测序数据重建矩阵~~

**保留的判断**：GSE217845 的公开表达矩阵是 **filtered**（`LPDAC_30_Tumor` 的 barcode 数为 9,246）。

**范围限定**：本更正**只针对 GSE217845**。G1 报告对原稿另外三个系列（GSE155698、GSE156405、GSE212966）的核查用的是同样的 `esummary` 路径，因此**同样存在漏检 SRA 的可能**，需用 `elink` 重新核验（见第五节）。

### 更正二：补充 GSE237846（同一研究的空间数据）

作者提示该研究另公开了 Molecular Cartography 数据。已核查，结论见第四节。

---

## 一、准入裁决（v2）

| 候选 | 类型 | 三基因覆盖 | 患者映射 | **SRA** | 裁决 |
|---|---|---|---|---|---|
| **GSE217845** | scRNA-seq | **3/3** | 样本名内含 | **有（16 个 experiment）** | **通过** |
| GSE310352 | CosMx | **1/3** | 需核 | — | **不通过** |
| GSE277782 | CosMx | **2/3** | 完整 | — | **不通过** |
| **GSE237846** | Molecular Cartography | **0/3** | 完整 | — | **不通过** |

**三个空间数据集全部落选，原因都是面板不含 SIGLEC15。**

---

## 二、GSE217845 —— 通过

**标题**：IL-1b+ tumor-associated macrophages fuel pathogenic inflammation in pancreatic cancer - scRNAseq
**公开日期**：2023-09-11　**样本数**：16（10 个肿瘤 + 6 个外周血）
**BioProject**：`PRJNA901009`　**SRA**：16 个 experiment

### 三基因覆盖

从 `GSM6727545_LPDAC_30_Tumor_features.tsv.gz` 读取：**33,514 个基因**。

```
[OK] SIGLEC15    [OK] ACP5    [OK] MMP9
```

### 患者映射（样本标题）

```
肿瘤样本（10）
  GSM6727542  Patient LPDAC_15  Tumor CD45+CD3-CD19- sorted scRNAseq   ← 分选
  GSM6727543  Patient LPDAC_25  Tumor scRNAseq
  GSM6727544  Patient LPDAC_26  Tumor scRNAseq
  GSM6727545  Patient LPDAC_30  Tumor scRNAseq
  GSM6727546  Patient PDAC_47   Tumor scRNAseq
  GSM6727547  Patient PDAC_48   Tumor scRNAseq
  GSM6727549  Patient PDAC_51   Tumor scRNAseq
  GSM6727550  Patient PDAC_55   Tumor scRNAseq
  （另有 2 个肿瘤样本，见 E1_ELIGIBILITY.json）

外周血（6）
  GSM6727552  Patient PDAC_47   Peripheral blood scRNAseq
  GSM6727553  Patient PDAC_48   Peripheral blood scRNAseq
  GSM6727555  Patient PDAC_51   Peripheral blood scRNAseq
  GSM6727556  Patient PDAC_55   Peripheral blood scRNAseq
  （另有 2 个外周血样本）
```

**样本标题显式写 `Patient <ID>`**，肿瘤与配对血样共用同一 ID。

### 关于 `PDAC_*` 与 `LPDAC_*`

**已核实的线索**：同一研究的空间数据 GSE237846 中，样本名为 `Patient LPDAC30`（原发灶）与 `Patient LIMET*`（肝转移）。**`LPDAC` 与肝转移并列出现，说明它指原发灶，而非转移灶。**

**仍未核实的**：`LPDAC` 是否意味着与 `PDAC_*` 不同的分期或治疗背景。**GEO 记录本身没有解释该前缀**，必须查原始论文及其补充材料。

**因此 E2 中的处理（分析规划，不是疾病认定）**：

| 类别 | E2 处理 |
|---|---|
| `PDAC_*` 肿瘤样本 | 进入主要候选池，先确认髓系细胞可估计性 |
| 未分选的 `LPDAC_*` 肿瘤样本 | 保留，单独标记，分层敏感性分析 |
| `LPDAC_15`（分选） | 单独标记；**不用于检出率与未分选样本的比较** |
| 外周血样本 | 不纳入肿瘤髓系对比量 |

---

## 三、两个 CosMx 数据集 —— 不通过

### GSE310352（990 基因面板，124,165 细胞）

```
[--] SIGLEC15      [--] ACP5      [OK] MMP9
```

SIGLEC15 与 ACP5 均缺失，**无法计算 ACP5/MMP9 程序分数，也无法定义 SIGLEC15 可检出细胞**。

### GSE277782（1000 基因面板，82,348 细胞）

```
[--] SIGLEC15      [OK] ACP5      [OK] MMP9
```

SIGLEC15 缺失，无法定义对比分组。

**该数据集在其他维度质量很高（患者 `Pt` 列完整、原发 + 肝转移配对、空间多边形与转录本坐标齐全、Macrophage 3,477 个），若将来找到含 SIGLEC15 的同平台数据可直接复用。**

---

## 四、GSE237846（Molecular Cartography）—— 不通过

**标题**：IL-1b+ tumor-associated macrophages fuel pathogenic inflammation in pancreatic cancer - Molecular Cartography
**公开日期**：2023-09-11　**样本数**：13　**平台**：98 基因定制面板

### 按作者指示，只读取基因面板

从 `GSM7655264_LPDAC30_B2_segmentation_counts.csv.gz`（158 KB）读取面板，**未下载任何 TIFF 图像**（单张 DAPI/KRT19 图像为 180–400 MB）。

**面板全部 98 个基因**：

```
ACTA2, ANPEP, APOC1, APOE, C1QC, CCR2, CCR7, CD163, CD1C, CD27, CD274, CD36,
CD3D, CD4, CD68, CD8A, CDK1, CHIT1, CLEC10A, CLEC4E, CLEC9A, CLPS, CPA3, CSF3R,
CTLA4, CXCL1, CXCL3, EOMES, F13A1, FABP5, FCGR3A, FCGR3B, FOLR2, FOXP3, FSCN1,
GCG, GPNMB, GPR34, GPR84, GZMB, HAVCR2, HLA-DQB2, HLA-DRA, IAPP, ICOS, IFNG,
IL1A, IL1B, IL1R1, IL2RA, IL6, IL7R, INHBA, KRT19, LAG3, LIPA, LPL, LST1,
LYVE1, MAF, MARCO, MCEMP1, MEF2A, MEF2B, MEF2C, MEF2D, MERTK, MKI67, MRC1,
MS4A1, MS4A7, MZB1, NEUROD1, NKG7, NLRP3, PDCD1, PDCD1LG2, PDGFB, PDGFRA,
PECAM1, PLP1, PTGER1, PTGER2, PTGER3, PTGER4, PTGS2, PTPRC, RGS5, SCTR,
SELL, SLCO2B1, STAB1, THBS1, TIGIT, TNF, TNFRSF9, TOP2A, VCAN, XCR1
```

**三基因比对结果**：

```
[--] SIGLEC15      [--] ACP5      [--] MMP9
```

**三个目标基因全部不在面板内。** 面板中含 SIGLEC/ACP/MMP 字样的仅 `CD68` 与 `CD163`，均为髓系身份标记，不是本研究的目标或程序基因。

**这是一个免疫细胞分型面板**（T 细胞检查点、髓系亚群、肿瘤/基质标记），**不含骨吸收与基质降解程序基因**。

**裁决：不通过。** 按作者指示立即停止，未下载图像文件。

### 该数据的其他属性（记录备查）

| 项 | 内容 |
|---|---|
| 样本 | `LPDAC30`（3 张切片：B2/C2/D2）、`LIMET10/49/139/162/163/165/325` 等 |
| 每样本文件 | `DAPI.tiff.gz`、`KRT19.tiff.gz`、`results.txt.gz`、`segmentation.csv.gz`、`segmentation_cell_stats.csv.gz`、`segmentation_counts.tsv.gz` |
| 细胞统计列 | `cell`、`x`、`y`、`cluster`、`n_transcripts`、`density`、`elongation`、`area`、`avg_confidence` |
| 空间坐标 | **有**（`x`、`y`） |
| 每细胞转录本数 | 约为 **36**（首行示例），偏低 |

---

## 五、由更正一引出的连带问题

**G1 报告对原稿三个系列的 SRA 判断需要重新核验。**

G1 用的是 ESummary 的 `relation` 字段；本次更正表明该字段**可能为空而实际存在 SRA 关联**。因此对 GSE155698、GSE156405、GSE212966 需要用 `elink` 重查。

**这一点直接关系到 SI §16 中"三个系列未关联 SRA"的表述。** 该表述目前**未经 `elink` 验证**，若重查后发现任一系列存在 SRA，则 SI 的这一句需要相应修改。

**建议在 E2 开始前先完成这项重查。**

---

## 六、对 E2 的影响

| E2 子项 | 可行性 | 依据 |
|---|---|---|
| **E2-A** 从 SRA 重建 | **本机不可执行** | 见下 |
| **E2-B** 独立细胞关联复核 | **可执行** | filtered 矩阵可下载（每样本 10–50 MB） |
| **E2-C** 样本构成敏感性 | **可执行** | 依赖 E2-B |

### E2-A 不可执行的具体依据

| 项 | 实测 |
|---|---|
| SRA experiment 数 | 16 |
| 单 Run 规模 | `total_spots` ≈ 3.1×10⁸，`total_bases` ≈ 4.0×10¹⁰ |
| 三个 Run 合计 bases | ≈ 1.22×10¹¹ |
| 16 个 experiment 估算 | ≈ **6.5×10¹¹ bases**，FASTQ 未压缩约 **200 GB** |
| 本地磁盘可用 | C: 63 GB｜D: 223 GB｜E: 499 GB（合计约 785 GB） |
| 本地内存 / CPU | 31.4 GB / 24 逻辑核 |
| **SRA Toolkit**（`fasterq-dump` / `prefetch`） | **未安装** |
| **Cell Ranger** | **未安装**（需商业许可） |
| **STARsolo / salmon / alevin-fry** | **均未安装** |

**结论**：E2-A 需要先安装 SRA Toolkit 与比对工具，再下载约 200 GB reads 并对 16 个样本逐一比对，属于数天量级的独立工作。**本轮不执行，如实记录。**

**这改变了 E2 的技术路线**：E2-B 将基于**公开的 filtered 矩阵**（与原稿三系列同一数据层级），因此**同样无法做数据驱动的环境 RNA 校正**。作者指出的"可以增加从 FASTQ 重建矩阵、检查未过滤 droplets 以及实施标准环境 RNA 校正"这一路线**成立但未被本轮执行**。

---

## 七、产出文件

```
QA\gate4\E1_ELIGIBILITY.json   四候选的系列元数据、样本清单、文件清单、面板探测
QA\gate4\E1_REPORT.md          本报告（v2）
Scripts\e1_eligibility.py      准入核查脚本
QA\gate4\_e1_cache\            实样缓存（面板、barcode、元数据；无图像文件）
```

**核查未修改任何稿件文件。E2、E3 尚未执行。**
