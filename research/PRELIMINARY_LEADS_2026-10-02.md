# 初步线索(非系统性,未按 PRISMA 检索)

**性质**:2026-10-02 由 AI 用通用网页搜索做的 5 次查询。数据库 API 被网络策略拒绝(403),因此**没有可复现检索式、没有命中数、没有双库覆盖**。下表条目来自**搜索结果的标题和摘要片段**,**我没有阅读任何一篇全文**,也没有按协议纳入标准逐篇判断。它只能用来决定下一步怎么查,不能当作"别人已做到哪一步"的结论。

## 5 次查询
1. 横断面肿瘤蛋白/磷酸化数据推断带反馈的信号网络
2. CARNIVAL / 因果推理 / OmniPath
3. 横断面稳态数据中有向边与反馈环的可识别性
4. 生物学约束神经网络,突变到表达,患者肿瘤
5. CPTAC 蛋白基因组网络推断与独立队列验证

## 线索(按与 Q1/Q2/Q3 的相关性分组;"未读"指未读全文)

### A. 先验网络加多组学的因果上下文化(Q1)
- CARNIVAL:从表达足迹推断因果通路,用 OmniPath 先验,整数线性规划。[Nature npj Syst Biol Appl](https://www.nature.com/articles/s41540-019-0118-z)(未读)
- Causal integration of multi-omics data with prior knowledge to generate mechanistic hypotheses。[Mol Syst Biol](https://link.springer.com/article/10.15252/msb.20209730)(未读)
- CORNETO:多样本网络推断,用 CPTAC 肺腺癌示例。[Nature Mach Intell](https://www.nature.com/articles/s42256-025-01069-9)(未读)
- Phoslink:CPTAC 上磷酸化与表达的因果调控链接。[ScienceDirect](https://www.sciencedirect.com/science/article/pii/S1535947625000039)(未读)
- DMPA:不依赖先验的多组学通路分析。[PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11259815/)(未读)
- Causal network analysis of omics data using prior knowledge databases。[Briefings in Bioinformatics](https://academic.oup.com/bib/article/26/6/bbaf654/8371896)(未读,可能是综述)

### B. 生物学约束神经网络(Q1)
- P-NET,前列腺癌。[Nature](https://www.nature.com/articles/s41586-021-03922-4)(未读)
- MPVNN。[Bioinformatics](https://academic.oup.com/bioinformatics/article/38/22/5026/6705227)(未读)
- PGLCN,胃癌 TMB。[PubMed](https://pubmed.ncbi.nlm.nih.gov/37810279/)(未读)
- **已有综述**:A systematic review of biologically-informed deep learning models for cancer。[BMC Bioinformatics](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-023-05262-8)(未读)。从标题看,这些模型多用于预测结局,而非推断反馈调控,但这一点**未核实**。

### C. 验证与基准(Q2)
- **Phosphoproteomics data-driven signalling network inference: Does it work?**[PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9798138/)(未读)。标题直接对应"这类推断能否可靠",优先阅读。
- Comprehensive evaluation of phosphoproteomic-based kinase activity inference。[PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12098709/)(未读)
- Reconstructing kinase network topologies from phosphoproteomics data(癌症相关重连)。[Nature Biotechnol](https://www.nature.com/articles/s41587-019-0391-9)(未读)

### D. 可识别性(Q3)
- Bicycle(含环的干预因果发现)。[PMLR](https://proceedings.mlr.press/v236/rohbeck24a/rohbeck24a.pdf)(未读)
- Directed Cyclic Graphs ... Instrumental Variables。[JMLR](https://www.jmlr.org/papers/volume26/23-0272/23-0272.pdf)(未读)
- Causal Discovery with Heterogeneous Observational Data。[PMLR](https://proceedings.mlr.press/v180/zhou22a/zhou22a.pdf)(未读)

搜索摘要中的说法:"仅凭观察数据,含环有向图一般只能识别到马尔可夫等价类"。这是搜索工具转述,**未在原文核实**。

## 这些线索暗示了什么(推断,不是结论)
1. "用先验通路加多组学做因果上下文化"这个方向**已有成熟方法**(A 组),所以"新"很可能不在这里。
2. 从标题看,A、B 组多数不明确处理反馈环,也未必在独立队列上验证边;这是 C、D 组要核实的问题。
3. 项目最可能的差异点,或最大的风险,都集中在"横断面数据下的反馈可识别性与验证",对应 Q2、Q3。

## 局限与下一步
- 非系统、不可复现、只读了标题与摘要片段;可能漏掉大量文献,也可能被搜索排序偏置。
- 下一步:取得学术数据库访问(放开网络,或由作者导出),按 `SEARCH_STRINGS_v1.0.md` 做正式检索。在此之前,不得把本文件当作 PRISMA 结果引用。

---

# 增补 2(同日):扩大搜索,并尝试逐篇读全文

## 全文读取结果:**0 篇成功**
对第一批 10 篇(PMC9798138、CARNIVAL PMC6848167、MSB PMC7838823、Phoslink PMC11889353、CausalPath PMC8633371、CORNETO Nature、COSMOS bioRxiv、k-hop GAT bioRxiv、综述 PMC12703490、BMC 系统综述)用 WebFetch 读取,全部返回 `EGRESS_BLOCKED`。随后用 curl 探测了 14 个其他主机(pmc.ncbi、europepmc、arxiv、PMLR、Springer、OUP、ScienceDirect、GitHub、Cell、Bioconductor、saezlab.github.io、Frontiers、MDPI、JMLR),均 403。被拒主机清单:`www.ncbi.nlm.nih.gov`、`pmc.ncbi.nlm.nih.gov`、`eutils.ncbi.nlm.nih.gov`、`europepmc.org`、`www.ebi.ac.uk`、`www.nature.com`、`www.biorxiv.org`、`link.springer.com`、`academic.oup.com`、`www.sciencedirect.com`、`www.cell.com`、`arxiv.org`、`export.arxiv.org`、`proceedings.mlr.press`、`www.jmlr.org`、`bmcbioinformatics.biomedcentral.com`、`www.frontiersin.org`、`www.mdpi.com`、`bioconductor.org`、`saezlab.github.io`、`github.com`、`api.openalex.org`、`api.crossref.org`、`api.semanticscholar.org`。
因此以下所有条目仍是**搜索片段层面的记录,未读全文**,不能填写数据提取表。

## 新增候选(第二轮 8 个查询;仍为"未读")
| 编号 | 条目 | 链接 | 与问题的相关性(仅据标题/片段) |
|---|---|---|---|
| N1 | CausalPath(Patterns 2021;含 CPTAC 11 癌种 1,110 患者的应用) | [Patterns](https://www.cell.com/patterns/fulltext/S2666-3899(21)00083-0) | Q1/Q2:先验加蛋白/磷酸化因果优先;片段称用于 CPTAC |
| N2 | COSMOS(含 11 例 ccRCC 患者的应用) | [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.07.15.603538.full.pdf) | Q1:跨层先验网络,样本量极小 |
| N3 | k-hop graph attention 的信号网络推断 | [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.09.16.508281.full.pdf) | Q1:GNN 加多组学加先验 |
| N4 | NEM-Tar:癌症调控网络概率图模型 | [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8100334/) | Q1/Q3:概率图模型 |
| N5 | pyPARAGON:疾病网络构建,整合磷酸化数据与互作组 | [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11334722/) | Q1 |
| N6 | 多层网络与有向随机游走做通路活性 | [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.07.22.163949.full.pdf) | Q1(结局预测取向) |
| N7 | Inferencing bulk tumor and single-cell multi-omics regulatory networks(综述性) | [MDPI](https://www.mdpi.com/2073-4409/12/1/101) | 片段称:缺少时间序列是核心障碍;贝叶斯网络多为无环 → 与 Q3 相关(**未核实**) |
| N8 | Patient-specific Boolean models(TCGA 前列腺,432-488 例) | [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9018074/) | Q1:布尔动力学,可含环 |
| N9 | Patient-specific logic models(活检筛选) | [MSB](https://link.springer.com/article/10.15252/msb.20188664) | Q1/Q2:需要药物筛选数据,可能触及排除标准 |
| N10 | Systematic analysis of somatic mutations impacting gene expression in 12 tumour types(trans 效应) | [Nat Commun](https://www.nature.com/articles/ncomms9554) | Q1:突变到表达,但可能无网络 |
| N11 | DriverNet:突变对转录网络的影响 | [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4056374/) | Q1 |
| N12 | Proteogenomic insights suggest druggable pathways in endometrial carcinoma(CPTAC-UCEC) | [Cancer Cell](https://www.sciencedirect.com/science/article/pii/S1535610823002477) | 与本项目数据直接相关;描述性,不一定为网络推断 |
| N13 | Pan-cancer proteogenomic investigations identify post-transcriptional kinase targets | [Commun Biol](https://www.nature.com/articles/s42003-021-02636-7) | Q1 |
| N14 | Deciphering the dark cancer phosphoproteome using machine-learned co-regulation | [Nat Commun](https://www.nature.com/articles/s41467-025-57993-2) | Q1/Q2 |
| N15 | Comprehensive evaluation of phosphoproteomic-based kinase activity inference | [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12098709/) | Q2:基准 |
| N16 | Causal Inference Methods to Integrate Omics and Complex Traits | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8091955/) | Q3:综述 |
| N17 | When Does GRN Inference Break?(单细胞诊断研究) | [arXiv](https://arxiv.org/html/2605.04930v1) | Q3:单细胞,可能超出纳入范围 |
| N18 | 干预型因果结构学习比较;CausalBench | [bioRxiv](https://www.biorxiv.org/content/10.64898/2025.12.05.692565.full.pdf) | Q3:依赖扰动,按协议属"附表登记",非主分析 |

按协议第 4 节,依赖扰动或时间序列的方法(Bicycle、CausalBench、N9、N18)不入主分析,但保留在附表。

## 目前能说和不能说的
- **能说**:搜索排名靠前的结果显示,"先验网络+多组学+因果推理"有许多成熟工具,至少包括 CARNIVAL、CausalPath、COSMOS、CORNETO、Phoslink;存在至少一篇关于生物学约束深度学习的系统综述和至少一篇关于该类因果分析的综述。
- **不能说**:它们是否允许反馈环、是否在独立队列验证、是否陈述可识别性假设。这三点正是 Q1–Q3 的核心,**全部需要全文,当前读不到**。
- 搜索片段里有两句相关的话,但都是搜索工具转述,未在原文核实:(a) 静态贝叶斯网络无法表示反馈环;(b) 含环有向图仅凭观察数据一般只能识别到马尔可夫等价类。

## 解除阻塞所需
任一:(1) 在环境设置里放开上面的被拒域名(至少 `pmc.ncbi.nlm.nih.gov`、`europepmc.org`、`www.biorxiv.org`、`arxiv.org`、`www.nature.com`);(2) 作者把关键论文的 PDF/全文文本放进 `research/papers/`;(3) 在有网络的新会话里继续,本文件和协议已在分支上。
