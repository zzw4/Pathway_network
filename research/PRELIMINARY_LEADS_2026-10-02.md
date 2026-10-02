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
