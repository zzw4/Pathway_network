# 检索式 v1.0(依据 PROTOCOL_v1.0.md 第 5 节)

由 AI 起草,未经 PRESS 同行核查;检索后在报告中披露。每个库检索后请记录:**检索日期、命中数、使用的过滤器**。

## 1. PubMed(https://pubmed.ncbi.nlm.nih.gov)

在首页检索框粘贴(整条一行):

```
(cancer[tiab] OR tumor[tiab] OR tumour[tiab] OR neoplasm*[tiab]) AND (multiomic*[tiab] OR multi-omic*[tiab] OR proteogenomic*[tiab] OR phosphoproteom*[tiab]) AND ("signaling network*"[tiab] OR "signalling network*"[tiab] OR "regulatory network*"[tiab] OR "causal network*"[tiab] OR pathway*[tiab]) AND (inference[tiab] OR infer*[tiab] OR reconstruct*[tiab] OR "network learning"[tiab]) AND ("prior knowledge"[tiab] OR "knowledge-guided"[tiab] OR "pathway-informed"[tiab] OR feedback[tiab])
```

过滤器:年份 2010-2026。导出:Save → Selection: All results → Format: CSV → Create file。

## 2. Web of Science(Core Collection,高级检索)

```
TS=((cancer OR tumor OR tumour OR neoplasm*) AND (multiomic* OR "multi-omic*" OR proteogenomic* OR phosphoproteom*) AND ("signaling network*" OR "signalling network*" OR "regulatory network*" OR "causal network*" OR pathway*) AND (inference OR infer* OR reconstruct* OR "network learning") AND ("prior knowledge" OR "knowledge-guided" OR "pathway-informed" OR feedback))
```

过滤器:2010-2026。导出:Export → Tab delimited / RIS,记录 Full Record。

## 3. Scopus

```
TITLE-ABS-KEY((cancer OR tumor OR tumour OR neoplasm*) AND (multiomic* OR "multi-omic*" OR proteogenomic* OR phosphoproteom*) AND ("signaling network*" OR "signalling network*" OR "regulatory network*" OR "causal network*" OR pathway*) AND (inference OR infer* OR reconstruct* OR "network learning") AND ("prior knowledge" OR "knowledge-guided" OR "pathway-informed" OR feedback))
```

过滤器:2010-2026。导出:Export → CSV(勾选 Citation + Abstract)。

## 4. 预印本(arXiv / bioRxiv)

在 arXiv 高级检索与 bioRxiv 高级搜索里,用关键词:`signaling network inference cancer multi-omics prior knowledge`;逐条人工记录候选(预印本无统一导出)。

## 记录表(检索后填)

| 数据库 | 检索日期 | 命中数 | 过滤器 | 导出文件名 |
|---|---|---|---|---|
| PubMed | | | 2010-2026 | |
| Web of Science | | | 2010-2026 | |
| Scopus | | | 2010-2026 | |
| arXiv/bioRxiv | | | | |

## 说明
- 如果某个库命中数过大(>3000)或过小(<50),把数字告诉我,检索式要调整,调整记入协议修订表。
- 没有机构访问权限的库(WoS、Scopus)可以跳过,只做 PubMed;这会在报告局限里写明。
