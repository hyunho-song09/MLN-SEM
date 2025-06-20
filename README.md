# MLN-SEM
Multi-Layer Network-SEM


# MLN-SEM for Microglia Metabolic Pathway Analysis

This repository contains a fully automated pipeline for identifying plausible causal pathways from **Aβ treatment (Ab)** to **EV or Media metabolites** via intracellular **microglial metabolic networks**. It leverages:

- **WGCNA** for module detection
- **Spearman correlation** for module–target association
- **igraph** for network path inference
- **lavaan** for Structural Equation Modeling (SEM)
- **Bootstrap-based mediation analysis** for causal inference

## 📂 Folder Structure



## ⚙️ Dependencies

Install the following R packages:

```r
install.packages(c(
  "lavaan", "igraph", "tidyverse", "progress", "KEGGREST", "xlsx",
  "WGCNA", "flashClust", "doSNOW", "foreach", "purrr", "parallel"
))
devtools::install_github("jeonghwanhyeon/LMSstat")  # for LMSstat

source("main.R")
main()


---

필요하다면 `functions.R`, `config.R`, `results_summary.csv` 등으로 스크립트를 더 나누는 것도 가능합니다. 원하시면 디렉토리화 및 분할 작업도 도와드릴게요.
