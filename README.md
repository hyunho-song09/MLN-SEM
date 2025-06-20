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
