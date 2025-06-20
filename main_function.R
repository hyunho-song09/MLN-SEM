# === 00. Load Required Packages ===
load_required_packages <- function() {
  pkgs <- c("lavaan", "igraph", "tidyverse", "progress", "KEGGREST", "xlsx", 
            "WGCNA", "flashClust", "doSNOW", "foreach", "purrr", "parallel", "LMSstat")
  lapply(pkgs, require, character.only = TRUE)
}

# === 01. Load and Filter Data ===
load_data <- function(file_path) {
  raw <- read.xlsx(file_path, sheetName = "Final")
  mapping <- read.xlsx(file_path, sheetName = "Mapping")
  list(raw = raw %>% filter(Ab %in% c(0, 1)), mapping = mapping)
}

# === 02. Identify Significant Metabolites ===
identify_targets <- function(data, sample_type, metabolite_cols, 
                             adjust_p = TRUE, method = "t", p_cutoff = 0.05) {
  subset <- data %>% filter(Sample_Type == sample_type)
  ab_0 <- subset %>% filter(Ab == 0) %>% select(all_of(metabolite_cols))
  ab_1 <- subset %>% filter(Ab == 1) %>% select(all_of(metabolite_cols))
  
  raw_p <- map_dbl(metabolite_cols, ~ {
    tryCatch({
      if (method == "t_test") t.test(ab_0[[.x]], ab_1[[.x]])$p.value
      else if (method == "u_test") wilcox.test(ab_0[[.x]], ab_1[[.x]], exact = FALSE)$p.value
      else stop("Invalid method")
    }, error = function(e) NA_real_)
  })
  
  p_val <- if (adjust_p) p.adjust(raw_p, method = "BH") else raw_p
  metabolite_cols[!is.na(p_val) & p_val <= p_cutoff]
}

# === 03. Statistical Boxplots ===
run_stats_and_plot <- function(df, target, color = NULL) {
  df$Group[df$Group == paste0(target, "_Ab")] <- 1
  df$Group[df$Group == paste0(target, "_Con")] <- 0
  stats <- All_stats(df, Adjust_p_value = TRUE, Adjust_method = "BH", parallel = FALSE)
  Boxplot(stats, asterisk = "u_test", color = color)
  return(stats)
}

# === 04. WGCNA Module Detection ===
run_wgcna <- function(dat, power = 6, network_type = "signed", cor_fun = "bicor",
                      cor_opts = "use = 'pairwise.complete.obs'", merge_thresh = 0.20,
                      min_size = 4, method = "average") {
  dat <- scale(dat, center = TRUE, scale = TRUE)
  A <- adjacency(dat, power = power, type = network_type, corFnc = cor_fun, corOptions = cor_opts)
  diss <- TOMdist(A, TOMType = network_type)
  tree <- flashClust(as.dist(diss), method = method)
  labels1 <- cutreeDynamic(tree, distM = diss, method = "hybrid", deepSplit = 4, 
                           pamRespectsDendro = TRUE, minClusterSize = min_size)
  labels1 <- labels2colors(labels1)
  merge <- mergeCloseModules(dat, labels1, corFnc = cor_fun, corOptions = list(use = 'pairwise.complete.obs'), cutHeight = merge_thresh)
  list(modules = merge$colors, eigengenes = orderMEs(merge$newMEs), dat = dat, A = A)
}

# === 05. Evaluate Module Differences by Ab ===
compare_modules_by_ab <- function(MEs, Ab_vector) {
  pvals <- map_dbl(colnames(MEs), ~ wilcox.test(MEs[, .x][Ab_vector == 0], MEs[, .x][Ab_vector == 1], exact = FALSE)$p.value)
  names(pvals) <- colnames(MEs)
  sig_modules <- names(pvals)[pvals < 0.05]
  list(significant = sig_modules, pvalues = pvals)
}

# === 06. Correlation with Media/EV ===
correlate_modules_targets <- function(MEs, data_list, targets, method = "spearman") {
  result <- data.frame()
  for (stype in names(data_list)) {
    for (mod in colnames(MEs)) {
      eig <- MEs[, mod]
      for (met in intersect(colnames(data_list[[stype]]), targets)) {
        valid <- complete.cases(eig, data_list[[stype]][[met]])
        if (sum(valid) >= 3) {
          ct <- cor.test(eig[valid], data_list[[stype]][[met]][valid], method = method)
          result <- rbind(result, data.frame(Sample_Type = stype, Module = mod, Target_Metabolite = met, Correlation = ct$estimate, P_value = ct$p.value))
        }
      }
    }
  }
  result %>% filter(abs(Correlation) > 0.5 & P_value < 0.05)
}

# === 07. SEM Modeling ===
generate_sem_model <- function(path) {
  if (length(path) < 2) return(NULL)
  eqs <- paste0(path[1], " ~ Ab")
  eqs <- c(eqs, paste0(path[-1], " ~ ", path[-length(path)]))
  paste(eqs, collapse = "\n")
}

fit_sem_models <- function(paths, data) {
  sem_models <- list()
  for (p in paths) {
    model <- generate_sem_model(p)
    if (!is.null(model)) {
      fit <- tryCatch(sem(model, data = data, estimator = "ML", std.ov = TRUE), error = function(e) NULL)
      if (!is.null(fit)) sem_models[[paste("Ab ->", paste(p, collapse = " -> "))]] <- fit
    }
  }
  sem_models
}

# === 08. Main Execution ===
main <- function() {
  load_required_packages()

  file_path <- "input/250503_final.xlsx"
  input <- load_data(file_path)
  data <- input$raw; mapping <- input$mapping
  
  media_tgts <- identify_targets(data, "Media", mapping$Mapping, FALSE, "u_test", 0.05)
  ev_tgts <- identify_targets(data, "EV", mapping$Mapping, FALSE, "u_test", 0.05)
  
  media_data <- data %>% filter(Sample_Type == "Media")
  ev_data <- data %>% filter(Sample_Type == "EV")
  micro_data <- data %>% filter(Sample_Type == "Microglia")
  
  media_df <- cbind(media_data[, c(1,4)], media_data[, media_tgts])
  ev_df <- cbind(ev_data[, c(1,4)], ev_data[, ev_tgts])
  micro_df <- cbind(micro_data[, c(1,4)], micro_data[, 5:ncol(micro_data)])
  
  run_stats_and_plot(micro_df, "Microglia", c("darkgreen", "darkred"))
  run_stats_and_plot(media_df, "Media")
  run_stats_and_plot(ev_df, "EV", c("darkgreen", "darkred"))
  
  net <- run_wgcna(micro_data[, 5:ncol(micro_data)])
  MEs <- net$eigengenes
  sig_modules <- compare_modules_by_ab(MEs, micro_data$Ab)
  targets <- unique(c(media_tgts, ev_tgts))
  
  cor_result <- correlate_modules_targets(MEs[, sig_modules$significant], list(Media = media_data, EV = ev_data), targets)
  filtered_modules <- unique(cor_result$Module)
  
  cat("Significant correlated modules:", filtered_modules, "\n")
  # Additional path finding and SEM fitting would follow here...
}

# Execute
main()
