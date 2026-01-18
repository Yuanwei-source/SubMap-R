# ==============================================================================
# Modern SubMap Implementation (2025 Edition)
# Optimized for: Speed (Vectorized/Parallel), Memory (No huge arrays), & Usability
# ==============================================================================

library(matrixStats)
library(parallel)
library(ComplexHeatmap)
library(circlize)

#' 内部辅助：快速计算 S2N (Signal-to-Noise Ratio)
#' @param mat 表达矩阵 (Gene x Sample)
#' @param labels 样本标签向量
#' @return 矩阵 (Gene x Group), 值为 S2N Score
calc_s2n_matrix <- function(mat, labels) {
  unique_labels <- sort(unique(labels))
  
  # 预先计算每组的 Mean 和 Sd，避免在比较时重复计算
  # 这是一个巨大的性能优化点
  res_list <- lapply(unique_labels, function(k) {
    idx <- labels == k
    # 组内 (Class k)
    m1 <- rowMeans2(mat, cols = which(idx))
    s1 <- rowSds(mat, cols = which(idx))
    
    # 组外 (Rest)
    m0 <- rowMeans2(mat, cols = which(!idx))
    s0 <- rowSds(mat, cols = which(!idx))
    
    # S2N Calculation
    # 加个极小值 1e-6 防止除以零
    (m1 - m0) / (s1 + s0 + 1e-6)
  })
  
  s2n_mat <- do.call(cbind, res_list)
  colnames(s2n_mat) <- unique_labels
  rownames(s2n_mat) <- rownames(mat)
  return(s2n_mat)
}

#' 内部辅助：计算 ES (Enrichment Score) - Kolmogorov-Smirnov 变体
#' @param rank_data 排序后的 S2N 值向量 (Desc sorted)
#' @param hit_indices Marker 基因在 rank_data 中的位置索引 (1-based)
#' @param weighted_score_type 权重 (1 = standard, 0 = unweighted)
calc_es_ks <- function(rank_data, hit_indices, weighted_score_type = 1) {
  N <- length(rank_data)
  Nh <- length(hit_indices)
  Nm <- N - Nh
  
  if (Nh == 0) return(0) # 容错
  
  # 构建 Hit 向量 (1 if hit, 0 if miss)
  # 使用 R 的高效向量操作
  
  # 计算权重
  if (weighted_score_type == 0) {
    prep_score <- rep(1, N)
  } else {
    prep_score <- abs(rank_data) ^ weighted_score_type
  }
  
  # 在 Hit 位置的值
  hit_scores <- prep_score[hit_indices]
  sum_hit_score <- sum(hit_scores)
  
  if (sum_hit_score == 0) return(0)
  
  # 计算步进 (Step)
  # Hit 时的上升步长
  step_up <- 1 / sum_hit_score
  # Miss 时的下降步长
  step_down <- 1 / Nm
  
  # 生成游走路径 (Running Sum)
  # 这是一个极度优化的技巧：我们不需要生成长度为 N 的完整向量，
  # 只需要计算在 Hit 点的变化即可推导最大值。
  # 但为了逻辑简单和兼容性，这里使用全向量逻辑，因为 N 通常不大 (<20000)
  
  indicator <- rep(-step_down, N)
  indicator[hit_indices] <- step_up * prep_score[hit_indices]
  
  running_sum <- cumsum(indicator)
  
  # SubMap 逻辑：比较绝对值最大的峰值，保留符号
  max_es <- max(running_sum)
  min_es <- min(running_sum)
  
  if (abs(max_es) > abs(min_es)) {
    return(max_es)
  } else {
    return(min_es)
  }
}

#' 执行 SubMap 单向映射 (Query -> Target)
#' @param s2n_target Target 数据的 S2N 矩阵
#' @param marker_indices_query Query 数据的 Marker 基因索引列表
#' @return 矩阵 (Query_Groups x Target_Groups) ES Score
run_unidirectional_es <- function(s2n_target, marker_indices_query, weighted_score_type = 1) {
  target_groups <- colnames(s2n_target)
  query_groups <- names(marker_indices_query)
  
  # 结果矩阵
  es_mat <- matrix(0, nrow = length(query_groups), ncol = length(target_groups),
                   dimnames = list(query_groups, target_groups))
  
  for (tg in target_groups) {
    # 对 Target 的某一列进行排序，得到基因顺序
    s2n_vals <- s2n_target[, tg]
    # 获取排序后的索引 (decreasing = TRUE)
    # sort.int 和 order 在这里是性能瓶颈，但必不可少
    ord <- order(s2n_vals, decreasing = TRUE)
    sorted_s2n <- s2n_vals[ord]
    
    # 映射回原始基因名 (为了查找 Marker 位置)
    # 优化：其实我们只需要知道 Marker 在 'ord' 中的位置
    # match(marker_gene_indices, ord) ? 
    # 更简单：在外部我们传递的是基因名，但在内部我们使用整数索引会更快。
    # 这里为了代码可读性，我在外部不做整数化，内部处理。
    
    current_gene_names <- rownames(s2n_target)[ord]
    
    for (qg in query_groups) {
      features <- marker_indices_query[[qg]] # 这是基因名
      
      # 找到 Marker 在排序列表中的位置 (Hit Indices)
      # 使用 fastmatch 或者 base match
      hits <- match(features, current_gene_names)
      hits <- hits[!is.na(hits)] # 移除未匹配的
      hits <- sort(hits)
      
      # 计算 ES
      es_mat[qg, tg] <- calc_es_ks(sorted_s2n, hits, weighted_score_type)
    }
  }
  return(es_mat)
}

#' 主函数：Modern SubMap
#' @param exp_a 表达矩阵 A (参考集/Ref), Rows=Genes, Cols=Samples
#' @param group_a 分组 A (Ref Classes)
#' @param exp_b 表达矩阵 B (目标集/Query), Rows=Genes, Cols=Samples
#' @param group_b 分组 B (Query Classes)
#' @param n_markers 用于富集分析的特征基因数 (Top N), 默认 50-100
#' @param n_perm 置换次数, 建议 1000
#' @param n_cores 并行核数 (Linux 下推荐 >1)
#' @param use_fisher 是否计算双向 Fisher 统计量 (强烈推荐 TRUE)
#' @param weighted_score_type 权重 (1 = standard, 0 = unweighted)
run_submap <- function(exp_a, group_a, exp_b, group_b, 
                       n_markers = 100, n_perm = 1000, 
                       n_cores = 8, use_fisher = TRUE, weighted_score_type = 1) {
  
  # --- 1. 数据对齐与预处理 ---
  message(">> [Preprocess] Aligning shared genes...")
  common_genes <- intersect(rownames(exp_a), rownames(exp_b))
  if (length(common_genes) < n_markers * 2) stop("Shared genes too few!")
  
  # 对齐基因
  dat_a <- as.matrix(exp_a[common_genes, ])
  dat_b <- as.matrix(exp_b[common_genes, ])
  
  # 预定义分组
  grps_a <- sort(unique(group_a))
  grps_b <- sort(unique(group_b))
  
  message(paste0("   Genes: ", length(common_genes), 
                 " | A Groups: ", length(grps_a), 
                 " | B Groups: ", length(grps_b)))

  # --- 2. 预计算 Observed S2N ---
  message(">> [Step 1] Calculating observed Signal-to-Noise Ratios (S2N)...")
  obs_s2n_a <- calc_s2n_matrix(dat_a, group_a) # Gene x GroupA
  obs_s2n_b <- calc_s2n_matrix(dat_b, group_b) # Gene x GroupB
  
  # --- 3. 提取特征基因 (Markers) ---
  # 策略：分别提取 A 和 B 中每个亚型的 Top N Markers (Names)
  get_top_markers <- function(s2n_mat, top_n) {
    res <- list()
    for (g in colnames(s2n_mat)) {
      # 排序并取 Top
      ord <- order(s2n_mat[, g], decreasing = TRUE)
      top_idx <- ord[1:top_n]
      res[[g]] <- rownames(s2n_mat)[top_idx]
    }
    return(res)
  }
  
  markers_a <- get_top_markers(obs_s2n_a, n_markers)
  markers_b <- get_top_markers(obs_s2n_b, n_markers)
  
  # --- 4. 计算 Observed ES (双向) ---
  message(">> [Step 2] Calculating Observed ES (Bi-directional)...")
  # A -> B (A 的 Marker 在 B 中富集)
  obs_es_A_on_B <- run_unidirectional_es(obs_s2n_b, markers_a, weighted_score_type)
  # B -> A (B 的 Marker 在 A 中富集)
  obs_es_B_on_A <- run_unidirectional_es(obs_s2n_a, markers_b, weighted_score_type) # 结果是 GroupB(rows) x GroupA(cols)
  
  # 为了后续 Fisher 合并，我们统一矩阵方向：Rows=A, Cols=B
  # obs_es_B_on_A 目前是 B x A，需要转置为 A x B
  obs_es_B_on_A_T <- t(obs_es_B_on_A) 
  
  # --- 5. 置换检验 (Pool Mode) ---
  message(paste0(">> [Step 3] Running Permutations (", n_perm, " times) on ", n_cores, " cores..."))
  
  # 只需要打乱 B 的标签即可同时模拟两个方向的 Null 分布：
  # 1. Null A->B: 真实的 Markers_A 映射到 打乱的 S2N_B
  # 2. Null B->A: 打乱的 Markers_B (由打乱的 S2N_B 产生) 映射到 真实的 S2N_A
  
  perm_func <- function(i) {
    # Shuffle B labels
    perm_group_b <- sample(group_b)
    
    # 1. 重新计算 B 的 S2N
    perm_s2n_b <- calc_s2n_matrix(dat_b, perm_group_b)
    
    # 2. Null A->B (A Markers are fixed, B S2N is permuted)
    perm_es_a_b <- run_unidirectional_es(perm_s2n_b, markers_a, weighted_score_type)
    
    # 3. Null B->A (B Markers must be re-derived from permuted S2N, A S2N is fixed)
    #    这里需要从 perm_s2n_b 重新提取 Markers
    perm_markers_b <- get_top_markers(perm_s2n_b, n_markers)
    perm_es_b_a <- run_unidirectional_es(obs_s2n_a, perm_markers_b, weighted_score_type) # B_perm x A_fixed
    
    # 返回向量化的结果以节省内存 (Pool 不需要知道哪个 cell 产生的 ES)
    list(
      es_a_b = as.vector(perm_es_a_b),
      es_b_a = as.vector(perm_es_b_a) # Flattened B x A
    )
  }
  
  # 并行执行
  perm_res_list <- mclapply(1:n_perm, perm_func, mc.cores = n_cores)
  
  # 收集 Null ES Pool
  pool_es_A_on_B <- unlist(lapply(perm_res_list, function(x) x$es_a_b))
  pool_es_B_on_A <- unlist(lapply(perm_res_list, function(x) x$es_b_a))
  
  # --- 6. 计算 Nominal P-values (Pool) ---
  message(">> [Step 4] Calculating Nominal P-values...")
  
  # 辅助函数：计算 P 值 (rank based)
  calc_p <- function(obs_val, null_dist) {
    # P = (sum(null >= obs) + 1) / (N_perm + 1)
    sum(null_dist >= obs_val) / length(null_dist)
  }
  
  p_mat_A_on_B <- obs_es_A_on_B
  p_mat_B_on_A_T <- obs_es_B_on_A_T
  
  for (i in seq_len(nrow(p_mat_A_on_B))) {
    for (j in seq_len(ncol(p_mat_A_on_B))) {
      p_mat_A_on_B[i, j] <- calc_p(obs_es_A_on_B[i, j], pool_es_A_on_B)
      # B->A (注意 obs_es_B_on_A_T 的索引是 [i,j] 对应 A_i, B_j)
      # 我们需要去 B_on_A 的 Pool 比较
      p_mat_B_on_A_T[i, j] <- calc_p(obs_es_B_on_A_T[i, j], pool_es_B_on_A)
    }
  }
  
  if (!use_fisher) {
    return(list(p_A_on_B = p_mat_A_on_B, p_B_on_A = p_mat_B_on_A_T))
  }
  
  # --- 7. Fisher Combination & Significance ---
  message(">> [Step 5] Computing Fisher Statistics & Final Significance...")
  
  # Fisher Stat = -2 * (ln(p1) + ln(p2))
  # 注意处理 p=0 的情况，设一个极小值
  p_min <- 1 / (length(pool_es_A_on_B) + 1)
  
  obs_p1 <- pmax(p_mat_A_on_B, p_min)
  obs_p2 <- pmax(p_mat_B_on_A_T, p_min) # 也就是 transposed B->A
  
  obs_fisher <- -2 * (log(obs_p1) + log(obs_p2))
  
  # 构建 Null Fisher Distribution (通过从两个 Pool 中随机抽样模拟)
  # 我们需要模拟 n_perm * (Rows * Cols) 那个量级的 null stats，或者直接采 n_perm * 10 次
  n_sim_fisher <- 20000 
  rand_p1 <- sample(pool_es_A_on_B, n_sim_fisher, replace = TRUE)
  rand_p2 <- sample(pool_es_B_on_A, n_sim_fisher, replace = TRUE)
  
  # 将 Null ES 转换为 P 值 (Self-rank in pool)
  # 为加速计算，直接在该分布中用 Rank/N 近似 P
  # 但更简单的方法是：直接利用 calculate 出来的 P-pool
  # 为了严谨，我们应该把 Pool ES 转换为 Pool P 
  null_p1 <- (rank(-rand_p1) / length(rand_p1))
  null_p2 <- (rank(-rand_p2) / length(rand_p2))
  
  null_fisher <- -2 * (log(pmax(null_p1, p_min)) + log(pmax(null_p2, p_min)))
  
  # 计算最终 SA P-value
  final_sa_p <- obs_fisher
  for (i in seq_len(nrow(final_sa_p))) {
    for (j in seq_len(ncol(final_sa_p))) {
      # Fisher Stat 越大越显著 (P 越小 -> -logP 越大)
      # P = (sum(Null_Stat >= Obs_Stat) + 1) 
      final_sa_p[i, j] <- sum(null_fisher >= obs_fisher[i, j]) / length(null_fisher)
    }
  }
  
  # FDR Correction
  final_sa_fdr <- matrix(p.adjust(as.vector(final_sa_p), method = "fdr"),
                         nrow = nrow(final_sa_p),
                         dimnames = dimnames(final_sa_p))
  
  # Bonferroni Correction
  final_sa_bonf <- matrix(p.adjust(as.vector(final_sa_p), method = "bonferroni"),
                          nrow = nrow(final_sa_p),
                          dimnames = dimnames(final_sa_p))
  
  message(">> Done.")
  
  return(list(
    sa_p_matrix = final_sa_p,
    sa_fdr_matrix = final_sa_fdr,
    sa_bonf_matrix = final_sa_bonf,
    fisher_stats = obs_fisher,
    raw_p_A_on_B = p_mat_A_on_B,
    raw_p_B_on_A = p_mat_B_on_A_T
  ))
}

#' 可视化 SubMap 结果 (发表级 Heatmap)
#' @param submap_res run_submap 的返回结果 List
#' @param threshold 显著性阈值 (默认 0.05)
#' @param row_title Y轴标题 (Dataset A)
#' @param col_title X轴标题 (Dataset B)
#' @param correction_method 显著性标注使用的校正方法: "Bonferroni", "FDR", "Raw" (默认 "Bonferroni")
plot_submap_heatmap <- function(submap_res, threshold = 0.05, 
                                row_title = "Dataset A", col_title = "Dataset B",
                                correction_method = "Bonferroni") {
  
  p_mat <- submap_res$sa_p_matrix
  
  # 绘图值： -log10 P-value
  # 为了防止 Inf，限制最大值
  plot_val <- -log10(pmax(p_mat, 1e-5))
  
  # 颜色映射： 0 -> 1.3 (-log10(0.05)) -> 3 (-log10(0.001))
  col_fun <- colorRamp2(c(0, -log10(threshold), 3), c("white", "orange", "firebrick"))
  
  ht <- Heatmap(plot_val,
          name = "-log10(P)",
          col = col_fun,
          cluster_rows = FALSE, 
          cluster_columns = FALSE,
          rect_gp = gpar(col = "gray90", lwd = 1),
          row_names_side = "left",
          column_names_side = "bottom",
          column_title = col_title,
          row_title = row_title,
          cell_fun = function(j, i, x, y, width, height, fill) {
            # 显著性标注
            val_to_check <- 1
            if (correction_method == "Bonferroni") {
              val_to_check <- submap_res$sa_bonf_matrix[i, j]
            } else if (correction_method == "FDR") {
              val_to_check <- submap_res$sa_fdr_matrix[i, j]
            } else {
              val_to_check <- submap_res$sa_p_matrix[i, j]
            }
            
            txt <- ""
            if (val_to_check < 0.01) {
              txt <- "**"
            } else if (val_to_check < threshold) {
              txt <- "*"
            }
            
            grid.text(txt, x, y, gp = gpar(fontsize = 16, fontface = "bold"))
          }
  )
  draw(ht)
}
