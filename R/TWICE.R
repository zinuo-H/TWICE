#' Calculate TWiCE Scores for Pathways
#'
#' This function calculates TWiCE scores, p-values, and secondary p-values for each pathway
#' based on correlation data and NS values. It supports three modes: 'mixed', 'directed', and 'separate'.
#'
#' @param all_NS_res A data frame containing pathway information with columns:
#'   Pathway_id, name, and NS.
#' @param cor_df A data frame containing correlation data with columns: name and cor.
#' @param mode Character string specifying the mode. One of 'mixed', 'directed', or 'separate'.
#'   Default is 'separate'.
#' @param universe Character vector of KO names to consider as the background universe.
#'   Default is cor_df$name.
#' @param p.adjust.method Character string specifying the method for p-value adjustment.
#'   Default is "BH".
#' @param spearate_res Logical. When mode is 'separate', determines whether to return
#'   a combined data frame with positive and negative results stacked (FALSE) or a single
#'   wide data frame (TRUE). Default is FALSE.
#' @param perm Integer. Number of permutations for null distribution generation.
#'   Default is 999.
#' @param parallel Logical. Whether to use parallel processing. Default is TRUE.
#' @param cores Integer. Number of CPU cores to use when parallel = TRUE. Default is 8.
#' @param min_exist_KO Integer. Minimum number of KOs in a pathway that must have
#'   corresponding correlation values to compute scores. Default is 1.
#' @param max_exist_KO Integer. Maximum number of KOs in a pathway that can have
#'   corresponding correlation values to compute scores. Default is 600.
#'
#' @return A data frame with columns:
#'   - Pathway_id: The ID of the pathway.
#'   - K_num: The total number of KOs in the pathway.
#'   - Exist_K_num: The number of KOs in the pathway that have corresponding correlation.
#'   - Exist_K: The names of the KOs in the pathway that have corresponding correlation,
#'     separated by "|".
#'   - S_mean: The mean TWiCE score for the pathway.
#'   - S_sum: The sum TWiCE score for the pathway (S_mean multiplied by Exist_K_num).
#'   - p_value: The p-value for the pathway based on permutation testing.
#'   - p_adjust: The adjusted p-value for the pathway.
#'   - Direction: (only for mode = 'separate' and spearate_res = FALSE) Indicates whether
#'     the result is for positive or negative correlations ("Positive" or "Negative").
#'
#' @seealso calculate_NS
#' @export
calculate_TWICE <- function(all_NS_res, cor_df,
                            mode = c("directed", "mixed", "separate")[3],
                            universe = cor_df$name,
                            p.adjust.method = "BH",
                            spearate_res = FALSE,
                            perm = 999,
                            parallel = TRUE,
                            cores = 4,
                            min_exist_KO = 5, max_exist_KO = 600) {
  # Validate input parameters
  stopifnot(mode %in% c("directed", "mixed", "separate"))
  stopifnot(is.logical(parallel))
  stopifnot(is.numeric(cores) && cores > 0)
  stopifnot(is.numeric(perm) && perm > 0)
  stopifnot(is.numeric(min_exist_KO) && min_exist_KO >= 0)
  stopifnot(is.numeric(max_exist_KO) && max_exist_KO > min_exist_KO)

  # Declare variable names to avoid R CMD check notes
  name <- Pathway_id <- NS <- p_value <- NULL

  # Filter pathways: keep only those with at least one gene in the universe
  # and at least one non-zero NS value
  filtered_NS_res <- all_NS_res %>%
    dplyr::filter(vapply(strsplit(name, " "), function(x) any(x %in% universe), logical(1))) %>%
    dplyr::group_by(Pathway_id) %>%
    dplyr::filter(!all(NS == 0)) %>% # Keep pathways with at least one NS != 0
    dplyr::ungroup()

  unique_paths <- unique(filtered_NS_res$Pathway_id)
  all_NS_res <- filtered_NS_res

  # Dispatch to appropriate mode-specific function
  if (mode == "separate") {
    result <- .calculate_TWICE_separate(
      all_NS_res = all_NS_res,
      cor_df = cor_df,
      unique_paths = unique_paths,
      p.adjust.method = p.adjust.method,
      spearate_res = spearate_res,
      perm = perm,
      parallel = parallel,
      cores = cores,
      min_exist_KO = min_exist_KO,
      max_exist_KO = max_exist_KO
    )
  } else if (mode == "directed") {
    result <- .calculate_TWICE_directed(
      all_NS_res = all_NS_res,
      cor_df = cor_df,
      unique_paths = unique_paths,
      p.adjust.method = p.adjust.method,
      perm = perm,
      parallel = parallel,
      cores = cores,
      min_exist_KO = min_exist_KO,
      max_exist_KO = max_exist_KO
    )
  } else if (mode == "mixed") {
    result <- .calculate_TWICE_mixed(
      all_NS_res = all_NS_res,
      cor_df = cor_df,
      unique_paths = unique_paths,
      p.adjust.method = p.adjust.method,
      perm = perm,
      parallel = parallel,
      cores = cores,
      min_exist_KO = min_exist_KO,
      max_exist_KO = max_exist_KO
    )
  } else {
    stop("Invalid mode specified. Choose from 'directed', 'mixed', or 'separate'.")
  }
  result <- dplyr::arrange(result, p_value)
  return(result)
}

#' Internal function for calculating TWiCE scores in 'separate' mode
#'
#' @return Data frame with TWiCE scores for separate mode
#' @keywords internal
.calculate_TWICE_separate <- function(all_NS_res, cor_df, unique_paths,
                                      p.adjust.method, spearate_res,
                                      perm, parallel, cores,
                                      min_exist_KO, max_exist_KO) { # 没有 min_coverage
  name <- Pathway_id <- NS <- W_pos <- W_neg <- cor <- NULL


  # 1. 自适应稳健修剪方差 (Adaptive Trimmed Variance via IQR)
  # 使用 Tukey 的 1.5*IQR 法则动态识别并剔除真实的生物学极值信号

  cor_pool_pos <- cor_df$cor[cor_df$cor > 0]
  cor_pool_neg <- cor_df$cor[cor_df$cor < 0]

  # --- 处理正向池 ---
  Q1_pos <- quantile(cor_pool_pos, 0.25, na.rm = TRUE)
  Q3_pos <- quantile(cor_pool_pos, 0.75, na.rm = TRUE)
  IQR_pos <- Q3_pos - Q1_pos
  # 设定正向信号上限阈值
  upper_bound_pos <- Q3_pos + 1.5 * IQR_pos

  # 只保留低于上限的纯背景噪音
  bg_pos <- cor_pool_pos[cor_pool_pos <= upper_bound_pos]
  mu_W_pos <- mean(bg_pos, na.rm = TRUE)
  var_W_pos <- var(bg_pos, na.rm = TRUE)

  # --- 处理负向池 ---
  Q1_neg <- quantile(cor_pool_neg, 0.25, na.rm = TRUE)
  Q3_neg <- quantile(cor_pool_neg, 0.75, na.rm = TRUE)
  IQR_neg <- Q3_neg - Q1_neg
  # 设定负向信号下限阈值 (因为全是负数，越小越极端)
  lower_bound_neg <- Q1_neg - 1.5 * IQR_neg

  # 只保留高于下限的纯背景噪音
  bg_neg <- cor_pool_neg[cor_pool_neg >= lower_bound_neg]
  mu_W_neg <- mean(bg_neg, na.rm = TRUE)
  var_W_neg <- var(bg_neg, na.rm = TRUE)

  compute_pathway <- function(path) {
    NS_res <- dplyr::filter(all_NS_res, Pathway_id == path) %>% dplyr::select(name, NS)
    tmp_K_num <- nrow(NS_res)

    path_name_vec <- unique(unlist(strsplit(NS_res$name, " ")))
    cor_name_vec_pos <- unique(cor_df$name[!is.na(cor_df$cor) & cor_df$cor > 0])
    cor_name_vec_neg <- unique(cor_df$name[!is.na(cor_df$cor) & cor_df$cor < 0])

    mapped_names_pos <- intersect(path_name_vec, cor_name_vec_pos)
    mapped_names_neg <- intersect(path_name_vec, cor_name_vec_neg)
    mapped_names_str_pos <- if (length(mapped_names_pos) > 0) paste(mapped_names_pos, collapse = "|") else NA_character_
    mapped_names_str_neg <- if (length(mapped_names_neg) > 0) paste(mapped_names_neg, collapse = "|") else NA_character_

    node_scores <- lapply(seq_len(nrow(NS_res)), function(j) {
      gene_list <- unlist(strsplit(NS_res$name[j], " "))
      cor_sub <- cor_df[cor_df$name %in% gene_list, , drop = FALSE]
      if (nrow(cor_sub) == 0) {
        return(data.frame(W_pos = NA, W_neg = NA, NS = NA))
      }

      W_pos <- if (any(cor_sub$cor > 0, na.rm = TRUE)) mean(cor_sub$cor[cor_sub$cor > 0], na.rm = TRUE) else NA
      W_neg <- if (any(cor_sub$cor < 0, na.rm = TRUE)) mean(cor_sub$cor[cor_sub$cor < 0], na.rm = TRUE) else NA

      data.frame(W_pos = W_pos, W_neg = W_neg, NS = NS_res$NS[j])
    }) %>%
      dplyr::bind_rows() %>%
      filter(!(is.na(W_pos) & is.na(W_neg)))

    tmp_exist_K_num_pos <- sum(!is.na(node_scores$W_pos))
    tmp_exist_K_num_neg <- sum(!is.na(node_scores$W_neg))

    S_mean_pos <- S_sum_pos <- p_mean_pos <- NA_real_
    S_mean_neg <- S_sum_neg <- p_mean_neg <- NA_real_

    # --- Z-score for Positive correlations ---
    # 仅使用 min_exist_KO 进行过滤
    if (tmp_exist_K_num_pos >= min_exist_KO & tmp_exist_K_num_pos <= max_exist_KO) {
      valid_nodes_pos <- node_scores[!is.na(node_scores$W_pos), ]
      S_mean_pos <- mean(valid_nodes_pos$W_pos * valid_nodes_pos$NS, na.rm = TRUE)
      S_sum_pos <- sum(valid_nodes_pos$W_pos * valid_nodes_pos$NS, na.rm = TRUE)

      mu_null_pos <- mu_W_pos * sum(valid_nodes_pos$NS)
      # 不除以 k_i，保留节点内基因的共表达特性
      var_null_pos <- var_W_pos * sum(valid_nodes_pos$NS^2)
      sigma_null_pos <- sqrt(var_null_pos)

      if (!is.na(sigma_null_pos) && sigma_null_pos > 0) {
        Z_pos <- (S_sum_pos - mu_null_pos) / sigma_null_pos
        p_mean_pos <- pnorm(Z_pos, lower.tail = FALSE)
      }
    }

    # --- Z-score for Negative correlations ---
    if (tmp_exist_K_num_neg >= min_exist_KO & tmp_exist_K_num_neg <= max_exist_KO) {
      valid_nodes_neg <- node_scores[!is.na(node_scores$W_neg), ]
      S_mean_neg <- mean(valid_nodes_neg$W_neg * valid_nodes_neg$NS, na.rm = TRUE)
      S_sum_neg <- sum(valid_nodes_neg$W_neg * valid_nodes_neg$NS, na.rm = TRUE)

      mu_null_neg <- mu_W_neg * sum(valid_nodes_neg$NS)
      var_null_neg <- var_W_neg * sum(valid_nodes_neg$NS^2)
      sigma_null_neg <- sqrt(var_null_neg)

      if (!is.na(sigma_null_neg) && sigma_null_neg > 0) {
        Z_neg <- (S_sum_neg - mu_null_neg) / sigma_null_neg
        p_mean_neg <- pnorm(Z_neg, lower.tail = TRUE)
      }
    }

    data.frame(
      Pathway_id = path, K_num = tmp_K_num,
      Exist_K_num_pos = tmp_exist_K_num_pos, Exist_K_pos = mapped_names_str_pos,
      S_mean_pos = S_mean_pos, S_sum_pos = S_sum_pos, p_pos = p_mean_pos, p_pos_adjust = NA_real_,
      Exist_K_num_neg = tmp_exist_K_num_neg, Exist_K_neg = mapped_names_str_neg,
      S_mean_neg = S_mean_neg, S_sum_neg = S_sum_neg, p_neg = p_mean_neg, p_neg_adjust = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  if (parallel) {
    pcutils::lib_ps("parallel")
    TWCE_scores_list <- parallel::mclapply(unique_paths, compute_pathway, mc.cores = cores)
  } else {
    TWCE_scores_list <- lapply(unique_paths, compute_pathway)
  }

  TWCE_scores <- dplyr::bind_rows(TWCE_scores_list)

  required_cols <- c(
    "Pathway_id", "K_num", "Exist_K_num_pos", "Exist_K_pos", "S_mean_pos", "S_sum_pos",
    "p_pos", "p_pos_adjust", "Exist_K_num_neg", "Exist_K_neg", "S_mean_neg", "S_sum_neg",
    "p_neg", "p_neg_adjust"
  )
  for (cn in required_cols) {
    if (!cn %in% names(TWCE_scores)) TWCE_scores[[cn]] <- NA_real_
  }

  # 2. 独立 FDR 矫正 (R 语言 p.adjust 会自动剔除 NA)

  TWCE_scores$p_pos_adjust <- p.adjust(TWCE_scores$p_pos, method = p.adjust.method)
  TWCE_scores$p_neg_adjust <- p.adjust(TWCE_scores$p_neg, method = p.adjust.method)

  if (spearate_res) {
    return(TWCE_scores)
  } else {
    Exist_K_num_pos <- Exist_K_pos <- S_mean_pos <- S_sum_pos <- p_pos <- p_pos_adjust <- S_mean <-
      Exist_K_num_neg <- Exist_K_neg <- S_mean_neg <- S_sum_neg <- p_neg <- p_neg_adjust <- NULL
    TWCE_scores_pos <- TWCE_scores %>%
      dplyr::select(dplyr::any_of(c("Pathway_id", "K_num", "Exist_K_num_pos", "Exist_K_pos", "S_mean_pos", "S_sum_pos", "p_pos", "p_pos_adjust"))) %>%
      dplyr::rename(Exist_K_num = Exist_K_num_pos, Exist_K = Exist_K_pos, S_mean = S_mean_pos, S_sum = S_sum_pos, p_value = p_pos, p_adjust = p_pos_adjust) %>%
      dplyr::mutate(Direction = "Positive") %>%
      dplyr::filter(!is.na(S_mean))

    TWCE_scores_neg <- TWCE_scores %>%
      dplyr::select(dplyr::any_of(c("Pathway_id", "K_num", "Exist_K_num_neg", "Exist_K_neg", "S_mean_neg", "S_sum_neg", "p_neg", "p_neg_adjust"))) %>%
      dplyr::rename(Exist_K_num = Exist_K_num_neg, Exist_K = Exist_K_neg, S_mean = S_mean_neg, S_sum = S_sum_neg, p_value = p_neg, p_adjust = p_neg_adjust) %>%
      dplyr::mutate(Direction = "Negative") %>%
      dplyr::filter(!is.na(S_mean))

    return(dplyr::bind_rows(TWCE_scores_pos, TWCE_scores_neg))
  }
}


#' Internal function for calculating TWiCE scores in 'directed' mode
#'
#' @return Data frame with TWiCE scores for directed mode
#' @keywords internal
#' Internal function for calculating TWiCE scores in 'directed' mode
#'
#' @return Data frame with TWiCE scores for directed mode
#' @keywords internal
.calculate_TWICE_directed <- function(all_NS_res, cor_df, unique_paths,
                                      p.adjust.method, perm, parallel, cores,
                                      min_exist_KO = 5, max_exist_KO = 600) { # 默认设为 5
  name <- Pathway_id <- NS <- W <- cor <- NULL

  # 1. 自适应稳健修剪方差 (Adaptive Trimmed Variance via IQR)
  # Directed 模式包含正负值，因此需要同时修剪上下限极值 (双尾)

  cor_pool <- cor_df$cor

  Q1 <- quantile(cor_pool, 0.25, na.rm = TRUE)
  Q3 <- quantile(cor_pool, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1

  lower_bound <- Q1 - 1.5 * IQR_val
  upper_bound <- Q3 + 1.5 * IQR_val

  # 仅提取落在上下限之间的稳健背景噪音
  bg_cor <- cor_pool[cor_pool >= lower_bound & cor_pool <= upper_bound]
  mu_W <- mean(bg_cor, na.rm = TRUE)
  var_W <- var(bg_cor, na.rm = TRUE)

  compute_pathway <- function(path) {
    NS_res <- dplyr::filter(all_NS_res, Pathway_id == path) %>% dplyr::select(name, NS)
    tmp_K_num <- nrow(NS_res)

    path_name_vec <- unique(unlist(strsplit(NS_res$name, " ")))
    cor_name_vec <- unique(cor_df$name[!is.na(cor_df$cor)])
    mapped_names <- intersect(path_name_vec, cor_name_vec)
    mapped_names_str <- if (length(mapped_names) > 0) paste(mapped_names, collapse = "|") else NA_character_

    node_scores <- lapply(seq_len(nrow(NS_res)), function(j) {
      gene_list <- unlist(strsplit(NS_res$name[j], " "))
      cor_sub <- cor_df[cor_df$name %in% gene_list, , drop = FALSE]
      if (nrow(cor_sub) == 0) {
        return(data.frame(W = NA, NS = NA))
      }

      data.frame(W = mean(cor_sub$cor, na.rm = TRUE), NS = NS_res$NS[j])
    }) %>%
      dplyr::bind_rows() %>%
      filter(!is.na(W))

    tmp_exist_K_num <- nrow(node_scores)
    S_mean <- S_sum <- p_mean <- NA_real_

    # --- Z-score 参数化检验核心 ---
    # 严格应用 min_exist_KO 过滤，过滤掉的通路 p_mean 保持为 NA
    if (tmp_exist_K_num >= min_exist_KO & tmp_exist_K_num <= max_exist_KO) {
      valid_nodes <- node_scores
      S_mean <- mean(valid_nodes$W * valid_nodes$NS, na.rm = TRUE)
      S_sum <- sum(valid_nodes$W * valid_nodes$NS, na.rm = TRUE)

      # 理论零分布的期望与方差
      mu_null <- mu_W * sum(valid_nodes$NS)
      # 不除以 k_i，保留节点内基因的共表达特征，保护特异性
      var_null <- var_W * sum(valid_nodes$NS^2)
      sigma_null <- sqrt(var_null)

      if (!is.na(sigma_null) && sigma_null > 0) {
        Z_score <- (S_sum - mu_null) / sigma_null
        # 根据得分方向计算单侧 p 值
        if (S_sum >= 0) {
          p_mean <- stats::pnorm(Z_score, lower.tail = FALSE) # 右尾
        } else {
          p_mean <- pnorm(Z_score, lower.tail = TRUE) # 左尾
        }
      }
    }

    data.frame(
      Pathway_id = path, K_num = tmp_K_num, Exist_K_num = tmp_exist_K_num,
      Exist_K = mapped_names_str, S_mean = S_mean, S_sum = S_sum,
      p_value = p_mean, p_adjust = NA_real_, stringsAsFactors = FALSE
    )
  }

  if (parallel) {
    pcutils::lib_ps("parallel")
    TWCE_scores_list <- parallel::mclapply(unique_paths, compute_pathway, mc.cores = cores)
  } else {
    TWCE_scores_list <- lapply(unique_paths, compute_pathway)
  }

  TWCE_scores <- dplyr::bind_rows(TWCE_scores_list)

  # 独立 FDR 校正，自动剔除 p_value 为 NA 的行
  TWCE_scores$p_adjust <- p.adjust(TWCE_scores$p_value, method = p.adjust.method)

  return(TWCE_scores)
}




#' Internal function for calculating TWiCE scores in 'mixed' mode
#'
#' @return Data frame with TWiCE scores for mixed mode
#' @keywords internal
#' Internal function for calculating TWiCE scores in 'mixed' mode
#'
#' @return Data frame with TWiCE scores for mixed mode
#' @keywords internal
.calculate_TWICE_mixed <- function(all_NS_res, cor_df, unique_paths,
                                   p.adjust.method, perm, parallel, cores,
                                   min_exist_KO = 5, max_exist_KO = 600) { # 默认设为 5
  name <- Pathway_id <- NS <- W <- cor <- NULL

  # 1. 自适应稳健修剪方差 (Adaptive Trimmed Variance via IQR)
  # Mixed 模式取了绝对值 (仅包含正数)，长尾仅在右侧，因此仅修剪上限 (单尾)
  cor_pool_abs <- abs(cor_df$cor)

  Q1 <- quantile(cor_pool_abs, 0.25, na.rm = TRUE)
  Q3 <- quantile(cor_pool_abs, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1

  upper_bound <- Q3 + 1.5 * IQR_val

  # 仅提取低于上限的稳健背景噪音
  bg_cor <- cor_pool_abs[cor_pool_abs <= upper_bound]
  mu_W <- mean(bg_cor, na.rm = TRUE)
  var_W <- var(bg_cor, na.rm = TRUE)

  compute_pathway <- function(path) {
    NS_res <- dplyr::filter(all_NS_res, Pathway_id == path) %>% dplyr::select(name, NS)
    tmp_K_num <- nrow(NS_res)

    path_name_vec <- unique(unlist(strsplit(NS_res$name, " ")))
    cor_name_vec <- unique(cor_df$name[!is.na(cor_df$cor)])
    mapped_names <- intersect(path_name_vec, cor_name_vec)
    mapped_names_str <- if (length(mapped_names) > 0) paste(mapped_names, collapse = "|") else NA_character_

    node_scores <- lapply(seq_len(nrow(NS_res)), function(j) {
      gene_list <- unlist(strsplit(NS_res$name[j], " "))
      cor_sub <- cor_df[cor_df$name %in% gene_list, , drop = FALSE]
      if (nrow(cor_sub) == 0) {
        return(data.frame(W = NA, NS = NA))
      }

      data.frame(W = mean(abs(cor_sub$cor), na.rm = TRUE), NS = NS_res$NS[j])
    }) %>%
      dplyr::bind_rows() %>%
      filter(!is.na(W))

    tmp_exist_K_num <- nrow(node_scores)
    S_mean <- S_sum <- p_mean <- NA_real_

    # --- Z-score 参数化检验核心 ---
    if (tmp_exist_K_num >= min_exist_KO & tmp_exist_K_num <= max_exist_KO) {
      valid_nodes <- node_scores
      S_mean <- mean(valid_nodes$W * valid_nodes$NS, na.rm = TRUE)
      S_sum <- sum(valid_nodes$W * valid_nodes$NS, na.rm = TRUE)

      mu_null <- mu_W * sum(valid_nodes$NS)
      # 同样不除以 k_i，保护特异性
      var_null <- var_W * sum(valid_nodes$NS^2)
      sigma_null <- sqrt(var_null)

      if (!is.na(sigma_null) && sigma_null > 0) {
        Z_score <- (S_sum - mu_null) / sigma_null
        # Mixed 模式由于取绝对值，得分必定为正，只计算极高值的概率（右尾）
        p_mean <- pnorm(Z_score, lower.tail = FALSE)
      }
    }

    data.frame(
      Pathway_id = path, K_num = tmp_K_num, Exist_K_num = tmp_exist_K_num,
      Exist_K = mapped_names_str, S_mean = S_mean, S_sum = S_sum,
      p_value = p_mean, p_adjust = NA_real_, stringsAsFactors = FALSE
    )
  }

  if (parallel) {
    pcutils::lib_ps("parallel")
    TWCE_scores_list <- parallel::mclapply(unique_paths, compute_pathway, mc.cores = cores)
  } else {
    TWCE_scores_list <- lapply(unique_paths, compute_pathway)
  }

  TWCE_scores <- dplyr::bind_rows(TWCE_scores_list)

  # 独立 FDR 校正
  TWCE_scores$p_adjust <- p.adjust(TWCE_scores$p_value, method = p.adjust.method)

  return(TWCE_scores)
}
