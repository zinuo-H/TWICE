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
                            cores = 8,
                            min_exist_KO = 1, max_exist_KO = 600) {
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
                                      min_exist_KO, max_exist_KO) {
  # Declare variable names to avoid R CMD check notes
  name <- Pathway_id <- NS <- W_pos <- W_neg <- cor <- NULL

  # Prepare correlation pools for positive and negative correlations
  cor_pool_pos <- cor_df$cor[cor_df$cor > 0]
  cor_pool_neg <- cor_df$cor[cor_df$cor < 0]
  pool_ns <- all_NS_res$NS

  # Main computation function for a single pathway
  compute_pathway <- function(path) {
    NS_res <- dplyr::filter(all_NS_res, Pathway_id == path) %>%
      dplyr::select(name, NS)
    tmp_K_num <- nrow(NS_res)

    # Identify mapped KOs for positive and negative correlations
    path_name_vec <- unique(unlist(strsplit(NS_res$name, " ")))
    cor_name_vec_pos <- unique(cor_df$name[!is.na(cor_df$cor) & cor_df$cor > 0])
    cor_name_vec_neg <- unique(cor_df$name[!is.na(cor_df$cor) & cor_df$cor < 0])
    mapped_names_pos <- intersect(path_name_vec, cor_name_vec_pos)
    mapped_names_neg <- intersect(path_name_vec, cor_name_vec_neg)
    mapped_names_str_pos <- if (length(mapped_names_pos) > 0) {
      paste(mapped_names_pos, collapse = "|")
    } else {
      NA_character_
    }
    mapped_names_str_neg <- if (length(mapped_names_neg) > 0) {
      paste(mapped_names_neg, collapse = "|")
    } else {
      NA_character_
    }

    ## --- Step A: Node aggregation ---
    node_scores <- lapply(seq_len(nrow(NS_res)), function(j) {
      node_name <- NS_res$name[j]
      node_ns <- NS_res$NS[j]
      gene_list <- unlist(strsplit(node_name, " "))
      cor_sub <- cor_df[cor_df$name %in% gene_list, , drop = FALSE]

      if (nrow(cor_sub) == 0) {
        return(data.frame(W_pos = NA, W_neg = NA, NS = NA))
      }

      W_pos <- if (any(cor_sub$cor > 0, na.rm = TRUE)) {
        mean(cor_sub$cor[cor_sub$cor > 0], na.rm = TRUE)
      } else {
        NA
      }
      W_neg <- if (any(cor_sub$cor < 0, na.rm = TRUE)) {
        mean(cor_sub$cor[cor_sub$cor < 0], na.rm = TRUE)
      } else {
        NA
      }

      data.frame(W_pos = W_pos, W_neg = W_neg, NS = node_ns)
    }) %>%
      dplyr::bind_rows() %>%
      filter(!(is.na(W_pos) & is.na(W_neg)))

    tmp_exist_K_num_pos <- sum(!is.na(node_scores$W_pos))
    tmp_exist_K_num_neg <- sum(!is.na(node_scores$W_neg))

    # Initialize result variables
    S_mean_pos <- S_sum_pos <- p_mean_pos <- NA_real_
    S_mean_neg <- S_sum_neg <- p_mean_neg <- NA_real_

    # --- Step B: Calculate positive pathway score ---
    if (tmp_exist_K_num_pos >= min_exist_KO & tmp_exist_K_num_pos <= max_exist_KO) {
      S_mean_pos <- mean(node_scores$W_pos * node_scores$NS, na.rm = TRUE)
      S_sum_pos <- S_mean_pos * tmp_exist_K_num_pos

      null_mean_pos <- replicate(perm, {
        perm_ns <- sample(pool_ns, tmp_exist_K_num_pos, replace = TRUE)
        perm_cor <- sample(cor_pool_pos, tmp_exist_K_num_pos, replace = TRUE)
        mean(perm_ns * perm_cor, na.rm = TRUE)
      })
      if (!is.na(sd(null_mean_pos)) && sd(null_mean_pos) > 0) {
        p_mean_pos <- (sum(null_mean_pos >= S_mean_pos, na.rm = TRUE) + 1) / (perm + 1)
      }
    }

    # --- Step C: Calculate negative pathway score ---
    if (tmp_exist_K_num_neg >= min_exist_KO & tmp_exist_K_num_neg <= max_exist_KO) {
      S_mean_neg <- mean(node_scores$W_neg * node_scores$NS, na.rm = TRUE)
      S_sum_neg <- S_mean_neg * tmp_exist_K_num_neg

      null_mean_neg <- replicate(perm, {
        perm_ns <- sample(pool_ns, tmp_exist_K_num_neg, replace = TRUE)
        perm_cor <- sample(cor_pool_neg, tmp_exist_K_num_neg, replace = TRUE)
        mean(perm_ns * perm_cor, na.rm = TRUE)
      })
      if (!is.na(sd(null_mean_neg)) && sd(null_mean_neg) > 0) {
        p_mean_neg <- (sum(null_mean_neg <= S_mean_neg, na.rm = TRUE) + 1) / (perm + 1)
      }
    }

    # --- Step D: Return results ---
    data.frame(
      Pathway_id = path,
      K_num = tmp_K_num,
      Exist_K_num_pos = tmp_exist_K_num_pos,
      Exist_K_pos = mapped_names_str_pos,
      S_mean_pos = S_mean_pos,
      S_sum_pos = S_sum_pos,
      p_pos = p_mean_pos,
      p_pos_adjust = NA_real_,
      Exist_K_num_neg = tmp_exist_K_num_neg,
      Exist_K_neg = mapped_names_str_neg,
      S_mean_neg = S_mean_neg,
      S_sum_neg = S_sum_neg,
      p_neg = p_mean_neg,
      p_neg_adjust = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  # Run computation (with optional parallelization)
  if (parallel) {
    pcutils::lib_ps("parallel")
    TWCE_scores_list <- parallel::mclapply(unique_paths, compute_pathway, mc.cores = cores)
  } else {
    TWCE_scores_list <- lapply(unique_paths, compute_pathway)
  }

  TWCE_scores <- dplyr::bind_rows(TWCE_scores_list)

  # Ensure all required columns exist before p-value adjustment
  required_cols <- c(
    "Pathway_id", "K_num", "Exist_K_num_pos", "Exist_K_pos", "S_mean_pos", "S_sum_pos",
    "p_pos", "p_pos_adjust", "Exist_K_num_neg", "Exist_K_neg", "S_mean_neg", "S_sum_neg",
    "p_neg", "p_neg_adjust"
  )

  for (cn in required_cols) {
    if (!cn %in% names(TWCE_scores)) {
      TWCE_scores[[cn]] <- NA_real_
    }
  }

  # Adjust p-values
  TWCE_scores$p_pos_adjust <- p.adjust(TWCE_scores$p_pos, method = p.adjust.method)
  TWCE_scores$p_neg_adjust <- p.adjust(TWCE_scores$p_neg, method = p.adjust.method)

  # Format output based on spearate_res parameter
  if (spearate_res) {
    return(TWCE_scores)
  } else {
    Exist_K_num_pos <- Exist_K_pos <- S_mean_pos <- S_sum_pos <- p_pos <- p_pos_adjust <-
      Exist_K_num_neg <- Exist_K_neg <- S_mean_neg <- S_sum_neg <- p_neg <- p_neg_adjust <- NULL
    TWCE_scores_pos <- TWCE_scores %>%
      dplyr::select(dplyr::any_of(c(
        "Pathway_id", "K_num", "Exist_K_num_pos", "Exist_K_pos",
        "S_mean_pos", "S_sum_pos", "p_pos", "p_pos_adjust"
      ))) %>%
      dplyr::rename(
        Exist_K_num = Exist_K_num_pos,
        Exist_K = Exist_K_pos,
        S_mean = S_mean_pos,
        S_sum = S_sum_pos,
        p_value = p_pos,
        p_adjust = p_pos_adjust
      ) %>%
      dplyr::mutate(Direction = "Positive")

    TWCE_scores_neg <- TWCE_scores %>%
      dplyr::select(dplyr::any_of(c(
        "Pathway_id", "K_num", "Exist_K_num_neg", "Exist_K_neg",
        "S_mean_neg", "S_sum_neg", "p_neg", "p_neg_adjust"
      ))) %>%
      dplyr::rename(
        Exist_K_num = Exist_K_num_neg,
        Exist_K = Exist_K_neg,
        S_mean = S_mean_neg,
        S_sum = S_sum_neg,
        p_value = p_neg,
        p_adjust = p_neg_adjust
      ) %>%
      dplyr::mutate(Direction = "Negative")

    return(dplyr::bind_rows(TWCE_scores_pos, TWCE_scores_neg))
  }
}


#' Internal function for calculating TWiCE scores in 'directed' mode
#'
#' @return Data frame with TWiCE scores for directed mode
#' @keywords internal
.calculate_TWICE_directed <- function(all_NS_res, cor_df, unique_paths,
                                      p.adjust.method, perm, parallel, cores,
                                      min_exist_KO, max_exist_KO) {
  # Declare variable names to avoid R CMD check notes
  name <- Pathway_id <- NS <- W <- cor <- NULL

  # Prepare correlation pool (preserving signs)
  cor_pool <- cor_df$cor
  pool_ns <- all_NS_res$NS

  # Main computation function for a single pathway
  compute_pathway <- function(path) {
    NS_res <- dplyr::filter(all_NS_res, Pathway_id == path) %>%
      dplyr::select(name, NS)
    tmp_K_num <- nrow(NS_res)

    # Identify mapped KOs
    path_name_vec <- unique(unlist(strsplit(NS_res$name, " ")))
    cor_name_vec <- unique(cor_df$name[!is.na(cor_df$cor)])
    mapped_names <- intersect(path_name_vec, cor_name_vec)
    mapped_names_str <- if (length(mapped_names) > 0) {
      paste(mapped_names, collapse = "|")
    } else {
      NA_character_
    }

    ## --- Step A: Node aggregation ---
    node_scores <- lapply(seq_len(nrow(NS_res)), function(j) {
      node_name <- NS_res$name[j]
      node_ns <- NS_res$NS[j]
      gene_list <- unlist(strsplit(node_name, " "))
      cor_sub <- cor_df[cor_df$name %in% gene_list, , drop = FALSE]

      if (nrow(cor_sub) == 0) {
        return(data.frame(W = NA, NS = NA))
      }

      W <- mean(cor_sub$cor, na.rm = TRUE)
      data.frame(W = W, NS = node_ns)
    }) %>%
      dplyr::bind_rows() %>%
      filter(!is.na(W))

    tmp_exist_K_num <- nrow(node_scores)

    # Initialize result variables
    S_mean <- S_sum <- p_mean <- NA_real_

    # --- Step B: Calculate pathway score ---
    if (tmp_exist_K_num >= min_exist_KO & tmp_exist_K_num <= max_exist_KO) {
      S_mean <- mean(node_scores$W * node_scores$NS, na.rm = TRUE)
      S_sum <- S_mean * tmp_exist_K_num

      # Null distribution generation
      null_mean <- replicate(perm, {
        perm_ns <- sample(pool_ns, tmp_exist_K_num, replace = TRUE)
        perm_cor <- sample(cor_pool, tmp_exist_K_num, replace = TRUE)
        mean(perm_ns * perm_cor, na.rm = TRUE)
      })

      # Calculate p-value based on direction of S_mean
      if (!is.na(sd(null_mean)) && sd(null_mean) > 0) {
        if (S_mean >= 0) {
          p_mean <- (sum(null_mean >= S_mean, na.rm = TRUE) + 1) / (perm + 1)
        } else {
          p_mean <- (sum(null_mean <= S_mean, na.rm = TRUE) + 1) / (perm + 1)
        }
      }
    }

    # --- Step D: Return results ---
    data.frame(
      Pathway_id = path,
      K_num = tmp_K_num,
      Exist_K_num = tmp_exist_K_num,
      Exist_K = mapped_names_str,
      S_mean = S_mean,
      S_sum = S_sum,
      p_value = p_mean,
      p_adjust = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  # Run computation (with optional parallelization)
  if (parallel) {
    pcutils::lib_ps("parallel")
    TWCE_scores_list <- parallel::mclapply(unique_paths, compute_pathway, mc.cores = cores)
  } else {
    TWCE_scores_list <- lapply(unique_paths, compute_pathway)
  }

  TWCE_scores <- dplyr::bind_rows(TWCE_scores_list)
  TWCE_scores$p_adjust <- p.adjust(TWCE_scores$p_value, method = p.adjust.method)

  return(TWCE_scores)
}


#' Internal function for calculating TWiCE scores in 'mixed' mode
#'
#' @return Data frame with TWiCE scores for mixed mode
#' @keywords internal
.calculate_TWICE_mixed <- function(all_NS_res, cor_df, unique_paths,
                                   p.adjust.method, perm, parallel, cores,
                                   min_exist_KO, max_exist_KO) {
  # Declare variable names to avoid R CMD check notes
  name <- Pathway_id <- NS <- W <- cor <- NULL

  # Prepare absolute correlation pool
  cor_pool_abs <- abs(cor_df$cor)
  pool_ns <- all_NS_res$NS

  # Main computation function for a single pathway
  compute_pathway <- function(path) {
    NS_res <- dplyr::filter(all_NS_res, Pathway_id == path) %>%
      dplyr::select(name, NS)
    tmp_K_num <- nrow(NS_res)

    # Identify mapped KOs
    path_name_vec <- unique(unlist(strsplit(NS_res$name, " ")))
    cor_name_vec <- unique(cor_df$name[!is.na(cor_df$cor)])
    mapped_names <- intersect(path_name_vec, cor_name_vec)
    mapped_names_str <- if (length(mapped_names) > 0) {
      paste(mapped_names, collapse = "|")
    } else {
      NA_character_
    }

    # Step A: Node aggregation (using absolute correlations)
    node_scores <- lapply(seq_len(nrow(NS_res)), function(j) {
      node_name <- NS_res$name[j]
      node_ns <- NS_res$NS[j]
      gene_list <- unlist(strsplit(node_name, " "))
      cor_sub <- cor_df[cor_df$name %in% gene_list, , drop = FALSE]

      if (nrow(cor_sub) == 0) {
        return(data.frame(W = NA, NS = NA))
      }

      W <- mean(abs(cor_sub$cor), na.rm = TRUE) # Take absolute values
      data.frame(W = W, NS = node_ns)
    }) %>%
      dplyr::bind_rows() %>%
      filter(!is.na(W))

    tmp_exist_K_num <- nrow(node_scores)
    S_mean <- S_sum <- p_mean <- NA_real_

    # Calculate scores if KO count is within bounds
    if (tmp_exist_K_num >= min_exist_KO & tmp_exist_K_num <= max_exist_KO) {
      S_mean <- mean(node_scores$W * node_scores$NS, na.rm = TRUE)
      S_sum <- S_mean * tmp_exist_K_num

      # Null distribution generation using absolute correlations
      null_mean <- replicate(perm, {
        perm_cor <- sample(cor_pool_abs, tmp_exist_K_num, replace = TRUE)
        perm_ns <- sample(pool_ns, tmp_exist_K_num, replace = TRUE)
        mean(perm_cor * perm_ns, na.rm = TRUE)
      })

      # Calculate one-sided p-value
      if (!is.na(sd(null_mean)) && sd(null_mean) > 0) {
        p_mean <- (sum(null_mean >= S_mean, na.rm = TRUE) + 1) / (perm + 1)
      }
    }

    # Return results
    data.frame(
      Pathway_id = path,
      K_num = tmp_K_num,
      Exist_K_num = tmp_exist_K_num,
      Exist_K = mapped_names_str,
      S_mean = S_mean,
      S_sum = S_sum,
      p_value = p_mean,
      p_adjust = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  # Run computation (with optional parallelization)
  if (parallel) {
    pcutils::lib_ps("parallel")
    TWCE_scores_list <- parallel::mclapply(unique_paths, compute_pathway, mc.cores = cores)
  } else {
    TWCE_scores_list <- lapply(unique_paths, compute_pathway)
  }

  TWCE_scores <- dplyr::bind_rows(TWCE_scores_list)
  TWCE_scores$p_adjust <- p.adjust(TWCE_scores$p_value, method = p.adjust.method)

  return(TWCE_scores)
}
