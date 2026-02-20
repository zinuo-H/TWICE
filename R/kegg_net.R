#' Calculate the NS (Node Score) for KEGG pathway networks
#'
#' @param path_index data.frame from `pathway_net_index`
#' @param lambda Numeric, a parameter for weighting the normalized degree and downstream node count. Default is `0.5`.
#'
#' @returns A data frame with additional columns for normalized degree, normalized downstream node count, and NS score.
#' @export
#'
calculate_NS <- function(path_index, lambda = 0.5) {
  Pathway_id <- Degree <- down_num <- down_num_degree <- SDegree <- Sdown_num <- Sdown_num_degree <- NS <- NULL
  # 计算下游节点数与度数的差值
  path_index$down_num_degree <- with(path_index, down_num - Out_degree)

  # 定义归一化函数
  normalize <- function(x) {
    if (sum(is.na(x)) == 0) {
      lx <- log(x + 2)
      s <- sum(lx, na.rm = TRUE)
      if (s > 0) {
        lx / s
      } else {
        rep(1 / length(x), length(x)) # 如果分母为0, 节点同等重要，平分1
      }
    } else {
      NA
    }
  }

  # 使用dplyr进行分组计算
  path_index <- path_index %>%
    dplyr::group_by(Pathway_id) %>%
    dplyr::mutate(
      Degree = as.numeric(Degree),
      SDegree = normalize(Degree),
      Sdown_num_degree = normalize(down_num_degree),
      NS = lambda * SDegree + (1 - lambda) * Sdown_num_degree
    ) %>%
    dplyr::ungroup()

  return(as.data.frame(path_index))
}

#' Update all KEGG pathway networks and calculate their NS scores
#'
#' @param pathway_xml_ls A list of KEGG pathway XML files, where each element is a `tbl_graph` or `igraph` object. If `NULL`, it will load the existing XML files.
#' @param org Character, the KEGG organism code (e.g., "hsa" for human). If `NULL`, uses "ko" as the default prefix for pathway IDs.
#' @inheritParams calculate_NS
#' @seealso calculate_NS
#' @export
#' @return No value. The function saves the NS results to a file.
update_all_pathway_NS <- function(pathway_xml_ls = NULL, org = NULL, lambda = 0.5) {
  all_pathway_index <- ReporterScore::get_all_pathway_net_index(pathway_xml_ls, org = org)
  all_NS_res <- calculate_NS(all_pathway_index, lambda = lambda)

  attributes(all_NS_res)$time <- Sys.time()
  # 保存结果
  pack_dir <- tools::R_user_dir("TWICE")
  if (!dir.exists(pack_dir)) {
    dir.create(pack_dir, recursive = TRUE)
  }
  if (is.null(org)) {
    save_path <- paste0(pack_dir, "/all_pathway_NS.rda")
  } else {
    save_path <- paste0(pack_dir, "/", org, "_all_pathway_NS.rda")
  }
  save(all_NS_res, file = save_path)
  message("All pathway NS results saved to: ", save_path)
}
#' Load all KEGG pathway NS results
#'
#' @param org Character, the KEGG organism code (e.g., "hsa" for human). If `NULL`, uses "ko" as the default prefix for pathway IDs.
#' @param verbose Logical, whether to print messages about the loading process. Default is `TRUE`.
#' @returns A data frame containing the NS scores for all pathways, with columns: Pathway_id, name, NS.
#' @export
#' @seealso calculate_NS
load_all_pathway_NS <- function(org = NULL, verbose = TRUE) {
  pack_dir <- tools::R_user_dir("TWICE")
  prefix <- "all_pathway_NS"
  if (is.null(org)) {
    file_path <- paste0(pack_dir, "/", prefix, ".rda")
  } else {
    file_path <- paste0(pack_dir, "/", org, "_", prefix, ".rda")
  }

  if (!file.exists(file_path)) {
    stop("All pathway NS results not found. Please run update_all_pathway_NS() first.")
  }

  envir <- environment()
  load(file_path, envir = envir)
  res <- get("all_NS_res", envir = envir)

  if (verbose) {
    pcutils::dabiao("load ", prefix)
    if (!is.null(attributes(res)$"time")) {
      pcutils::dabiao(paste0(prefix, " time: ", attributes(res)$"time"))
      message("If you want to update ", prefix, ", use `update_all_pathway_NS()`")
    }
  }
  return(res)
}
