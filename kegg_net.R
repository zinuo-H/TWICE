#' Download KEGG pathway XML files and create networks
#'
#' @param download_dir Directory to save the downloaded XML files.
#' @param org kegg organism, listed in https://www.genome.jp/kegg/catalog/org_list.html such as "hsa", default NULL means ko.
#'
#' @returns No value
#' @export
#'
update_pathway_xml_ls <- function(download_dir, org = NULL) {
  # 检查依赖包
  if (!requireNamespace("ggkegg", quietly = TRUE)) {
    stop("Package 'ggkegg' is required but not installed. Run: install.packages('ggkegg')")
  }

  # 创建目录
  if (!dir.exists(download_dir)) {
    dir.create(download_dir, recursive = TRUE)
  }
  pack_dir <- tools::R_user_dir("ReporterScore")
  if (!dir.exists(pack_dir)) dir.create(pack_dir, recursive = TRUE)

  # 加载通路ID表
  if (is.null(org)) {
    Pathway_htable <- load_Pathway_htable()
    pathway_ids <- gsub("map", "ko", Pathway_htable$Pathway_id)
  } else {
    org_pathway <- KEGGREST::keggList("pathway", org) # 获取'KEGG'数据库中所有人类通路
    if (is.null(org_pathway)) {
      stop("No pathways found for organism: ", org)
    }
    pathway_ids <- names(org_pathway)
  }

  # 初始化结果列表和失败记录
  pathway_xml_ls <- list()
  failed_paths <- character()

  # 遍历所有通路ID
  for (path in pathway_ids) {
    message("Processing pathway: ", path)

    # 直接尝试，报错则跳过
    result <- tryCatch(
      {
        pathway_xml_ls[[path]] <- ggkegg::pathway(path, directory = download_dir)
        message("Success: ", path)
      },
      error = function(e) {
        failed_paths <<- c(failed_paths, path)
        message("Failed to download ", path, ": ", e$message)
        return(NULL)
      }
    )
  }

  # 保存结果（含失败记录）
  attributes(pathway_xml_ls)$download_time <- Sys.time()
  attributes(pathway_xml_ls)$failed_paths <- failed_paths
  if (is.null(org)) {
    save_path <- paste0(pack_dir, "/pathway_xml_ls.rda")
  } else {
    save_path <- paste0(pack_dir, "/", org, "_pathway_xml_ls.rda")
  }

  save(pathway_xml_ls, file = save_path)
  message("\nResults saved to: ", save_path)

  # 打印失败汇总
  if (length(failed_paths) > 0) {
    warning(
      "Failed to download ", length(failed_paths), " pathways:\n",
      paste(failed_paths, collapse = ", ")
    )
  }
}

#' Load KEGG pathway XML file list
#'
#' @inheritParams update_pathway_xml_ls
#' @param verbose Logical, whether to print messages about the loading process. Default is `TRUE`.
#'
#' @returns A list of KEGG pathway XML files, where each element is a `tbl_graph` or `igraph` object.
#' @export
load_pathway_xml_ls <- function(org = NULL, verbose = TRUE) {
  # 加载KEGG通路XML文件列表
  pack_dir <- tools::R_user_dir("ReporterScore")
  prefix <- "pathway_xml_ls"

  if (is.null(org)) {
    file_path <- paste0(pack_dir, "/", prefix, ".rda")
  } else {
    file_path <- paste0(pack_dir, "/", org, "_", prefix, ".rda")
  }

  if (!file.exists(file_path)) {
    stop("Pathway XML list not found. Please run update_pathway_xml_ls() first.")
  }

  envir <- environment()
  if (file.exists(file_path)) {
    load(file_path, envir = envir)
  } else {
    message("Not find ", prefix, ", please run `update_pathway_xml_ls()` first!")
    return(invisible())
  }
  res <- get(prefix, envir = envir)

  if (verbose) {
    pcutils::dabiao("load ", prefix)
    if (!is.null(attributes(res)$"download_time")) {
      pcutils::dabiao(paste0(prefix, " download time: ", attributes(res)$"download_time"))
      message("If you want to update ", prefix, ", use `update_pathway_xml_ls()`")
    }
  }
  return(res)
}

#' Create a network from KEGG pathway XML files
#'
#' @param pathway_xml A `tbl_graph` or `igraph` object, or a file path to a KEGG XML file.
#'
#' @returns A `metanet` object representing the pathway network.
#' @export
#' @seealso plot_pathway_net
c_net_from_pathway_xml <- function(pathway_xml) {
  lib_ps("MetaNet", "igraph", "ggkegg", library = FALSE)
  name <- NULL
  if (is.null(pathway_xml)) {
    return(NULL)
  }
  if (identical(class(pathway_xml), c("tbl_graph", "igraph"))) {
    g <- pathway_xml
  } else if (is.character(pathway_xml) && file.exists(pathway_xml)) {
    pid <- basename(pathway_xml)
    pid <- gsub("\\.xml$", "", pid) # Remove .xml extension if present
    g <- ggkegg::pathway(pid, directory = dirname(pathway_xml))
  } else {
    stop("pathway_xml must be a tbl_graph, igraph object or a file path to a KEGG XML file.")
  }

  node1 <- as.data.frame(g)
  pathway_id <- unique(node1$pathway_id)[1]
  edge1 <- igraph::as_data_frame(g, what = "edges")

  path_gene_compound <- dplyr::filter(node1, type %in% c("ortholog", "gene", "compound")) %>%
    dplyr::distinct(name, .keep_all = TRUE)
  if (nrow(edge1) > 0) {
    path_net <- dplyr::filter(
      edge1, from %in% path_gene_compound$name,
      to %in% path_gene_compound$name
    )
  } else {
    path_net <- edge1
    path_net$subtype_name <- character()
  }
  MetaNet::c_net_from_edgelist(path_net, vertex_df = path_gene_compound, direct = TRUE) -> path_net_c
  # igraph::V(path_net_c)$label <- gsub(" .*", "", igraph::V(path_net_c)$name)
  igraph::V(path_net_c)$label <- igraph::V(path_net_c)$graphics_name
  MetaNet::c_net_set(path_net_c,
    vertex_group = "type", vertex_class = "type",
    edge_type = "subtype_name"
  ) -> path_net_c

  attributes(path_net_c)$pathway_id <- pathway_id
  path_net_c
}

#' Plot a KEGG pathway network
#'
#' @param path_net_c A `metanet` object representing the pathway network.
#' @param ... Additional arguments passed to `MetaNet::c_net_plot`.
#' @param simplify Logical, whether to simplify the network by removing multiple edges and loops. Default is `FALSE`.
#' @param plot_depth Logical, whether to plot the network as a tree layout. Default is `FALSE`.
#'
#' @returns A plot of the pathway network.
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("MetaNet") && requireNamespace("ggkegg")) {
#'   tmp_dir <- tempdir()
#'   pcutils::download2("https://rest.kegg.jp/get/ko01521/kgml", file.path(tmp_dir, "ko01521.xml"))
#'   path_net_c <- c_net_from_pathway_xml(file.path(tmp_dir, "ko01521.xml"))
#'   plot_pathway_net(path_net_c)
#'   pathway_net_index(path_net_c)
#' }
#' }
plot_pathway_net <- function(path_net_c, simplify = FALSE, plot_depth = FALSE, ...) {
  lib_ps("MetaNet", "igraph", library = FALSE)
  name <- NULL

  if (simplify) igraph::simplify(path_net_c, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "first") -> path_net_c

  default_arg <- list(
    labels_num = "all", main = attributes(path_net_c)$pathway_id,
    vertex.shape = c("rectangle", "circle"), vertex_size_range = list(c(12, 12), c(8, 8))
  )

  if (plot_depth) {
    default_arg <- append(default_arg, list(
      coors = igraph::as_tree(), rescale = TRUE
    ))
    do.call(MetaNet::c_net_plot, append(
      list(go = path_net_c),
      pcutils::update_param(default_arg, list(...))
    ))
  }
  MetaNet::get_v(path_net_c) -> tmp_v
  if (any(is.na(tmp_v$x)) || any(is.na(tmp_v$y))) {
    warning("Some node coordinates are NA. This may affect the plot layout.")
  }
  dplyr::filter(tmp_v, !is.na(x) & !is.na(y)) -> tmp_v
  MetaNet::c_net_filter(path_net_c, name %in% tmp_v$name) -> path_net_c
  do.call(MetaNet::c_net_plot, append(
    list(go = path_net_c),
    pcutils::update_param(default_arg, list(...))
  ))
}

#' Create a pathway network index data.frame
#'
#' @param path_net_c A `tbl_graph` or `igraph` object representing a KEGG pathway network.
#'
#' @returns A data frame containing vertex indices and their attributes, including in-degree and out-degree.
#' @export
#' @seealso plot_pathway_net
pathway_net_index <- function(path_net_c) {
  lib_ps("MetaNet", "igraph", library = FALSE)
  message("Removing multiple edges and loops from the network.")
  igraph::simplify(path_net_c, remove.multiple = TRUE, remove.loops = TRUE) -> path_net_c2

  igraph::vertex.attributes(path_net_c2) <- igraph::vertex.attributes(path_net_c2)[c("name", "type")]
  MetaNet::net_par(path_net_c2, mode = "v")$v_index -> path_net_index
  path_net_index$In_degree <- igraph::degree(path_net_c2, mode = "in")
  path_net_index$Out_degree <- igraph::degree(path_net_c2, mode = "out")
  attributes(path_net_index)$pathway_id <- attributes(path_net_c)$pathway_id

  # 计算下游节点数
  for (node in path_net_index$name) {
    neighbors <- MetaNet::c_net_neighbors(path_net_c2, node, order = 1000, mode = "out")
    path_net_index$down_num[path_net_index$name == node] <- length(neighbors) - 1
  }

  # 计算深度
  igraph::layout_as_tree(path_net_c) -> layout_tree

  data.frame(
    row.names = NULL,
    Pathway_id = attributes(path_net_c)$pathway_id,
    path_net_index,
    Depth = layout_tree[, 2]
  )
}


#' Create an index for all KEGG pathway networks
#'
#' @param pathway_xml_ls A list of KEGG pathway XML files, where each element is a `tbl_graph` or `igraph` object.
#' @param org Character, the KEGG organism code (e.g., "hsa" for human). If `NULL`, uses "ko" as the default prefix for pathway IDs.
#'
#' @returns A data frame containing indices for all pathways, including pathway IDs and their attributes.
#' @export
#'
get_all_pathway_net_index <- function(pathway_xml_ls = NULL, org = NULL) {
  lib_ps("MetaNet", "igraph", library = FALSE)
  if (is.null(pathway_xml_ls)) pathway_xml_ls <- load_pathway_xml_ls(org = org, verbose = FALSE)
  level2_name <- NULL

  load_Pathway_htable(verbose = FALSE) -> Pathway_htable
  tmp_prefix <- "ko"
  if (!is.null(org)) {
    tmp_prefix <- org
  }
  Pathway_htable$Pathway_id <- gsub("map", tmp_prefix, Pathway_htable$Pathway_id)
  Pathway_htable2 <- filter(Pathway_htable, !level2_name %in% c("Global and overview maps"))

  intersect(Pathway_htable2$Pathway_id, names(pathway_xml_ls)) -> pathway_ids

  pathway_net <- list()
  pathway_index_list <- list()
  for (path in pathway_ids) {
    message("Processing pathway: ", path)
    pathway_net[[path]] <- c_net_from_pathway_xml(pathway_xml_ls[[path]])
    if (is.null(pathway_net[[path]])) {
      message("Failed to create network for pathway: ", path)
      next
    }
    pathway_index_list[[path]] <- pathway_net_index(pathway_net[[path]])
  }
  do.call(rbind, pathway_index_list) -> all_pathway_index
  all_pathway_index
}


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
      lx <- log(x+2)
      s <- sum(lx, na.rm = TRUE) 
      if  (s > 0) {
        lx / s
      } else {
        rep(1/length(x), length(x)) # 如果分母为0, 节点同等重要，平分1
      }
    } else {
      NA
    }
  }

  # 使用dplyr进行分组计算
  path_index <- path_index %>%
    group_by(Pathway_id) %>%
    mutate(
      Degree = as.numeric(Degree),
      SDegree = normalize(Degree),
      Sdown_num_degree = normalize(down_num_degree),
      NS = lambda * SDegree + (1 - lambda) * Sdown_num_degree
    ) %>%
    ungroup()

  return(as.data.frame(path_index))
}

#' Update all KEGG pathway networks and calculate their NS scores
#'
#' @param pathway_xml_ls A list of KEGG pathway XML files, where each element is a `tbl_graph` or `igraph` object. If `NULL`, it will load the existing XML files.
#' @param org Character, the KEGG organism code (e.g., "hsa" for human). If `NULL`, uses "ko" as the default prefix for pathway IDs.
#' @inheritParams calculate_NS
#' @seealso calculate_NS
#' @export
update_all_pathway_NS <- function(pathway_xml_ls = NULL, org = NULL, lambda = 0.5) {
  all_pathway_index <- get_all_pathway_net_index(pathway_xml_ls, org = org)
  all_NS_res <- calculate_NS(all_pathway_index, lambda = lambda)

  attributes(all_NS_res)$time <- Sys.time()
  # 保存结果
  pack_dir <- tools::R_user_dir("ReporterScore")
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
  pack_dir <- tools::R_user_dir("ReporterScore")
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

#' Calculate TWiCE Scores for Pathways
#'
#' This function calculates TWiCE scores, p-values, and secondary p-values for each pathway
#' based on correlation data and NS values.
#'
#' @param all_NS_res A data frame containing pathway information with columns:
#'   Pathway_id, name, and NS.
#' @param cor_df A data frame containing correlation data with columns: name and cor.
#' @param mode 'mixed', 'directed' or 'separate'. Default is 'separate'.
#' @param p.adjust.method The method used for p-value adjustment (default: "BH").
#' @param spearate_res Logical, whether to return separate positive and negative results when mode is 'separate'. Default is FALSE.
#' @param min_exist_KO min exist KO number in a pathway (default, 1, when a pathway contains KOs less than 1, there will be no RS)
#' @param max_exist_KO max exist KO number in a pathway (default, 600, when a pathway contains KOs more than 600, there will be no RS)
#' @return A data frame with columns:.
#' @seealso calculate_NS
calculate_TWICE <- function(all_NS_res, cor_df, 
                            mode = c("directed", "mixed", "separate")[3],
                            universe = cor_df$name,
                            p.adjust.method = "BH",
                            spearate_res = FALSE,
                            perm = 999,
                            parallel = TRUE,
                            cores = 8, 
                            min_exist_KO = 1, max_exist_KO = 600
) {
  name <- Pathway_id <- NS <- W_pos <- W_neg <- W <- cor <- NULL
  
  filtered_NS_res <- all_NS_res %>%
    filter(sapply(strsplit(name, " "), function(x) any(x %in% universe))) %>%
    group_by(Pathway_id) %>%
    filter(!all(NS == 0)) %>% # 保留那些至少有一个 NS ≠ 0 的 Pathway
    ungroup()
  
  unique_paths <- unique(filtered_NS_res$Pathway_id)
  n_paths <- length(unique_paths)
  all_NS_res <- filtered_NS_res
  
  # 模式	含义	处理方式
  # "directed"	保留正负号，允许相互抵消	使用原始 correlation 值（正为正，负为负）
  # "mixed"	忽略方向性，只看绝对相关性	对 correlation 取绝对值
  # "separate"	分别计算正相关和负相关通路分数	正负方向分开算、分别置换
  
  # ============SEPARATE==============
  
  if (mode == "separate") {
    cor_pool_pos <- cor_df$cor[cor_df$cor > 0] 
    cor_pool_neg <- cor_df$cor[cor_df$cor < 0]
    pool_ns <- all_NS_res$NS
    
    # Step 1. 主体计算函数（单个 Pathway）
    compute_pathway <- function(path) {
      NS_res <- dplyr::filter(all_NS_res, Pathway_id == path) %>%
        dplyr::select(name, NS)
      tmp_K_num <- nrow(NS_res)
      
      path_name_vec  <- unique(unlist(strsplit(NS_res$name, " ")))
      cor_name_vec_pos  <- unique(cor_df$name[!is.na(cor_df$cor) & cor_df$cor > 0])
      cor_name_vec_neg  <- unique(cor_df$name[!is.na(cor_df$cor) & cor_df$cor < 0])
      mapped_names_pos   <- intersect(path_name_vec, cor_name_vec_pos)
      mapped_names_neg   <- intersect(path_name_vec, cor_name_vec_neg)
      mapped_names_str_pos <- if (length(mapped_names_pos) > 0) {
        paste(mapped_names_pos, collapse = "｜")
      } else {
        NA_character_
      }
      mapped_names_str_neg <- if (length(mapped_names_neg) > 0) {
        paste(mapped_names_neg, collapse = "｜")
      } else {
        NA_character_
      }
      
      ## --- Step A: 节点聚合 ---
      node_scores <- lapply(seq_len(nrow(NS_res)), function(j) {
        node_name <- NS_res$name[j]
        node_ns   <- NS_res$NS[j]
        gene_list <- unlist(strsplit(node_name, " "))
        cor_sub <- cor_df[cor_df$name %in% gene_list, , drop = FALSE] # 取这些基因在 cor_df 的相关性
        if (nrow(cor_sub) == 0)
          return(data.frame(W_pos = NA, W_neg = NA, NS = NA))
        
        W_pos <- if (any(cor_sub$cor > 0, na.rm = TRUE))
          mean(cor_sub$cor[cor_sub$cor > 0], na.rm = TRUE) else NA
        W_neg <- if (any(cor_sub$cor < 0, na.rm = TRUE))
          mean(cor_sub$cor[cor_sub$cor < 0], na.rm = TRUE) else NA
        
        data.frame(W_pos = W_pos, W_neg = W_neg, NS = node_ns)
      }) %>% bind_rows() %>%
        filter(!(is.na(W_pos) & is.na(W_neg)))
      
      tmp_exist_K_num_pos <- sum(!is.na(node_scores$W_pos))
      tmp_exist_K_num_neg <- sum(!is.na(node_scores$W_neg))
      
      # 初始化结果变量
      S_mean_pos <- S_sum_pos <- p_mean_pos <- NA_real_
      S_mean_neg <- S_sum_neg <- p_mean_neg <- NA_real_
      
      # --- Step B: 计算正向通路得分 ---
      if (tmp_exist_K_num_pos >= min_exist_KO & tmp_exist_K_num_pos <= max_exist_KO) {
        S_mean_pos <- mean(node_scores$W_pos * node_scores$NS, na.rm = TRUE)
        S_sum_pos  <- S_mean_pos * tmp_exist_K_num_pos
        
        null_mean_pos <- replicate(perm, {
          perm_ns=sample(pool_ns, tmp_exist_K_num_pos, replace = TRUE)
          perm_cor <- sample(cor_pool_pos, tmp_exist_K_num_pos, replace = TRUE)
          mean(perm_ns * perm_cor ,na.rm = TRUE)
        })
        if (!is.na(sd(null_mean_pos)) && sd(null_mean_pos) > 0) {
          p_mean_pos <- (sum(null_mean_pos >= S_mean_pos, na.rm = TRUE) + 1) / (perm + 1)
          Z_mean_pos <- (S_mean_pos - mean(null_mean_pos, na.rm = TRUE)) / sd(null_mean_pos, na.rm = TRUE)
        }
      }
      
      # --- Step C:  计算负向通路得分 ---
      if (tmp_exist_K_num_neg >= min_exist_KO & tmp_exist_K_num_neg <= max_exist_KO) {
        S_mean_neg <- mean(node_scores$W_neg * node_scores$NS, na.rm = TRUE)
        S_sum_neg  <- S_mean_neg * tmp_exist_K_num_neg
        
        null_mean_neg <- replicate(perm, {
          perm_ns=sample(pool_ns, tmp_exist_K_num_neg, replace = TRUE)
          perm_cor <- sample(cor_pool_neg, tmp_exist_K_num_neg, replace = TRUE)
          mean(perm_ns * perm_cor ,na.rm = TRUE)
        })
        if (!is.na(sd(null_mean_neg)) && sd(null_mean_neg) > 0) {
          p_mean_neg <- (sum(null_mean_neg <= S_mean_neg, na.rm = TRUE) + 1) / (perm + 1)
        }
      }
      
      # --- Step D: 返回结果 ---
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
    
    # Step 4: 运行（可选并行）
    if (parallel) {
      library(parallel)
      TWCE_scores_list <- mclapply(unique_paths, compute_pathway, mc.cores = cores)
    } else {
      TWCE_scores_list <- lapply(unique_paths, compute_pathway)
    }
    
    TWCE_scores <- bind_rows(TWCE_scores_list)
    
    # Step 5: 调整 p 值（先兜底缺列，避免 p.adjust 报错）
    # 保证这些列一定存在
    for (cn in c( "Pathway_id" ,"K_num",  "Exist_K_num_pos", "Exist_K_pos",  "S_mean_pos" , "S_sum_pos",    
                  "p_pos",    "p_pos_adjust" ,"Exist_K_num_neg" ,"Exist_K_neg",  "S_mean_neg",  "S_sum_neg",   
                 "p_neg" ,  "p_neg_adjust")) {
      if (!cn %in% names(TWCE_scores)) {
        # p_* 和 Z_*、S_* 是数值列
        TWCE_scores[[cn]] <- NA_real_
      }
    }
    # 再做 p.adjust（对全 NA 列 p.adjust 也能返回 NA）
    TWCE_scores$p_pos_adjust <- p.adjust(TWCE_scores$p_pos, method = p.adjust.method)
    TWCE_scores$p_neg_adjust <- p.adjust(TWCE_scores$p_neg, method = p.adjust.method)
    
    # Step 6: 输出（用 any_of 防止列缺失时报错）
    if (spearate_res) {
      return(TWCE_scores)
    } else {
      TWCE_scores_pos <- TWCE_scores %>%
        dplyr::select(dplyr::any_of(c("Pathway_id","K_num","Exist_K_num_pos","Exist_K_pos",
                                      "S_mean_pos","S_sum_pos","p_pos","p_pos_adjust"))) %>%
        dplyr::rename(Exist_K_num = Exist_K_num_pos,
                      Exist_K = Exist_K_pos,
                      S_mean = S_mean_pos,
                      S_sum  = S_sum_pos,
                      p_value= p_pos,
                      p_adjust = p_pos_adjust) %>%
        dplyr::mutate(Direction = "Positive")
      
      TWCE_scores_neg <- TWCE_scores %>%
        dplyr::select(dplyr::any_of(c("Pathway_id","K_num","Exist_K_num_neg","Exist_K_neg",
                                      "S_mean_neg","S_sum_neg","p_neg","p_neg_adjust"))) %>%
        dplyr::rename(Exist_K_num = Exist_K_num_neg,
                      Exist_K = Exist_K_neg,
                      S_mean = S_mean_neg,
                      S_sum  = S_sum_neg,
                      p_value= p_neg,
                      p_adjust = p_neg_adjust) %>%
        dplyr::mutate(Direction = "Negative")
      
      return(dplyr::bind_rows(TWCE_scores_pos, TWCE_scores_neg))
    }
  }
  
  # ==========DIRECTED============
  if (mode == "directed") {
    cor_pool <- cor_df$cor
    pool_ns <- all_NS_res$NS
    
    # Step 1. 主体计算函数（单个 Pathway）
    compute_pathway <- function(path) {
      NS_res <- dplyr::filter(all_NS_res, Pathway_id == path) %>%
        dplyr::select(name, NS)
      tmp_K_num <- nrow(NS_res)
      
      path_name_vec  <- unique(unlist(strsplit(NS_res$name, " ")))
      cor_name_vec  <- unique(cor_df$name[!is.na(cor_df$cor)])
      mapped_names <- intersect(path_name_vec, cor_name_vec)
      mapped_names_str <- if (length(mapped_names) > 0) {
        paste(mapped_names, collapse = "｜")
      } else {
        NA_character_
      }
      
      ## --- Step A: 节点聚合 ---
      node_scores <- lapply(seq_len(nrow(NS_res)), function(j) {
        node_name <- NS_res$name[j]
        node_ns   <- NS_res$NS[j]
        gene_list <- unlist(strsplit(node_name, " "))
        
        cor_sub <- cor_df[cor_df$name %in% gene_list, , drop = FALSE] 
        
        if (nrow(cor_sub) == 0)
          return(data.frame(W = NA, NS = NA))
        
        W <- mean(cor_sub$cor, na.rm = TRUE)
        
        data.frame(W = W, NS = node_ns)
      }) %>% bind_rows() %>%
        filter(!is.na(W))
      
      tmp_exist_K_num <- nrow(node_scores)
      
      # 初始化结果变量
      S_mean <- S_sum  <- p_mean <- NA_real_
      
      # --- Step B: 计算通路得分 ---
      if (tmp_exist_K_num >= min_exist_KO & tmp_exist_K_num <= max_exist_KO) {
        S_mean <- mean(node_scores$W * node_scores$NS, na.rm = TRUE)
        S_sum  <- S_mean * tmp_exist_K_num
        
        # Null: 从全局cor中重采样（保留正负）
        null_mean <- replicate(perm, {
          perm_ns=sample(pool_ns, tmp_exist_K_num, replace = TRUE)
          perm_cor <- sample(cor_pool, tmp_exist_K_num, replace = TRUE)
          mean(perm_ns * perm_cor ,na.rm = TRUE)
        })
        if (!is.na(sd(null_mean)) && sd(null_mean) > 0) {
          if (S_mean >= 0) {
            p_mean <- (sum(null_mean >= S_mean, na.rm = TRUE) + 1) / (perm + 1)
          } else {
            p_mean <- (sum(null_mean <= S_mean, na.rm = TRUE) + 1) / (perm + 1)
          }
        }
      }
      # --- Step D: 返回结果 ---
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
    
    # Step 4: 运行（可选并行）
    if (parallel) {
      library(parallel)
      TWCE_scores_list <- mclapply(unique_paths, compute_pathway, mc.cores = cores)
    } else {
      TWCE_scores_list <- lapply(unique_paths, compute_pathway)
    }
    
    TWCE_scores <- bind_rows(TWCE_scores_list)
    TWCE_scores$p_adjust <- p.adjust(TWCE_scores$p_value, method = p.adjust.method)
    
    return(TWCE_scores)
  }
  # ==========MIXED============
  if (mode == "mixed") {
    cor_pool_abs <- abs(cor_df$cor)
    pool_ns <- all_NS_res$NS
    
    compute_pathway <- function(path) {
      NS_res <- dplyr::filter(all_NS_res, Pathway_id == path) %>%
        dplyr::select(name, NS)
      tmp_K_num <- nrow(NS_res)
      
      path_name_vec  <- unique(unlist(strsplit(NS_res$name, " ")))
      cor_name_vec  <- unique(cor_df$name[!is.na(cor_df$cor)])
      mapped_names <- intersect(path_name_vec, cor_name_vec)
      mapped_names_str <- if (length(mapped_names) > 0) {
        paste(mapped_names, collapse = "｜")
      } else {
        NA_character_
      }
      
      # Step A: 节点聚合
      node_scores <- lapply(seq_len(nrow(NS_res)), function(j) {
        node_name <- NS_res$name[j]
        node_ns   <- NS_res$NS[j]
        gene_list <- unlist(strsplit(node_name, " "))
        cor_sub <- cor_df[cor_df$name %in% gene_list, , drop = FALSE]
        if (nrow(cor_sub) == 0)
          return(data.frame(W = NA, NS = NA))
        W <- mean(abs(cor_sub$cor), na.rm = TRUE)  # 取绝对值
        data.frame(W = W, NS = node_ns)
      }) %>% bind_rows() %>% filter(!is.na(W))
      
      tmp_exist_K_num <- nrow(node_scores)
      S_mean <- S_sum <- p_mean <- NA_real_
      
      if (tmp_exist_K_num >= min_exist_KO & tmp_exist_K_num <= max_exist_KO) {
        S_mean <- mean(node_scores$W * node_scores$NS, na.rm = TRUE)
        S_sum  <- S_mean * tmp_exist_K_num
        
        # Null: 从绝对相关池中重采样
        null_mean <- replicate(perm, {
          perm_cor <- sample(cor_pool_abs, tmp_exist_K_num, replace = TRUE)
          perm_ns <- sample(pool_ns, tmp_exist_K_num, replace = TRUE)
          mean(perm_cor * perm_ns, na.rm = TRUE)
        })
        
        if (!is.na(sd(null_mean)) && sd(null_mean) > 0) {
          p_mean <- (sum(null_mean >= S_mean, na.rm = TRUE) + 1) / (perm + 1)
        }
      }
      
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
    
    if (parallel) {
      library(parallel)
      TWCE_scores_list <- mclapply(unique_paths, compute_pathway, mc.cores = cores)
    } else {
      TWCE_scores_list <- lapply(unique_paths, compute_pathway)
    }
    
    TWCE_scores <- bind_rows(TWCE_scores_list)
    TWCE_scores$p_adjust <- p.adjust(TWCE_scores$p_value, method = p.adjust.method)
    return(TWCE_scores)
  }
}




