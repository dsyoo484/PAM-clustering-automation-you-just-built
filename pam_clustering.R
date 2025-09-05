#!/usr/bin/env Rscript
# ==============================================================================
# PAM Clustering 자동화 스크립트 (Linux Rscript)
#   - 입력 : ASV/Taxon count table (행=ASV/Taxon, 열=Sample)  ← count 권장
#   - 단계 : 1) Argument  2) Import & Preprocess(TSS)
#           3) 최적 클러스터(CH + Silhouette)
#           4) PAM 결과 테이블/PCoA 저장
#           5) 바플롯 저장(개별 + 통합)
#           6) Excel report
#   - 출력 : PNG(지표/PCoA/Barplot), TSV(지표/결과), XLSX(report)
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(reshape2)
  library(cluster)
  library(clusterSim)
  library(ade4)
  library(ggplot2)
  library(openxlsx)
  # grid는 base에 포함, 함수에서 grid:: 로 사용
})

# --------------------------- 1) Argument importing -----------------------------
# -i / --input     : 입력 테이블 (TSV/CSV, 첫열=feature/ASV, 첫행=sample)
# -o / --outdir    : 출력 폴더
# -c / --cluster_n : 시험할 최대 k (>=2)
# --count          : TRUE이면 column별 합으로 TSS scaling 적용
# --stage          : 1..6 (부분 실행용; 기본 6=모두)
option_list <- list(
  make_option(c("-i","--input"),     type="character", help="Input count table (TSV/CSV)"),
  make_option(c("-o","--outdir"),    type="character", default="result_PAM", help="Output directory"),
  make_option(c("-c","--cluster_n"), type="integer",   default=10, help="Max k to test (>=2)"),
  make_option(c("--count"),          type="logical",   default=TRUE, help="TSS scaling (TRUE/FALSE)"),
  make_option(c("--stage"),          type="integer",   default=6, help="Run stage (1..6)")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input) || !file.exists(opt$input)) stop("❌ 입력 파일을 찾을 수 없습니다: ", opt$input)
if (opt$cluster_n < 2) stop("❌ --cluster_n 은 2 이상이어야 합니다.")
dir.create(file.path(opt$outdir, "plot"), recursive = TRUE, showWarnings = FALSE)
stage_at_least <- function(n) isTRUE(opt$stage >= n)
msg <- function(...) cat(sprintf(...), "\n")

# 팔레트
top10_col   <- c("#496989","#58A399","#A8CD9F","#cf685d","#652475",
                 "#41C9E2","#db74b9","#9b9e9e","#E8C872","#69595a")
cluster_col <- c("red","#7b54bf","blue","#e38b07","cyan","#de04ab",
                 "black","#3da302","yellow","#91422c")

# 표 저장 시 첫 컬럼 이름 지정
save_table <- function(data, col = "Sample") {
  new_df <- data.frame(col = rownames(data))
  new_df <- cbind(new_df, data)
  colnames(new_df)[1] <- col
  rownames(new_df) <- NULL
  new_df
}

# 여러 ggplot을 한 장 PNG로 (ggpubr 없이 grid만 사용)
combine_ggplots_grid <- function(plots, file, ncol, nrow, width, height, dpi = 300) {
  png(filename = file, width = width, height = height, units = "in", res = dpi)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow, ncol)))
  k <- 1
  for (r in seq_len(nrow)) {
    for (c in seq_len(ncol)) {
      if (k <= length(plots)) {
        print(plots[[k]], vp = grid::viewport(layout.pos.row = r, layout.pos.col = c))
        k <- k + 1
      }
    }
  }
  dev.off()
}

# --------------------------- 2) Data importing & preprocessing -----------------
if (stage_at_least(1)) {
  # biom 헤더(# Constructed from biom file) 한 줄 skip
  first  <- readLines(opt$input, n = 1, warn = FALSE)
  skip_n <- ifelse(grepl("^#\\s*Constructed from biom file", first), 1, 0)

  data <- tryCatch(
    read.csv(opt$input, header = TRUE, row.names = 1,
             sep = ifelse(grepl("\\.csv$", opt$input, ignore.case = TRUE), ",", "\t"),
             dec = ".", check.names = FALSE, skip = skip_n),
    error = function(e) stop("입력 파일 파싱 실패: ", conditionMessage(e))
  )

  # taxonomy 마지막 레벨만 남기기 (k__..;..;s__X → X)
  rn <- rownames(data)
  if (any(grepl("__", rn, fixed = TRUE))) rownames(data) <- sub(".*__", "", rn)

  writeLines(sprintf("dim=%d x %d", nrow(data), ncol(data)),
             file.path(opt$outdir, "_read_ok.txt"))
  msg("Stage1 OK: read (%d features × %d samples)", nrow(data), ncol(data))
}

if (stage_at_least(2)) {
  # TSS scaling (각 열 합=1)
  if (isTRUE(opt$count)) {
    data <- apply(data, 2, function(x) if (sum(x) == 0) x else x / sum(x))
  }
  writeLines("TSS scaling applied.", file.path(opt$outdir, "_preprocess_ok.txt"))
  msg("Stage2 OK: preprocessing (TSS=%s)", as.character(opt$count))
}

# --------------------------- 공용 함수 (JSD, PAM) ------------------------------
dist.JSD <- function(inMatrix, pseudocount = 1e-6) {
  KLD <- function(x, y) sum(x * log(x / y))
  JSD <- function(x, y) sqrt(0.5 * KLD(x, (x + y) / 2) + 0.5 * KLD(y, (x + y) / 2))
  inMatrix <- apply(inMatrix, 1:2, function(x) ifelse(x == 0, pseudocount, x))
  p <- ncol(inMatrix)
  m <- matrix(0, p, p)
  for (i in 1:p) {
    xi <- inMatrix[, i]
    for (j in i:p) {
      v <- JSD(xi, inMatrix[, j])
      m[i, j] <- v; m[j, i] <- v
    }
  }
  colnames(m) <- rownames(m) <- colnames(inMatrix)
  as.dist(m)
}
pam.clustering <- function(x, k) pam(as.dist(x), k, diss = TRUE)$clustering

# --------------------------- 3) 최적 클러스터 결정 -----------------------------
if (stage_at_least(3)) {
  data.dist <- dist.JSD(data)
  kmax <- opt$cluster_n

  ch_index         <- data.frame(cluster = 1:kmax, CH = NA_real_)
  silhouette_index <- data.frame(cluster = 1:kmax, Silhouette = NA_real_)
  sil_detail_list  <- list()

  for (k in 2:kmax) {
    cl <- pam.clustering(data.dist, k)

    # CH index
    ch_index$CH[k] <- tryCatch(
      index.G1(t(data), cl, d = data.dist, centrotypes = "medoids"),
      error = function(e) NA_real_
    )

    # Silhouette (존재할 때만 기록)
    sil <- tryCatch(silhouette(cl, data.dist), error = function(e) NULL)
    if (!is.null(sil)) {
      sil_mat   <- as.matrix(sil)
      sil_names <- rownames(sil_mat)
      if (is.null(sil_names) || length(sil_names) != nrow(sil_mat)) sil_names <- colnames(data)
      silhouette_index$Silhouette[k] <- mean(sil_mat[, 3], na.rm = TRUE)
      sil_detail_list[[as.character(k)]] <-
        data.frame(k = k, sample = sil_names, sil_width = sil_mat[, 3], row.names = NULL)
    } else {
      silhouette_index$Silhouette[k] <- NA_real_
    }
  }

  # 저장
  optimal_df <- merge(ch_index, silhouette_index, by = "cluster")
  write.table(optimal_df, file = file.path(opt$outdir, "0_cluster_index.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  if (length(sil_detail_list)) {
    sil_detail <- do.call(rbind, sil_detail_list)
    write.table(sil_detail, file = file.path(opt$outdir, "0_silhouette_detail.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }

  # 시각화
  p_ch <- ggplot(ch_index, aes(x = factor(cluster), y = CH)) +
    geom_bar(stat = "identity", fill = "grey50") + theme_bw() +
    labs(x = "k cluster", y = "CH index")
  ggsave(file.path(opt$outdir, "plot/1.CH_index.png"), p_ch, dpi = 300, width = 5, height = 3)

  p_sil <- ggplot(silhouette_index, aes(x = factor(cluster), y = Silhouette)) +
    geom_bar(stat = "identity", fill = "grey50") + theme_bw() +
    labs(x = "k cluster", y = "Silhouette index")
  ggsave(file.path(opt$outdir, "plot/2.Silhouette.png"), p_sil, dpi = 300, width = 5, height = 3)

  if (exists("sil_detail")) {
    p_sil_box <- ggplot(sil_detail, aes(x = factor(k), y = sil_width)) +
      geom_boxplot(outlier.size = 0.6) + theme_bw() +
      labs(x = "k cluster", y = "Silhouette width")
    ggsave(file.path(opt$outdir, "plot/2b.Silhouette_boxplot.png"),
           p_sil_box, dpi = 300, width = 5.5, height = 3.2)
  }

  # 최적 k (CH 최대, 동률 시 Silhouette 최대)
  k_by_ch <- ch_index$cluster[which(ch_index$CH == max(ch_index$CH, na.rm = TRUE))]
  if (length(k_by_ch) > 1) {
    sub <- silhouette_index[silhouette_index$cluster %in% k_by_ch, ]
    ocluster <- sub$cluster[which.max(sub$Silhouette)]
  } else {
    ocluster <- k_by_ch
  }
  assign("data.dist", data.dist, inherits = TRUE)
  msg("Stage3 OK: indices computed; best k = %d", ocluster)
}

# --------------------------- 4) 결과 테이블 & PCoA ------------------------------
if (stage_at_least(4)) {
  final_cluster <- pam.clustering(data.dist, ocluster)

  # 결과 테이블
  cluster_sample <- data.frame(SampleID = colnames(data),
                               Cluster  = final_cluster,
                               stringsAsFactors = FALSE)
  write.table(cluster_sample, file = file.path(opt$outdir, "1_clustering_results.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # PCoA
  obs.pcoa <- dudi.pco(data.dist, scannf = FALSE, nf = 2)
  pcoa_df  <- data.frame(PC1 = obs.pcoa$li[, 1], PC2 = obs.pcoa$li[, 2],
                         Cluster = factor(final_cluster))
  p_pcoa <- ggplot(pcoa_df, aes(PC1, PC2, color = Cluster)) +
    geom_point(size = 3) + stat_ellipse() + theme_bw() +
    scale_color_manual(values = cluster_col[seq_len(length(unique(pcoa_df$Cluster)))]) +
    labs(title = paste("PCoA (k =", ocluster, ")"))
  ggsave(file.path(opt$outdir, "plot/3.PCoA.png"), p_pcoa, dpi = 300, width = 5.5, height = 4)

  assign("cluster_sample", cluster_sample, inherits = TRUE)
  msg("Stage4 OK: PAM & PCoA")
}

# --------------------------- 5) 바플롯(개별 + 통합) ----------------------------
if (stage_at_least(5)) {
  # 클러스터 평균 → Top10 바플롯
  clusters_index <- cluster_sample$Cluster[match(colnames(data), cluster_sample$SampleID)]
  new_data <- as.data.frame(cbind(Cluster = clusters_index, t(data) * 100))  # %
  cluster_means <- new_data %>% group_by(Cluster) %>%
    summarise(across(everything(), mean), .groups = "drop")

  plots <- list()  # 통합용 컨테이너

  for (k in sort(unique(cluster_means$Cluster))) {
    row_vec <- as.numeric(cluster_means[cluster_means$Cluster == k, -1])
    feats   <- colnames(cluster_means)[-1]
    df_top  <- data.frame(ASV = feats, Abundance = row_vec)
    df_top  <- df_top[order(-df_top$Abundance), , drop = FALSE]
    df_top  <- df_top[seq_len(min(10, nrow(df_top))), , drop = FALSE]
    df_top$ASV <- factor(df_top$ASV, levels = rev(df_top$ASV))

    p_bar <- ggplot(df_top, aes(x = ASV, y = Abundance, fill = ASV)) +
      geom_bar(stat = "identity") + coord_flip() + theme_bw() +
      scale_fill_manual(values = rep(top10_col, length.out = nrow(df_top))) +
      labs(title = paste0("Top10 - Cluster ", k),
           x = "Feature (ASV/Taxon)", y = "Mean Relative Abundance (%)") +
      theme(legend.position = "none")

    # 개별 저장
    ggsave(file.path(opt$outdir, sprintf("plot/Barplot_Cluster%d.png", k)),
           p_bar, dpi = 300, width = 6, height = 4)

    plots[[as.character(k)]] <- p_bar
  }

  # 모든 클러스터 통합 PNG 저장 (grid 기반)
  if (length(plots) > 0) {
    ord <- as.character(sort(as.numeric(names(plots))))
    plots <- plots[ord]
    n <- length(plots)
    ncol <- if (n >= 3) 3 else n
    nrow <- ceiling(n / ncol)
    combine_ggplots_grid(
      plots,
      file   = file.path(opt$outdir, "plot/Barplot_AllClusters.png"),
      ncol   = ncol,
      nrow   = nrow,
      width  = 5.4 * ncol,
      height = 3.6 * nrow,
      dpi    = 300
    )
  }

  msg("Stage5 OK: barplots (individual + combined)")
}

# --------------------------- 6) Excel report -----------------------------------
if (stage_at_least(6)) {
  wb <- createWorkbook()
  addWorksheet(wb, "OptimalCluster")
  writeData(wb, "OptimalCluster",
            read.delim(file.path(opt$outdir, "0_cluster_index.tsv"), sep = "\t"))
  addWorksheet(wb, "ClusteringResult")
  writeData(wb, "ClusteringResult",
            read.delim(file.path(opt$outdir, "1_clustering_results.tsv"), sep = "\t"))
  if (file.exists(file.path(opt$outdir, "0_silhouette_detail.tsv"))) {
    addWorksheet(wb, "SilhouetteDetail")
    writeData(wb, "SilhouetteDetail",
              read.delim(file.path(opt$outdir, "0_silhouette_detail.tsv"), sep = "\t"))
  }
  saveWorkbook(wb, file = file.path(opt$outdir, "PAM_clustering_report.xlsx"), overwrite = TRUE)
  msg("Stage6 OK: excel report")
}

msg("🎉 All done → %s", normalizePath(opt$outdir))

