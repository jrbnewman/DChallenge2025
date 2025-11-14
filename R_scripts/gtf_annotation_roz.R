# Required packages (Bioconductor + CRAN)
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
pkgs <- c("GenomicRanges","GenomicFeatures","rtracklayer","AnnotationHub","AnnotationDbi",
          "S4Vectors","IRanges","GenomeInfoDb","dplyr","tibble")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p)

library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(AnnotationHub)
library(AnnotationDbi)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(dplyr)
library(tibble)

# Main function
annotateRangesFromGTF <- function(
  query_ranges,           # GRanges or data.frame with chr, start, end, (optional) strand
  gtf_path = NULL,        # local GTF/GFF path (if provided, used)
  download_ensembl = FALSE, # if TRUE, attempt to fetch GRCh38 Ensembl GTF via AnnotationHub
  ensembl_release = NULL,   # optional release number (if NULL use latest available in AH)
  genome = "GRCh38",        # only hg38/GRCh38 supported for download path option
  splice_window = 6L,       # bp considered "near" splice site (default 6)
  upstream_dist = 2000L,    # bp upstream of gene (default 2kb)
  downstream_dist = 2000L,  # bp downstream of gene (default 2kb)
  tss_window = 10000L,      # bp window to report distances to TSS (default 10kb)
  tes_window = 10000L       # bp window to report distances to TES (default 10kb)
) {
  # ---- input sanitization ----
  if (is.data.frame(query_ranges)) {
    # expect columns: seqnames (or chr), start, end, optionally strand
    df <- query_ranges
    colnames(df) <- gsub("^chrom$","seqnames",colnames(df),ignore.case=TRUE)
    if (!any(c("seqnames","chr","chrom") %in% colnames(query_ranges))) {
      stop("If providing data.frame, please include 'seqnames' (or 'chr'/'chrom'), 'start', 'end'.")
    }
    if (!"seqnames" %in% colnames(df)) {
      if ("chr" %in% colnames(df)) df$seqnames <- df$chr else df$seqnames <- df$chrom
    }
    if (!all(c("start","end") %in% colnames(df))) stop("data.frame must contain start and end columns")
    gr <- GRanges(seqnames = df$seqnames, ranges = IRanges(df$start, df$end),
                  strand = if ("strand" %in% colnames(df)) df$strand else "*")
  } else if (inherits(query_ranges, "GRanges")) {
    gr <- query_ranges
  } else {
    stop("query_ranges must be a GRanges or a data.frame")
  }

  # ---- load / create TxDb ----
  if (!is.null(gtf_path)) {
    message("Importing GTF from local path: ", gtf_path)
    txdb <- makeTxDbFromGFF(gtf_path, format = "gtf")
  } else if (download_ensembl) {
    message("Querying AnnotationHub for an Ensembl GRCh38 GTF (may download files)...")
    ah <- AnnotationHub::AnnotationHub()
    # Prefer resources that contain 'Ensembl' + 'GTF' + 'Homo sapiens' + genome GRCh38
    q <- query(ah, c("Homo sapiens", "Ensembl", "GTF", genome))
    if (!is.null(ensembl_release)) {
      q <- q[grep(as.character(ensembl_release), names(q) %||% metadata(q)$resourceName)]
    }
    if (length(q) == 0) stop("Could not find a suitable Ensembl GRCh38 GTF in AnnotationHub.")
    # Choose top hit (most relevant)
    ah_id <- names(q)[1]
    gtf_res <- q[[ah_id]]
    # gtf_res might already be a TxDb or a GRanges; try to coerce to TxDb if needed
    if (inherits(gtf_res, "TxDb")) {
      txdb <- gtf_res
    } else if (inherits(gtf_res, "GRanges")) {
      # Make TxDb from the GRanges (requires writing to a temporary GFF and using makeTxDbFromGFF)
      tmp <- tempfile(fileext = ".gtf")
      rtracklayer::export(gtf_res, con = tmp, format = "gtf")
      txdb <- makeTxDbFromGFF(tmp, format="gtf")
      file.remove(tmp)
    } else {
      # try to import as GTF object and then makeTxDb
      tmp <- tempfile(fileext = ".gtf")
      rtracklayer::export(gtf_res, con = tmp, format = "gtf")
      txdb <- makeTxDbFromGFF(tmp, format="gtf")
      file.remove(tmp)
    }
  } else {
    stop("Provide either gtf_path or set download_ensembl = TRUE")
  }

  # ---- extract annotation objects ----
  genes_gr <- genes(txdb)                    # gene ranges
  transcripts_gr <- transcripts(txdb)        # transcript ranges
  exons_tx <- exonsBy(txdb, by = "tx", use.names = TRUE)  # list of exons by transcript
  exons_gr <- unlist(exons_tx, use.names = TRUE)          # all exons
  cds_tx <- cdsBy(txdb, by = "tx", use.names = TRUE)      # CDS by transcript
  cds_gr <- unlist(cds_tx, use.names = TRUE)
  five_utrs_tx <- fiveUTRsByTranscript(txdb, use.names = TRUE)
  three_utrs_tx <- threeUTRsByTranscript(txdb, use.names = TRUE)
  introns_tx <- intronsByTranscript(txdb, use.names = TRUE)

  # Attempt to read gene symbols / gene_id if available in txdb metadata (optional)
  # We'll try to use 'org' from the GTF if present: use rtracklayer import? But TxDb may not contain gene_symbol metadata.
  # We will attempt to extract gene_id/gene_name from the 'genes_gr' metadata columns if present
  gene_mcols <- mcols(genes_gr)
  has_symbol <- any(c("gene_id","gene_name","symbol") %in% colnames(gene_mcols))
  # For safety, create a mapping table (gene_id -> symbol if available)
  gene_id <- if ("gene_id" %in% colnames(gene_mcols)) gene_mcols$gene_id else names(genes_gr)
  gene_symbol <- if ("gene_name" %in% colnames(gene_mcols)) gene_mcols$gene_name else if ("symbol" %in% colnames(gene_mcols)) gene_mcols$symbol else NA
  gene_map <- DataFrame(gene_id = gene_id, symbol = gene_symbol, row.names = names(genes_gr))

  # ---- helper: transcript <-> gene mapping ----
  tx2gene <- select(txdb, keys=names(transcripts_gr), keytype="TXNAME", columns=c("TXNAME","GENEID"))
  tx2gene_df <- as.data.frame(tx2gene)
  colnames(tx2gene_df) <- c("txname","gene_id")
  # create mapping vector
  tx2gene_map <- setNames(tx2gene_df$gene_id, tx2gene_df$txname)

  # ---- build splice site GRanges (donors and acceptors) per transcript ----
  # We'll generate splice junction single-base positions (the exon boundary positions)
  splice_sites_list <- list()
  for (tx in names(exons_tx)) {
    exs <- exons_tx[[tx]]
    if (length(exs) <= 1) next  # no splice sites in single-exon transcript
    # Ensure exons are ordered by transcript strand
    exs_ord <- if (as.character(unique(strand(exs))) == "-") sort(exs, decreasing = TRUE) else sort(exs)
    # donors and acceptors definition:
    # For '+' strand: donor = end(exon_n) (splice out), acceptor = start(exon_{n+1})
    # For '-' strand: donor/acceptor are reversed relative to coordinates, but we'll mark type considering transcript strand
    # We'll produce a GRanges of single positions at exon boundaries
    starts <- start(exs_ord)
    ends <- end(exs_ord)
    n <- length(exs_ord)
    # donors: boundary after exon i (i from 1..n-1) -> use end(exon i) as donor site coordinate
    donors_pos <- ends[1:(n-1)]
    donors_gr <- GRanges(seqnames = seqnames(exs_ord)[1:(n-1)],
                         ranges = IRanges(donors_pos, donors_pos),
                         strand = strand(exs_ord)[1:(n-1)])
    mcols(donors_gr)$tx <- tx
    mcols(donors_gr)$type <- "donor"
    # acceptors: boundary before exon i (i from 2..n) -> start(exon i) as acceptor coordinate
    acceptors_pos <- starts[2:n]
    acceptors_gr <- GRanges(seqnames = seqnames(exs_ord)[2:n],
                            ranges = IRanges(acceptors_pos, acceptors_pos),
                            strand = strand(exs_ord)[2:n])
    mcols(acceptors_gr)$tx <- tx
    mcols(acceptors_gr)$type <- "acceptor"
    splice_sites_list[[tx]] <- c(donors_gr, acceptors_gr)
  }
  if (length(splice_sites_list) > 0) {
    splice_sites_gr <- unlist(GRangesList(splice_sites_list), use.names = FALSE)
  } else {
    splice_sites_gr <- GRanges()
  }

  # ---- function to convert hits to unique collapsed strings for summary ----
  collapse_unique <- function(x) {
    if (length(x) == 0) return(NA_character_)
    x <- unique(as.character(x))
    paste(x, collapse = ";")
  }

  # ---- compute overlaps and annotations ----
  N <- length(gr)
  # initialize result lists
  out_gene_ids <- vector("list", N)
  out_gene_symbols <- vector("list", N)
  out_tx_ids <- vector("list", N)
  out_feature_types <- vector("list", N)   # exon/intron/5UTR/3UTR/CDS
  out_coding_tx <- vector("list", N)      # transcript IDs where overlaps CDS
  out_splice_site <- logical(N)
  out_near_splice <- vector("list", N)    # transcripts/types within splice_window with distances
  out_upstream_genes <- vector("list", N)
  out_downstream_genes <- vector("list", N)
  out_tss_distances <- vector("list", N)
  out_tes_distances <- vector("list", N)

  # Precompute flanking regions for upstream/downstream for each gene (strand aware)
  # For each gene, define upstream region coordinates based on strand
  genes_df <- as.data.frame(genes_gr)
  g_up <- promoters(genes_gr, upstream = upstream_dist, downstream = 0) # works strand-aware; for '-' it gives region after end
  g_down <- promoters(genes_gr, upstream = 0, downstream = downstream_dist) # downstream relative to TSS, not gene end; safer compute explicit
  # Safer approach to get downstream region relative to gene body:
  # For '+' strand: downstream region = end+1 .. end + downstream_dist
  # For '-' strand: downstream region = start - downstream_dist .. start -1
  g_down2 <- GRanges(seqnames=seqnames(genes_gr),
                     ranges = IRanges(ifelse(as.character(strand(genes_gr))=="+", end(genes_gr)+1, start(genes_gr)-downstream_dist),
                                      ifelse(as.character(strand(genes_gr))=="+", end(genes_gr)+downstream_dist, start(genes_gr)-1)),
                     strand = strand(genes_gr))
  g_up2 <- GRanges(seqnames=seqnames(genes_gr),
                   ranges = IRanges(ifelse(as.character(strand(genes_gr))=="+", start(genes_gr)-upstream_dist, end(genes_gr)+1),
                                    ifelse(as.character(strand(genes_gr))=="+", start(genes_gr)-1, end(genes_gr)+upstream_dist)),
                   strand = strand(genes_gr))
  # Adjust negative coordinates
  start(g_up2)[start(g_up2) < 1] <- 1
  start(g_down2)[start(g_down2) < 1] <- 1

  # compute TSS/TES positions for transcripts
  tss_pos <- resize(transcripts_gr, width = 1, fix = "start") # takes strand into account for start
  # But for '-' transcripts, start is actually transcript end, so ensure TSS is start() when using transcripts()
  # GenomicFeatures::transcripts uses correct start/end w.r.t. strand; start() is TSS for plus, end() for minus? To avoid confusion:
  # define tss per strand explicitly:
  tx_str <- as.character(strand(transcripts_gr))
  tss_coords <- ifelse(tx_str == "+", start(transcripts_gr), end(transcripts_gr))
  tes_coords <- ifelse(tx_str == "+", end(transcripts_gr), start(transcripts_gr))
  tss_gr <- GRanges(seqnames = seqnames(transcripts_gr),
                    ranges = IRanges(tss_coords, tss_coords),
                    strand = strand(transcripts_gr))
  mcols(tss_gr)$tx <- names(transcripts_gr)
  tes_gr <- GRanges(seqnames = seqnames(transcripts_gr),
                    ranges = IRanges(tes_coords, tes_coords),
                    strand = strand(transcripts_gr))
  mcols(tes_gr)$tx <- names(transcripts_gr)

  # iterate queries (vectorized approach where possible)
  # Overlaps: genes, transcripts, exons, introns, UTRs, CDS
  ov_genes <- findOverlaps(gr, genes_gr, ignore.strand = FALSE)
  ov_txs <- findOverlaps(gr, transcripts_gr, ignore.strand = FALSE)
  ov_exons <- findOverlaps(gr, exons_gr, ignore.strand = FALSE)
  ov_cds <- findOverlaps(gr, cds_gr, ignore.strand = FALSE)
  # intronsByTranscript is a list; combine for overlaps:
  introns_all <- unlist(introns_tx, use.names = TRUE)
  ov_introns <- findOverlaps(gr, introns_all, ignore.strand = FALSE)
  # UTR overlaps possible
  utr5_all <- unlist(five_utrs_tx, use.names = TRUE)
  utr3_all <- unlist(three_utrs_tx, use.names = TRUE)
  ov_5utr <- if (length(utr5_all)>0) findOverlaps(gr, utr5_all, ignore.strand = FALSE) else IRanges::IntegerList()
  ov_3utr <- if (length(utr3_all)>0) findOverlaps(gr, utr3_all, ignore.strand = FALSE) else IRanges::IntegerList()

  # splice site distances: compute distance to splice_sites_gr (strand-aware)
  if (length(splice_sites_gr) > 0) {
    # we'll use distanceToNearest to get minimal distances
    d_splice <- distanceToNearest(gr, splice_sites_gr, ignore.strand = FALSE)
    # d_splice is Hits; we can get distances via mcols(d_splice)$distance
    d_splice_df <- as.data.frame(d_splice)
    # distances vector aligned to queries: initialize with NA
    min_splice_dist <- rep(NA_integer_, N)
    splice_hit_idx <- rep(NA_integer_, N)
    if (nrow(d_splice_df) > 0) {
      min_splice_dist[d_splice_df$queryHits] <- d_splice_df$distance
      splice_hit_idx[d_splice_df$queryHits] <- d_splice_df$subjectHits
    }
  } else {
    min_splice_dist <- rep(NA_integer_, N)
    splice_hit_idx <- rep(NA_integer_, N)
  }

  # For TSS/TES distances findNearest (within windows)
  d_tss <- distanceToNearest(gr, tss_gr, ignore.strand = FALSE)
  df_tss <- if (length(d_tss)>0) as.data.frame(d_tss) else data.frame(queryHits=integer(0), subjectHits=integer(0), distance=integer(0))
  d_tes <- distanceToNearest(gr, tes_gr, ignore.strand = FALSE)
  df_tes <- if (length(d_tes)>0) as.data.frame(d_tes) else data.frame(queryHits=integer(0), subjectHits=integer(0), distance=integer(0))

  # Now fill per-query results
  for (i in seq_len(N)) {
    q <- gr[i]

    # genes overlapping
    gi <- subjectHits(ov_genes[queryHits(ov_genes) == i])
    if (length(gi) > 0) {
      gids <- names(genes_gr)[gi]
      out_gene_ids[[i]] <- gids
      # gene symbols if available
      syms <- if (!all(is.na(gene_map$symbol))) {
        sapply(gids, function(g) {
          if (g %in% rownames(gene_map)) gene_map[g,"symbol"] else NA_character_
        })
      } else rep(NA_character_, length(gids))
      out_gene_symbols[[i]] <- syms
    } else {
      out_gene_ids[[i]] <- character(0)
      out_gene_symbols[[i]] <- character(0)
    }

    # transcripts overlapping
    txi <- subjectHits(ov_txs[queryHits(ov_txs) == i])
    if (length(txi) > 0) {
      txs <- names(transcripts_gr)[txi]
      out_tx_ids[[i]] <- txs
    } else {
      out_tx_ids[[i]] <- character(0)
    }

    # features: exon/intron/5'UTR/3'UTR/CDS
    feats <- character(0)
    if (length(subjectHits(ov_exons[queryHits(ov_exons) == i])) > 0) feats <- c(feats, "exon")
    if (length(subjectHits(ov_introns[queryHits(ov_introns) == i])) > 0) feats <- c(feats, "intron")
    if (length(subjectHits(ov_5utr[queryHits(ov_5utr) == i])) > 0) feats <- c(feats, "5UTR")
    if (length(subjectHits(ov_3utr[queryHits(ov_3utr) == i])) > 0) feats <- c(feats, "3UTR")
    if (length(subjectHits(ov_cds[queryHits(ov_cds) == i])) > 0) feats <- c(feats, "CDS")
    out_feature_types[[i]] <- unique(feats)

    # coding transcripts (overlap with CDS)
    cds_hits <- subjectHits(ov_cds[queryHits(ov_cds) == i])
    if (length(cds_hits) > 0) {
      # cds_gr may have names like tx:exon; attempt to map to tx names
      txnames_cds <- names(cds_gr)[cds_hits]
      # if cds_gr names are like "txname" or include TXNAME, use them; otherwise map by overlap to transcripts
      # fallback: find which transcripts have CDS that overlap this query
      txs_with_cds <- unique(sapply(txnames_cds, function(x) {
        if (!is.null(x) && x != "") {
          # sometimes names(cds_gr) are "txname"
          x
        } else NA_character_
      }))
      txs_with_cds <- txs_with_cds[!is.na(txs_with_cds)]
      out_coding_tx[[i]] <- txs_with_cds
    } else {
      out_coding_tx[[i]] <- character(0)
    }

    # splice site exact hit?
    if (!is.na(min_splice_dist[i]) && min_splice_dist[i] == 0) {
      out_splice_site[i] <- TRUE
    } else {
      out_splice_site[i] <- FALSE
    }

    # near splice: collect all splice sites within splice_window distance (and report transcript and type and distance)
    if (length(splice_sites_gr) > 0) {
      # compute distances to all splice sites on same seqnames/strand quickly: use distanceToNearest already provided minimal; but we may want all within splice_window
      # We'll find splice sites on same seq and strand with distance <= splice_window
      candidate_idx <- which(seqnames(splice_sites_gr) == seqnames(q) & (as.character(strand(splice_sites_gr)) %in% c(as.character(strand(q), "*"))))
      if (length(candidate_idx) > 0) {
        dists <- distanceToNearest(q, splice_sites_gr[candidate_idx], ignore.strand = FALSE)
        if (length(dists) > 0) {
          ddf <- as.data.frame(dists)
          if (nrow(ddf) > 0) {
            within_idx <- which(ddf$distance <= splice_window)
            if (length(within_idx) > 0) {
              hits <- ddf[within_idx, , drop=FALSE]
              recs <- lapply(seq_len(nrow(hits)), function(r) {
                sidx <- candidate_idx[hits$subjectHits[r]]
                list(tx = mcols(splice_sites_gr)$tx[sidx],
                     type = mcols(splice_sites_gr)$type[sidx],
                     distance = hits$distance[r])
              })
              out_near_splice[[i]] <- recs
            } else {
              out_near_splice[[i]] <- list()
            }
          } else {
            out_near_splice[[i]] <- list()
          }
        } else {
          out_near_splice[[i]] <- list()
        }
      } else {
        out_near_splice[[i]] <- list()
      }
    } else {
      out_near_splice[[i]] <- list()
    }

    # upstream genes overlapping our upstream region
    up_ov <- findOverlaps(q, g_up2, ignore.strand = FALSE)
    if (length(up_ov) > 0) {
      gidx <- subjectHits(up_ov)
      out_upstream_genes[[i]] <- names(genes_gr)[gidx]
    } else out_upstream_genes[[i]] <- character(0)

    # downstream genes
    down_ov <- findOverlaps(q, g_down2, ignore.strand = FALSE)
    if (length(down_ov) > 0) {
      gidx <- subjectHits(down_ov)
      out_downstream_genes[[i]] <- names(genes_gr)[gidx]
    } else out_downstream_genes[[i]] <- character(0)

    # TSS distances within tss_window
    tss_near <- df_tss[df_tss$queryHits == i & df_tss$distance <= tss_window, , drop=FALSE]
    if (nrow(tss_near) > 0) {
      recs <- lapply(seq_len(nrow(tss_near)), function(r) {
        tx <- mcols(tss_gr)$tx[tss_near$subjectHits[r]]
        list(tx = tx, distance = tss_near$distance[r])
      })
      out_tss_distances[[i]] <- recs
    } else out_tss_distances[[i]] <- list()

    # TES distances within tes_window
    tes_near <- df_tes[df_tes$queryHits == i & df_tes$distance <= tes_window, , drop=FALSE]
    if (nrow(tes_near) > 0) {
      recs <- lapply(seq_len(nrow(tes_near)), function(r) {
        tx <- mcols(tes_gr)$tx[tes_near$subjectHits[r]]
        list(tx = tx, distance = tes_near$distance[r])
      })
      out_tes_distances[[i]] <- recs
    } else out_tes_distances[[i]] <- list()
  }

  # ---- assemble results into a tibble ----
  result <- tibble(
    query_index = seq_len(N),
    seqnames = as.character(seqnames(gr)),
    start = start(gr),
    end = end(gr),
    strand = as.character(strand(gr)),
    gene_ids = out_gene_ids,
    gene_symbols = out_gene_symbols,
    transcript_ids = out_tx_ids,
    feature_overlap = out_feature_types,
    coding_transcripts = out_coding_tx,
    is_splice_site = out_splice_site,
    splice_near = out_near_splice,
    upstream_genes = out_upstream_genes,
    downstream_genes = out_downstream_genes,
    tss_within_window = out_tss_distances,
    tes_within_window = out_tes_distances
  )

  return(result)
}

# ---- Example usage ----
# Example 1: Use a local GTF and annotate sample GRanges
# gtf_file <- "~/data/Homo_sapiens.GRCh38.109.gtf"
# my_gr <- GRanges("chr1", IRanges(c(1000000,1000200), width=1), strand = c("+","-"))
# res <- annotateRangesFromGTF(query_ranges = my_gr, gtf_path = gtf_file)
# View(res)

# Example 2: Download Ensembl GRCh38 and annotate a data.frame of ranges
# df <- data.frame(seqnames = c("chr1"), start = c(11868), end = c(11868), strand = c("+"))
# res2 <- annotateRangesFromGTF(query_ranges = df, download_ensembl = TRUE)
# print(res2)













# annotateRangesFromGTF_with_flatten_and_CDS_frame.R
# Enhanced R function that:
#  - returns either list-column tibble (default) or a flattened table (flatten = TRUE)
#  - reports exon number for overlaps (exon index within transcript)
#  - reports coding information including reading-frame (0/1/2) and predicted amino-acid positions
#
# Dependencies: GenomicRanges, GenomicFeatures, rtracklayer, AnnotationHub (optional), dplyr, tibble, tidyr
# (install Bioconductor packages with BiocManager::install() as needed)

# Example usage (after loading this file):
# res <- annotateRangesFromGTF(query_ranges = my_gr, gtf_path = gtf_file, flatten = TRUE)
# head(res)

annotateRangesFromGTF <- function(
  query_ranges,
  gtf_path = NULL,
  download_ensembl = FALSE,
  ensembl_release = NULL,
  genome = "GRCh38",
  splice_window = 6L,
  upstream_dist = 2000L,
  downstream_dist = 2000L,
  tss_window = 10000L,
  tes_window = 10000L,
  flatten = FALSE,                # If TRUE return a flattened data.frame (one row per annotation hit)
  include_exon_number = TRUE,    # If TRUE compute exon number (if overlapping exon)
  include_coding_frame = TRUE    # If TRUE compute CDS frame and aa positions for CDS overlaps
) {
  # ---- load libs ----
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) stop("Install GenomicRanges")
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) stop("Install GenomicFeatures")
  if (!requireNamespace("rtracklayer", quietly = TRUE)) stop("Install rtracklayer")
  if (!requireNamespace("S4Vectors", quietly = TRUE)) stop("Install S4Vectors")
  if (!requireNamespace("IRanges", quietly = TRUE)) stop("Install IRanges")
  if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) stop("Install GenomeInfoDb")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Install dplyr")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Install tibble")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Install tidyr")

  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
  library(S4Vectors)
  library(IRanges)
  library(GenomeInfoDb)
  library(dplyr)
  library(tibble)
  library(tidyr)

  # ---- input -> GRanges ----
  if (is.data.frame(query_ranges)) {
    df <- query_ranges
    colnames(df) <- gsub("^chrom$","seqnames",colnames(df),ignore.case=TRUE)
    if (!"seqnames" %in% colnames(df)) {
      if ("chr" %in% colnames(df)) df$seqnames <- df$chr else if ("chrom" %in% colnames(df)) df$seqnames <- df$chrom
    }
    if (!all(c("start","end") %in% colnames(df))) stop("data.frame must contain start and end")
    gr <- GRanges(seqnames = df$seqnames, ranges = IRanges(df$start, df$end), strand = if ("strand" %in% colnames(df)) df$strand else "*")
  } else if (inherits(query_ranges, "GRanges")) {
    gr <- query_ranges
  } else stop("query_ranges must be a GRanges or data.frame")

  # ---- load TxDb ----
  if (!is.null(gtf_path)) {
    txdb <- makeTxDbFromGFF(gtf_path, format = "gtf")
  } else if (download_ensembl) {
    ah <- AnnotationHub::AnnotationHub()
    q <- query(ah, c("Homo sapiens", "Ensembl", "GTF", genome))
    if (!is.null(ensembl_release)) q <- q[grep(as.character(ensembl_release), names(q) %||% metadata(q)$resourceName)]
    if (length(q) == 0) stop("No Ensembl GRCh38 GTF found in AnnotationHub")
    res <- q[[names(q)[1]]]
    if (inherits(res, "TxDb")) txdb <- res else {
      tmp <- tempfile(fileext = ".gtf"); rtracklayer::export(res, tmp, format = "gtf"); txdb <- makeTxDbFromGFF(tmp, format = "gtf"); file.remove(tmp)
    }
  } else stop("Provide gtf_path or set download_ensembl = TRUE")

  # ---- extract features ----
  genes_gr <- genes(txdb)
  transcripts_gr <- transcripts(txdb)
  exons_by_tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
  exons_all <- unlist(exons_by_tx, use.names = TRUE)
  cds_by_tx <- cdsBy(txdb, by = "tx", use.names = TRUE)
  cds_all <- unlist(cds_by_tx, use.names = TRUE)
  five_by_tx <- tryCatch(fiveUTRsByTranscript(txdb, use.names = TRUE), error = function(e) GRangesList())
  three_by_tx <- tryCatch(threeUTRsByTranscript(txdb, use.names = TRUE), error = function(e) GRangesList())
  introns_by_tx <- tryCatch(intronsByTranscript(txdb, use.names = TRUE), error = function(e) GRangesList())
  introns_all <- if (length(introns_by_tx) > 0) unlist(introns_by_tx, use.names = TRUE) else GRanges()

  # gene symbol mapping attempt
  gene_meta <- mcols(genes_gr)
  gene_symbol_map <- NULL
  if ("gene_name" %in% colnames(gene_meta)) gene_symbol_map <- setNames(gene_meta$gene_name, names(genes_gr))
  if (is.null(gene_symbol_map) & "symbol" %in% colnames(gene_meta)) gene_symbol_map <- setNames(gene_meta$symbol, names(genes_gr))

  # ---- helper: exon number mapping per transcript ----
  # create a data.frame mapping exon GRanges indices -> exon_number and tx
  exon_map <- NULL
  if (length(exons_by_tx) > 0) {
    exon_map_list <- lapply(names(exons_by_tx), function(tx) {
      exs <- exons_by_tx[[tx]]
      # order by transcript strand
      if (as.character(unique(strand(exs))) == "-") exs <- sort(exs, decreasing = TRUE) else exs <- sort(exs)
      data.frame(tx = tx, exon_idx = seq_along(exs), seqnames = as.character(seqnames(exs)), start = start(exs), end = end(exs), stringsAsFactors = FALSE)
    })
    exon_map <- do.call(rbind, exon_map_list)
  }

  # ---- CDS: compute CDS-start position per transcript for frame calculations ----
  cds_start_by_tx <- list()
  if (length(cds_by_tx) > 0) {
    for (tx in names(cds_by_tx)) {
      cds_parts <- cds_by_tx[[tx]]
      # order CDS fragments according to transcript strand
      if (as.character(unique(strand(cds_parts))) == "-") cds_parts <- sort(cds_parts, decreasing = TRUE) else cds_parts <- sort(cds_parts)
      # CDS start coordinate (genomic) is the first base of CDS in transcript coordinate
      cds_start_by_tx[[tx]] <- start(cds_parts)[1]
      # also store cumulative widths for computing offsets across multi-exon CDS if necessary
      # store GRanges for easier offset computation later
    }
  }

  # ---- compute overlaps: genes, transcripts, exons, introns, UTRs, CDS ----
  ov_genes <- findOverlaps(gr, genes_gr, ignore.strand = FALSE)
  ov_txs <- findOverlaps(gr, transcripts_gr, ignore.strand = FALSE)
  ov_exons <- findOverlaps(gr, exons_all, ignore.strand = FALSE)
  ov_introns <- if (length(introns_all) > 0) findOverlaps(gr, introns_all, ignore.strand = FALSE) else IntegerList()
  ov_5utr <- if (length(unlist(five_by_tx))>0) findOverlaps(gr, unlist(five_by_tx, use.names=TRUE), ignore.strand = FALSE) else IntegerList()
  ov_3utr <- if (length(unlist(three_by_tx))>0) findOverlaps(gr, unlist(three_by_tx, use.names=TRUE), ignore.strand = FALSE) else IntegerList()
  ov_cds <- if (length(cds_all) > 0) findOverlaps(gr, cds_all, ignore.strand = FALSE) else IntegerList()

  # Build splice sites (single-base donors/acceptors) from exon boundaries (like previous implementation)
  splice_list <- list()
  if (length(exons_by_tx) > 0) {
    for (tx in names(exons_by_tx)) {
      exs <- exons_by_tx[[tx]]
      if (length(exs) <= 1) next
      exs_ord <- if (as.character(unique(strand(exs))) == "-") sort(exs, decreasing = TRUE) else sort(exs)
      starts <- start(exs_ord); ends <- end(exs_ord); n <- length(exs_ord)
      donors <- GRanges(seqnames = seqnames(exs_ord)[1:(n-1)], ranges = IRanges(ends[1:(n-1)], ends[1:(n-1)]), strand = strand(exs_ord)[1:(n-1)])
      mcols(donors)$tx <- tx; mcols(donors)$type <- "donor"
      accepts <- GRanges(seqnames = seqnames(exs_ord)[2:n], ranges = IRanges(starts[2:n], starts[2:n]), strand = strand(exs_ord)[2:n])
      mcols(accepts)$tx <- tx; mcols(accepts)$type <- "acceptor"
      splice_list[[tx]] <- c(donors, accepts)
    }
  }
  splice_sites <- if (length(splice_list)>0) unlist(GRangesList(splice_list), use.names = FALSE) else GRanges()

  # distances to splice sites (min), TSS/TES distance
  # prepare transcript TSS/TES
  tx_strand <- as.character(strand(transcripts_gr))
  tss_coords <- ifelse(tx_strand == "+", start(transcripts_gr), end(transcripts_gr))
  tes_coords <- ifelse(tx_strand == "+", end(transcripts_gr), start(transcripts_gr))
  tss_gr <- GRanges(seqnames = seqnames(transcripts_gr), ranges = IRanges(tss_coords, tss_coords), strand = strand(transcripts_gr)); mcols(tss_gr)$tx <- names(transcripts_gr)
  tes_gr <- GRanges(seqnames = seqnames(transcripts_gr), ranges = IRanges(tes_coords, tes_coords), strand = strand(transcripts_gr)); mcols(tes_gr)$tx <- names(transcripts_gr)

  # helper to compute frame and AA positions for a given overlap with CDS
  compute_frame_info <- function(query_idx, cds_idx) {
    # cds_idx indexes into cds_all
    # cds_all has names like "tx.." if use.names=TRUE; but we will find transcript by mapping cds_by_tx
    gr_q <- gr[query_idx]
    cds_range <- cds_all[cds_idx]
    # determine transcript name: cds_all may have names like "tx:cds-range" or be anonymous; try to find parent tx by overlap in cds_by_tx
    tx_candidates <- names(cds_by_tx)[sapply(cds_by_tx, function(x) any(overlapsAny(x, cds_range)))]
    tx <- if (length(tx_candidates)>0) tx_candidates[1] else NA_character_
    if (is.na(tx) | is.null(tx)) return(list(is_coding = TRUE, frame = NA_integer_, aa_start = NA_integer_, aa_end = NA_integer_, tx = NA_character_))

    # compute offset from CDS start in transcript coordinates
    # We'll compute number of coding bases from the beginning of the CDS to the start of overlap
    # Approach: gather ordered CDS parts for this tx, compute cumulative widths, and find how many bases before the overlap start.
    cds_parts <- cds_by_tx[[tx]]
    # ensure order by transcript orientation
    if (as.character(unique(strand(cds_parts))) == "-") cds_parts <- sort(cds_parts, decreasing = TRUE) else cds_parts <- sort(cds_parts)
    # compute cumulative widths
    widths <- width(cds_parts)
    cum_width <- cumsum(widths)
    # find which CDS part the cds_range corresponds to by matching start and end
    # match by any overlap
    part_idx <- which(sapply(cds_parts, function(x) any(overlapsAny(x, cds_range))))
    if (length(part_idx) == 0) part_idx <- 1

    # compute bases before overlap start: sum widths of prior parts + offset within current part
    # compute query overlap relative to cds_parts[part_idx]
    cur_part <- cds_parts[part_idx]
    if (as.character(strand(cur_part)) == "+") {
      offset_in_part <- start(cds_range) - start(cur_part)
    } else {
      # for negative strand, coordinates reverse
      offset_in_part <- end(cur_part) - end(cds_range)
    }
    bases_before <- if (part_idx > 1) sum(widths[1:(part_idx-1)]) + offset_in_part else offset_in_part
    # compute frame: modulo 3
    frame0 <- (bases_before %% 3)
    # compute aa position of the first affected codon
    aa_start <- floor(bases_before / 3) + 1
    # compute length of overlap in coding bases
    overlap_width <- width(intersect(gr_q, cds_range))
    aa_end <- floor((bases_before + overlap_width - 1) / 3) + 1
    return(list(is_coding = TRUE, frame = as.integer(frame0), aa_start = aa_start, aa_end = aa_end, tx = tx))
  }

  # ---- assemble flattened records ----
  records <- list()
  qN <- length(gr)
  for (i in seq_len(qN)) {
    q <- gr[i]
    # find overlapping genes
    ghits <- subjectHits(ov_genes[queryHits(ov_genes) == i])
    gh_names <- if (length(ghits) > 0) names(genes_gr)[ghits] else character(0)
    g_symbols <- if (!is.null(gene_symbol_map) && length(gh_names)>0) sapply(gh_names, function(g) gene_symbol_map[g]) else rep(NA_character_, length(gh_names))

    # transcripts overlapping
    thits <- subjectHits(ov_txs[queryHits(ov_txs) == i])
    tx_names <- if (length(thits) > 0) names(transcripts_gr)[thits] else character(0)

    # if flatten and there are multiple hits, create one row per transcript/gene/feature; otherwise create a single row with list-columns
    # We'll generate per-feature rows: exon, CDS, intron, 5UTR, 3UTR. Also if none, mark intergenic and optionally upstream/downstream

    # helper to add a record
    add_rec <- function(feature, tx = NA_character_, gene = NA_character_, gene_symbol = NA_character_, exon_idx = NA_integer_, is_coding = FALSE, frame = NA_integer_, aa_start = NA_integer_, aa_end = NA_integer_) {
      rec <- list(
        query_index = i,
        query_seqnames = as.character(seqnames(q)),
        query_start = start(q),
        query_end = end(q),
        query_strand = as.character(strand(q)),
        gene_id = gene,
        gene_symbol = gene_symbol,
        tx_id = tx,
        feature = feature,
        exon_number = exon_idx,
        is_coding = is_coding,
        coding_frame = frame,
        aa_start = aa_start,
        aa_end = aa_end
      )
      records[[length(records) + 1]] <<- rec
    }

    any_hit <- FALSE

    # exons
    exon_hits <- subjectHits(ov_exons[queryHits(ov_exons) == i])
    if (length(exon_hits) > 0) {
      any_hit <- TRUE
      for (eh in exon_hits) {
        # find transcript containing this exon by matching coordinates in exon_map
        ex_row <- exon_map[exon_map$seqnames == as.character(seqnames(exons_all[eh])) & exon_map$start == start(exons_all[eh]) & exon_map$end == end(exons_all[eh]), , drop = FALSE]
        if (nrow(ex_row) >= 1) {
          tx <- ex_row$tx[1]; exon_idx <- ex_row$exon_idx[1]
        } else {
          tx <- NA_character_; exon_idx <- NA_integer_
        }
        # map gene via tx -> gene (use select)
        gene <- NA_character_
        if (!is.na(tx)) {
          sel <- tryCatch(GenomicFeatures::select(txdb, keys = tx, keytype = "TXNAME", columns = c("TXNAME","GENEID")), error = function(e) NULL)
          if (!is.null(sel) && nrow(sel)>0) gene <- sel$GENEID[1]
        }
        gsym <- if (!is.na(gene) && !is.null(gene_symbol_map) && gene %in% names(gene_symbol_map)) gene_symbol_map[gene] else NA_character_
        add_rec("exon", tx = tx, gene = gene, gene_symbol = gsym, exon_idx = exon_idx, is_coding = FALSE)
      }
    }

    # CDS
    cds_hits <- subjectHits(ov_cds[queryHits(ov_cds) == i])
    if (length(cds_hits) > 0) {
      any_hit <- TRUE
      for (ch in cds_hits) {
        # attempt compute frame/aa
        frame_info <- if (include_coding_frame) compute_frame_info(i, ch) else list(is_coding = TRUE, frame = NA_integer_, aa_start = NA_integer_, aa_end = NA_integer_, tx = NA_character_)
        # figure gene id
        gene <- NA_character_
        if (!is.null(frame_info$tx) && !is.na(frame_info$tx)) {
          sel <- tryCatch(GenomicFeatures::select(txdb, keys = frame_info$tx, keytype = "TXNAME", columns = c("TXNAME","GENEID")), error = function(e) NULL)
          if (!is.null(sel) && nrow(sel)>0) gene <- sel$GENEID[1]
        }
        gsym <- if (!is.na(gene) && !is.null(gene_symbol_map) && gene %in% names(gene_symbol_map)) gene_symbol_map[gene] else NA_character_
        add_rec("CDS", tx = frame_info$tx, gene = gene, gene_symbol = gsym, exon_idx = NA_integer_, is_coding = TRUE, frame = frame_info$frame, aa_start = frame_info$aa_start, aa_end = frame_info$aa_end)
      }
    }

    # UTRs
    utr5_hits <- if (length(unlist(five_by_tx))>0) subjectHits(ov_5utr[queryHits(ov_5utr) == i]) else integer(0)
    if (length(utr5_hits) > 0) { any_hit <- TRUE; for (uh in utr5_hits) add_rec("5UTR") }
    utr3_hits <- if (length(unlist(three_by_tx))>0) subjectHits(ov_3utr[queryHits(ov_3utr) == i]) else integer(0)
    if (length(utr3_hits) > 0) { any_hit <- TRUE; for (uh in utr3_hits) add_rec("3UTR") }

    # introns
    if (length(introns_all) > 0) {
      intr_hits <- subjectHits(ov_introns[queryHits(ov_introns) == i])
      if (length(intr_hits) > 0) { any_hit <- TRUE; for (ih in intr_hits) add_rec("intron") }
    }

    # if no transcript/gene feature hit: consider upstream/downstream or intergenic
    if (!any_hit) {
      # upstream/downstream
      # compute nearest gene and whether within upstream_dist/downstream_dist
      dist_genes <- distanceToNearest(q, genes_gr, ignore.strand = FALSE)
      if (length(dist_genes) > 0) {
        ddf <- as.data.frame(dist_genes)
        if (ddf$queryHits == i) {
          gidx <- ddf$subjectHits
          gene <- names(genes_gr)[gidx]
          # determine relative position w.r.t. gene strand
          gobj <- genes_gr[gidx]
          # compute if upstream or downstream
          rel <- if (as.character(strand(gobj)) == "+") {
            if (end(q) < start(gobj)) "upstream" else if (start(q) > end(gobj)) "downstream" else "overlapping"
          } else {
            if (start(q) > end(gobj)) "upstream" else if (end(q) < start(gobj)) "downstream" else "overlapping"
          }
          if (rel %in% c("upstream","downstream")) add_rec(rel, gene = gene, gene_symbol = if (!is.null(gene_symbol_map)) gene_symbol_map[gene] else NA_character_)
          else add_rec("intergenic")
        } else add_rec("intergenic")
      } else add_rec("intergenic")
    }
  }

  # convert records into tibble
  if (length(records) == 0) {
    # no annotations found -> return empty tibble matching schema
    res <- tibble(query_index = integer(0), query_seqnames = character(0), query_start = integer(0), query_end = integer(0), query_strand = character(0), gene_id = character(0), gene_symbol = character(0), tx_id = character(0), feature = character(0), exon_number = integer(0), is_coding = logical(0), coding_frame = integer(0), aa_start = integer(0), aa_end = integer(0))
  } else {
    res <- bind_rows(lapply(records, as_tibble))
  }

  # If flatten = FALSE, return list-column summary similar to earlier function
  if (!flatten) {
    # build summary per query
    qN <- length(gr)
    summary_list <- vector("list", qN)
    for (i in seq_len(qN)) {
      sub <- res %>% filter(query_index == i)
      if (nrow(sub) == 0) {
        summary_list[[i]] <- tibble(gene_id = NA_character_, gene_symbol = NA_character_, tx_id = NA_character_, feature = NA_character_)[0,]
      } else {
        summary_list[[i]] <- sub %>% select(gene_id, gene_symbol, tx_id, feature, exon_number, is_coding, coding_frame, aa_start, aa_end)
      }
    }
    out <- tibble(query_index = seq_len(qN), seqnames = as.character(seqnames(gr)), start = start(gr), end = end(gr), strand = as.character(strand(gr)), annotations = summary_list)
    return(out)
  } else {
    return(res)
  }
}

# End of file


