# Required packages (Bioconductor + CRAN)
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
pkgs <- c("GenomicRanges","GenomicFeatures","rtracklayer","AnnotationHub","AnnotationDbi",
          "S4Vectors","IRanges","GenomeInfoDb","dplyr","tibble","txdbmaker")
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
library(txdbmaker)

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


query_ranges = proxy.ranges
gtf_path = NULL
download_ensembl = TRUE
ensembl_release = NULL
genome = "GRCh38"
splice_window = 6
upstream_dist = 2000
downstream_dist = 2000
tss_window = 10000
tes_window = 10000



  if (is.data.frame(query_ranges)) {
    # expect columns: seqnames (or chr), start, end, optionally strand
    df <- query_ranges
    colnames(df) <- gsub("^chrom$","seqnames",colnames(df),ignore.case=TRUE)
    if (!any(c("seqnames","chr","chrom") %in% colnames(query_ranges))) {
      stop("If providing data.frame, please include 'seqnames' (or 'chr'/'chrom'), 'start', 'end'.")
    }
    if (!"seqnames" %in% colnames(df)) {
      if ("chr" %in% colnames(df)) df$seqnames <- gsub("chr","",df$chr) else df$seqnames <- gsub("chr","",df$chrom)
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
    txdb <- txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf")
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
      txdb <- txdbmaker::makeTxDbFromGFF(tmp, format="gtf")
      file.remove(tmp)
    } else {
      # try to import as GTF object and then makeTxDb
      tmp <- tempfile(fileext = ".gtf")
      rtracklayer::export(gtf_res, con = tmp, format = "gtf")
      txdb <- txdbmaker::makeTxDbFromGFF(tmp, format="gtf")
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
  tx2gene <- select(txdb, keys=transcripts_gr$tx_name, keytype="TXNAME", columns=c("TXNAME","GENEID"))
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
      txs <- transcripts_gr$tx_name[txi]
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
      candidate_idx <- which(as.character(seqnames(splice_sites_gr)) == as.character(seqnames(q)))
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

annot.df <- as.data.frame(result)
annot.df$gene_ids <- sapply(annot.df$gene_ids, paste, collapse = "|")
annot.df$transcript_ids <- sapply(annot.df$transcript_ids, paste, collapse = "|")
annot.df$feature_overlap <- sapply(annot.df$feature_overlap, paste, collapse = "|")
annot.df$coding_transcripts <- sapply(annot.df$coding_transcripts, paste, collapse = "|")
annot.df$upstream_genes <- sapply(annot.df$upstream_genes, paste, collapse = "|")
annot.df$downstream_genes <- sapply(annot.df$downstream_genes, paste, collapse = "|")
annot.df$gene_symbols <- NULL
annot.df$tss_within_window <- sapply(annot.df$tss_within_window, paste, collapse = "|")
annot.df$tes_within_window <- sapply(annot.df$tes_within_window, paste, collapse = "|")
annot.df$splice_near <- sapply(annot.df$splice_near, paste, collapse = "|")



write.csv(annot.df, "/orange/brusko/jrbnewman/full_SNP_annotation_example.csv")

exon.snps <- annot.df["exon" %in% annot.df$feature_overlap,]