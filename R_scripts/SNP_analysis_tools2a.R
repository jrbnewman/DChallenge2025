#' Load or install the most recent versions of SNPloc.Hsapiens.dbSNP and XtraSNPloc.Hsapiens.dbSNP packages
#' 
#' Detects installed SNPloc.Hsapiens and XtraSNPloc.Hsapiens, installs the newest ones from
#' Bioconductor if needed, loads it, and records the version in options(SNPloc.current)
#' and options(XtraSNPloc.current).
#' @param verbose Logical; print progress messages.
#' @return The name(s) of loaded packages (invisibly)
#' @export

verbose=TRUE


## Helper function for splitting a vector or dataframe column by a delimiter
## and then returning the n-th element
## Surprised there isn't a straightforward way like this in base R already

split_string_delim <- function(x, delimiter, n) {
    # Ensure input is character
    x <- as.character(x)

    # Split and extract n-th element
    sapply(strsplit(x, delimiter, fixed=TRUE), function(parts) {
        if(length(parts) >= n) {
	    parts[[n]]
	} else {
	   NA  # return NA if element doesn't exist
	   }
	})
}



# Load or install the most recent dbSNP.VERSION package
load_latest_dbSNP <- function(verbose = TRUE) {
  # Ensure BiocManager is available
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  # List installed package names
  pkgs <- installed.packages()[, "Package"]

  # Match dbSNP.NUMERIC pattern
  dbsnp_pkgs <- grep("^SNPlocs.Hsapiens.dbSNP[0-9]+\\.GRCh38$", pkgs, value = TRUE)  
  xtradbsnp_pkgs <- grep("^XtraSNPlocs.Hsapiens.dbSNP[0-9]+\\.GRCh38$", pkgs, value = TRUE)
  
  # Check for SNPlocs first
  if (length(dbsnp_pkgs) == 0) {
    if (verbose) cat("No local SNPlocs.Hsapiens.dbSNP packages found. Checking Bioconductor...\n")

    # Check Bioconductor for available dbSNP packages
    bioc_pkgs <- BiocManager::available()
    bioc_dbsnp <- grep("^SNPlocs.Hsapiens.dbSNP[0-9]+\\.GRCh38$", bioc_pkgs, value = TRUE)

    if (length(bioc_dbsnp) == 0) {
      stop("No SNPlocs.Hsapiens.dbSNP packages found on Bioconductor.")
    }

    # Determine most recent Bioconductor dbSNP version
    bioc_versions <- as.numeric(sub("SNPlocs.Hsapiens.dbSNP([0-9]+)\\.GRCh38", "\\1", bioc_dbsnp))
    latest_pkg <- bioc_dbsnp[which.max(bioc_versions)]

    if (verbose) cat("Installing latest SNPlocs.Hsapiens.dbSNP package from Bioconductor:", latest_pkg, "\n")
    BiocManager::install(latest_pkg, ask = FALSE, update = FALSE)
  } else {
    # Find the most recent installed version
    versions <- as.numeric(sub("SNPlocs.Hsapiens.dbSNP([0-9]+)\\.GRCh38", "\\1", dbsnp_pkgs))
    latest_pkg <- dbsnp_pkgs[which.max(versions)]
  }

  # Now check for XtraSNPlocs first
  if (length(xtradbsnp_pkgs) == 0) {
    if (verbose) cat("No local XtraSNPlocs.Hsapiens.dbSNP packages found. Checking Bioconductor...\n")

    # Check Bioconductor for available dbSNP packages
    bioc_pkgs.xtra <- BiocManager::available()
    bioc_dbsnp.xtra <- grep("^XtraSNPlocs.Hsapiens.dbSNP[0-9]+\\.GRCh38$", bioc_pkgs.xtra, value = TRUE)

    if (length(bioc_dbsnp.xtra) == 0) {
      stop("No XtraSNPlocs.Hsapiens.dbSNP packages found on Bioconductor.")
    }

    # Determine most recent Bioconductor dbSNP version
    bioc_versions.xtra <- as.numeric(sub("XtraSNPlocs.Hsapiens.dbSNP([0-9]+)\\.GRCh38", "\\1", bioc_dbsnp.xtra))
    latest_pkg.xtra <- bioc_dbsnp.xtra[which.max(bioc_versions.xtra)]

    if (verbose) cat("Installing latest XtraSNPlocs.Hsapiens.dbSNP package from Bioconductor:", latest_pkg.xtra, "\n")
    BiocManager::install(latest_pkg.xtra, ask = FALSE, update = FALSE)
  } else {
    # Find the most recent installed version
    versions.xtra <- as.numeric(sub("XtraSNPlocs.Hsapiens.dbSNP([0-9]+)\\.GRCh38", "\\1", xtradbsnp_pkgs))
    latest_pkg.xtra <- xtradbsnp_pkgs[which.max(versions.xtra)]
  }



  # Load the chosen packages
  if (verbose) cat("Loading package:", latest_pkg, "\n")
  suppressPackageStartupMessages(library(latest_pkg, character.only = TRUE))
  if (verbose) cat("Loading package:", latest_pkg.xtra, "\n")
  suppressPackageStartupMessages(library(latest_pkg.xtra, character.only = TRUE))

  # Store in global options for easy recall
  options(SNPloc.current = latest_pkg)
  options(XtraSNPloc.current = latest_pkg.xtra)

  if (verbose) cat("Current SNPloc.Hsapiens.dbSNP.GRCh38 version set to:", getOption("SNPloc.current"), "\n")
  if (verbose) cat("Current XtraSNPloc.Hsapiens.dbSNP.GRCh38 version set to:", getOption("XtraSNPloc.current"), "\n")

  invisible(latest_pkg)
  invisible(latest_pkg.xtra)
}


# Get the currently active dbSNP package version
get_dbSNP_version <- function() {
  current <- getOption("SNPloc.current", NULL)
  if (is.null(current)) {
    warning("No SNPloc.Hsapiens.dbSNP version currently loaded. Run load_latest_dbSNP() first.")
    return(NULL)
  }
  # Check if still loaded
  if (!(current %in% loadedNamespaces())) {
    warning("dbSNP package not loaded in this session. Reloading...")
    suppressPackageStartupMessages(library(current, character.only = TRUE))
  }
  return(current)
}

get_XtradbSNP_version <- function() {
  current <- getOption("XtraSNPloc.current", NULL)
  if (is.null(current)) {
    warning("No XtraSNPloc.Hsapiens.dbSNP version currently loaded. Run load_latest_dbSNP() first.")
    return(NULL)
  }
  # Check if still loaded
  if (!(current %in% loadedNamespaces())) {
    warning("XtraSNPloc.Hsapiens.dbSNP package not loaded in this session. Reloading...")
    suppressPackageStartupMessages(library(current, character.only = TRUE))
  }
  return(current)
}


# Switch to a specific dbSNP version (e.g., switch_dbSNP_version(155))
switch_dbSNP_version <- function(version, verbose = TRUE) {
  pkg_name <- if (is.numeric(version)) sprintf("SNPlocs.Hsapiens.dbSNP%s.GRCh38", version) else version

  # Check if installed
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop(sprintf("Package '%s' is not installed.", pkg_name))
  }

  # Load package
  if (verbose) cat("Switching to:", pkg_name, "\n")
  suppressPackageStartupMessages(library(pkg_name, character.only = TRUE))

  # Store in options
  options(SNPloc.current = pkg_name)

  if (verbose) cat("Current dbSNP version set to:", getOption("SNPloc.current"), "\n")

  invisible(pkg_name)
}


# Switch to a specific XtradbSNP version (e.g., switch_dbSNP_version(155))
switch_XtradbSNP_version <- function(version, verbose = TRUE) {
  pkg_name <- if (is.numeric(version)) sprintf("XtraSNPlocs.Hsapiens.dbSNP%s.GRCh38", version) else version

  # Check if installed
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop(sprintf("Package '%s' is not installed.", pkg_name))
  }

  # Load package
  if (verbose) cat("Switching to:", pkg_name, "\n")
  suppressPackageStartupMessages(library(pkg_name, character.only = TRUE))

  # Store in options
  options(XtraSNPloc.current = pkg_name)

  if (verbose) cat("Current dbSNP version set to:", getOption("XtraSNPloc.current"), "\n")

  invisible(pkg_name)
}


# Example usage:
# load_latest_dbSNP()          # Detects and loads latest
# get_dbSNP_version()          # Returns currently active version
# switch_dbSNP_version(155)    # Switch to a specific version



#################################################

#snp_file = "/orange/brusko/jrbnewman/software/dev/snp_list2.txt"
snp_file = "/orange/brusko/jrbnewman/software/dev/mahajan.txt"

# Match SNPs to coordinates via rsIDs or match positions to rsIDs
# Input is a text file (i.e. path to a text file), where each
# line is a SNP id in the formats:
#        rs#
#	 chromosome:position
#	 chromosome:position:ref_allele:alt_allele
#	 chromosome:start_position-end_position
#	 chromosome:start_position-end_position:ref_allele:alt:allele
#
# Output is a dataframe with input ID, rsID, chromsome, position/
# start and end positions, alleles (if possible)
# Output can then be used in subsequent functions

match_SNPs_to_coordinates <- function(snp_file) {

	require(BSgenome)
	require(BSgenome.Hsapiens.UCSC.hg38)
	require(biomaRt)
	load_latest_dbSNP()
	# Read SNP file
	message("Reading SNP file...")
	snps <- readLines(snp_file)
	snps <- trimws(snps)
	snps <- snps[nzchar(snps)]  # remove blank lines

	message("Parsing SNP IDs...")
	# Parse input identifiers: rsIDs, coordinates, coordinates with alleles
	parse_snp <- function(x) {
		parts <- unlist(strsplit(x, ":", fixed = TRUE))
		if (length(parts) == 1 && grepl("^rs", x, ignore.case = TRUE)) {
			return(data.frame(input = x, type = "rsid", rsid = x, chr = NA, start = NA, end = NA, ref = NA, alt = NA, stringsAsFactors = FALSE))
		} else if (length(parts) == 2) {
			return(data.frame(input = x, type = "pos", rsid = NA, chr = parts[1], start = unlist(strsplit(parts[2],"-",fixed=TRUE))[1],
			                  end=unlist(strsplit(parts[2],"-", fixed=TRUE))[2], ref = NA, alt = NA, stringsAsFactors = FALSE))
		} else if (length(parts) == 4) {
			return(data.frame(input = x, type = "allelic", rsid = NA, chr = parts[1], start = unlist(strsplit(parts[2],"-",fixed=TRUE))[1],
			                  end=unlist(strsplit(parts[2],"-", fixed=TRUE))[2],
			                  ref = parts[3], alt = parts[4], stringsAsFactors = FALSE))
		} else {
			warning(sprintf("Unrecognized SNP format: %s", x))
			return(data.frame(input = x, type = NA, rsid = NA, chr = NA, start=NA, end = NA, ref = NA, alt = NA, stringsAsFactors = FALSE))
		}
	}
	snp_df <- do.call(rbind, lapply(snps, parse_snp))
	
	# First query rsIDs and positions in loaded dbSNP databases
	# Anything not found will be later passed to BioMart (minimizes Ensembl server load)
	
	# Load SNP databases
	message("Loading SNP databases...")
	snp.db <- get(get_dbSNP_version())
	xtra.db <- get(get_XtradbSNP_version())
	
	# Query dbSNP databases with rsIDs to get coordinates
	message("Querying SNP databases on rsID...")
	rsids <- snp_df[snp_df$type == "rsid",]$input
	rsids.snv <- as.data.frame(snpsById(snp.db, rsids, ifnotfound="drop", genome="hg38"))
	rsids.xtra <- as.data.frame(snpsById(xtra.db, rsids, ifnotfound="drop"))

	
	# Query dbSNP databases with chr:position coordinates to get rsIDs
	## Set up query table for ease of transferring metadata over to a GRanges object
	message("Querying SNP databases on coordinates...")
	positions <- snp_df[snp_df$type != "rsid",]
	positions$seqnames <- gsub("chr","",positions$chr)
	
	## dbSNP and dbSNP Xtra have slightly different chromosome forms so account for that here
	positions$range <- ifelse(is.na(positions$end) == TRUE, paste0(positions$seqnames,":",positions$start,"-",positions$start),
	                                                        paste0(positions$seqnames,":",positions$start,"-",positions$end))
	positions$range.xtra <- ifelse(is.na(positions$end) == TRUE, paste0("ch",positions$seqnames,":",positions$start,"-",positions$start),
	                                                        paste0("ch",positions$seqnames,":",positions$start,"-",positions$end))
	range.snv <- GRanges(positions$range)
	mcols(range.snv)$inputID <- positions$input
	mcols(range.snv)$q.chr <- positions$chr
	mcols(range.snv)$q.start <- positions$start
	mcols(range.snv)$q.end <- positions$end
	range.xtra <- GRanges(positions$range.xtra)
	mcols(range.xtra)$inputID <- positions$input
	mcols(range.xtra)$q.chr <- positions$chr
	mcols(range.xtra)$q.start <- positions$start
	mcols(range.xtra)$q.end <- positions$end
	
	## Get SNPs within overlapping ranges
	pos.snv <- snpsByOverlaps(snp.db, range.snv)
	pos.xtra <- snpsByOverlaps(xtra.db, range.xtra)
	
	## relate hits back to input ranges
	pos.snv.overlaps <- findOverlaps(GRanges(positions$range), pos.snv)
	pos.xtra.overlaps <- findOverlaps(GRanges(positions$range.xtra), pos.xtra)
	pos.snv.overlaps.df <- as.data.frame(GRanges(seqnames = seqnames(pos.snv[subjectHits(pos.snv.overlaps)]),
	                               ranges = ranges(pos.snv[subjectHits(pos.snv.overlaps)]),
				       strand = strand(pos.snv[subjectHits(pos.snv.overlaps)]),
				       mcols = cbind(mcols(pos.snv[subjectHits(pos.snv.overlaps)]),
				                     mcols(range.snv[queryHits(pos.snv.overlaps)]))))
	pos.xtra.overlaps.df <- as.data.frame(GRanges(seqnames = seqnames(pos.xtra[subjectHits(pos.xtra.overlaps)]),
	                               ranges = ranges(pos.xtra[subjectHits(pos.xtra.overlaps)]),
				       strand = strand(pos.xtra[subjectHits(pos.xtra.overlaps)]),
				       mcols = cbind(mcols(pos.xtra[subjectHits(pos.xtra.overlaps)]),
				                     mcols(range.xtra[queryHits(pos.xtra.overlaps)]))))
						     
	## filter positional matches to match the input
	message("Filter positional matches on expected ID")
	pos.snv.overlaps.df$mcols.q.end <- ifelse(is.na(pos.snv.overlaps.df$mcols.q.end) == TRUE,
	                                        pos.snv.overlaps.df$mcols.q.start,
						pos.snv.overlaps.df$mcols.q.end)
	pos.snv.list <- unique(pos.snv.overlaps.df$mcols.inputID)
	pos.snv.overlaps.df.filtered <- pos.snv.overlaps.df[(pos.snv.overlaps.df$start == pos.snv.overlaps.df$mcols.q.start) &
	                                                    (pos.snv.overlaps.df$end == pos.snv.overlaps.df$mcols.q.end),]
	pos.snv.list.olap <- unique(pos.snv.overlaps.df.filtered$mcols.inputID)

	# Get the list of SNV positions without a match
 	pos.snv.missing <- setdiff(pos.snv.list, pos.snv.list.olap)
	pos.xtra.overlaps.df$mcols.q.end <- ifelse(is.na(pos.xtra.overlaps.df$mcols.q.end) == TRUE,
	                                        pos.xtra.overlaps.df$mcols.q.start,
						pos.xtra.overlaps.df$mcols.q.end)
	pos.xtra.list <- unique(pos.xtra.overlaps.df$mcols.inputID)
	pos.xtra.overlaps.df.filtered <- pos.xtra.overlaps.df[(pos.xtra.overlaps.df$start == pos.xtra.overlaps.df$mcols.q.start) &
	                                                    (pos.xtra.overlaps.df$end == pos.xtra.overlaps.df$mcols.q.end),]
	pos.xtra.list.olap <- unique(pos.xtra.overlaps.df.filtered$mcols.inputID)
	
	## Get the positional inputs without perfect coordinate matches, and match on either start OR end
	## This is in case the recorded coordinates are slightly off, or we're extracting overlapping variants
	## that have yet to be merged
	message("Filtrering positional matches on either start or end coordinate...")
	pos.xtra.missing <- setdiff(pos.xtra.list, pos.xtra.list.olap)
	pos.xtra.overlaps.df.missing <- pos.xtra.overlaps.df[pos.xtra.overlaps.df$mcols.inputID %in% pos.xtra.missing,]
	pos.xtra.missing.df <- pos.xtra.overlaps.df.missing[(pos.xtra.overlaps.df.missing$mcols.q.start == pos.xtra.overlaps.df.missing$start) |
	                                                    (pos.xtra.overlaps.df.missing$mcols.q.end == pos.xtra.overlaps.df.missing$end),]
	
	## Get the final list of positional inputs without a perfect or imperfect match
	pos.xtra.missing2 <- setdiff(pos.xtra.list,c(pos.xtra.missing,pos.xtra.list.olap))
	## Merge with list of unmatched SNV positional inputs
	pos.any.missing <- setdiff(pos.snv.missing, c(pos.xtra.missing,pos.xtra.list.olap))
	
	## Report which positional inputs have no match
	if (length(pos.any.missing) > 0) { message("The following inputs have no matching SNP IDs:", paste(pos.any.missing, collaspe=" ")) }
	
	# Combine all outputs and merge with original SNP dataframe
	message("Combining positional hits")
	pos.xtra.overlaps.df.filtered <- pos.xtra.overlaps.df.filtered[,c("seqnames","start","end","mcols.RefSNP_id","mcols.inputID")]
	pos.xtra.missing.df <- pos.xtra.missing.df[,c("seqnames","start","end","mcols.RefSNP_id","mcols.inputID")]
	pos.snv.overlaps.df.filtered <- pos.snv.overlaps.df.filtered[,c("seqnames","start","end","mcols.RefSNP_id","mcols.inputID")]
	

        if (nrow(rsids.snv) > 0) {
	    rsids.snv$seqnames <- as.character(rsids.snv$seqnames)
  	    rsids.snv$start <- rsids.snv$pos
	    rsids.snv$end <- rsids.snv$pos
	    rsids.snv$inputID <- rsids.snv$RefSNP_id
	    rsids.snv <- rsids.snv[,c("seqnames","start","end","RefSNP_id","inputID","ref_allele","alt_alleles")]
	    names(rsids.snv) <- c("match.chr","match.start","match.end","match.RefSNP_id","input","ref_allele","alt_alleles")
	    } else {
	    rsids.snv <- data.frame(match.chr = character(),
		                                           match.start = integer(),
							   match.end = integer(),
							   match.RefSNP_id = character(),
							   input = character(),
							   ref_allele = character(),
							   alt_alleles = character())
	    }
	if (nrow(rsids.xtra) > 0) {
		rsids.xtra$seqnames <- as.character(rsids.xtra$seqnames)
		rsids.xtra$inputID <- rsids.xtra$RefSNP_id
		rsids.xtra <- rsids.xtra[,c("seqnames","start","end","RefSNP_id","inputID")]
		rsids.xtra$ref_allele <- NA
		rsids.xtra$alt_alleles <- NA
		names(rsids.xtra) <-  c("match.chr","match.start","match.end","match.RefSNP_id","input","ref_allele","alt_alleles")
		}  else {
	   	rsids.xtra <-  data.frame(match.chr = character(),
		                                           match.start = integer(),
							   match.end = integer(),
							   match.RefSNP_id = character(),
							   input = character(),
							   ref_allele = character(),
							   alt_alleles = character())
	    }
	if (nrow(pos.xtra.overlaps.df.filtered) > 0) {
	   pos.xtra.overlaps.df.filtered$seqnames <- as.character(pos.xtra.overlaps.df.filtered$seqnames)
	   pos.xtra.overlaps.df.filtered$ref_allele <- NA
	   pos.xtra.overlaps.df.filtered$alt_alleles <- NA
	   names(pos.xtra.overlaps.df.filtered) <-   c("match.chr","match.start","match.end","match.RefSNP_id","input","ref_allele","alt_alleles")
	   } else {
	   pos.xtra.overlaps.df.filtered   <- data.frame(match.chr = character(),
		                                           match.start = integer(),
							   match.end = integer(),
							   match.RefSNP_id = character(),
							   input = character(),
							   ref_allele = character(),
							   alt_alleles = character())
	   }
	if (nrow(pos.xtra.missing.df) > 0) {
		pos.xtra.missing.df$seqnames <- as.character(pos.xtra.missing.df$seqnames)
		pos.xtra.missing.df$ref_allele <- NA
		pos.xtra.missing.df$alt_alleles <- NA
		names(pos.xtra.missing.df) <-   c("match.chr","match.start","match.end","match.RefSNP_id","input","ref_allele","alt_alleles")
		} else {
		 pos.xtra.missing.df <- data.frame(match.chr = character(),
		                                           match.start = integer(),
							   match.end = integer(),
							   match.RefSNP_id = character(),
							   input = character(),
							   ref_allele = character(),
							   alt_alleles = character())
	        }
	if (nrow(pos.snv.overlaps.df.filtered) > 0) {
	   pos.snv.overlaps.df.filtered$seqnames <- as.character(pos.snv.overlaps.df.filtered$seqnames)
   	   pos.snv.overlaps.df.filtered$ref_allele <- NA
	   pos.snv.overlaps.df.filtered$alt_alleles <- NA
	   names(pos.snv.overlaps.df.filtered) <-  c("match.chr","match.start","match.end","match.RefSNP_id","input","ref_allele","alt_alleles")
	   } else {
	       	pos.snv.overlaps.df.filtered <- data.frame(match.chr = character(),
		                                           match.start = integer(),
							   match.end = integer(),
							   match.RefSNP_id = character(),
							   input = character(),
							   ref_allele = character(),
							   alt_alleles = character())
	   }
	   
		
	
	## And combine
	all.matched.snps <- rbind(rsids.snv,
	                          rsids.xtra,
				  pos.xtra.overlaps.df.filtered,
				  pos.xtra.missing.df,
				  pos.snv.overlaps.df.filtered)

	# Merge with original dataframe
	merged.snp.coord <- merge(snp_df, all.matched.snps, by="input", all.x=TRUE)
	# Fix up chromosome names
	merged.snp.coord$match.chr <- paste0("chr",gsub("ch","",
	                 gsub("Ch","",
	                 gsub("Chr","",
	                 gsub("chr","",merged.snp.coord$match.chr)))))

	# Add alleles from positional inputs if possible/necessary
	merged.snp.coord$ref_allele <- ifelse(is.na(merged.snp.coord$ref_allele)==TRUE,
				              ifelse(is.na(merged.snp.coord$ref) == TRUE, NA, merged.snp.coord$ref),
					      merged.snp.coord$ref_allele)
					      
	merged.snp.coord$alt_alleles <- ifelse(is.na(merged.snp.coord$alt_alleles)==TRUE,
				              ifelse(is.na(merged.snp.coord$alt) == TRUE, NA, merged.snp.coord$alt),
					      merged.snp.coord$alt_alleles)

        # Now take the matched rsIDs from the merged output and query THOSE for ref and alt alleles
	new.rsids <- unique(merged.snp.coord[is.na(merged.snp.coord$match.RefSNP_id)==FALSE,]$match.RefSNP_id)
	new.rsids.snv <- as.data.frame(snpsById(snp.db, new.rsids, ifnotfound="drop", genome="hg38"))
	new.rsids.snv <- new.rsids.snv[,c("RefSNP_id","ref_allele","alt_alleles")]
	names(new.rsids.snv) <- c("match.RefSNP_id","ref_allele2","alt_alleles2")
	merged.snp.coord <- merge(merged.snp.coord, new.rsids.snv, by="match.RefSNP_id", all.x=TRUE)

	merged.snp.coord$ref_allele <- ifelse(is.na(merged.snp.coord$ref_allele) == TRUE,
	                                      merged.snp.coord$ref_allele2,
					      merged.snp.coord$ref_allele)
	merged.snp.coord$alt_alleles <- ifelse(is.na(merged.snp.coord$alt_alleles) == TRUE,
	                                      merged.snp.coord$alt_alleles2,
					      merged.snp.coord$alt_alleles)
	merged.snp.coord$ref_allele2 <- NULL
	merged.snp.coord$alt_alleles2 <- NULL
	
	# Report inputs with no rsIDs					      					      
        no.rsid <- merged.snp.coord[is.na(merged.snp.coord$match.RefSNP_id) == TRUE,]
	if (dim(no.rsid)[1] > 0 ) { message("The following inputs have no associated rsIDs: ", paste(unique(no.rsid$input), collapse=" "))
	                             } else { message("All input SNPs and/or coordinates have annotated rsIDs!") }


        # For missing entries, query BioMart
	## connect to biomart
	message("Fetching annotations for remaining SNPs from Biomart...")
	mart <- biomaRt::useEnsembl(biomart = "snp", dataset = "hsapiens_snp")

	rs_results <- data.frame()
	if (nrow(no.rsid) > 0 ) {
	    rs_results <- biomaRt::getBM(attributes = c("refsnp_id","chr_name","chrom_start","chrom_end","allele", "minor_allele","validated"),
	                                 filters = "snp_filter",
					 values = no.rsid$input,
					 mart=mart)
	    rs_results$match_level <- "rsid"
	    }
	# Drop scaffold SNPs
	allowed.chr.list <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16",
	                    "17","18","19","20","21","22","X","Y","MT","chr1","chr2","chr3","chr4",
			    "chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13",
			    "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22",
			    "chrX","chrY","chrM")
	
	rs_results <- rs_results[rs_results$chr_name %in% allowed.chr.list,]
	
	# Format and add back to merged.snp.coord
	if (nrow(rs_results) > 0) {
 	    rs_results <- rs_results[,c("chr_name","chrom_start","chrom_end","refsnp_id","allele")]
	    rs_results$input <- rs_results$refsnp_id
	    rs_results$ref_allele <- split_string_delim(rs_results$allele, "/", 1)
	    rs_results$alt_alleles <- split_string_delim(rs_results$allele, "/", 2)
	    rs_results$allele <- NULL
	    colnames(rs_results) <- c("match.chr","match.start","match.end","match.RefSNP_id","input","ref_allele","alt_alleles")
	    no.rsid2 <- no.rsid[,c("input","type","chr","rsid", "start","end", "ref","alt")]
	    no.rsid2 <- merge(x=no.rsid2, y=rs_results, by="input", all.x=TRUE)
	    no.rsid2 <- no.rsid2[,c("match.RefSNP_id","input","type","rsid","chr", "start","end","ref","alt",
	                        "match.chr","match.start","match.end","ref_allele","alt_alleles")]
				
	    merged.snp.coord2 <- merged.snp.coord[is.na(merged.snp.coord$match.RefSNP_id) == FALSE,]
	    merged.snp.coord2 <- rbind(merged.snp.coord2, no.rsid2)
	    } else {
	    merged.snp.coord2 <- merged.snp.coord # Nothing in Biomart, just pass the previous results
	    }
	    
	    # Report unmatched inputs
	    missing.all <- setdiff(unique(merged.snp.coord2$input), unique(snp_df$input))

	if (length(missing.all) > 0) { message("The following inputs have no matches in dbSNP:", paste(missing.all, collapse=" "))
	                             } else { message("All input SNPs and/or coordinates annotated!") }
				     
        # Return SNP-to-coordinate index
        return(merged.snp.coord2)

}

snp_index <- match_SNPs_to_coordinates(snp_file)

# topLD LD proxies
## Function to see what topMed populations there are

show_topmed_populations <- function() {
			       message(
			           paste0(c("Available TOPMed populations:",
					           "AFR - African ancestry panel",
						   "EAS - East Asian ancestry panel",
						   "EUR - European ancestry panel",
						   "SAS - South Asian ancestry panel")
						 , sep = "\n"))
				}

show_topmed_populations()



## Function to use topLD
## Ensure you have downloaded the topld_api from https://github.com/linnabrown/topld_api
## and either copied the topld_api binary to wherever your binaries are located (e.g. /usr/local/bin/ or /path/to/your/conda_environment/bin)
## or added the path to the topld_api binary in your $PATH environmental variable (e.g. export PATH=$PATH:/path/to/topld_api)


find_proxies_topmed <- function(snp_list = snp_list,
		                populations = c("AFR","EAS","EUR","SAS"),
				min_r2 = 0.8,
				min_maf = 0.01,
				topld_bin_path = Sys.which("topld_api")) {

			require(dplyr)
		       #require(data.table)
		       try(suppressWarnings(Sys.chmod(topld_bin_path, mode = "0755")), silent=TRUE)
		       
		       # set up final list objects to store data from each population
		       proxies.all <- list()
		       annot.all <- list()

		       # Iterate over selected populations
		       for (p in populations) {
		       	   ## Set up temporary list objects for storing SNP proxies and annotations
			   ## for the current population
			   proxies.list <- list()
			   annot.list <- list()
			      		       
		       	   # Iterate over list of SNPs
			   # use tempfile() to create temporary input and output files
			   message("Running TopLD...")
			   snp_list <- snp_list[is.na(snp_list) == FALSE]
			   for (s in snp_list) {
		               if (is.na(s) == TRUE) { next } else {
			       # Set up temporary I/O
			       inTmp <- tempfile(pattern = "topld_in", fileext=".txt")
			       ldTmp <- tempfile(pattern = "topld_outLD", fileext=".txt")
			       annotTmp <- tempfile(pattern = "topld_outAnnot", fileext=".txt")
			       writeLines(s, inTmp)
	      
			       # Run topLD
			       args <- c("-pop", p,
			                 "-thres", as.character(min_r2),
					 "-maf", as.character(min_maf),
					 "-inFile", inTmp,
					 "-outputLD", ldTmp,
					 "-outputInfo", annotTmp)
			       
			       message(paste0("[TOPLD] Fetching proxy SNPs for ", s, " in the ", p, " population..."))
			         stdout_stderr <- system2(topld_bin_path, args = args, stdout = TRUE, stderr = TRUE)

			       # Read in temporary output files and store in list objects
			       if(!file.exists(ldTmp) || file.size(ldTmp) == 0) {
			           message(paste0("No fetchable proxy SNPs for ", s, " in the ", p, " population. Variant may be missing from reference panel"))
				   next
		 	       }
			       s2 <- gsub(":","_",s)
			       read.delim(ldTmp, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE) %>%
			           dplyr::mutate(across(everything(), as.character)) -> proxies.list[[s2]]
			       read.delim(annotTmp, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE) %>%
			           dplyr::mutate(across(everything(), as.character)) -> annot.list[[s2]]
			       }
			       }
			       proxies.all[[p]] <- dplyr::bind_rows(proxies.list, .id = "input_snp")
			       annot.all[[p]] <- dplyr::bind_rows(annot.list, .id = "input_snp")
			  }
			  # report which input SNPs in which populations have no proxies
			  for (p in populations) {
			      proxies.all[[p]]$input_snp <- gsub("_",":",proxies.all[[p]]$input_snp)
			      no.proxies <- setdiff(snp_list, unique(proxies.all[[p]]$input_snp))
    		       	      if (length(no.proxies) > 0) {
			           message(paste0("The following inputs had no proxies in the ", p, " population matching your criteria or were not found in the panel selected:"))
				   message(paste0(c(no.proxies), sep="\n"))
				   }
				   }

			  proxies.all.merge <- dplyr::bind_rows(proxies.all, .id="population")
			  annot.all.merge <- dplyr::bind_rows(annot.all, .id="population")
			  proxy2input <- proxies.all.merge[,c("population","input_snp","rsID2")]
			  input2input <- data.frame(population = rep(populations, each = length(snp_list)),
			  	      	 	    input_snp = rep(snp_list, rep=length(populations)),
						    rsID2 = rep(snp_list, rep = length(populations)))
			  proxy.index <- rbind(proxy2input, input2input)
			  names(proxy.index) <- c("population","input_snp","proxy_snp")
			  proxy.index$source <- "TOPMed"
			  ld.source <- "topLD"
			  topld_out <- list(proxy_index = proxy.index, proxy_info = proxies.all.merge, input_info = annot.all.merge, ld_source = ld.source)
			  message("Done!")
			  return(topld_out)
		       }


#test.list <- unique(snp_index$match.RefSNP_id)[1:10]
test.list <- unique(c(snp_index$match.RefSNP_id, snp_index$input, snp_index$rsid))


topld.test <- find_proxies_topmed(snp_list = test.list2,
                                  populations = "EUR",
             		    	  min_r2 = 0.8,
				  min_maf = 0.05)
				  
# LDlink LD proxies

## Function to print 1000 Genomes populations
## This is just a wrapper for LDlinkR::list_pop()

show_1000G_populations <-
show_LDlink_populations <- function() {
            require("LDlinkR")

	    LDlinkR::list_pop()
	    }

show_1000G_populations()

show_LDlink_populations()


# Function to take a list of variants and find proxies with LDlink

find_proxies_LDlink <- function(snp_list = snp_list,
                                populations = LDlinkR::list_pop()$pop_code, # List of populations. By default it will run EVERY population, so be cafeful
				r2d = c("r2","d"),			   # r2 or D' as LD statistic
				min_r2 = 0,				   # Minimum r2. 0 = keep all (default)
				min_d = 0,				   # Minimum D'. 0 = keep all (default)
				token = NULL,				   # LDlink token. Use your own.
				max_distance = 1000000,			   # Minimum distance from input SNP (1000000 default)
				genome_build = "grch38") {		   # Genome. GRCh38 is default
                           require("LDlinkR")
			   require("dplyr")
			   message("Fetching proxies via LDLink")
			   proxy.ldlink <- list()
			   for (p in populations) {
			       message(paste0("[LDLink] Fetching proxies for population: ",p))
			       proxy.list <- list()
			       snp_list <- snp_list[is.na(snp_list) == FALSE]
 			       for (i in snp_list) {
			           message(paste0("[LDLink] Fetching proxies for ",i," in the ", p," population..."))
			       	   proxy.list[[i]] <- tryCatch(LDlinkR::LDproxy(snp = i, pop = p, r2d = r2d,
                               	   token = token,
                               	   genome_build=genome_build), error = function(e) { NULL })
			       }
			       proxy.ldlink[[p]] <- dplyr::bind_rows(proxy.list, .id="input_snp")
			       }
			   ldlink.merge <- dplyr::bind_rows(proxy.ldlink, .id="population")

			   message(paste0("Filtering proxies (r2 >= ",
			                  as.character(min_r2),
					  ", D-prime >= ",
					  as.character(min_d),
					  ") ..."))
			   ldlink.merge <- ldlink.merge[(ldlink.merge$R2 >= min_r2) &
			                                (ldlink.merge$Dprime >= min_d) &
							(abs(ldlink.merge$Distance) <= abs(max_distance)),]
			   for (p in names(proxy.ldlink)) {
			       no.proxies <- setdiff(snp_list, unique(proxy.ldlink[[p]]$input_snp))
			       if (length(no.proxies) > 0) {
			           message(paste0("The following inputs had no proxies in the ", p, " population matching your criteria or were not found in the panel selected:"))
				   message(paste0(c(no.proxies), sep="\n"))
				   }
				   }
			   # remove blank/NA rows and replace empty proxy IDs with coord ID
			   ldlink.merge <- ldlink.merge[!is.na(ldlink.merge$input_snp) == TRUE,]
			   ldlink.merge$RS_Number <- ifelse(ldlink.merge$RS_Number == ".", paste0(ldlink.merge$Coord,":",
												  gsub("\\(","",split_string_delim(ldlink.merge$Allele,"/",1)), ":",
												  gsub("\\)","",split_string_delim(ldlink.merge$Allele,"/",2))),
											   ldlink.merge$RS_Number)
												  
			   proxy2input <- ldlink.merge[c("population","input_snp","RS_Number")]
			   input2input <- data.frame(population = rep(populations, each = length(snp_list)),
			  	      	 	    input_snp = rep(snp_list, rep=length(populations)),
						    RS_Number = rep(snp_list, rep = length(populations)))
			   proxy.index <- rbind(proxy2input, input2input)
			   proxy.index <- unique(proxy.index)
 			   names(proxy.index) <- c("population","input_snp","proxy_snp")
			   proxy.index$source <- "1000G"
			   ld.source <- "LDLink"
			   message("Done!")
			   ldlink.list <- list(proxy_index = proxy.index, proxy_info = ldlink.merge, ld_source = ld.source)
			   return(ldlink.list)
			   }

test.list2 <- test.list[c(695,1:357)]


r2.test <- find_proxies_LDlink(snp_list = test.list2,
	            population = "EUR",
		    r2d = "r2",
		    min_r2 = 0.8,
		    token = "62510a7ba674",
		    genome_build="grch38")
		    

# Merge proxies and create an input-to-proxy index with inputs also as proxies
# Need to add a method for carrying forward coordinates and alleles too...
# Add in distance to index and r2 stats too

merge_snp_proxies <- function(x, populations=NULL) {
merged_index <- data.frame(input_snp = c(),
                           proxy_snp = c(),
			   distance_to_input_snp = c(),
			   r2_input_snp = c(),
			   Dprime_input_snp = c(),
			   population = c())
				    
         if (length(populations) == 0) { message("Combining all proxies...")
	   } else { message(paste0(c("Combining proxies from populations:", populations), sep=" "))   }
	    
         for (i in 1:length(x)) {
	# Process topLD and LDLink outputs differently
	     if (x[[i]]$ld_source == "topLD") {
	         if (length(populations) == 0) {
	           # if no populations are selected then combine all
	             roz1 <- x[[i]]$proxy_index[,c("input_snp","proxy_snp","population")]
		     roz2 <- x[[i]]$proxy_info[,c("input_snp","POS1(gb38)","POS2(gb38)","R2","Dprime","rsID2","population")]
		   } else {
  	             roz1 <- x[[i]]$proxy_index[x[[i]]$proxy_index$population %in% populations,][,c("input_snp","proxy_snp","population")]
	             roz2 <- x[[i]]$proxy_info[x[[i]]$proxy_index$population %in% populations,][,c("input_snp","POS1(gb38)","POS2(gb38)","R2","Dprime","rsID2","population")]
		   }
	       	   roz2$distance <- as.numeric(roz2$`POS2(gb38)`) - as.numeric(roz2$`POS1(gb38)`)
		   roz2 <- roz2[,c("input_snp","rsID2","distance","R2","Dprime","population")]
		   names(roz2) <- c("input_snp","proxy_snp","distance_to_input_snp","r2_input_snp","Dprime_input_snp","population")
		   roz <- merge(roz1, roz2, by=c("input_snp","proxy_snp","population"), all.x=TRUE)
		   roz$distance_to_input_snp <- ifelse(is.na(roz$distance_to_input_snp) == TRUE, 0, roz$distance_to_input_snp)
		   roz$r2_input_snp <- ifelse(is.na(roz$r2_input_snp) == TRUE, 1, roz$r2_input_snp)
		   roz$Dprime_input_snp <- ifelse(is.na(roz$Dprime_input_snp) == TRUE, 1, roz$Dprime_input_snp)
	       } else if (x[[i]]$ld_source == "LDLink") {
	         if (length(populations) == 0) {
	           # if no populations are selected then combine all
	             roz1 <- x[[i]]$proxy_index[,c("input_snp","proxy_snp","population")]
		     roz2 <- x[[i]]$proxy_info[,c("input_snp","RS_Number","Distance","R2","Dprime","population")]
		   } else {
		     roz1 <- x[[i]]$proxy_index[x[[i]]$proxy_index$population %in% populations,][,c("input_snp","proxy_snp","population")]
		     roz2 <- x[[i]]$proxy_info[x[[i]]$proxy_index$population %in% populations,][,c("input_snp","RS_Number","Distance","R2","Dprime","population")]
		   }
		   names(roz2) <- c("input_snp","proxy_snp","distance_to_input_snp","r2_input_snp","Dprime_input_snp","population")
		   roz <- merge(roz1, roz2, by=c("input_snp","proxy_snp","population"), all.x=TRUE)
		   roz$distance_to_input_snp <- ifelse(is.na(roz$distance_to_input_snp) == TRUE, 0, roz$distance_to_input_snp)
		   roz$r2_input_snp <- ifelse(is.na(roz$r2_input_snp) == TRUE, 1, roz$r2_input_snp)
		   roz$Dprime_input_snp <- ifelse(is.na(roz$Dprime_input_snp) == TRUE, 1, roz$Dprime_input_snp)
	       }
               merged_index <- rbind(merged_index, roz)
	    	      }
	  # For each input-proxy pair, REGARDLESS OF POPULATION OR LD SOURCE, the best available R2/Dprime stats are used
	    	         merged_index2 <- merged_index %>%
	                        dplyr::group_by(input_snp, proxy_snp) %>%
				arrange(desc(r2_input_snp)) %>%
				slice_head(n = 1)
		merged_index2$population <- NULL
          merged_index2 <- as.data.frame(merged_index2)
          message("Proxies combined! Merged proxy-to-input list is ready for further analysis.")
	 return(merged_index2)
	 }
	 
	       
#merged.index2 <- merge_snp_proxies(list(topld.test, r2.test))	      
merged.index <- merge_snp_proxies(list(topld.test, r2.test), populations="EUR")
#merged.index3 <- merge_snp_proxies(list(topld.test, r2.test, d.test), populations=c("EUR","SAS","AFR"))

##### Create an LD network to group proxies into LD regions

make_LD_regions <- function(x) {
    require(igraph)
    # Create graph from edges
    g <- igraph::graph_from_data_frame(x, directed = FALSE, vertices = unique(c(x$index_snp, x$proxy_snp)))
    # Get connected components
    comp <- components(g)
    # output region to SNP index
    groups.tall <- data.frame(snp_id = names(comp$membership), LD_group_id = comp$membership)

    return(groups.tall)
}


ld_network <- make_LD_regions(merged.index)

############################

source("/orange/brusko/jrbnewman/software/dev/snptools_mahajan_auto.R")


#### Proxy optimizer

vcf.path <- "/orange/brusko/UFDIchip/pipeline_2.0/TOPMed_imputation/QC87_89_90_91_92/chr_"
#chr.list <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
chr.list <- c(1,10,11,3,12,17)
vcf.list <- paste0(vcf.path, chr.list, "/chr", chr.list, ".dose.vcf.gz")

# Create a proxy summary tables to aid in proxy selection

summarize_proxies <- function(proxy_index = proxy_index,
                                      vcfs = c(),
				      min_r2 = 0.9,
				      min_d = 0.9,
				      max_distance = 100000,
				      min_imputed_r2 = 1,
				      min_avg_cs = 1) {
				      
    # packages
    library(data.table)
    require(R.utils)

    # load imputed vcfs and prep info index
    message("Loading VCF(s) and fetching imputation statistics...")
    vcf.all <- list()
    for (v in vcfs) {
        message(paste0("Loading VCF at ", v))
        vcf.in <- fread(v, sep = "\t", header = TRUE, skip = "#CHROM", select=c(1,2,3,4,5,8)) # we only need the first 5 columns and the info column
	vcf.in$INFO <- gsub("TYPED;IMPUTED", "TYPED_IMPUTED", vcf.in$INFO)
	vcf.in$status <- split_string_delim(vcf.in$INFO, ";", 1)
	vcf.in$impute_avg_cs <- gsub("AVG_CS=", "", split_string_delim(vcf.in$INFO, ";", 4))
	vcf.in$impute_r2 <- gsub("R2=", "", split_string_delim(vcf.in$INFO, ";", 5))
	vcf.in <- vcf.in[,c("#CHROM","POS","ID","REF","ALT","status","impute_avg_cs","impute_r2")]
	names(vcf.in) <- c("chr","pos","snp_id","ref_allele","alt_allele","status","impute_avg_cs","impute_r2")
	vcf.all[[v]] <- vcf.in
	rm(vcf.in)
	gc(verbose=FALSE)
     }	
     vcf.all <- dplyr::bind_rows(vcf.all, .id="vcf.chrom")
     message("Replacing missing rsIDs with coordinate ID [chr:pos:ref:alt]...")
     vcf.all$snp_id <- ifelse(vcf.all$snp_id == ".",
                              paste0(vcf.all$chr,":",
			             vcf.all$pos,":",
			             vcf.all$ref_allele,":",
				     vcf.all$alt_allele), vcf.all$snp_id)
     gc(verbose=FALSE)
     
     # (1) Identify which input SNPs are present in user VCF and which are missing
     # NOTE: this is solely on the ID column in your input VCF
     # It is expecting an rsID, but if it is missing, then we've changed it to a
     # chr:pos:ref:alt format
     message("Identifying input SNPs present in imputed VCF(s)...")
     inputs.in.vcf <- proxy_index$input_snp
     inputs.in.vcf <- intersect(inputs.in.vcf, vcf.all$snp_id)
     
     # (2) Report input SNPs without proxies meeting filitering criteria
     # (a) First, merge with proxy index and take the common hits
     # We merge on rsID. I might provide functionality for a coordinate-based backup selection later
     # but not now
     message("Identifying proxy SNPs in imputed VCF(s)...")
     vcf2proxy <- merge(x=proxy_index, y=vcf.all, by.x="proxy_snp", by.y="snp_id")

     # (b) Filter proxies based on distance, min LD statistics, min imputation statistics
     message("Filtering proxy SNPs based on LD, imputation statistics, and distance to input SNP ...")
     vcf2proxy.filtered <- vcf2proxy[(vcf2proxy$distance_to_input_snp <= max_distance) &
                                     (vcf2proxy$impute_avg_cs >= min_avg_cs) &
				     (vcf2proxy$impute_r2 >= min_imputed_r2) &
				     (vcf2proxy$r2_input_snp >= min_r2) &
				     (vcf2proxy$Dprime_input_snp >= min_d),]

     # (c) Now get list of input SNPs in filtered vcf2proxy dataset
     message("Identifying input and proxy SNPs in filtered imputed VCF SNP list")
     inputs.in.filtered <- unique(vcf2proxy.filtered$input_snp)
     proxies.in.filtered <- unique(vcf2proxy.filtered$proxy_snp)


     # Get first best proxy that is not itself by grouping
     # proxy SNPs by input SNP and arrange in order:
     # 	     	   (1) r2 to input SNP, descending
     #		   (2) Dprime to input SNP, descending (redundant with r2)
     #		   (3) Average imputation confidence score, descending
     #		   (4) Average imputation r2, descending 
     #		   (5) Distance to input SNP, ascending
     #		   Then take the first one on the list
     #
     #		   Basically we prioritize LD and imputation quality over distance
     #		   LD has distance indirectly baked-in (as SNPs are more likely to be correlated
     #		   the closer they are together), but we don't want to let chromosomal distance
     #		   define this because LD can work over a larger distance and there might be a
     #		   better variant a bit further away.
     #		   Sorting by Dprime is somewhat redundant, as it will typically have a higher
     #		   value than r2, but it's in here just in case
     #		   Add later functionality: user-defined order of sort
     message("Identifying optimal proxy SNPs  ...")
     vcf2proxy.filtered$abs_distance <- abs(vcf2proxy.filtered$distance_to_input_snp)
     
     vcf2proxy.top <- vcf2proxy.filtered[vcf2proxy.filtered$input_snp != vcf2proxy.filtered$proxy_snp,] %>%
                      group_by(input_snp) %>%
		      arrange(desc(r2_input_snp),
		              desc(Dprime_input_snp),
			      desc(impute_avg_cs),
			      desc(impute_r2),
			      abs_distance) %>%
		      slice_head(n=1) %>%
		      as.data.frame()

     # Annotate proxy index
     message("Annotating proxy index ...")
     proxy_index$input_snp_in_vcf <- ifelse(proxy_index$input_snp %in%  inputs.in.vcf, TRUE, FALSE)
     proxy_index$input_snp_pass_filter <- ifelse(proxy_index$input_snp %in%  inputs.in.filtered, TRUE, FALSE)
     proxy_index$proxy_snp_pass_filter <- ifelse(proxy_index$input_snp %in%  proxies.in.filtered, TRUE, FALSE)
     
     proxy_index$proxy_status <- ifelse(proxy_index$input_snp == proxy_index$proxy_snp, "Same as input",
                                        ifelse(proxy_index$proxy_snp %in%  proxies.in.filtered,
					       ifelse(proxy_index$proxy_snp %in% vcf2proxy.top$proxy_snp,
					              "Best proxy",
						      "Other proxy"),
					       "Proxy excluded"))
					
     
     # Add imputation stats to table as well so decisions can be tracked
     impute.stats <- vcf2proxy[,c("proxy_snp","input_snp","ref_allele","alt_allele",
                                  "status","impute_avg_cs","impute_r2")]
				  
     proxy_index2 <- merge(proxy_index, impute.stats, by=c("input_snp","proxy_snp"), all.x=TRUE)
     names(proxy_index2) <- c("input_snp","proxy_snp","distance_to_input_snp","r2_input_snp","Dprime_input_snp",
			      "input_snp_in_vcf","input_snp_pass_filter","proxy_snp_pass_filter","proxy_status",
			      "vcf_ref_allele","vcf_alt_allele","imputation_status",
			      "imputation_avg_cs","imputation_r2")
			      
     message("Done summarizing proxy selection!")
     return(proxy_index2)
     }




proxy_summary <- summarize_proxies(proxy_index = merged.index,
                                      vcfs = vcf.list,
				      min_r2 = 0.5,
				      min_d = 0.5,
				      max_distance = 1000000,
				      min_imputed_r2 = 0,
				      min_avg_cs = 0)
				      



as.matrix(proxy_summary[proxy_summary$input_snp_in_vcf == FALSE,])

write.csv(proxy_summary, "/orange/brusko/jrbnewman/Mahajan_missing_snps_proxies.csv")



####
# Function to output the input SNPs and the chosen proxies (if any)
# this can then be used as input for updating GRS weight files and such

get_final_proxies <- function(proxy_summary = proxy_summary,
                              snps_to_update = c(),
			      all = FALSE)
     # By default, if there are no SNPs selected to update then update all
     # input SNPs that don't pass VCF filter otherwise, why would you be
     # running this function?

     if (all) {
     	# if all input SNPs are to be replaced then run this
	inputs <- as.data.frame(unique(proxy_summary$input_snp))
	input2proxy <- proxy_summary[proxy_summary$proxy_status == "Best proxy", c("input_snp","proxy_snp")]
	names(inputs) <- "input_snp"
	inputs2 <- merge(inputs, input2proxy, by="input_snp", all.x=TRUE)
	inputs2$proxy_snp <- ifelse(is.na(inputs2$proxy_snp) == TRUE, inputs2$input_snp, inputs2$proxy_snp)
	inputs2$note <- ifelse(inputs2$input_snp == inputs2$proxy_snp,
	                       "No proxies found in VCF meeting criteria. Default to input SNP",
			       "Best proxy SNP in VCF")
	final_proxy_index <- inputs2
	
     } else {
         # otherwise decide whether to do all that don't pass or just a specific set to update
         if (length(snps_to_update) > 0) {
       	     to.keep <- proxy_summary[proxy_summary$input_snp %in% snps_to_update,]

	

     } else {
	 # Update all inputs not passing VCF filter
     }

      




#############################

# Annotate SNP: genes, transcripts, exon/intron/utr, up/downstream, splice site, intergenic, etc.
# Intended input is your proxy SNPs but you could feed any vector of SNP IDs/coordinates
# Will use COORDINATES not rsID




## Annotate coordinates from GTF



## Annotate coordinates from BioMart


## Annotate rsIDs from BioMart










####################################


# Annotate proxy lists: alleles (if don't already have), genes, transcriptsion features

######################################

# Filter proxies by feature: exon, utr, intron, promoter, distance from gene
# If multiple per region (input snp), then output all
# Output regions with no proxies that match criteria

####################################

# LD network generation and plots

################################

# Functional enrichment (all annotated genes, genic only)

# Motif enrichment - by region/proxy set, all proxies
# Just make this based on a set of input SNPs and/or coordinates, then call as needed
# Input SNP list, define sequence flank (default +/- 30 nt), motif enrichment
# For a given motif hit, show an MSA plot of the sequences contributing towards that hit

# Motif discovery? Ref vs alt allele

# Motif GSEA?



















######################################










# FUTURE FUNCTION: Import SNP list (as rsID and/or coordinate)
# parse and create a full data table of:
# input ID, chr, pos, ref, alt, rsID
# We'll be using this table as the basis for all future analysis


load_latest_dbSNP() 

# Load Packages
	
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

# something like library(paste("XtraSNPlocs.Hsapiens.dbSNP",dbsnpXtraver,".GRCh38"))
snps <- readLines(snp_file)
snps <- trimws(snps)
snps <- snps[nzchar(snps)]  # remove blank lines

# ---- Parse identifiers ----
parse_snp <- function(x) {
parts <- unlist(strsplit(x, ":", fixed = TRUE))
if (length(parts) == 1 && grepl("^rs", x, ignore.case = TRUE)) {
return(data.frame(input = x, type = "rsid", rsid = x, chr = NA, start = NA, end = NA, ref = NA, alt = NA, stringsAsFactors = FALSE))
} else if (length(parts) == 2) {
return(data.frame(input = x, type = "pos", rsid = NA, chr = parts[1], start = unlist(strsplit(parts[2],"-",fixed=TRUE))[1],
                  end=unlist(strsplit(parts[2],"-", fixed=TRUE))[2], ref = NA, alt = NA, stringsAsFactors = FALSE))
} else if (length(parts) == 4) {
return(data.frame(input = x, type = "allelic", rsid = NA, chr = parts[1], start = unlist(strsplit(parts[2],"-",fixed=TRUE))[1],
                  end=unlist(strsplit(parts[2],"-", fixed=TRUE))[2],
                  ref = parts[3], alt = parts[4], stringsAsFactors = FALSE))
} else {
warning(sprintf("Unrecognized SNP format: %s", x))
return(data.frame(input = x, type = NA, rsid = NA, chr = NA, start=NA, end = NA, ref = NA, alt = NA, stringsAsFactors = FALSE))
}
}

snp_df <- do.call(rbind, lapply(snps, parse_snp))

# Load SNP databases
## since the SNPloc packages have an eponymously named function, call the output of get_dbSNP_version()
## as a function name using get()
snp.db <- get(get_dbSNP_version())
xtra.db <- get(get_XtradbSNP_version())



# If rsid, get coordinates and any alternative rsID
rsids <- snp_df[snp_df$type == "rsid",]$input

# First try SNVs only, then try other variant types
rsids.snv <- as.data.frame(snpsById(snp.db, rsids, ifnotfound="warning", genome="hg38"))
rsids.xtra <- as.data.frame(snpsById(xtra.db, rsids, ifnotfound="warning"))

# If allelic or positional, get rsIDs
positions <- snp_df[snp_df$type != "rsid",]

## parse positions into a vector
positions$seqnames <- gsub("chr","",positions$chr)

positions$range <- ifelse(is.na(positions$end) == TRUE, paste0(positions$seqnames,":",positions$start,"-",positions$start),
                                                        paste0(positions$seqnames,":",positions$start,"-",positions$end))
positions$range.xtra <- ifelse(is.na(positions$end) == TRUE, paste0("ch",positions$seqnames,":",positions$start,"-",positions$start),
                                                        paste0("ch",positions$seqnames,":",positions$start,"-",positions$end))

range.snv <- GRanges(positions$range)
mcols(range.snv)$inputID <- positions$input
mcols(range.snv)$q.chr <- positions$chr
mcols(range.snv)$q.start <- positions$start
mcols(range.snv)$q.end <- positions$end

range.xtra <- GRanges(positions$range.xtra)
mcols(range.xtra)$inputID <- positions$input
mcols(range.xtra)$q.chr <- positions$chr
mcols(range.xtra)$q.start <- positions$start
mcols(range.xtra)$q.end <- positions$end

pos.snv <- snpsByOverlaps(snp.db, range.snv)
pos.xtra <- snpsByOverlaps(xtra.db, range.xtra)

## relate hits back to input ranges
pos.snv.overlaps <- findOverlaps(GRanges(positions$range), pos.snv)
pos.xtra.overlaps <- findOverlaps(GRanges(positions$range.xtra), pos.xtra)

pos.snv.overlaps.df <- as.data.frame(GRanges(seqnames = seqnames(pos.snv[subjectHits(pos.snv.overlaps)]),
                               ranges = ranges(pos.snv[subjectHits(pos.snv.overlaps)]),
			       strand = strand(pos.snv[subjectHits(pos.snv.overlaps)]),
			       mcols = cbind(mcols(pos.snv[subjectHits(pos.snv.overlaps)]),
			                     mcols(range.snv[queryHits(pos.snv.overlaps)]))))

pos.xtra.overlaps.df <- as.data.frame(GRanges(seqnames = seqnames(pos.xtra[subjectHits(pos.xtra.overlaps)]),
                               ranges = ranges(pos.xtra[subjectHits(pos.xtra.overlaps)]),
			       strand = strand(pos.xtra[subjectHits(pos.xtra.overlaps)]),
			       mcols = cbind(mcols(pos.xtra[subjectHits(pos.xtra.overlaps)]),
			                     mcols(range.xtra[queryHits(pos.xtra.overlaps)]))))


					     
# filter positional matches to match the input
pos.snv.overlaps.df$mcols.q.end <- ifelse(is.na(pos.snv.overlaps.df$mcols.q.end) == TRUE,
                                        pos.snv.overlaps.df$mcols.q.start,
					pos.snv.overlaps.df$mcols.q.end)
pos.snv.list <- unique(pos.snv.overlaps.df$mcols.inputID)
pos.snv.overlaps.df.filtered <- pos.snv.overlaps.df[(pos.snv.overlaps.df$start == pos.snv.overlaps.df$mcols.q.start) &
                                                    (pos.snv.overlaps.df$end == pos.snv.overlaps.df$mcols.q.end),]
pos.snv.list.olap <- unique(pos.snv.overlaps.df.filtered$mcols.inputID)
pos.snv.missing <- setdiff(pos.snv.list, pos.snv.list.olap)



pos.xtra.overlaps.df$mcols.q.end <- ifelse(is.na(pos.xtra.overlaps.df$mcols.q.end) == TRUE,
                                        pos.xtra.overlaps.df$mcols.q.start,
					pos.xtra.overlaps.df$mcols.q.end)
pos.xtra.list <- unique(pos.xtra.overlaps.df$mcols.inputID)
pos.xtra.overlaps.df.filtered <- pos.xtra.overlaps.df[(pos.xtra.overlaps.df$start == pos.xtra.overlaps.df$mcols.q.start) &
                                                    (pos.xtra.overlaps.df$end == pos.xtra.overlaps.df$mcols.q.end),]
pos.xtra.list.olap <- unique(pos.xtra.overlaps.df.filtered$mcols.inputID)
pos.xtra.missing <- setdiff(pos.xtra.list, pos.xtra.list.olap)
pos.xtra.overlaps.df.missing <- pos.xtra.overlaps.df[pos.xtra.overlaps.df$mcols.inputID %in% pos.xtra.missing,]
pos.xtra.missing.df <- pos.xtra.overlaps.df.missing[(pos.xtra.overlaps.df.missing$mcols.q.start == pos.xtra.overlaps.df.missing$start) |
                                                    (pos.xtra.overlaps.df.missing$mcols.q.end == pos.xtra.overlaps.df.missing$end),]


pos.xtra.missing2 <- setdiff(pos.xtra.list,c(pos.xtra.missing,pos.xtra.list.olap))
pos.any.missing <- setdiff(pos.snv.missing, c(pos.xtra.missing,pos.xtra.list.olap))

if (length(pos.any.missing) > 0) { message("The following inputs have no matching SNP IDs:", paste(pos.any.missing, collaspe=" ")) }

# Combine all outputs and merge with original SNP dataframe

pos.xtra.overlaps.df.filtered <- pos.xtra.overlaps.df.filtered[,c("seqnames","start","end","mcols.RefSNP_id","mcols.inputID")]
pos.xtra.missing.df <- pos.xtra.missing.df[,c("seqnames","start","end","mcols.RefSNP_id","mcols.inputID")]
pos.snv.overlaps.df.filtered <- pos.snv.overlaps.df.filtered[,c("seqnames","start","end","mcols.RefSNP_id","mcols.inputID")]

pos.xtra.overlaps.df.filtered$ref_allele <- NA
pos.xtra.missing.df$ref_allele <- NA
pos.snv.overlaps.df.filtered$ref_allele <- NA
pos.xtra.overlaps.df.filtered$alt_alleles <- NA
pos.xtra.missing.df$alt_alleles <- NA
pos.snv.overlaps.df.filtered$alt_alleles <- NA


rsids.snv$start <- rsids.snv$pos
rsids.snv$end <- rsids.snv$pos
rsids.snv$inputID <- rsids.snv$RefSNP_id
rsids.xtra$inputID <- rsids.xtra$RefSNP_id
rsids.snv <- rsids.snv[,c("seqnames","start","end","RefSNP_id","inputID","ref_allele","alt_alleles")]
rsids.xtra <- rsids.xtra[,c("seqnames","start","end","RefSNP_id","inputID")]
rsids.xtra$ref_allele <- NA
rsids.xtra$alt_alleles <- NA

names(pos.xtra.overlaps.df.filtered) <-
    names(pos.xtra.missing.df) <-
    names(pos.snv.overlaps.df.filtered) <-
    names(rsids.snv) <-
    names(rsids.xtra) <-
    c("match.chr","match.start","match.end","match.RefSNP_id","input","ref_allele","alt_alleles")

all.matched.snps <- rbind(rsids.snv,
                          rsids.xtra,
			  pos.xtra.overlaps.df.filtered,
			  pos.xtra.missing.df,
			  pos.snv.overlaps.df.filtered)

merged.snp.coord <- merge(snp_df, all.matched.snps, by="input")
merged.snp.coord$match.chr <- paste0("chr",gsub("ch","",
                 gsub("Ch","",
                 gsub("Chr","",
                 gsub("chr","",merged.snp.coord$match.chr)))))

merged.snp.coord$ref_allele <- ifelse(is.na(merged.snp.coord$ref_allele)==TRUE,
			              ifelse(is.na(merged.snp.coord$ref) == TRUE, NA, merged.snp.coord$ref),
				      merged.snp.coord$ref_allele)
				      
merged.snp.coord$alt_alleles <- ifelse(is.na(merged.snp.coord$alt_alleles)==TRUE,
			              ifelse(is.na(merged.snp.coord$alt) == TRUE, NA, merged.snp.coord$alt),
				      merged.snp.coord$alt_alleles)
				      

missing.all <- setdiff(unique(merged.snp.coord$input), unique(snp_df$input))

if (length(missing.all) > 0) { message("The following inputs have no matching SNP IDs:", paste(missing.all, collaspe=" "))
                             } else { message("All input SNPs and/or coordinates annotated!")


#### I have my FIRST output table, from which I can probably do everything else






## reformat chromosome column as below

# Now merge everything back onto the original list of SNPs, so we have input ID, chr, pos, 
# Also maybe re-add ref and alt alleles if possible (re-query SNP database? Might need to be SNVs only)
# Then it should be ready for: LD proxies, gene/transcript annotation, motif enrichment (SNV only)
# if ref and alt alleles provided in input, then check against rsID alleles
# if no alleles reported (e.g. structural variants) then leave blank and skip motif enrichment/discovery
# If positional has alleles, then might need to also filter on alleles (or report all and have the disclaimer
# that some variants that don't have rsIDs are queried on chromosomal range/position and so may
# include additional variants not anticipated. It's an unfortunate limitation of the approach.

### Keep some variables and make others
rsids.snv$start <- rsids.snv$pos
rsids.snv$end <- rsids.snv$pos
pos.snv$start <- pos.snv$pos
pos.snv$end <- pos.snv$pos

rsids.snv$chr <- paste0("chr",gsub("ch","",
gsub("Ch","",
                 gsub("Chr","",
                 gsub("chr","",rsids.snv$seqnames)))))
rsids.xtra$chr <- paste0("chr",gsub("ch","",
                 gsub("Ch","",
                 gsub("Chr","",
                 gsub("chr","",rsids.xtra$seqnames)))))
pos.snv$chr <- paste0("chr",gsub("ch","",
                 gsub("Ch","",
                 gsub("Chr","",
                 gsub("chr","",pos.snv$seqnames)))))
pos.xtra$chr <- paste0("chr",gsub("ch","",
                 gsub("Ch","",
                 gsub("Chr","",
                 gsub("chr","",pos.xtra$seqnames)))))



rsids.snv <- rsids.snv[,c("chr","start","end","RefSNP_id")]
rsids.xtra <- rsids.xtra[,c("chr","start","end","RefSNP_id")]
pos.snv <- pos.snv[,c("chr","start","end","RefSNP_id")]
pos.xtra <- pos.xtra[,c("chr","start","end","RefSNP_id")]

