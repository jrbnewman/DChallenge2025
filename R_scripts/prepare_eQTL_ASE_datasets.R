### Prepare TIGER and DICE eQTL and ASE summary statistics into
### queryable datasets


# Helper function: streamlined syntax for splitting a vector of strings
# on a delimiter and taking the n-th element
# literal syntax is "take string X, split it on given delimiter, take element N"

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





# Install dependencies
  for (pkg in c("vcfR")) {
      if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos="https://cloud.r-project.org/")
  }

# Load libraries
  library(vcfR)
  library(dplyr)
  library(biomaRt)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)

# DICE and TIGER use hg19 coordinates, so we're going to lift those
# over to hg38. First, we set up the chain file
  message("Downloading hg19-to-hg38 chain file and loading into R...")
  chain_url <- "https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz"
  chain_file_local <- file.path(tempdir(), "hg19ToHg38.over.chain.gz")
  download.file(chain_url, destfile = chain_file_local)

  # Uncompress the chain file
  system(paste("gzip -d -f", chain_file_local))
  chain_file_uncompressed <- gsub("\\.gz$", "", chain_file_local)
  
  # import chain
  chain <- rtracklayer::import.chain(chain_file_uncompressed)
  message("Chain file loaded!")


# DICE data -- this is eQTL data from mutliple immune cell types
# These are formatted as VCFs
  cell.types <- c("B_CELL_NAIVE", "CD4_NAIVE", "CD4_STIM", "CD8_NAIVE",
  	     	  "CD8_STIM", "M2", "MONOCYTES", "NK", "TFH", "TH17",
		  "TH1", "TH2", "THSTAR", "TREG_MEM", "TREG_NAIVE")
  in.path = "/orange/brusko/jrbnewman/software/dev/datasets"
  
  message("Processing DICE eQTL data...")
# Set up a list object for appending data
  dice.sumstats <- list()
 
  for (c in cell.types) {
      message(paste0("Processing eQTLs for cell type: ", c))
      sumstats <- vcfR::read.vcfR(paste0(in.path,"/",c,".vcf"))
      # Get variant metadata
      fixstats <- as.data.frame(vcfR::getFIX(sumstats, getINFO = TRUE))

      # Parse info field into separate columns
      message("Parsing info field...")
      fixstats$gene_id <- gsub("Gene=","", split_string_delim(fixstats$INFO, ";", 1))
      fixstats$gene_symbol <- gsub("GeneSymbol=","", split_string_delim(fixstats$INFO, ";", 2))
      fixstats$p.value <- gsub("Pvalue=","", split_string_delim(fixstats$INFO, ";", 3))
      fixstats$beta <- gsub("Beta=","", split_string_delim(fixstats$INFO, ";", 4))
      fixstats<- fixstats[,c("CHROM","POS","ID","REF","ALT","FILTER",
                                        "gene_id","gene_symbol","p.value","beta")]
      names(fixstats) <- c("chr.hg19","pos.hg19","rsid","ref_allele","alt_allele","filter",
                                     "gene_id","gene_symbol","p.value","beta")

      # liftover to hg38, flag/drop any that do not liftover
      # Set up GRanges object
      message("Defining genomic ranges to convert...")
      hg19.coord <- data.frame(chr = fixstats$chr.hg19,
                               start = fixstats$pos.hg19,
			       end = fixstats$pos.hg19,
			       metadata = paste0(rownames(fixstats),"|",fixstats$chr.hg19,":",fixstats$pos.hg19))
      hg19.ranges <- GenomicRanges::makeGRangesFromDataFrame(hg19.coord, keep.extra.columns = TRUE)
      GenomeInfoDb::seqlevelsStyle(hg19.ranges) <- "UCSC"

      message("Lifting over hg19 genomic ranges to hg38 coordinates...")
      hg38.ranges <- liftOver(hg19.ranges, chain)
      hg38.ranges.df <- lapply(hg38.ranges, as.data.frame)
      hg38.coord <- dplyr::bind_rows(hg38.ranges.df, .id="rowid")

      # parse metadata column back to hg19 coordinates and merge with fixstats
      message("Annotating hg38 coordinates to eQTL data...")
      hg38.coord$pos.hg19 <- split_string_delim(hg38.coord$metadata, ":", 2)
      hg38.coord$chr.hg19 <- split_string_delim(split_string_delim(hg38.coord$metadata, "|", 2), ":", 1)
      hg38.coord <- hg38.coord[,c("rowid","seqnames","start","end","pos.hg19","chr.hg19")]

      names(hg38.coord) <- c("range_index","chr.hg38","start.hg38","stop.hg38","pos.hg19","chr.hg19")
      fixstats.hg38 <- merge(fixstats, hg38.coord, by=c("chr.hg19","pos.hg19"), all.x=TRUE)
      dice.sumstats[[c]] <- fixstats.hg38
      message("Done!")
  }
  dice.all <- dplyr::bind_rows(dice.sumstats, .id="cell_type")

# save flatfile
save(dice.all, file="/orange/brusko/jrbnewman/software/dev/DICE_eqtl.rda")
message("All DICE data processed and saved to an RDA file!")


# TIGER eQTL
   chrom.list <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15",
                   "16", "17", "18", "19", "20", "21", "22", "23_female", "23_male")

   tiger.eqtl.list <- list()

   message("Processing TIGER eQTL data...")

for (c in chrom.list) {
       message(paste0("Reading TIGER eQTL data for chromosome ", c))
       tiger <- read.table(gzfile(paste0(in.path,"/tiger_summary_statistics/tiger_eqtl_stats_chr",c,".tsv.gz")),
                           sep="\t", header=TRUE)
       tiger <- tiger[tiger$p.value < 0.1,]
       tiger$Effect_allele <- toupper(tiger$Effect_allele)
       tiger$Non.effect_allele <- toupper(tiger$Non.effect_allele)

       tiger.eqtl.list[[c]] <- tiger
       }
    
tiger.merge <- dplyr::bind_rows(tiger.eqtl.list, .id="chrom_dataset")

save(tiger.merge, file="/orange/brusko/jrbnewman/software/dev/TIGER_eqtl.hg19.rda")




#####







       # use BioMart to get SNP IDs
       # TIGER has data as hg19 coordinates, so we'll lift them over to hg38 like with the DICE data
       message("Lifting over coordinates from hg19 to hg38")
       names(tiger) <- c("gene_id","chr.hg19","pos.hg19","effect_allele","noneffect_allele",
                         "weight", "Z_score","p_value","Direction","HetISq","HetChiSq","HetDF",
			 "HetPVal")
       tiger$chr.hg19 <- paste0("chr",tiger$chr.hg19)

       
       message("Defining genomic ranges to convert...")
       hg19.coord <- data.frame(chr = tiger$chr.hg19,
                               start = tiger$pos.hg19,
			       end = tiger$pos.hg19,
			       metadata = paste0(rownames(tiger),"|",tiger$chr.hg19,":",tiger$pos.hg19))
      hg19.ranges <- GenomicRanges::makeGRangesFromDataFrame(hg19.coord, keep.extra.columns = TRUE)
      GenomeInfoDb::seqlevelsStyle(hg19.ranges) <- "UCSC"

      message("Lifting over coordinates...")
      hg38.ranges <- liftOver(hg19.ranges, chain)
      message("Parsing output...")
      hg38.ranges.df <- lapply(hg38.ranges, as.data.frame)
      message("Binding rows...")
      hg38.coord <- dplyr::bind_rows(hg38.ranges.df, .id="rowid")

      # parse metadata column back to hg19 coordinates and merge with fixstats
      hg38.coord$pos.hg19 <- split_string_delim(hg38.coord$metadata, ":", 2)
      hg38.coord$chr.hg19 <- split_string_delim(split_string_delim(hg38.coord$metadata, "|", 2), ":", 1)
      hg38.coord <- hg38.coord[,c("rowid","seqnames","start","end","pos.hg19","chr.hg19")]

      tiger.hg38 <- merge(tiger, hg38.coord, by=c("chr.hg19","pos.hg19"), all.x=TRUE)

      message("Fetching possible rsIDs from BioMart...")
      tiger$range <- paste0(tiger$chr.hg19,":",tiger$pos.hg19,"-",tiger$pos.hg19)
      snp_data <- biomaRt::getBM(attributes = c("refsnp_id","chr_name","chrom_start","chrom_end","allele"),
                                    filters = c("chr_name","start","end"),
				    values = list(tiger.hg38$seqnames, tiger.hg38$start, tiger.hg38$end),
				    mart = snpMart)
				    
	    }



M2.vcf
MONOCYTES.vcf
NK.vcf
summary_statistics_tiger.tar.gz
TFH.vcf
TH17.vcf
TH1.vcf
TH2.vcf
THSTAR.vcf
tiger_case_stats.txt
TREG_MEM.vcf
TREG_NAIVE.vcf
