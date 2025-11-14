## Query the PRS catalogue using the "quincunx" R package
## May want to incorporate this into toolbox

#install.packages("quincunx")

#library("quincunx")

# Set trait
#trait <- "type 2 diabetes"

# Get trait EFO ID
#trait_results <- quincunx::get_traits(trait_term = trait)

# Check if traits found

#if (nrow(trait_results@traits) == 0) {
#   stop("No traits found for search term.")
#   }

# Get EFO ID
#efo_id <- trait_results@traits$efo_id[1]

# get all PGS scores associated with EFO ID
#all_scores <- quincunx::get_scores(efo_id = efo_id)

#variants.per.prs <- all_scores@scores[c("pgs_id","n_variants")]

#prs.sub.10k <- variants.per.prs[variants.per.prs$n_variants <= 100000,]

#pgs.list <- unique(prs.sub.10k$pgs_id)


#scores.list <- list()
#for (p in pgs.list) {
#roz <- quincunx::read_scoring_file(source = p, harmonized = TRUE, assembly="GRCh38")
#scores.list[[p]] <- roz[[p]]$data
#}

#repeats <- c()

#for (p in pgs.list) {
#   if(nrow(scores.list[[p]]) == 0) {
#       repeats <- c(repeats, p)
#       }
#       }
    
#for (p in repeats) {
#p="PGS003095"
#roz <- quincunx::read_scoring_file(source = p, harmonized = TRUE, assembly="GRCh38")
#scores.list[[p]] <- roz[[p]]$data
#}


# Import list of TOPMed SNPs
snp.by.chr <- list()

for (chr in (1:22)) {
    snp.by.chr[[chr]] <- read.table(paste0("/orange/brusko/UFDIchip/pipeline_2.0/TOPMed_imputation/QC87_95_96_98_99/chr_",chr,"/chr",chr,".bim"),
                      sep = "\t", header=F)
    snp.by.chr[[chr]]$V2 <- ifelse(snp.by.chr[[chr]]$V2 == ".", paste0(c(snp.by.chr[[chr]]$V1,
                                                                         snp.by.chr[[chr]]$V4,
								         snp.by.chr[[chr]]$V5,
								         snp.by.chr[[chr]]$V6), sep=":"),
							        snp.by.chr[[chr]]$V2)
    names(snp.by.chr[[chr]]) <- c("chromosome","rsID","cM","position","A1","A2")
    }
    
library(dplyr)

topmed.snp <- dplyr::bind_rows(snp.by.chr, .id = "chr")

rm(snp.by.chr)
save.image(file = "/orange/brusko/jrbnewman/prs_workspace.RData")

load(file = "/orange/brusko/jrbnewman/prs_workspace.RData")

# iterate over SNP lists and for each output a vector of:
# total SNPs in PRS
# SNPs match TOPMed output
# Percent match

match.stats <- list()

topmed.snp$coord <- paste0(topmed.snp$chromosome, ":", topmed.snp$position)
topmed.snp$allelic1 <- paste0(topmed.snp$coord,":",topmed.snp$A1,":",topmed.snp$A2)
topmed.snp$allelic2 <- paste0(topmed.snp$coord,":",topmed.snp$A2,":",topmed.snp$A1)

gc()

for (p in names(scores.list)) {
     # set up a coordinate ID and an allelic ID
     message(p)
     if (nrow(scores.list[[p]] > 0)) {
     message("Creating coordinate and allelic IDs")
     coord.score <- paste0(scores.list[[p]]$hm_chr, ":", scores.list[[p]]$hm_pos)
     allelic.score <- paste0(coord.score,":",scores.list[[p]]$other_allele,":",scores.list[[p]]$effect_allele)

     # Reset variables
     message("Reset some variables just in case")
     n_snps <- NULL
     n_match_rsID <- NULL
     perc_match_rsID <- NULL
     n_match_coord <- NULL
     perc_match_coord <- NULL
     n_match_allele1 <- NULL
     perc_match_allele1 <- NULL
     n_match_allele2 <- NULL
     perc_match_allele2 <- NULL
     n_match_alleleTotal <- NULL
     perc_match_alleleTotal <- NULL
     out.df <- NULL

     # Count total SNPs in PRS
     message("Counting total SNPs in PRS")
     n_snps <- nrow(scores.list[[p]])

     # Matches based on rsID
     message("Matching SNPs based on rsID")
     n_match_rsID <- length(intersect(scores.list[[p]]$hm_rsID, topmed.snp$rsID))
     perc_match_rsID <- n_match_rsID / n_snps

     # Matches based on coordinate only (alleles not considered)
     message("Matching SNPs based on chromosome and position")
     n_match_coord <- length(intersect(coord.score,topmed.snp$coord))
     perc_match_coord <- n_match_coord / n_snps

     # Matches based on coordinate and allele (relies on alleles being coded the same)
     message("Matching SNPs based on chromosome, position, and alleles")
     message("Allele 1")
     allele1.match <- intersect(allelic.score, topmed.snp$allelic1)
     message("Allele 2")
     allele2.match <- intersect(allelic.score, topmed.snp$allelic2)
     n_match_allele1 <- length(allele1.match)
     n_match_allele2 <- length(allele2.match)
     n_match_alleleTotal <- length(unique(c(allele1.match, allele2.match)))

     perc_match_allele1 <- n_match_allele1 / n_snps
     perc_match_allele2 <- n_match_allele1 / n_snps
     perc_match_alleleTotal <- n_match_alleleTotal / n_snps

     message("Output dataframe")
     out.df <- data.frame(PRS_ID = as.character(p),
                          n_snps = as.integer(n_snps),
			  n_match_rsID = as.integer(n_match_rsID),
			  perc_match_rsID = as.numeric(perc_match_rsID),
			  n_match_coord = as.integer(n_match_coord),
			  perc_match_coord = as.numeric(perc_match_coord),
			  n_match_allele1 = as.integer(n_match_allele1),
			  perc_match_allele1 = as.numeric(perc_match_allele1),
			  n_match_allele2 = as.integer(n_match_allele2),
			  perc_match_allele2 = as.numeric(perc_match_allele2),
			  n_match_alleleTotal = as.integer(n_match_alleleTotal),
			  perc_match_alleleTotal = as.numeric(perc_match_alleleTotal)
			  			  )
     match.stats[[p]] <- out.df
     } else { next }
     }

matched.df <- dplyr::bind_rows(match.stats, .id = "EFO_ID")

saveRDS(matched.df, "/orange/brusko/jrbnewman/T2D_grs_matches_harmonized.rds")

############################





matched.df <- readRDS( "/orange/brusko/jrbnewman/T2D_grs_matches_harmonized.rds")


## then, try matching rsIDs to coordinates and using a unified identifier using SNPtools functions

## CISTRALATIN
# 10X data processing/QC
# SNP analysis functions (annotation, LD, coming soon: motif analysis)
# DiffExp analysis + GSEA with plots, wrapped into a single function
# Flow cytometry functions
# DNA and protein-based multiple sequence alignment packages with rich annotation
# Genetic database queries: eQTL, GWAS, PRS catalogs, GTEx
# SNP/gene/cell type-based queries of DICE and TIGER eQTL/sQTL/ASE summary statistics
# More coming!

#Considering rolling this in

##########################################


# Load packages
library(DBI)
library(RSQLite)

 db_connection <- dbConnect(RSQLite::SQLite(), "/orange/brusko/jrbnewman/software/dev/variants.db")

 table_names <- dbListTables(db_connection)

# read each table into a list

prsedm.list <- list()

for (table_name in table_names) {
  # Read table into R data frame
  query <- paste0("SELECT * FROM ", table_name)
  data_frame <- dbGetQuery(db_connection, query) # For DBI/RSQLite
  prsedm.list[[table_name]] <- data_frame
}

dbDisconnect(db_connection) # For DBI/RSQLite

for (p in names(prsedm.list)) {
   print(colnames(prsedm.list[[p]]))
   }

# keep only linear models

prsedm.list[["cdgrs_hlainteraction"]] <- NULL
prsedm.list[["hla_ranking_klitz"]] <- NULL
prsedm.list[["pbcgrs_hlainteraction"]] <- NULL
prsedm.list[["t1dgrs2_hlainteraction"]] <- NULL


for (p in names(prsedm.list)) {
   prsedm.list[[p]]$rsID <- prsedm.list[[p]]$rsid
     message(p)
     if (nrow(prsedm.list[[p]] > 0)) {
     message("Creating coordinate and allelic IDs")
     coord.score<- paste0(prsedm.list[[p]]$contig_id, ":", prsedm.list[[p]]$position_hg38)
     allelic.score <- paste0(coord.score,":",prsedm.list[[p]]$effect_allele,":",prsedm.list[[p]]$effect_allele)

     # Reset variables
     message("Reset some variables just in case")
     n_snps <- NULL
     n_match_rsID <- NULL
     perc_match_rsID <- NULL
     n_match_coord <- NULL
     perc_match_coord <- NULL
     n_match_allele1 <- NULL
     perc_match_allele1 <- NULL
     n_match_allele2 <- NULL
     perc_match_allele2 <- NULL
     n_match_alleleTotal <- NULL
     perc_match_alleleTotal <- NULL
     out.df <- NULL

     # Count total SNPs in PRS
     message("Counting total SNPs in PRS")
     n_snps <- nrow(prsedm.list[[p]])

     # Matches based on rsID
     message("Matching SNPs based on rsID")
     n_match_rsID <- length(intersect(prsedm.list[[p]]$rsid, topmed.snp$rsID))
     perc_match_rsID <- n_match_rsID / n_snps

     # Matches based on coordinate only (alleles not considered)
     message("Matching SNPs based on chromosome and position")
     n_match_coord <- length(intersect(coord.score,topmed.snp$coord))
     perc_match_coord <- n_match_coord / n_snps

     # Matches based on coordinate and allele (relies on alleles being coded the same)
     message("Matching SNPs based on chromosome, position, and alleles")
     message("Allele 1")
     allele1.match <- intersect(allelic.score, topmed.snp$allelic1)
     message("Allele 2")
     allele2.match <- intersect(allelic.score, topmed.snp$allelic2)
     n_match_allele1 <- length(allele1.match)
     n_match_allele2 <- length(allele2.match)
     n_match_alleleTotal <- length(unique(c(allele1.match, allele2.match)))

     perc_match_allele1 <- n_match_allele1 / n_snps
     perc_match_allele2 <- n_match_allele1 / n_snps
     perc_match_alleleTotal <- n_match_alleleTotal / n_snps

     message("Output dataframe")
     out.df <- data.frame(PRS_ID = as.character(p),
                          n_snps = as.integer(n_snps),
			  n_match_rsID = as.integer(n_match_rsID),
			  perc_match_rsID = as.numeric(perc_match_rsID),
			  n_match_coord = as.integer(n_match_coord),
			  perc_match_coord = as.numeric(perc_match_coord),
			  n_match_allele1 = as.integer(n_match_allele1),
			  perc_match_allele1 = as.numeric(perc_match_allele1),
			  n_match_allele2 = as.integer(n_match_allele2),
			  perc_match_allele2 = as.numeric(perc_match_allele2),
			  n_match_alleleTotal = as.integer(n_match_alleleTotal),
			  perc_match_alleleTotal = as.numeric(perc_match_alleleTotal)
			  			  )
     match.stats[[p]] <- out.df
     } else { next }
     }

matched.df <- dplyr::bind_rows(match.stats, .id = "EFO_ID")

saveRDS(matched.df, "/orange/brusko/jrbnewman/T2D_grs_matches_to_PRSedm_lists.rds")




}

# Mahajan SNPs from Gaulton: what do these better match?


match.mahajan <- list()
gc()

for (p in names(scores.list)) {
     # set up a coordinate ID and an allelic ID
     message(p)
     if (nrow(scores.list[[p]] > 0)) {
     message("Creating coordinate and allelic IDs")
     coord.score <- paste0(scores.list[[p]]$hm_chr, ":", scores.list[[p]]$hm_pos)
     allelic.score <- paste0(coord.score,":",scores.list[[p]]$other_allele,":",scores.list[[p]]$effect_allele)

     # Reset variables
     message("Reset some variables just in case")
     n_snps <- NULL
     n_match_rsID <- NULL
     perc_match_rsID <- NULL
     n_match_coord <- NULL
     perc_match_coord <- NULL
     n_match_allele1 <- NULL
     perc_match_allele1 <- NULL
     n_match_allele2 <- NULL
     perc_match_allele2 <- NULL
     n_match_alleleTotal <- NULL
     perc_match_alleleTotal <- NULL
     out.df <- NULL

     # Count total SNPs in PRS
     message("Counting total SNPs in PRS")
     n_snps <- nrow(scores.list[[p]])

     # Matches based on rsID
     message("Matching SNPs based on rsID")
     n_match_rsID <- length(intersect(scores.list[[p]]$hm_rsID, prsedm.list$t2d_mahajan22_ma$rsid))
     perc_match_rsID <- n_match_rsID / n_snps

     message("Output dataframe")
     out.df <- data.frame(PRS_ID = as.character(p),
                          n_snps = as.integer(n_snps),
			  n_match_rsID = as.integer(n_match_rsID),
			  perc_match_rsID = as.numeric(perc_match_rsID)  )
     match.mahajan[[p]] <- out.df
     } else { next }
     }
