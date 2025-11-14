

load(file = "/orange/brusko/jrbnewman/ASE_T1D_variants.Rdata")

library(biomaRt)




ensembl_snp = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")

snp.features <- c("refsnp_id","chr_name","chrom_start","chrom_end","allele",
		  "allele_1","minor_allele","clinical_significance",
		  "associated_gene","ensembl_gene_stable_id",
		  "ensembl_gene_name","ensembl_transcript_stable_id",
		  "ensembl_type","consequence_type_tv","cdna_start",
		  "cdna_end","cds_start","cds_end","distance_to_transcript",
		  "polyphen_prediction","polyphen_score","sift_prediction",
		  "sift_score","reg_feature_stable_id","reg_allele_string",
		  "reg_consequence_types")

snp.ids <- unique(ld_network.annot$match.RefSNP_id)

	rs_results <- data.frame()
        rs_results <- biomaRt::getBM(attributes = snp.features,
	                                 filters = "snp_filter",
					 values = snp.ids,
					 mart=ensembl_snp)

	# Drop scaffold SNPs
	allowed.chr.list <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16",
	                    "17","18","19","20","21","22","X","Y","MT","chr1","chr2","chr3","chr4",
			    "chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13",
			    "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22",
			    "chrX","chrY","chrM")
	
	rs_results <- rs_results[rs_results$chr_name %in% allowed.chr.list,]
        rs_results$match.RefSNP_id <- rs_results$refsnp_id
	rs_results$ref_allele <- ifelse(rs_results$allele_1 == "",
	                            split_string_delim(rs_results$allele, "/", 1),
				    rs_results$allele_1)
	rs_results$minor_allele <- ifelse(rs_results$minor_allele == "",
	                            split_string_delim(rs_results$allele, "/", 2),
				    rs_results$minor_allele)
	rs_results$alt_alleles <- sapply(strsplit(rs_results$allele, "/"), function(x) paste(x[-1], collapse=","))

	ld_network.bm <- merge(ld_network.annot, rs_results, by="match.RefSNP_id", all.x=TRUE)
	

# subset transcribed SNPs

xscript.snp <- ld_network.bm[ld_network.bm$consequence_type_tv %in%
               c("synonymous_variant","5_prime_UTR_variant","3_prime_UTR_variant","missense_variant","splice_acceptor_variant","stop_gained","frameshift_variant","splice_donor_variant","start_lost"),]

xscript.snp2 <- merge(xscript.snp, ld_network, by.x="match.RefSNP_id", by.y="snp_id", all.x=T)


saveRDS(ld_network.bm, "/orange/brusko/jrbnewman/T1D_ASE_SNP_proxies.rds")

xscript.snp2 <- xscript.snp2[is.na(xscript.snp2$associated_gene) == FALSE,]
xscript.snp2 <- xscript.snp2[xscript.snp2$associated_gene != "",]
xscript.snp2 <- xscript.snp2[xscript.snp2$associated_gene != "intergenic",]
xscript.snp2 <- xscript.snp2[xscript.snp2$associated_gene != "intergenicN1L",]

snp.queries <- unique(xscript.snp2$match.RefSNP_id)



library("tidyverse")
library("httr")
library("glue")
library("dplyr")
library("coloc")
library("jsonlite")
library(ggplot2)
library("ggrepel")
library("utils")


## Get list of available datasets to query
	
	# Change parameters
	max_pulled_rows = 1000 #All datasets will be pulled if this parameter is bigger than the actual number of datasets
	
	URL = glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size={max_pulled_rows}")
	
	# Make a request
	r <- GET(utils::URLencode(URL), accept_json())
	# Check status
	status_code(r)
	
	# Extract content
	cont <- content(r, "text", encoding = "UTF-8")
	# Convert content to dataframe
	datasets <- fromJSON(cont)
	knitr::kable(head(datasets, n = 20), format="markdown")
	
	request_datasets_from_api <- function(study_id = "",
	                                      quant_method = "",
	                                      sample_group = "",
	                                      tissue_id = "",
	                                      study_label = "",
	                                      tissue_label = "",
	                                      condition_label = "",
					      size = 1000,
					      start = 0) {

	  parameter_values = c(study_id,quant_method,sample_group,tissue_id,study_label, 
	                       tissue_label,condition_label)
	  parameter_names = c('study_id','quant_method','sample_group','tissue_id',
	                      'study_label','tissue_label','condition_label')
	  
	  while (T) {
	    URL = glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size={size}&start={start}")
	    
	    #Adding defined parameters to the request
	    for (i in 1:length(parameter_values)) {
	      par = parameter_values[i]
	      par_name = parameter_names[i]
	      if (par != "")
	        URL = glue("{URL}&{par_name}={par}")
	    }
	    
	    r <- GET(utils::URLencode(URL), accept_json())
	    cont <- content(r, "text", encoding = "UTF-8")
	    
	    # If the request was unsuccessful
	    if (status_code(r) != 200) {
	      #If we get no results at all, print error
	      if (start == 0) {
	        print(glue("Error {status_code(r)}"))
	        print(cont)
	        return ()
	      }
	      #else just break
	      break
	    }
	    
	    cont_df <- fromJSON(cont)
	    
	    if (start == 0) {
	      responses <- cont_df
	    }
	    else{
	      responses <- rbind(responses, cont_df)
	    }
	    start <- start + size
	  }
	  return(responses)
	}


datasets_ge = request_datasets_from_api(quant_method = "ge")


get_assoc_over_datasets <- function(datasets, variant, variant_type = c("coordinate","rsid")) {
  size = 1000
  first = T
  final_df <- list()
  for (i in rownames(datasets)) {
    row = datasets[i, ]
    dataset_id = row$dataset_id
    
        if (variant_type == "coordinate") {
    URL = glue(
          "https://www.ebi.ac.uk/eqtl/api/v2/datasets/{dataset_id}/associations?size={size}&variant={variant}"
	  )
	  } else if (variant_type == "rsid") {
    URL = glue(
          "https://www.ebi.ac.uk/eqtl/api/v2/datasets/{dataset_id}/associations?size={size}&rsid={variant}"
	  ) } else {
	  message("Please specify if variant is in [chr_position_ref_alt] coordinate or rsid!")
	  break }
    
        
    r <- GET(utils::URLencode(URL), accept_json())
    
    if (status_code(r) != 200) {
      next
    }
    
    cont <- content(r, "text", encoding = "UTF-8")
    cont_df <- fromJSON(cont)
    cont_with_metadata <- cbind(cont_df, row)
    final_df[[i]] <- cont_with_metadata
    }
    final_df <- dplyr::bind_rows(final_df, .id="dataset_num")
    return(final_df)
}


eQTL.list <- list()
for (s in snp.queries) {
    eQTL.list[[s]] <- get_assoc_over_datasets(datasets_ge, s, variant_type = "rsid")
    }
    
    
saveRDS(eQTL.list, "/orange/brusko/jrbnewman/eQTLs_for_T1D_ASE.rds")


eQTL.plots <- list()

for (s in names(eQTL.list)) {
   if (nrow(eQTL.list[[s]]) > 0) {
    eQTL.plots[[s]] <- eQTL.list[[s]] %>%
   unite(study, study_label, tissue_label, sep = " ", remove = FALSE) %>%
   ggplot(aes(x=beta, y=nlog10p, label=study)) +
   geom_point() +
   geom_text_repel() +
   ggtitle(paste0("Allelic effect of ", s, " on expression"))
  }
  }

eQTL.snps <- names(eQTL.plots)
   
xscript.snp2$has_eQTL <- ifelse(xscript.snp2$match.RefSNP_id %in% eQTL.snps, TRUE, FALSE)

# Query DICE data

load("/orange/brusko/jrbnewman/software/dev/DICE_eqtl.rda")

# subset DICE data for queried SNPs

dice.xs.snp <- dice.all[dice.all$rsid %in% snp.queries,]
dice.xs.eQTL <- dice.all[dice.all$rsid %in% eQTL.snps,]


dice.proxies <- dice.all[dice.all$rsid %in% unique(ld_network.annot$match.RefSNP_id),]

### import TIGER data

tiger.ase <- read.table("/orange/brusko/jrbnewman/software/dev/datasets/tiger_case_stats.txt", header=T)



# Load libraries
  library(dplyr)
  library(biomaRt)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)

# DICE and TIGER use hg19 coordinates, so we're going to lift those
# over to hg38. First, we set up the chain file
  message("Downloading hg19-to-hg38 chain file and loading into R...")
  chain_url <- "https://hgdownload.soe.ucsc.edu/gbdb/hg38/liftOver/hg38ToHg19.over.chain.gz"
  chain_file_local <- file.path(tempdir(), "hg38ToHg19.over.chain.gz")
  download.file(chain_url, destfile = chain_file_local)

  # Uncompress the chain file
  system(paste("gzip -d -f", chain_file_local))
  chain_file_uncompressed <- gsub("\\.gz$", "", chain_file_local)
  
  # import chain
  chain <- rtracklayer::import.chain(chain_file_uncompressed)
  message("Chain file loaded!")

ranges <- data.frame(CHROM = ld_network.annot$chr,
                     START=as.numeric(ld_network.annot$start),
		     END=as.numeric(ld_network.annot$end),
		     METADATA=ld_network.annot$match.RefSNP_id)
gr.hg38 <- GRanges(seqnames = ranges$CHROM,
                   ranges = IRanges(start = ranges$START, end = ranges$END),
		   strand = "*",
		   metadata = ranges$METADATA)

gr.hg19 <- liftOver(gr.hg38, chain)
gr.hg19.u <- unlist(gr.hg19)

snp.query.hg19 <- data.frame(chr=as.character(seqnames(gr.hg19.u)),
                             start = start(gr.hg19.u),
			     end = end(gr.hg19.u),
			     snp_id = mcols(gr.hg19.u)$metadata)

ase2keep <- list()

for (s in 1:nrow(snp.query.hg19)) {
    chr.roz <- snp.query.hg19$chr[s]
    pos.roz <- snp.query.hg19$start[s]
    snp.roz <- snp.query.hg19$snp_id[s]
    keeper <- tiger.ase[(tiger.ase$Chrom == chr.roz) & (tiger.ase$Position == pos.roz),]
    if (nrow(keeper) > 0) {
    ase2keep[[s]] <- keeper
    ase2keep[[s]]$snp_id <- snp.roz
    } else { ase2keep[[s]] <- data.frame(ID=c(),
					 Chrom=c(),
					 Position=c(),
					 Gene=c(),
					 Effect_allele=c(),
					 Non.effect_allele=c(),
					 FDR=c(),
					 p.value=c(),
					 Z.score=c(),
					 AR.sigma=c(),
					 N_significant=c(),
					 N_heterozygous = c(),
					 snp_id=c())
					 }
}

tiger.ase.t1d.risk <- dplyr::bind_rows(ase2keep, .id="index")

write.table(snp.query.hg19$start, "/orange/brusko/jrbnewman/software/dev/datasets/hg19_pos.txt")



gunzip -c tiger_eqtl_stats_chr10.tsv.gz | head -n1 > ../tiger_eqtl_olap.txt
for FILE in *.tsv.gz; do echo $FILE; gunzip -c $FILE | grep -f ../hg19_pos.txt >> ../tiger_eqtl_olap.txt; done


tiger.eqtl <- read.table("/orange/brusko/jrbnewman/software/dev/datasets/tiger_eqtl_olap.txt", header=T)
tiger.eqtl$coord <- paste0("chr",tiger.eqtl$Chrom,":",tiger.eqtl$Position)
snp.query.hg19$coord <- paste0(snp.query.hg19$chr,":",snp.query.hg19$start)
tiger.eqtl.t1d.risk <- merge(tiger.eqtl, snp.query.hg19, by="coord")


library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_symbols <- getBM(attributes = c("hgnc_symbol","ensembl_gene_id"), filters="ensembl_gene_id",
values = unique(tiger.eqtl.t1d.risk$Gene), mart=ensembl)

gene_symbols$hgnc_symbol <- ifelse(gene_symbols$hgnc_symbol == "", gene_symbols$ensembl_gene_id,
                                   gene_symbols$hgnc_symbol)

tiger.eqtl.t1d.risk2 <- merge(tiger.eqtl.t1d.risk, gene_symbols, by.x="Gene", by.y="ensembl_gene_id", all.x=T)
tiger.ase.t1d.risk2 <- merge(tiger.ase.t1d.risk, gene_symbols, by.x="Gene", by.y="ensembl_gene_id", all.x=T)



save.image("/orange/brusko/jrbnewman/T1D_ASE_working.Rdata")


###############################################################

load("/orange/brusko/jrbnewman/T1D_ASE_working.Rdata")


# I want a table with the following info:
# Proxy/SNP info (genes, features, transcripts)
# SNP LD region info
# If SNP is DICE eQTL
# If SNP is TIGER eQTL
# If SNP is TIGER ASE


ld_network.bm2 <- merge(x=ld_network.bm, y=ld_network, by.x="match.RefSNP_id", by.y="snp_id")
annot.all <- read.csv("/orange/brusko/jrbnewman/full_SNP_annotation_example.csv")


tiger.ase.snp <- unique(tiger.ase.t1d.risk2[as.numeric(tiger.ase.t1d.risk2$p.value) <0.05,]$snp_id)
tiger.eqtl.snp <- unique(tiger.eqtl.t1d.risk2[as.numeric(tiger.eqtl.t1d.risk2$p.value) < 0.05,]$snp_id)
dice.snp <- unique(dice.proxies[as.numeric(dice.proxies$p.value)<0.05,]$rsid)

tiger.eqtl.gene <- unique(tiger.eqtl.t1d.risk2[as.numeric(tiger.eqtl.t1d.risk2$p.value) < 0.05,]$Gene)
dice.gene <- unique(dice.proxies[as.numeric(dice.proxies$p.value)<0.05,]$gene_id)

ld_network.bm2$coordinate <- paste0("chr",gsub("chr","",ld_network.bm2$chr),":",ld_network.bm2$start,":",ld_network.bm2$end)
annot.all$coordinate <- paste0("chr",annot.all$seqnames,":",annot.all$start,":",annot.all$end)

snps.full.annot <- merge(x=ld_network.bm2, y=annot.all, by="coordinate")

snps.full.annot$flag.TIGER.ASE <- ifelse(snps.full.annot$match.RefSNP_id %in% tiger.ase.snp, TRUE, FALSE)
snps.full.annot$flag.TIGER.eQTL <- ifelse(snps.full.annot$match.RefSNP_id %in% tiger.eqtl.snp, TRUE, FALSE)
snps.full.annot$flag.DICE.eQTL <- ifelse(snps.full.annot$match.RefSNP_id %in% dice.snp, TRUE, FALSE)

snps.full.annot$flag.TIGER.eQTL.gene <- ifelse(any(tiger.eqtl.gene %in% snps.full.annot$match.RefSNP_id), TRUE, FALSE)
snps.full.annot$flag.DICE.eQTL.gene <- ifelse(any(dice.gene %in% snps.full.annot$match.RefSNP_id), TRUE, FALSE)





# install.packages("VennDiagram") # Only run if the package is not installed


library(VennDiagram)


    list_of_sets <- list(
      Set1 = tiger.ase.snp[!is.na(tiger.ase.snp)],
      Set2 = tiger.eqtl.snp[!is.na(tiger.eqtl.snp)],
      Set3 = dice.snp[!is.na(dice.snp)]
    )

    venn.diagram(
      x = list_of_sets,
      category.names = c("TIGER islet ASE", "TIGER islet eQTL", "DICE immune eQTL"),
      filename = "/orange/brusko/jrbnewman/T1D_ASE_eQTL_Venn.png",
      output = TRUE # Set to FALSE to display directly in RStudio plot pane
    )

  # devtools::install_github("gaospecial/ggVennDiagram")
    library(ggVennDiagram)


   ggVennDiagram(list_of_sets,
                  category.names = c("TIGER islet ASE", "TIGER islet eQTL", "DICE immune eQTL"),
                  label_alpha = 0) +
      ggplot2::scale_fill_gradient(low = "blue", high = "red") # Customize fill colors


# Then summarize by region

tiger.ase.region <- unique(snps.full.annot[snps.full.annot$flag.TIGER.ASE==TRUE,]$LD_group_id)
tiger.eqtl.region <- unique(snps.full.annot[snps.full.annot$flag.TIGER.eQTL==TRUE,]$LD_group_id)
dice.eqtl.region <- unique(snps.full.annot[snps.full.annot$flag.DICE.eQTL==TRUE,]$LD_group_id)



    list_of_sets <- list(
      Set1 = tiger.ase.region[!is.na(tiger.ase.region)],
      Set2 = tiger.eqtl.region[!is.na(tiger.eqtl.region)],
      Set3 = dice.eqtl.region[!is.na(dice.eqtl.region)]
    )

    venn.diagram(
      x = list_of_sets,
      category.names = c("TIGER islet ASE", "TIGER islet eQTL", "DICE immune eQTL"),
      filename = "/orange/brusko/jrbnewman/T1D_ASE_eQTL_Venn_by_region.png",
      output = TRUE # Set to FALSE to display directly in RStudio plot pane
    )


# THEN repeat on ONLY SNPs in transcribed regions

# subset exonic, UTR, and splice site SNPs


snps.exon <- unique(snps.full.annot[grepl("exon",snps.full.annot$feature_overlap) & !is.na(snps.full.annot$match.RefSNP_id),]$match.RefSNP_id)
snps.utr <- unique(snps.full.annot[grepl("UTR",snps.full.annot$feature_overlap)& !is.na(snps.full.annot$match.RefSNP_id),]$match.RefSNP_id)
snps.splice <- unique(snps.full.annot[(snps.full.annot$is_splice_site == TRUE) & !is.na(snps.full.annot$match.RefSNP_id),]$match.RefSNP_id)


    list_of_sets <- list(
      Set1 = snps.exon,
      Set2 = snps.utr,
      Set3 = snps.splice
    )

    venn.diagram(
      x = list_of_sets,
      category.names = c("Exon SNPs", "UTR SNPs", "Splice site SNPs"),
      filename = "/orange/brusko/jrbnewman/T1D_ASE_eQTL_Venn_transcribed_SNP_overlap.png",
      output = TRUE # Set to FALSE to display directly in RStudio plot pane
    )

snps.exon <- unique(snps.full.annot[grepl("exon",snps.full.annot$feature_overlap) & !is.na(snps.full.annot$match.RefSNP_id),])
snps.utr <- unique(snps.full.annot[grepl("UTR",snps.full.annot$feature_overlap)& !is.na(snps.full.annot$match.RefSNP_id),])
snps.splice <- unique(snps.full.annot[(snps.full.annot$is_splice_site == TRUE) & !is.na(snps.full.annot$match.RefSNP_id),])

snps.xs.annot <- unique(rbind(snps.exon, snps.utr, snps.splice))



tiger.ase.snp.xs <- unique(snps.xs.annot[(snps.xs.annot$flag.TIGER.ASE==TRUE) & !is.na(snps.xs.annot$match.RefSNP_id),]$match.RefSNP_id)
tiger.eqtl.snp.xs <- unique(snps.xs.annot[(snps.xs.annot$flag.TIGER.eQTL==TRUE) & !is.na(snps.xs.annot$match.RefSNP_id),]$match.RefSNP_id)
dice.eqtl.snp.xs <- unique(snps.xs.annot[(snps.xs.annot$flag.DICE.eQTL==TRUE) & !is.na(snps.xs.annot$match.RefSNP_id),]$match.RefSNP_id)





    list_of_sets <- list(
      Set1 = tiger.ase.snp.xs,
      Set2 = tiger.eqtl.snp.xs,
      Set3 = dice.eqtl.snp.xs
    )

    venn.diagram(
      x = list_of_sets,
      category.names = c("TIGER islet ASE", "TIGER islet eQTL", "DICE immune eQTL"),
      filename = "/orange/brusko/jrbnewman/T1D_ASE_eQTL_Venn_transcribed_SNP_eQTL_ASE.png",
      output = TRUE # Set to FALSE to display directly in RStudio plot pane
    )



# Then summarize by region

tiger.ase.region.xs <- unique(snps.xs.annot[snps.xs.annot$flag.TIGER.ASE==TRUE,]$LD_group_id)
tiger.eqtl.region.xs <- unique(snps.xs.annot[snps.xs.annot$flag.TIGER.eQTL==TRUE,]$LD_group_id)
dice.eqtl.region.xs <- unique(snps.xs.annot[snps.xs.annot$flag.DICE.eQTL==TRUE,]$LD_group_id)


    list_of_sets <- list(
      Set1 = tiger.ase.region.xs[!is.na(tiger.ase.region.xs)],
      Set2 = tiger.eqtl.region.xs[!is.na(tiger.eqtl.region.xs)],
      Set3 = dice.eqtl.region.xs[!is.na(dice.eqtl.region.xs)]
    )

    venn.diagram(
      x = list_of_sets,
      category.names = c("TIGER islet ASE", "TIGER islet eQTL", "DICE immune eQTL"),
      filename = "/orange/brusko/jrbnewman/T1D_ASE_eQTL_Venn_transcribed_SNPs_eqtl_ase_by_region.png",
      output = TRUE # Set to FALSE to display directly in RStudio plot pane
    )


### 


tiger.ase.t1d.risk2$flag.dice <- ifelse(tiger.ase.t1d.risk2$snp_id %in% dice.snp, "DICE eQTL", "Non-DICE SNP")
tiger.ase.t1d.risk2$ASE.P <-  ifelse(tiger.ase.t1d.risk2$p.value < 0.05, "TIGER ASE P<0.05", "TIGER ASE N.S.")
tiger.ase.t1d.risk2$DICE.ASE <- paste0(tiger.ase.t1d.risk2$flag.dice, ", ",tiger.ase.t1d.risk2$ASE.P)


library(ggplot2)
library(ggrepel)




png("/orange/brusko/jrbnewman/TIGER_ASE_colored_by_DICE_eqtl_status.png", width=500, height=500)
print(ggplot(data = tiger.ase.t1d.risk2, aes(x = Z.score, y = -log10(p.value), col = DICE.ASE, label = snp_id)) + geom_point(size = 2) +
  labs(color = 'Significance', x = expression("Z.score"), y = expression("-log"[10]*"p-value")) + # Add axis labels
  theme_classic() )
  dev.off()


tiger.ase.t1d.risk2.xs <- tiger.ase.t1d.risk2[tiger.ase.t1d.risk2$snp_id %in% tiger.ase.snp.xs,]

png("/orange/brusko/jrbnewman/TIGER_ASE_colored_by_DICE_eqtl_status_xs_SNP_only.png", width=500, height=500)
print(ggplot(data = tiger.ase.t1d.risk2.xs, aes(x = Z.score, y = -log10(p.value), col = DICE.ASE, label = snp_id)) + geom_point(size = 2) +
  labs(color = 'Significance', x = expression("Z.score"), y = expression("-log"[10]*"p-value")) + # Add axis labels
  theme_classic() )
  dev.off()


# Correlate ASE effect and eQTL effect


tiger.ase.eff <- tiger.ase.t1d.risk2[,c("snp_id","Z.score","p.value")]
tiger.eqtl.eff <- tiger.eqtl.t1d.risk2[,c("snp_id","Z.score","p.value")]
dice.eqtl.eff <- dice.proxies[,c("rsid","beta","p.value")]

names(tiger.ase.eff) <- c("snp_id","TIGER_ASE_Zscore","TIGER_ASE_P")
names(tiger.eqtl.eff) <- c("snp_id","TIGER_eQTL_Zscore","TIGER_eQTL_P")
names(dice.eqtl.eff) <- c("snp_id","DICE_beta","DICE_P")

tiger.ase2eqtl <- merge(tiger.ase.eff, tiger.eqtl.eff, by="snp_id")
tiger2dice <- merge(tiger.ase2eqtl, dice.eqtl.eff, by="snp_id")


tiger2dice2 <- tiger2dice[(tiger2dice$DICE_P < 0.05) &
                          (tiger2dice$TIGER_ASE_P  < 0.05) &
			   (tiger2dice$TIGER_eQTL_P  < 0.05),]
dim(tiger2dice2)			   

tiger.ase2eqtl2 <- tiger.ase2eqtl[(tiger.ase2eqtl$TIGER_ASE_P <0.05) &
		                  (tiger.ase2eqtl$TIGER_eQTL_P <0.05),]


ggplot(tiger.ase2eqtl, aes(x = TIGER_ASE_Zscore, y = TIGER_eQTL_Zscore)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, col = "blue") + # Customize color and remove CI
  labs(title = "TIGER eQTL vs ASE effect size",
       x = "TIGER ASE Z-score", y = "TIGER eQTL Z-score")


tiger2dice2 <- tiger2dice[(tiger2dice$DICE_P < 0.05) &
                          (tiger2dice$TIGER_ASE_P  < 0.05),]
dim(tiger2dice2)

ggplot(tiger2dice, aes(x = TIGER_ASE_Zscore, y = DICE_beta)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, col = "blue") + # Customize color and remove CI
  labs(title = "TIGER ASE vs DICE eQTL effect size",
       x = "TIGER ASE Z-score", y = "DICE eQTL beta")




tiger2dice2 <- tiger2dice[(tiger2dice$DICE_P < 0.05) &
			   (tiger2dice$TIGER_eQTL_P  < 0.05),]
dim(tiger2dice2)			   
ggplot(tiger2dice, aes(x = TIGER_eQTL_Zscore, y = DICE_beta)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, col = "blue") + # Customize color and remove CI
  labs(title = "TIGER vs DICE eQTL effect size",
       x = "TIGER eQTL Z-score", y = "DICE eQTL beta")




##### Extract SNPs for the 16 regions that have TIGER islet ASE


ase.region.snps <- snps.full.annot[snps.full.annot$LD_group_id %in% tiger.ase.region.xs,]
#3548 SNPs. This will take some time, so let's remove intronic and non-genic SNPs

ase.region.snps2 <- ase.region.snps[grepl("exon", ase.region.snps$feature_overlap),]

# 362. Still alot. Go back to all ASE region SNPs and subset just the sig ones

ase.region.snps2 <- ase.region.snps[ase.region.snps$flag.TIGER.ASE == TRUE,]

#54 SNPs is better

ase.snps <- unique(ase.region.snps2$match.RefSNP_id)




datasets_ge = request_datasets_from_api(quant_method = "ge")
datasets_ex = request_datasets_from_api(quant_method = "exon")
datasets_tx = request_datasets_from_api(quant_method = "tx")


dim(datasets_ge)
dim(datasets_ex)
dim(datasets_tx)

unique(datasets_ge$tissue_label)
unique(datasets_ex$tissue_label)
unique(datasets_tx$tissue_label)

query.tissue <- c("B cell","blood","CD16+ monocyte","CD4+ CTL cell","CD4+ memory T cell","CD4+ T cell","CD4+ TCM cell","CD4+ TEM cell","CD56+ NK cell","CD8+ T cell","CD8+ TCM cell","CD8+ TEM cell","dendritic cell","dnT cell","gdT cell","hematopoietic precursor cell","LCL","macrophage","MAIT cell","memory B cell","monocyte","neutrophil","NK cell","pancreas","pancreatic islet","plasmablast","plasmacytoid dendritic cell","spleen","T cell","Tfh cell","Th1 cell","Th17 cell","Th2 cell","Treg memory","Treg naive")

datasets_ge2 <- datasets_ge[datasets_ge$tissue_label %in% query.tissue,]
datasets_ex2 <- datasets_ex[datasets_ex$tissue_label %in% query.tissue,]
datasets_tx2 <- datasets_tx[datasets_tx$tissue_label %in% query.tissue,]


dim(datasets_ge2)
dim(datasets_ex2)
dim(datasets_tx2)




ASE.eqtl_ge.list <- list()
ASE.eqtl_tx.list <- list()
ASE.eqtl_ex.list <- list()

for (s in ase.snps) {
print(s)
    ASE.eqtl_ge.list[[s]] <- get_assoc_over_datasets(datasets_ge, s, variant_type = "rsid")
    ASE.eqtl_tx.list[[s]] <- get_assoc_over_datasets(datasets_tx, s, variant_type = "rsid")
    ASE.eqtl_ex.list[[s]] <- get_assoc_over_datasets(datasets_ex, s, variant_type = "rsid")
    }


############








