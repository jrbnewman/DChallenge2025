###SEE: https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tutorials/API_v2/eQTL_API_tutorial.md


# If you do not already have all packages installed you can use this syntax to install them:
#install.packages(c("tidyverse", "httr", "jsonlite", "dplyr", "coloc", "ggrepel", "glue"))


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
	

## Functon to request datasets meeting criteria

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
	
	request_datasets_from_api(quant_method = "ge", tissue_label = "T cell", size=1000, start=0)


### to get all datasets: datasets_ge = request_datasets_from_api(quant_method = "ge")
	
	
## See dataset metadata (useful for querying)
	# Change parameters
	set_dataset_id = 'QTD000105'
	  
	URL = glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/{set_dataset_id}")
	
	r <- GET(utils::URLencode(URL), accept_json())
	cont <- content(r, "text", encoding = "UTF-8")
	dataset_metadata <- fromJSON(cont)
	dataset_metadata
	
	
### Function to extract eQTLs based on position

request_associations_around_position <- function(dataset_id,
				                 position,
						 chromosome_id,
						 gene_id = NULL,
						 offset = 500000,
						 start = 0,
						 size = 1000) {

range_start = position - offset
  range_end = position + offset
  
  while (TRUE){
    if(!is.null(gene_id) == TRUE) {
    URL = glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/{dataset_id}/associations?size={size}&start={start}&pos={chromosome_id}:{range_start}-{range_end}&gene_id={gene_id}")
    } else {
    URL = glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/{dataset_id}/associations?size={size}&start={start}&pos={chromosome_id}:{range_start}-{range_end}")
    }
    message(URL)

    r <- GET(utils::URLencode(URL), accept_json())
    cont <- content(r, "text", encoding = "UTF-8")
    
    if (status_code(r) != 200) {
      # Loop will break if the request was unsuccessful
      if(start==0) {
        print(glue("Error {status_code(r)}"))
        print(cont)
        return()}
      break
    }
    
   
    cont_df <- fromJSON(cont)
    
    if (start == 0){
      responses <- cont_df
    }
    else{
      responses <- rbind(responses, cont_df)
    }
    start <- start + size
  }
  return(responses)
}


sirpg <- request_associations_around_position(dataset_id = "QTD000105",
				                 position = 1635560 ,
						 chromosome_id = 20)

## Function to extract associations based on rsID


request_associations_from_api <- function(
    dataset_id, 
    pos="",
    variant="", 
    rsid="",
    molecular_trait_id="",
    gene_id="",
    nlog10p=""){
  
  size = 1000
  start = 0
  
  parameter_values = c(dataset_id,pos,variant,rsid,molecular_trait_id, 
                       gene_id,nlog10p)
  parameter_names = c('dataset_id','pos','variant','rsid','molecular_trait_id', 
                       'gene_id','nlog10p')
  
  while (T) {
    URL = glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/{dataset_id}/associations?size={size}&start={start}")
    
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

sirpg2 <- request_associations_from_api(dataset_id = "QTD000105",  gene_id="ENSG00000089012")

sirpg3<- request_associations_from_api(dataset_id = "QTD000105", rsid="rs6043409")


#### Get associations across all gene expression eQTL datasets

get_assoc_over_datasets <- function(datasets, variant, variant_type = c("coordinate","rsid"), gene_id) {
  size = 1000
  first = T
  final_df <- list()
  for (i in rownames(datasets)) {
    row = datasets[i, ]
    dataset_id = row$dataset_id
    
        if (variant_type == "coordinate") {
    URL = glue(
          "https://www.ebi.ac.uk/eqtl/api/v2/datasets/{dataset_id}/associations?size={size}&variant={variant}&gene_id={gene_id}"
	  )
	  } else if (variant_type == "rsid") {
    URL = glue(
          "https://www.ebi.ac.uk/eqtl/api/v2/datasets/{dataset_id}/associations?size={size}&rsid={variant}&gene_id={gene_id}"
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





variant = "chr1_109274570_A_G"
gene_id = "ENSG00000134243"

datasets_ge = request_datasets_from_api(quant_method = "ge")

assoc.coord <- get_assoc_over_datasets(datasets=datasets_ge,
                                       variant = variant,
				       variant_type = "coordinate",
				       gene_id = gene_id)


assoc.rsid <- get_assoc_over_datasets(datasets=datasets_ge,
                                       variant = rsid,
				       variant_type = "rsid",
				       gene_id = gene_id)
				       
# visualize

assoc.coord %>%
   unite(study, study_label, tissue_label, sep = " ", remove = FALSE) %>%
   ggplot(aes(x=beta, y=nlog10p, label=study)) +
   geom_point() +
   geom_text_repel() +
   ggtitle("SNP effect on SORT1 expression")
   