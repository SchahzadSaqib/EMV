old_path <- Sys.getenv("PATH")
old_path
base::Sys.setenv(PATH = paste(old_path, "~/opt/anaconda3/bin", sep = ":"))

##### Load packages and custom scripts ##-----
source(here::here("Scripts", "load_packages.R"))
devtools::install_github("SchahzadSaqib/taxminer")

## Set NCBI key 
rentrez::set_entrez_key("Place.NCBI.key.here")
base::Sys.getenv("ENTREZ_KEY")

seqtab_10_13 <- base::data.frame(readRDS(here::here("Data", 
                                              "MiSeq_10_13_seqtab_nochim.rds")))
seqtab_20 <- base::data.frame(readRDS(here::here("Data", 
                                           "MiSeq_20_seqtab_nochim.rds")))
seqtab_24 <- base::data.frame(readRDS(here::here("Data", 
                                           "MiSeq_24_seqtab_nochim.rds")))

## Read in metadata
metadata <- utils::read.delim(here::here("Data", 
                                  "meta_all_1st_Rversion.txt"))

## compile and filter sequence table
seqtab_Pa_Ga_324 <- seqtab_10_13 %>%
  dplyr::bind_rows(seqtab_20, seqtab_24) %>%
  dplyr::mutate(across(.cols = everything(), 
                       .fns = ~ifelse(is.na(.x), 0, .x))) %>%
  tibble::rownames_to_column(var = "Sample") %>%
  dplyr::filter(Sample %in% metadata$Sample) %>%
  dplyr::arrange(Sample) %>%
  tibble::column_to_rownames(var = "Sample") %>%
  dplyr::select(where(~any(. > 0))) %>%
  base::as.matrix()


##### Taxonomic annotations ##-----
## write ASV table
base::saveRDS(seqtab_Pa_Ga_324, file = here::here("Data", 
                                            "Pa_Ga_324_seqtab_nochim.rds"))


## Align against nt database, restricted by accession IDs 
## (no environmental samples)

## BLAST based alignment
BlastAnnot_Pa_Ga_324_noUC <- taxminer::txm_align(
  input = seqtab_Pa_Ga_324, 
  database_path = here::here("~",
                             "Documents",
                             "NCBI_databases",
                             "nt"), 
  database_name = "nt", 
  output_path = here::here("Results"), 
  output_name = "Pa_Ga_324_align_noUC", 
  accession_list = "Bacteria_noUC.seq", 
  accession_path = here::here("Data"), 
  do_acc_check = F,
  Run_Blast = F,
  threads = 5, 
  qcvg = 98, 
  pctidt = 98, 
  max_out = 5000) %>%
  
  ## Text-mining based filtration of hits
  taxminer::txm_ecosrc(
    filter_host = "human", 
    filter_site = c("vagina+FRS", "gut+skin+oral+clinical+human"), 
    filter_negate = "non_human", 
    Precomp_tbl = here::here("Data", "Bacteria_eco.rds")
  ) %>%
  
  ## Clean data
  dplyr::group_by(ID) %>%
  dplyr::mutate(TaxID = as.numeric(TaxID)) %>%
  dplyr::filter(bitscore == max(bitscore)) %>%
  dplyr::filter(Pct == max(Pct)) %>%
  dplyr::filter(qcovs == max(qcovs)) %>%
  dplyr::filter(Evalue == min(Evalue)) %>%
  dplyr::distinct(Species, .keep_all = T) %>%
  dplyr::mutate(across(where(is.list), as.character)) %>% 
  dplyr::distinct(ID, .keep_all = T) %>%
  dplyr::ungroup() %>%
  base::as.data.frame()


## Align against nt database, restricted by accession IDs 
## (only environmental samples)

## BLAST based alignment
BlastAnnot_Pa_Ga_324_UC <- taxminer::txm_align(
  input = seqtab_Pa_Ga_324, 
  database_path = here::here("~", 
                             "Documents", 
                             "NCBI_databases", 
                             "nt"), 
  database_name = "nt", 
  output_path = here::here("Results"), 
  output_name = "Pa_Ga_324_align_UC", 
  accession_list = "Bacteria_UC.seq", 
  accession_path = here::here("Data"), 
  do_acc_check = F,
  Run_Blast = F,
  threads = 5, 
  qcvg = 98, 
  pctidt = 98, 
  max_out = 5000
) %>%
  
  ## Text-mining based filtration of hits
  taxminer::txm_ecosrc(
    filter_host = "human", 
    filter_site = c("vagina+FRS", "gut+skin+oral+clinical+human"), 
    filter_negate = "non_human", 
    Precomp_tbl = here::here("Data", "Bacteria_eco.rds")
  ) %>%
  
  ## Clean data
  dplyr::group_by(ID) %>%
  dplyr::mutate(TaxID = as.numeric(TaxID)) %>%
  dplyr::filter(bitscore == max(bitscore)) %>%
  dplyr::filter(Pct == max(Pct)) %>%
  dplyr::filter(qcovs == max(qcovs)) %>%
  dplyr::filter(Evalue == min(Evalue)) %>%
  dplyr::distinct(Species, .keep_all = T) %>%
  dplyr::mutate(across(where(is.list), as.character)) %>% 
  dplyr::distinct(ID, .keep_all = T) %>%
  dplyr::ungroup() %>%
  base::as.data.frame() %>%
  dplyr::filter(!ID %in% BlastAnnot_Pa_Ga_324_noUC$ID)


## Combining annotations and adding lineage
Annot_Pa_Ga_324 <- BlastAnnot_Pa_Ga_324_noUC %>%
  dplyr::bind_rows(BlastAnnot_Pa_Ga_324_UC) %>%
  dplyr::arrange(ID) %>%
  taxminer::txm_lineage(Precomp_tbl = here::here("Data", 
                                                 "Bacteria_lineage.rds"))


## Create taxonomy table
taxonomy <- Annot_Pa_Ga_324 %>%
  dplyr::select(superkingdom, phylum, class, order, family, genus, species) %>%
  dplyr::rename_with(tolower)


## Filter out sequences that were removed
seqtab_Pa_Ga_324_filt <- seqtab_Pa_Ga_324 %>%
  base::as.data.frame() %>%
  base::t() %>%
  base::as.data.frame() %>%
  dplyr::mutate(ids = seq(1:nrow(.))) %>%
  dplyr::filter(ids %in% Annot_Pa_Ga_324$ID) %>%
  dplyr::select(-ids) %>%
  base::t()

## Assign sequences to rownames of taxonomy  
base::rownames(taxonomy) <- base::colnames(seqtab_Pa_Ga_324_filt)
taxonomy <- as.matrix(taxonomy)

## save final ASV and taxonomy tables
base::saveRDS(taxonomy, file = here::here("Data", 
                                    "Pa_Ga_324_taxonomy.rds"))
base::saveRDS(seqtab_Pa_Ga_324_filt, file = here::here("Data", 
                                                 "seqtab_Pa_Ga_324_filt.rds"))


## track reads through the annotation pipeline
track_list <- readxl::read_xlsx(here::here("Results", 
                                           "MiSeq_10_13_readtrack.xlsx")) %>%
  dplyr::bind_rows(readxl::read_xlsx(here::here("Results", 
                                                "MiSeq_20_readtrack.xlsx"))) %>%
  dplyr::bind_rows(readxl::read_xlsx(here::here("Results", 
                                                "MiSeq_24_readtrack.xlsx"))) %>%
  dplyr::filter(samples %in% metadata$Sample)

## Add annotated read percentage to tracking
track_cmpl <- data.frame(seqtab_Pa_Ga_324_filt) %>%
  dplyr::mutate(reads_annotated = rowSums(.)) %>%
  tibble::rownames_to_column("samples") %>%
  dplyr::select(samples, reads_annotated) %>%
  dplyr::inner_join(track_list) %>%
  dplyr::relocate(reads_annotated, .after = last_col()) %>%
  dplyr::mutate(Perc_reads_annotated = reads_annotated/nonchim*100) %>%
  dplyr::arrange(str_detect(samples, "\\+"), desc(Perc_retained_afterfilt))

## write final read track excel file
writexl::write_xlsx(track_cmpl, 
                    path = here::here("Results", 
                                      "Pa_Ga_324_readtrack.xlsx"))
