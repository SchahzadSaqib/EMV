##### Load packages and custom scripts ##-----
source(here::here("Scripts", "load_packages.R"))

##### Read in and clean data ##-----
## metadata
metadata <- readr::read_delim(
  here::here("Data",
             "meta_all_1st_Rversion.txt")) %>%
  dplyr::rename("Primipara" = Primipara...46) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(var = "Sample") %>%
  dplyr::mutate(GP3 = case_when(
    (GP3 == "Multipara" ~ "Multiparous"),
    (GP3 == "Nulliparous primigravida" ~ "Primigravida"),
    (GP3 == "Nulliparous non-primigravida" ~ "Nulliparous multigravida"))) %>%
  dplyr::mutate(Primipara = if_else(Primipara == "Yes", 
                                    "Nulliparous", 
                                    "Multiparous")) %>%
  
  ## convert all numeric columns to factor
  dplyr::mutate(across(where(is.numeric), .fns = as.factor)) %>%
  
  ## convert back to numeric if more than 10 levels
  ## factor to numeric needs an intermediate step to convert to character
  dplyr::mutate(across(where(
    ~all(nlevels(.x) > 10)), .fns = ~as.numeric(as.character(.x)))) %>%
  
  ## convert all character columns to factor
  dplyr::mutate(across(where(is.character), .fns = as.factor)) %>%
  
  ## convert back to character if more than 10 levels
  dplyr::mutate(across(where(
    ~all(nlevels(.x) > 10)), .fns = as.character)) %>%
  dplyr::rename_with(~paste("B.Breve"), starts_with("Bacteria_Actinobacteria"))


## load ASVs
seqtab_sub <- readr::read_rds(
  here::here("Data",
             "seqtab_Pa_Ga_324_filt.rds")
) %>%
  base::as.data.frame() %>%
  tibble::rownames_to_column(var = "SampleID") %>%
  dplyr::mutate(SampleID = str_remove(SampleID, pattern = "_.*")) %>%
  dplyr::filter(SampleID %in% rownames(metadata)) %>%
  tibble::column_to_rownames(var = "SampleID")

## taxonomy 
taxonomy_sub <- readr::read_rds(
  here::here("Data",
             "Pa_Ga_324_taxonomy.rds")
) %>%
  base::as.data.frame() %>%
  dplyr::mutate(species = paste(species, " ", sep = "")) %>%
  dplyr::mutate(species = ifelse(str_ends(species, ".1"), 
                                 "uncultured bacteria", 
                                 species)) %>%
  dplyr::mutate(species = str_extract(species, "^(?:[^ ]+ ){2}")) %>%
  dplyr::mutate(species = str_replace_all(species, " ", ".")) %>%
  dplyr::mutate(species = str_replace(species, "\\.", " ")) %>%
  dplyr::mutate(species = str_remove_all(species, "\\.")) %>%
  base::as.matrix()

## Calculate richness and diversity
metadata <- metadata %>%
  dplyr::mutate(
    Richness = vegan::specnumber(seqtab_sub)
  ) %>%
  dplyr::mutate(
    Diversity = vegan::diversity(seqtab_sub, index = "simpson")
  )


##### Compile phyloseq object ##-----
## Phyloseq object 
EMMI_phy324 <- phyloseq::phyloseq(otu_table(seqtab_sub, taxa_are_rows = F), 
                                  sample_data(metadata),
                                  tax_table(taxonomy_sub)) %>%
  
  ## agglomerate ASVs based on species
  phyloseq::tax_glom(taxrank = "species")


## Constructing a phylogenetic tree
#seqs <- phyloseq::taxa_names(EMMI_phy324)
#names(seqs) <- phyloseq::taxa_names(EMMI_phy324)
#alignment <- DECIPHER::AlignSeqs(DNAStringSet(seqs), anchor = NA, verbose = T)
#phangAlign <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
#dm <- phangorn::dist.ml(phangAlign)
#treeNJ <- phangorn::NJ(dm)
#fit = phangorn::pml(treeNJ, data = phangAlign)
#fitGTR <- stats::update(fit, k = 4, inv = 0.2)
#fitGTR <- phangorn::optim.pml(fitGTR, 
#                              model = "GTR", 
#                              optInv = T, 
#                              optGamma = T, 
#                              rearrangement = "stochastic",
#                              control = phangorn::pml.control(trace = 0))
#
#saveRDS(fitGTR, file = here::here("Data", "Phy_tree.rds"))
fitGTR <- readr::read_rds(here::here("Data", "Phy_tree.rds"))


## Add phy tree to phyloseq object
EMMI_phy324 <- phyloseq::merge_phyloseq(EMMI_phy324, phy_tree(fitGTR$tree)) %>%
  
  ## filter out taxa that are present in less than 1% of the data
  phyloseq::filter_taxa(function(x){sum(x > 0) > (0.01*length(x))}, 
                        prune = TRUE) %>%
  
  ## filter out taxa that have less than 500 reads
  phyloseq::filter_taxa(function(x){sum(x) > 500}, prune = TRUE)


## Change ASV sequences to species names
phyloseq::taxa_names(EMMI_phy324) <- as.character(
  EMMI_phy324@tax_table@.Data[,7]
  )


## Adjust reads
reads <- data.frame(EMMI_phy324@otu_table) %>%
  dplyr::summarise(rowSums(.)) %>%
  purrr::set_names("ReadCount")

## Save read counts
phyloseq::sample_data(EMMI_phy324)$ReadCount <- reads$ReadCount

## Filter samples below 500 reads
EMMI_phy324 <- phyloseq::subset_samples(EMMI_phy324, ReadCount > 500)

## Save phyloseq object
base::saveRDS(EMMI_phy324, file = here::here("Data", "EMMI_phy324.rds"))

