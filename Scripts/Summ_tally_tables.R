load.pkgs <- list(here::here("Scripts", "load_packages.R"),
                  here::here("Scripts", "plot_theme.R")) %>%
  purrr::map(source)

## Create directory for tables 
if (!dir.exists(here::here("Results", "tables"))) {
  dir.create(here::here("Results", "tables"))
}

##### Read data ##----------------------------------------------------------------------------------------------
EMMI_phy324 <- readr::read_rds(here::here("Data", "EMMI_phy324.rds"))

## Convert raw counts to relative abundance
EMMI_phy324.trf <- transform_sample_counts(EMMI_phy324, function(OTU) OTU/sum(OTU)) 
EMMI_phy324.trf <- transform_sample_counts(EMMI_phy324.trf, function(OTU) ifelse(is.na(OTU), 0, OTU))


## Extract L.crispatus abundance in each sample 
L.crispatus_prc <- data.frame(EMMI_phy324.trf@otu_table) %>%
  dplyr::select("Lactobacillus.crispatus")

## Add update L.crispatus abundance to metadata
sample_data(EMMI_phy324.trf)$L.crispatus_prc <- L.crispatus_prc$Lactobacillus.crispatus



## Extract L.iners abundance in each sample 
L.iners_prc <- data.frame(EMMI_phy324.trf@otu_table) %>%
  dplyr::select("Lactobacillus.iners")

## Add update L.crispatus abundance to metadata
sample_data(EMMI_phy324.trf)$L.iners_prc <- L.iners_prc$Lactobacillus.iners


## Extract L.crispatus abundance in each sample 
G.vag_prc <- data.frame(EMMI_phy324.trf@otu_table) %>%
  dplyr::select("Gardnerella.vaginalis")

## Add update L.crispatus abundance to metadata
sample_data(EMMI_phy324.trf)$G.vag_prc <- G.vag_prc$Gardnerella.vaginalis





##### Rearrange data for plotting ##-------------------------------------------------------------------------
## Extract ASV table and elongate
ASV_table <- data.frame(EMMI_phy324.trf@otu_table) %>%
  tibble::rownames_to_column(var = "ID") %>%
  tidyr::pivot_longer(cols = !ID, names_to = "taxa", values_to = "Abundance")

## Extract meta table
meta_plot <- data.frame(EMMI_phy324.trf@sam_data) %>%
  tibble::rownames_to_column(var = "ID") %>%
  dplyr::mutate(GP3 := as.character(GP3)) %>%
  dplyr::mutate(GP3 := case_when((GP3 == "Nulliparous non-primigravida") ~ paste("Nulliparous", 
                                                                                 "non-primigravida", 
                                                                                 sep = "\n"),
                                 (GP3 == "Nulliparous primigravida") ~ paste("Nulliparous", 
                                                                             "primigravida", 
                                                                             sep = "\n"),
                                 TRUE ~ GP3)) %>%
  dplyr::mutate(across(.cols = everything(), .fns = ~as.factor(.x)))

## Combine
plot_data <- ASV_table %>%
  dplyr::inner_join(meta_plot)


taxa_table <- data.frame(EMMI_phy324@tax_table) %>%
  dplyr::mutate("Total reads" = phyloseq::taxa_sums(EMMI_phy324)) %>%
  dplyr::arrange(desc(`Total reads`))

writexl::write_xlsx(taxa_table, 
                    path = here::here("Results", 
                                      "tables",
                                      "taxa.table.xlsx"))

summ_genus_count <- data.frame(EMMI_phy324@tax_table) %>%
  dplyr::filter(!is.na(genus)) %>%
  dplyr::distinct(genus)

summ_genus_count_NA <- data.frame(EMMI_phy324@tax_table) %>%
  dplyr::filter(is.na(genus))


## genus tally 5% (presence) ##
summ_genus_5pct <- plot_data %>%
  dplyr::select(ID, taxa, Abundance, Ga_at_sample_weeks, Primipara, GP3, ends_with("_prc")) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  tidyr::separate(col = taxa, 
                  into = c("genus", "species"),
                  sep = " ") %>%
  dplyr::filter(Abundance > 0.05) %>%
  dplyr::group_by(ID, genus) %>%
  dplyr::summarise(Abundance = sum(Abundance)) %>%
  dplyr::group_by(genus) %>%
  dplyr::add_tally(name = "genus_tally") %>%
  dplyr::mutate(genus_tally_prc = paste(round(genus_tally/324*100, digits = 1),
                                        "%", 
                                        sep = "")) %>%
  dplyr::distinct(genus, .keep_all = T) %>%
  dplyr::select(genus, contains("tally")) %>%
  dplyr::arrange(desc(genus_tally))

# write table 
writexl::write_xlsx(summ_genus_5pct, 
                    path = here::here("Results", 
                                      "tables", 
                                      "Tally_genus_5.xlsx"))

## genus tally 50% (dominance) ##
summ_genus_50pct <- plot_data %>%
  dplyr::select(ID, taxa, Abundance, Ga_at_sample_weeks, Primipara, GP3, ends_with("_prc")) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  tidyr::separate(col = taxa, 
                  into = c("genus", "species"),
                  sep = " ") %>%
  dplyr::filter(Abundance > 0.50) %>%
  dplyr::group_by(ID, genus) %>%
  dplyr::summarise(Abundance = sum(Abundance)) %>%
  dplyr::group_by(genus) %>%
  dplyr::add_tally(name = "genus_tally") %>%
  dplyr::mutate(genus_tally_prc = paste(round(genus_tally/324*100, digits = 1), 
                                        "%", 
                                        sep = "")) %>%
  dplyr::distinct(genus, .keep_all = T) %>%
  dplyr::select(genus, contains("tally")) %>%
  dplyr::arrange(desc(genus_tally))

# write table 
writexl::write_xlsx(summ_genus_50pct, 
                    path = here::here("Results", 
                                      "tables", 
                                      "Tally_genus_50.xlsx"))


summ_data <- plot_data %>%
  dplyr::select(ID, taxa, Abundance, Ga_at_sample_weeks, Primipara, GP3, ends_with("_prc")) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::group_by(taxa) %>%
  dplyr::mutate(Pooled = "Pooled") %>%
  tidyr::pivot_longer(cols = c(Primipara, GP3, Pooled), 
                      names_to = "sub", 
                      values_to = "split") %>%
  dplyr::filter(if_else(sub == "GP3", !split == "Multiparous", TRUE)) %>%
  dplyr::ungroup()


sample_tally <- summ_data %>%
  dplyr::select(ID, split) %>%
  dplyr::group_by(split) %>%
  dplyr::distinct(ID, .keep_all = T) %>%
  dplyr::tally(name = "sample_tally") %>%
  dplyr::ungroup()

## species tally 5% (presence) ##
summ_species_5pct <- summ_data %>%
  dplyr::filter(Abundance > 0.05) %>%
  dplyr::group_by(taxa, split) %>%
  dplyr::add_tally(name = "split_tally") %>%
  dplyr::distinct(taxa, .keep_all = T) %>%
  dplyr::select(taxa, ends_with("tally")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(split_tally)) %>%
  dplyr::mutate(num_fct = as.numeric(factor(taxa, levels = unique(taxa)))) %>%
  dplyr::arrange(num_fct) %>%
  dplyr::select(-num_fct) %>%
  dplyr::inner_join(sample_tally) %>%
  dplyr::mutate(split_prc = paste(round(split_tally/sample_tally*100, digits = 1),
                                  "%",
                                  sep = ""))

# write table 
writexl::write_xlsx(summ_species_5pct, 
                    path = here::here("Results", 
                                      "tables", 
                                      "Tally_species_5.xlsx"))

## species tally >50% (dominance) ##
summ_species_50pct <- summ_data %>%
  dplyr::filter(Abundance > 0.50) %>%
  dplyr::group_by(taxa, split) %>%
  dplyr::add_tally(name = "split_tally") %>%
  dplyr::distinct(taxa, .keep_all = T) %>%
  dplyr::select(taxa, ends_with("tally")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(split_tally)) %>%
  dplyr::mutate(num_fct = as.numeric(factor(taxa, levels = unique(taxa)))) %>%
  dplyr::arrange(num_fct) %>%
  dplyr::select(-num_fct) %>%
  dplyr::inner_join(sample_tally) %>%
  dplyr::mutate(split_prc = paste(round(split_tally/sample_tally*100, digits = 1),
                                  "%",
                                  sep = ""))


# write table 
writexl::write_xlsx(summ_species_50pct, 
                    path = here::here("Results", 
                                      "tables", 
                                      "Tally_species_50.xlsx"))                
