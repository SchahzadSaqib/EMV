##### load packages and custom scripts ##-----
source(here::here("Scripts", 
                  "load_packages.R"))

##### Read data ##-----
EMMI_phy324 <- readr::read_rds(here::here("Data", 
                                          "EMMI_phy324.rds"))

## Convert raw counts to relative abundance
EMMI_phy324.trf <- transform_sample_counts(EMMI_phy324, 
                                           function(OTU) OTU/sum(OTU)) 
EMMI_phy324.trf <- transform_sample_counts(EMMI_phy324.trf, 
                                           function(OTU) 
                                             ifelse(is.na(OTU), 
                                                    0, 
                                                    OTU))


##### Adonis helpers ##-----
## permutations
permt_to_use <- 99999

## cleaner
adonis_clean <- function(x) {
  x <- x %>%
    tibble::rownames_to_column(var = "Var") %>%
    dplyr::rename_with(~paste("Pr_F_stat"), starts_with("Pr")) %>%
    dplyr::rename_with(~paste("F_statistic"), starts_with("F")) %>%
    dplyr::mutate(
      Sigif = case_when((Pr_F_stat > 0 & Pr_F_stat <= 0.001) ~ "***",
                        (Pr_F_stat > 0.001 & Pr_F_stat <= 0.01) ~ "**",
                        (Pr_F_stat > 0.01 & Pr_F_stat <= 0.05) ~ "*",
                        (Pr_F_stat > 0.05 & Pr_F_stat <= 0.1) ~ "."),
      Sigif = tidyr::replace_na(.data$Sigif, replace = ""),
      Pr_F_stat = as.numeric(.data$Pr_F_stat)) %>%
    dplyr::filter(Pr_F_stat > 0) %>%
    dplyr::arrange(Pr_F_stat)
}


## write excel file 
if (!dir.exists(here::here("Results", "adonis"))) {
  dir.create(here::here("Results", "adonis"), recursive = T)
}

to.test <- c("Pooled")
for (lev in 1:length(to.test)) {
  
  ## asv table
  ASV_tab <- data.frame(EMMI_phy324.trf@otu_table) %>%
    dplyr::filter(sum(.) != 0) 
  
  ## metadata
  meta <- data.frame(EMMI_phy324.trf@sam_data) %>%
    droplevels() %>%
    dplyr::select(where(
      ~all(nlevels(.x) > 1) | is.numeric(.x) | is.character(.x)
    ))
  
  ## Run adonis
  adonis_out <- data.frame()
  
  for (i in 1:length(meta)) {
    set.seed(26)
    print(colnames(meta[i]))
    ASV_tab_ad <- ASV_tab %>%
      filter(rownames(.) %in% rownames(na.omit(meta[i])))
    if (nrow(ASV_tab_ad) == 0) next
    adonis_hold <- adonis2(formula = ASV_tab_ad ~ na.omit(meta[,i]), 
                           data = meta, 
                           method = "bray", 
                           permutations = permt_to_use)
    rownames(adonis_hold)[1] <- colnames(meta[i])
    adonis_out <- rbind(adonis_out, adonis_hold)
  }
  
  ## clean output
  adonis_bac <- adonis_clean(adonis_out) %>%
    as.data.frame()
  
  assign(paste("adonis_", to.test[lev], sep = ""), adonis_bac, .GlobalEnv)
}

## write excel file 
writexl::write_xlsx(adonis_Pooled, 
                    path = here::here("Results", 
                                      "adonis", 
                                      "adonis_pooled.xlsx"))

## save .RData 
save(adonis_Pooled,
     file = here::here("Results", 
                       "adonis", 
                       "adonis_Pooled.RData"))
