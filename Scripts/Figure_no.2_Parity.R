##### Figure 2: Parity 


##### load packages and custom scripts ##-----
load.pkgs <- list(here::here("Scripts", "load_packages.R"),
                  here::here("Scripts", "plot_theme.R"),
                  here::here("Scripts", "GroupTest_mod.R"))
purrr::map(load.pkgs, source)

##### Read data ##-----
EMMI_phy324 <- readr::read_rds(here::here("Data", "EMMI_phy324.rds"))

# Convert raw counts to relative abundance 
EMMI_phy324.trf <- phyloseq::transform_sample_counts(EMMI_phy324, 
                                                     function(OTU) OTU/sum(OTU)) 
EMMI_phy324.trf <- phyloseq::transform_sample_counts(EMMI_phy324.trf, 
                                                     function(OTU) 
                                                       ifelse(is.na(OTU), 
                                                              0, 
                                                              OTU))


##### Rearrange data for plotting ##-----
# Extract ASV table and elongate
ASV_table <- data.frame(EMMI_phy324.trf@otu_table) %>%
  tibble::rownames_to_column(var = "ID") %>%
  tidyr::pivot_longer(cols = !ID, 
                      names_to = "taxa", 
                      values_to = "Abundance") %>%
  as.data.frame()

# Extract and clean meta table
meta_plot <- data.frame(EMMI_phy324.trf@sam_data) %>%
  tibble::rownames_to_column(var = "ID") %>%
  dplyr::mutate(GP3 = as.character(GP3)) %>%
  dplyr::mutate(GP3 = ifelse(GP3 == "Nulliparous multigravida", 
                             paste("Nulliparous", 
                                   "multigravida", 
                                   sep = "\n"), GP3)) %>%
  dplyr::mutate(across(.cols = everything(), .fns = ~as.factor(.x)))

# Compile data for plotting
plot_data <- ASV_table %>%
  dplyr::inner_join(meta_plot) %>%
  dplyr::select(ID, 
                taxa, 
                Abundance, 
                Ga_at_sample_weeks, 
                Primipara, 
                GP3, 
                ends_with("_prc")) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::group_by(taxa) %>%
  tidyr::pivot_longer(cols = c(Primipara, GP3), 
                      names_to = "sub", 
                      values_to = "split") %>%
  dplyr::filter(if_else(sub == "GP3", !split == "Multiparous", TRUE)) %>%
  dplyr::ungroup() %>%
  as.data.frame()




##### Define plot variables ##-----
# Path to figures
path_to_fig <- here::here("Results", "Mscpt_plots")

# Global variables
tidyvar <- sym("split")
var <- "split"
gg_grobs <- list() 

# Check directories
if (!dir.exists(here::here(path_to_fig))) {
  dir.create(here::here(path_to_fig), recursive = T)
}

# Plot aesthetics
plot_colour <- "grey50"
text_size <- 14
plot_theme(text_size = text_size)


# ---------------------------------------------------------------------------- #
                           #### Compile plot ####
# ---------------------------------------------------------------------------- #
##### Bar plot ##-----
## Sample summaries
sample_tally <- plot_data %>%
  dplyr::select(ID, !!tidyvar) %>%
  dplyr::group_by(!!tidyvar) %>%
  dplyr::distinct(ID, .keep_all = T) %>%
  dplyr::tally() %>%
  dplyr::ungroup()

## Clean and rearrange data
p1 <- plot_data %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::group_by(taxa, !!tidyvar) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  dplyr::mutate(taxa = if_else(Abundance < 0.002, "Others", taxa)) %>%
  dplyr::inner_join(sample_tally) %>%
  dplyr::arrange(desc(Abundance)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
  as.data.frame() %>%
  dplyr::mutate(split = factor(split, levels = c("Multiparous",
                                                 "Nulliparous\nmultigravida",
                                                 "Primigravida",
                                                 "Nulliparous")))

## Define species colours
mycols <- viridis::viridis(5, begin = 0.30, 
                           option = "mako") %>%
  base::append(viridis::viridis(5, option = "rocket", 
                                begin = 0.3, 
                                end = 0.90, 
                                direction = -1)) %>%
  base::append(PNWColors::pnw_palette(n = 5, name = "Sailboat")) %>%
  base::append(PNWColors::pnw_palette(n = nlevels(p1$taxa) - 10, 
                                      name = "Shuksan")[6:(nlevels(p1$taxa) - 10)]) %>%
  purrr::set_names(levels((p1$taxa))) %>%
  base::replace("Others", scico::scico(1, 
                                       palette = "grayC", 
                                       begin = 0.65)) %>%
  as.data.frame() %>%
  purrr::set_names("colours") %>%
  tibble::rownames_to_column(var = "taxa") %>%
  ## appropriate contrast and colour for text
  dplyr::mutate(lightness = shades::lightness(colours)) %>%
  dplyr::mutate(shades = if_else(lightness > 65, "light", "dark"))


p1 <- p1 %>%
  dplyr::inner_join(mycols) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(Abundance)) %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
  dplyr::mutate(shades = factor(shades, levels = c("dark", "light"))) %>%
  as.data.frame() 


bar_plot <- ggplot2::ggplot(p1, 
                            aes(x = !!tidyvar, 
                                y = Abundance, 
                                fill = taxa)) +
  
  ## add bar plot
  ggplot2::geom_bar(stat = "identity", 
                    position = position_stack(reverse = T)) +
  
  ## Add label for sample tally  
  ggplot2::geom_text(data = sample_tally, 
                     aes(x = !!tidyvar,
                         y = 1.02,
                         label = paste("n =", n, sep = " ")), 
                     angle = 0, 
                     hjust = "left",
                     size = text_size*0.42,
                     colour = plot_colour,
                     inherit.aes = F) +
  
  ## flip coordinates
  ggplot2::coord_flip() +
  
  ## Add mean relative abundance label for each taxa
  ggplot2::geom_text(aes(label = if_else(Abundance > 0.01, 
                                         paste(round(Abundance*100, digits = 2), 
                                               "%"), 
                                         ""), 
                         size = Abundance,
                         colour = shades,
                         group = taxa, 
                         angle = if_else(Abundance > 0.10, 0, 270)), 
                     position = position_stack(0.5, reverse = T), 
                     show.legend = F) +
  
  ## assign manual set of colours
  ggplot2::scale_fill_manual(values = mycols$colours, 
                             name = "Species") +
  
  ## scale colour for labels
  ggplot2::scale_colour_manual(values = c("grey85", "grey35"), 
                               name = "Species") +
  
  ## convert y-axis scale to percentage
  ggplot2::scale_y_continuous(labels = scales::percent, 
                              limits = c(0, 1.1), 
                              expand = c(0,0)) +
  
  ## adjust labels and theme
  ggplot2::xlab("") +
  ggplot2::ylab("") +
  ggplot2::theme(legend.position = "bottom", 
                 axis.text.x = element_text(angle = 270)) +
  ggplot2::guides(fill = guide_legend(ncol = 5))


## append object as grob
gg_grobs <- append(gg_grobs, list(bar_plot))


## assign variables to be tested
grp <- c("Primipara", "GP3", "GP3")
assign_names <- list("Nulliparous", 
                     "Nulliparous multigravida", 
                     "Primigravida")
table_comb <- c()

## creat directory for results
if (!dir.exists(here::here("Results", "table"))) {
  dir.create(here::here("Results", "tables"), recursive = T)
}


for (fct in 1:length(assign_names)) {
  ##### Grouptest ##-----
  gt_res <- try({GroupTest_mod(species.table = EMMI_phy324,
                               meta = EMMI_phy324,
                               group = grp[[fct]], 
                               group_name = grp[[fct]],
                               compare.to = assign_names[[fct]], 
                               dir_for_res = here::here("Results", "mare", var), 
                               min.prevalence = 0.05, 
                               min.abundance = 0.05,
                               outlier.cutoff = 3, 
                               p.cutoff = 0.05, 
                               keep.result = T, 
                               nonzero = F,
                               pdf = T, 
                               show_quartz = F)
  }) 
  
  if (class(gt_res) == "try-error") {
    gt_res <- c()
  }
  
  ## add mean relative abundances
  pavg <- plot_data %>%
    dplyr::group_by(taxa, !!tidyvar) %>%
    dplyr::summarise(Abundance = mean(Abundance))
  
  ## extract FDR values and clean the results
  gt_summ_cmp <- gt_res %>%
    dplyr::select(taxon, model, contains("FDR")) %>%
    tidyr::pivot_longer(cols = contains("FDR"), 
                        names_to = "comp",
                        values_to = "pvals") %>%   
    dplyr::group_by(taxon) %>%
    dplyr::filter(pvals <= 0.05) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(comp = str_replace(comp, " ", "\n")) %>%
    dplyr::mutate(variable = str_extract(comp, 
                                         pattern = paste(
                                           unique(plot_data[,var]), 
                                           collapse = "|"))) %>%
    dplyr::group_by(taxon, variable) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sig = case_when((pvals > 0 & pvals <= 0.001) ~ "***",
                                  (pvals > 0.001 & pvals <= 0.01) ~ "**",
                                  (pvals > 0.01 & pvals <= 0.05) ~ "*",
                                  (pvals > 0.05 & pvals <= 0.1) ~ ".")) %>%
    dplyr::rename(!!tidyvar := variable) %>%
    dplyr::rename("taxa" = taxon) %>%
    dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
    dplyr::select(-comp) %>%
    dplyr::inner_join(pavg) %>%
    dplyr::mutate(Group = assign_names[[fct]])
  
  ## append to table
  table_comb <- table_comb %>%
    dplyr::bind_rows(gt_summ_cmp) %>%
    dplyr::arrange(taxa)
  
  ## aundances for reference group
  gt_summ_ref <- pavg %>%
    dplyr::filter(!!tidyvar == assign_names[[fct]]) %>%
    dplyr::select(taxa, Abundance, !!tidyvar) %>%
    dplyr::filter(taxa %in% gt_summ_cmp$taxa)
  
  ## combine tables to write
  gt_summ <- gt_summ_cmp %>%
    dplyr::bind_rows(gt_summ_ref) %>%
    dplyr::arrange(taxa) %>%
    dplyr::relocate(c(model, pvals, sig), .after = last_col())
  
  ## write excel of final results
  writexl::write_xlsx(gt_summ, 
                      path = here::here(
                        "Results", 
                        "tables", 
                        paste("Fig2_Ref_", 
                              assign_names[[fct]], 
                              ".xlsx", sep = "")))
}


## Table 1: NM vs multi
table_edit_NnP_Mult <- table_comb %>% 
  dplyr::rename("comp" = all_of(var)) %>%
  dplyr::relocate(comp, .before = Group) %>%
  dplyr::filter(Group == "Nulliparous multigravida") %>%
  dplyr::mutate(Group = str_replace(Group, " ", "\n")) %>%
  dplyr::filter(comp == "Multiparous") %>%
  dplyr::select(-c(comp, Abundance)) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::mutate(pvals = signif(pvals, digits = 4)) %>%
  tidyr::unite("FDR (q.value)", c(pvals,sig), sep = "") %>%
  dplyr::relocate("FDR (q.value)", .before = model) %>%
  dplyr::mutate(model = str_replace(model, "::", " - ")) %>%
  dplyr::mutate(stats = paste("FDR (q.value): ", 
                              `FDR (q.value)`, 
                              "\n",
                              "Model: ",
                              model, sep = "")) %>%
  dplyr::select(-c(`FDR (q.value)`, model)) %>%
  dplyr::rename("Species" = taxa) %>%
  tidyr::pivot_longer(!c(Species, Group)) %>%
  tidyr::pivot_wider(names_from = Species, values_from = c(value)) %>%
  dplyr::select(-name)

## Table 2: P vs multi
table_edit_NP_Mult <- table_comb %>% 
  dplyr::rename("comp" = all_of(var)) %>%
  dplyr::relocate(comp, .before = Group) %>%
  dplyr::filter(Group == "Primigravida") %>%
  dplyr::mutate(Group = str_replace(Group, " ", "\n")) %>%
  dplyr::filter(comp == "Multiparous") %>%
  dplyr::select(-c(comp, Abundance)) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::mutate(pvals = signif(pvals, digits = 4)) %>%
  tidyr::unite("FDR (q.value)", c(pvals,sig), sep = "") %>%
  dplyr::relocate("FDR (q.value)", .before = model) %>%
  dplyr::mutate(model = str_replace(model, "::", " - ")) %>%
  dplyr::mutate(stats = paste("FDR (q.value): ", 
                              `FDR (q.value)`, 
                              "\n",
                              "Model: ",
                              model, sep = "")) %>%
  dplyr::select(-c(`FDR (q.value)`, model)) %>%
  dplyr::rename("Species" = taxa) %>%
  tidyr::pivot_longer(!c(Species, Group)) %>%
  tidyr::pivot_wider(names_from = Species, values_from = c(value)) %>%
  dplyr::select(-name)

## Table 3: nulli vs multi
table_edit_NL_Mult <- table_comb %>% 
  dplyr::rename("comp" = all_of(var)) %>%
  dplyr::relocate(comp, .before = Group) %>%
  dplyr::filter(Group == "Nulliparous") %>%
  dplyr::mutate(Group = str_replace(Group, " ", "\n")) %>%
  dplyr::filter(comp == "Multiparous") %>%
  dplyr::select(-c(comp, Abundance)) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::mutate(pvals = signif(pvals, digits = 4)) %>%
  tidyr::unite("FDR (q.value)", c(pvals,sig), sep = "") %>%
  dplyr::relocate("FDR (q.value)", .before = model) %>%
  dplyr::mutate(model = str_replace(model, "::", " - ")) %>%
  dplyr::mutate(stats = paste("FDR (q.value): ", 
                              `FDR (q.value)`, 
                              "\n",
                              "Model: ",
                              model, sep = "")) %>%
  dplyr::select(-c(`FDR (q.value)`, model)) %>%
  dplyr::rename("Species" = taxa) %>%
  tidyr::pivot_longer(!c(Species, Group)) %>%
  tidyr::pivot_wider(names_from = Species, values_from = c(value)) %>%
  dplyr::select(-name) %>%
  dplyr::bind_rows(table_edit_NP_Mult, table_edit_NnP_Mult) %>%
  dplyr::mutate(across(.cols = everything(), .fns = ~ifelse(is.na(.x), 
                                                            " - ", 
                                                            .x)))



## write excel table 
writexl::write_xlsx(table_edit_NL_Mult, 
                    path = here::here("Results", 
                                      "tables", 
                                      "Grouptest_Parity_sig_tab.xlsx"))




##### Violin + boxplots for significant taxa ##-----
## Clean and rearrange data
p2 <- plot_data %>%
  dplyr::select(ID, taxa, Abundance, !!tidyvar) %>%
  dplyr::filter(taxa %in% table_comb$taxa) %>%
  dplyr::mutate(split = as.character(split)) %>%
  dplyr::full_join(y = table_comb %>%
                     dplyr::ungroup() %>%
                     dplyr::select(-Abundance) %>%
                     dplyr::mutate(split = Group) %>%
                     dplyr::mutate(split = str_replace(split, " ", "\n")) %>%
                     as.data.frame()) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::mutate(split = factor(split, levels = c("Nulliparous",
                                                 "Primigravida",
                                                 "Nulliparous\nmultigravida",
                                                 "Multiparous")))
## assign colours for plot
pal <- PNWColors::pnw_palette("Sunset", 7) %>%
  .[c(1, 3, 5, 7)] %>%
  purrr::set_names(levels(p2$split))

## significance labels
p2_lab <- p2 %>%
  dplyr::filter(!is.na(sig)) %>%
  dplyr::group_by(taxa, split) %>%
  dplyr::distinct(split, .keep_all = T)


## construct viobox
gt_viobox <- ggplot2::ggplot(data = p2 %>%
                               dplyr::filter(Abundance > 0), 
                             aes(x = !!tidyvar,
                                 y = Abundance,
                                 fill = !!tidyvar)) +
  
  ## Add jitter plot
  ggplot2::geom_jitter(aes(colour = !!tidyvar),
                       show.legend = F, 
                       size = 3,
                       alpha = 0.8) +
  
  ## half violin on the left
  gghalves::geom_half_violin(mapping = aes(colour = !!tidyvar),
                             trim = F, 
                             side = "l",
                             draw_quantiles = c(0.25, 0.5, 0.75), 
                             scale = "area",
                             colour = "grey70",
                             size = .3,
                             alpha = 0.8) +
  
  ## half box on the right
  gghalves::geom_half_boxplot(mapping = aes(colour = !!tidyvar),
                              side = "r", 
                              nudge = .05,
                              size = .3, 
                              colour = "grey70",
                              outlier.shape = NA, 
                              alpha = 0.8,
                              show.legend = F) +
  
  ## facet by taxa
  ggplot2::facet_grid(.~taxa, scales = "free_y") +
  
  ## add asterisks labels for significant taxa
  ggplot2::geom_text(data = p2_lab, 
                     aes(x = !!tidyvar,
                         y = 4,
                         label = sig),
                     colour = "red3",
                     size = 8,
                     inherit.aes = F) +
  
  ## add colours
  ggplot2::scale_fill_manual(values = pal, name = "") +
  ggplot2::scale_colour_manual(values = pal)  +
  
  ## start new scale
  ggnewscale::new_scale_colour() +
  
  ## log10 transform y-axis
  ggplot2::scale_y_log10() +

  ## add point for median
  ggplot2::stat_summary(fun = "median", 
                        geom = "point", 
                        colour = "red3", 
                        size = 3.5, 
                        show.legend = F) +
  
  ## adjust labels and theme
  ggplot2::xlab("") +
  ggplot2::ylab("Relative abundance (log10)") +
  ggplot2::theme(legend.position = "bottom",
                 axis.text.x = element_blank(),
                 legend.key.size = unit(1, "cm")) +
  ggplot2::guides(fill = guide_legend(ncol = 4))


## append object as grob
gg_grobs <- base::append(gg_grobs, list(gt_viobox))


##### Arrange grobs and create illustration -----
## define layout design
layout_des <- rbind(c(1,1,1),
                    c(1,1,1),
                    c(2,2,2),
                    c(2,2,2))

## arrange grobs
g <- gridExtra::arrangeGrob(grobs = gg_grobs,
                            layout_matrix = layout_des,
                            padding = unit(2.5, "cm"),
                            
                            ## add margins around the plot
                            left = "\n",
                            right = "\n")

## add panel labels
p <- as_ggplot(g) + 
  cowplot::draw_plot_label(label = c("a", "b"), 
                           x = c(0.02, 0.02), 
                           y = c(0.98, 0.45), 
                           size = 25,
                           fontface = "bold",
                           colour = plot_colour) 

## save illustration
ggplot2::ggsave(path = here::here(path_to_fig),
                plot = p,
                filename = paste("Fig2_Parity.pdf", 
                                 sep = ""), 
                device = "pdf", 
                dpi = 320, 
                height = 25, 
                width = 43,
                units = "cm", 
                limitsize = F)
