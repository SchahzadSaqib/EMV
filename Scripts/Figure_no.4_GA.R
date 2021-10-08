##### load packages and custom scripts ##-----
load.pkgs <- list(here::here("Scripts", "load_packages.R"),
                  here::here("Scripts", "plot_theme.R"),
                  here::here("Scripts", "GroupTest_mod.R"))
purrr::map(load.pkgs, source)



##### Read data ##-----
EMMI_phy324 <- readr::read_rds(here::here("Data", "EMMI_phy324.rds"))

## Modify sample data
TL_mod <- data.frame(EMMI_phy324@sam_data) %>%
  dplyr::select(Term_vs_late) %>%
  dplyr::mutate(Term_vs_late = as.character(Term_vs_late)) %>%
  dplyr::mutate(Term_vs_late = if_else(Term_vs_late == "Term", 
                                       Term_vs_late, "Late term")) %>%
  dplyr::mutate(Term_vs_late = factor(Term_vs_late, 
                                      levels = unique(Term_vs_late)))

## update phyloseq object
phyloseq::sample_data(EMMI_phy324)$Term_vs_late <- TL_mod$Term_vs_late

## Convert raw counts to relative abundance
EMMI_phy324.trf <- transform_sample_counts(EMMI_phy324, 
                                           function(OTU) OTU/sum(OTU))
EMMI_phy324.trf <- transform_sample_counts(EMMI_phy324.trf, 
                                           function(OTU) 
                                             ifelse(is.na(OTU), 
                                                    0, 
                                                    OTU))



##### Rearrange data for plotting ##-----
## Extract ASV table and elongate
ASV_table <- data.frame(EMMI_phy324.trf@otu_table) %>%
  tibble::rownames_to_column(var = "ID") %>%  
  tidyr::pivot_longer(cols = !ID, 
                      names_to = "taxa", 
                      values_to = "Abundance") %>%
  as.data.frame()

## Extract and clean meta table
meta_plot <- data.frame(EMMI_phy324.trf@sam_data) %>%
  tibble::rownames_to_column(var = "ID") %>%
  dplyr::mutate(across(.cols = everything(), .fns = ~as.factor(.x)))

## Combine
plot_data <- ASV_table %>%
  dplyr::inner_join(meta_plot)


##### Define plot variables ##-----
tidyvar <- sym("Primipara")
var <- "Primipara"
tidy_var_grp <- sym("Term_vs_late")
var_grp <- "Term_vs_late"
gg_grobs <- list() 

## Define path to figures
path_to_fig <- here::here("Results", "Mscpt_plots")

if (!dir.exists(here::here(path_to_fig))) {
  dir.create(here::here(path_to_fig), recursive = T)
}

plot_colour <- "grey50"
text_size <- 14
plot_theme(text_size = text_size)



##### Bar plot ##-----
## Sample summaries
sample_tally <- plot_data %>%
  dplyr::select(ID, !!tidyvar, !!tidy_var_grp) %>%
  dplyr::distinct(ID, .keep_all = T) %>%
  dplyr::count(!!tidyvar, !!tidy_var_grp) %>%
  dplyr::ungroup()

## Check for sample totals
if (sum(sample_tally$n) != 324) stop(print("Samples do not add up"))

## Clean and rearrange data
p1 <- plot_data %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::group_by(taxa, !!tidyvar, !!tidy_var_grp) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  dplyr::mutate(taxa = if_else(Abundance < 0.002, "Others", taxa)) %>%
  dplyr::arrange(desc(Abundance)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
  base::as.data.frame() 

## define species colours
mycols <- viridis::viridis(5, begin = 0.30, 
                           option = "mako") %>%
  base::append(viridis::viridis(5, option = "rocket", 
                                begin = 0.3, 
                                end = 0.90, 
                                direction = -1)) %>%
  base::append(PNWColors::pnw_palette(n = 5, name = "Sailboat")) %>%
  base::append(PNWColors::pnw_palette(
    n = nlevels(p1$taxa) - 10, 
    name = "Shuksan")[6:(nlevels(p1$taxa) - 10)]) %>%
  purrr::set_names(levels((p1$taxa))) %>%
  base::replace("Others", scico::scico(1, 
                                       palette = "grayC", 
                                       begin = 0.65))  %>%
  base::as.data.frame() %>%
  purrr::set_names("colours") %>%
  tibble::rownames_to_column(var = "taxa") %>%
  dplyr::mutate(lightness = shades::lightness(colours)) %>%
  dplyr::mutate(shades = if_else(lightness > 65, "light", "dark"))

## add colours and shading to dataframe
p1 <- p1 %>%
  dplyr::inner_join(mycols) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(Abundance)) %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
  dplyr::mutate(shades = factor(shades, levels = c("dark", "light"))) %>%
  base::as.data.frame()


## Plot
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
                         y = 1.01, 
                         label = paste("n =", n, sep = " ")), 
                     angle = 0,
                     hjust = "left",
                     size = text_size*0.32,
                     colour = plot_colour,
                     inherit.aes = F) +
  
  ## facet by grouping var
  ggplot2::facet_wrap(as.formula(paste(var_grp, "~.", sep = "")), 
                      scales = "free_x") +
  
  ## flip coordinates
  ggplot2::coord_flip() +
  
  ## Add mean relative abundance label for each taxa
  ggplot2::geom_text(aes(label = if_else(Abundance > 0.01, 
                                         paste(round(Abundance*100, 
                                                     digits = 2), 
                                               "%"), 
                                         ""), 
                         size = Abundance,
                         colour = shades,
                         group = taxa, 
                         angle = if_else(Abundance > 0.10, 0, 270)), 
                     position = position_stack(0.5, reverse = T), 
                     show.legend = F) +
  
  ## assign manual set of colours
  ggplot2::scale_fill_manual(values = mycols$colours, name = "Species") +
  
  ## scale colour for labels
  ggplot2::scale_colour_manual(values = c("grey85", "grey35"), 
                               name = "Species") +
  
  ## convert y-axis scale to percentage
  ggplot2::scale_y_continuous(labels = scales::percent, 
                              limits = c(0,1.1), 
                              expand = c(0,0)) +
  
  ## adjust labels and plot theme
  ggplot2::xlab("") +
  ggplot2::ylab("") +
  ggplot2::theme(legend.position = "bottom", 
                 axis.text.x = element_text(angle = 270)) +
  ggplot2::guides(fill = guide_legend(ncol = 6))


gg_grobs <- base::append(gg_grobs, list(bar_plot))

if (!dir.exists(here::here("Results", "table"))) {
  dir.create(here::here("Results", "tables"), recursive = T)
}

## List variables to be tested in GroupTest
grp <- c("Primipara", "Primipara")
sel <- c("Term_vs_late", "Term_vs_late")
assign_names <- c("Nulliparous", "Nulliparous")
assign_select <- c(levels(plot_data$Term_vs_late))
table_comb <- c()


for (fct in 1:length(assign_names)) {
  print(grp[[fct]])
  print(assign_names[[fct]])
  print(sel[[fct]])
  print(assign_select[[fct]])
  
  tidyvar <- sym(grp[[fct]])
  tidysel <- sym(sel[[fct]])
  var <- grp[[fct]]
  
  ##### Grouptest ##-----
  gt_res <- try({GroupTest_mod(species.table = EMMI_phy324,
                               meta = EMMI_phy324,
                               group = grp[[fct]], 
                               group_name = grp[[fct]],
                               compare.to = assign_names[[fct]],
                               select.by = sel[[fct]],
                               select = assign_select[[fct]],
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
    dplyr::group_by(taxa, !!tidyvar, !!tidysel) %>%
    dplyr::summarise(Abundance = mean(Abundance)) %>%
    dplyr::filter(!!tidysel == assign_select[[fct]])
  
  ## extract FDR values and clean the results
  gt_summ_cmp <- gt_res %>%
    dplyr::select(taxon, model, contains("FDR")) %>%
    tidyr::gather("comp", "pvals", 3:length(.)) %>%
    dplyr::group_by(taxon) %>%
    dplyr::filter(pvals <= 0.05) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(variable = str_extract(comp, 
                                         pattern = paste(
                                           levels(plot_data[,var]),
                                           collapse = "|"))) %>%
    dplyr::group_by(taxon, variable) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sig = case_when((pvals > 0 & pvals <= 0.001) ~ "***",
                                  (pvals > 0.001 & pvals <= 0.01) ~ "**",
                                  (pvals > 0.01 & pvals <= 0.05) ~ "*",
                                  (pvals > 0.05 & pvals <= 0.1) ~ ".")) %>%
    dplyr::rename(!!tidyvar := variable) %>%
    dplyr::rename("taxa" = taxon) %>%
    dplyr::select(-comp) %>%
    dplyr::inner_join(pavg)
  
  ## append to table
  table_comb <- table_comb %>%
    dplyr::bind_rows(gt_summ_cmp) %>%
    dplyr::arrange(taxa)
  
  gt_summ_ref <- pavg %>%
    dplyr::filter(!!tidyvar == assign_names[[fct]]) %>%
    dplyr::filter(!!tidysel == assign_select[[fct]]) %>%
    dplyr::select(taxa, Abundance, !!tidyvar, !!tidysel) %>%
    dplyr::filter(taxa %in% gt_summ_cmp$taxa)
  
  gt_summ <- gt_summ_cmp %>%
    dplyr::bind_rows(gt_summ_ref) %>%
    dplyr::arrange(taxa) %>%
    dplyr::relocate(c(model, pvals, sig), .after = last_col())
  
  writexl::write_xlsx(gt_summ, 
                      path = here::here(
                        "Results", 
                        "tables", 
                        paste("Fig4_Ref_", 
                              grp[[fct]],
                              assign_names[[fct]],
                              "_",
                              assign_select[[fct]],
                              ".xlsx", sep = "")))
}

tidyvar <- sym("Primipara")
var <- "Primipara"
tidy_var_grp <- sym("Term_vs_late")
var_grp <- "Term_vs_late"

## Table 1: late term selection
table_edit_lt <- table_comb %>%
  dplyr::select(-Abundance) %>%
  dplyr::filter(!!tidy_var_grp == "Late term") %>%
  dplyr::select(-!!tidyvar) %>%
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
  tidyr::pivot_longer(!c(Species, !!tidy_var_grp)) %>%
  tidyr::pivot_wider(names_from = Species, values_from = c(value)) %>%
  dplyr::select(-name) %>%
  dplyr::rename("Group" = !!tidy_var_grp)

## Table 2: term
table_edit_tr <- table_comb %>%
  dplyr::select(-Abundance) %>%
  dplyr::filter(!!tidy_var_grp == "Term") %>%
  dplyr::select(-!!tidyvar) %>%
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
  tidyr::pivot_longer(!c(Species, !!tidy_var_grp)) %>%
  tidyr::pivot_wider(names_from = Species, values_from = c(value)) %>%
  dplyr::select(-name) %>%
  dplyr::rename("Group" = !!tidy_var_grp) %>%
  dplyr::bind_rows(table_edit_lt) %>%
  dplyr::mutate(across(.cols = everything(), 
                       .fns = ~ifelse(is.na(.x), " - ", .x)))

## wite excel table 
writexl::write_xlsx(table_edit_tr, 
                    path = here::here("Results", 
                                      "tables", 
                                      "Grouptest_Ga_sig_tab.xlsx"))



## viobox for significant taxa
p2 <- plot_data %>%
  dplyr::select(ID, taxa, Abundance, !!tidyvar, !!tidy_var_grp) %>%
  dplyr::filter(taxa %in% table_comb$taxa) %>%
  dplyr::mutate(!!tidy_var_grp := as.character(!!tidy_var_grp)) %>%
  dplyr::full_join(y = table_comb %>%
                     dplyr::ungroup() %>%
                     dplyr::select(-Abundance) %>%
                     as.data.frame()) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::mutate(!!tidyvar := factor(!!tidyvar, levels = c("Nulliparous",
                                                          "Multiparous"))) %>%
  dplyr::mutate(!!tidy_var_grp := factor(!!tidy_var_grp, 
                                         levels = c("Term", 
                                                    "Late term")))


pal <- PNWColors::pnw_palette("Sunset", 7) %>%
  .[c(1, 3)] %>%
  purrr::set_names(levels(p2$Primipara))

p2_lab <- p2 %>%
  dplyr::filter(!is.na(sig)) %>%
  dplyr::group_by(taxa, !!tidy_var_grp) %>%
  dplyr::distinct(!!tidy_var_grp, .keep_all = T) %>%
  dplyr::mutate(!!tidyvar := "Nulliparous")

## construct viobox
gt_viobox <- ggplot2::ggplot(p2, aes(x = !!tidyvar, 
                                     y = Abundance, 
                                     fill = !!tidyvar)) +
  
  ## Add jitter plot
  ggplot2::geom_jitter(aes(colour = !!tidyvar),
                       fill = "grey70",
                       show.legend = F, 
                       size = 3,
                       alpha = 0.8) +
  
  ## half violin in the left
  gghalves::geom_half_violin(mapping = aes(colour = !!tidyvar),
                             trim = F, 
                             side = "l",
                             draw_quantiles = c(0.25, 0.5, 0.75), 
                             scale = "area",
                             colour = "grey70",
                             size = .2,
                             alpha = 0.8) +
  
  ## half box on the right
  gghalves::geom_half_boxplot(mapping = aes(colour = !!tidyvar),
                              side = "r", 
                              nudge = .05,
                              size = .2, 
                              colour = "grey70",
                              outlier.shape = NA, 
                              alpha = 0.8, 
                              show.legend = F) +
  
  ## facet by taxa
  ggplot2::facet_grid(as.formula(paste( var_grp, "~taxa", sep = "")), 
                      scales = "free_y") +
  
  ggplot2::geom_text(data = p2_lab, 
                     aes(x = !!tidyvar,
                         y = 9,
                         label = sig),
                     colour = "red3",
                     size = 8,
                     inherit.aes = F) +
  
  ## add discrete fill and colour
  ggplot2::scale_fill_manual(values = pal, name = "") +
  ggplot2::scale_colour_manual(values = pal) +
  
  ggnewscale::new_scale_colour() +
  
  ## log10 transform y-axis
  ggplot2::scale_y_log10() + 
  ggplot2::theme(legend.position = "bottom",
                 axis.text.x = element_blank(),
                 legend.key.size = unit(1, "cm")) +
  
  ## add point for median
  ggplot2::stat_summary(fun = "median", 
                        geom = "point", 
                        colour = "red3", 
                        size = 3.5, 
                        show.legend = F) +
  
  ## adjust labels and theme
  ggplot2::xlab("") +
  ggplot2::ylab("Relative abundance (log10)") +
  ggplot2::guides(fill = guide_legend(ncol = 2))

## append object as grob
gg_grobs <- base::append(gg_grobs, list(gt_viobox))


##### Arrange grobs and create illustration ##-----
## define layout design
layout_des <- rbind(c(1,1,1),
                    c(1,1,1),
                    c(2,2,2),
                    c(2,2,2),
                    c(2,2,2))

## arrange grobs
g <- gridExtra::arrangeGrob(grobs = gg_grobs,
                            layout_matrix = layout_des,
                            padding = unit(1.5, "cm"), 
                            
                            ## add margins around the plot
                            right = "\n",
                            left = "\n")

## add panel labels
p <- as_ggplot(g) + 
  cowplot::draw_plot_label(label = c("a", "b"), 
                           x = c(0.01, 0.01), 
                           y = c(0.96, 0.57), 
                           size = 28,
                           colour = "grey50") 

## save illustration
ggplot2::ggsave(path = here::here(path_to_fig),
                plot = p,
                filename = paste("Fig4_GA.pdf", 
                                 sep = ""), 
                device = "pdf", 
                dpi = 320, 
                height = 30, 
                width = 50,
                units = "cm", 
                limitsize = F)

