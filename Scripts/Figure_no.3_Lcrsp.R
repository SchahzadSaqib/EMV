##### Figure 3: scatter infographic L. crispatus 


##### Load packages and custom scripts ##-----
load.pkgs <- list(here::here("Scripts", "load_packages.R"),
                  here::here("Scripts", "plot_theme.R"))
purrr::map(load.pkgs, source)


##### Read data ##-----
EMMI_phy324 <- readr::read_rds(here::here("Data", "EMMI_phy324.rds"))

## Convert raw counts to relative abundance
EMMI_phy324.trf <- phyloseq::transform_sample_counts(EMMI_phy324, 
                                                     function(OTU) OTU/sum(OTU)) 
EMMI_phy324.trf <- phyloseq::transform_sample_counts(EMMI_phy324.trf, 
                                                     function(OTU) 
                                                       ifelse(is.na(OTU), 
                                                              0, 
                                                              OTU))


## Extract L.crispatus abundance in each sample 
L.crispatus_prc <- data.frame(EMMI_phy324.trf@otu_table) %>%
  dplyr::select("Lactobacillus.crispatus")

## Add update L.crispatus abundance to metadata
phyloseq::sample_data(
  EMMI_phy324.trf)$L.crispatus_prc <- L.crispatus_prc$Lactobacillus.crispatus

## Extract L.iners abundance in each sample 
L.iners_prc <- data.frame(EMMI_phy324.trf@otu_table) %>%
  dplyr::select("Lactobacillus.iners")

## Add update L.crispatus abundance to metadata
phyloseq::sample_data(
  EMMI_phy324.trf)$L.iners_prc <- L.iners_prc$Lactobacillus.iners


## Extract L.crispatus abundance in each sample 
G.vag_prc <- data.frame(EMMI_phy324.trf@otu_table) %>%
  dplyr::select("Gardnerella.vaginalis")

## Add update L.crispatus abundance to metadata
phyloseq::sample_data(
  EMMI_phy324.trf)$G.vag_prc <- G.vag_prc$Gardnerella.vaginalis


##### Rearrange data for plotting ##-----
## Extract ASV table and elongate
ASV_table <- data.frame(EMMI_phy324.trf@otu_table) %>%
  tibble::rownames_to_column(var = "ID") %>%
  tidyr::pivot_longer(cols = !ID, names_to = "taxa", values_to = "Abundance")

## Extract and clean meta table
meta_plot <- data.frame(EMMI_phy324.trf@sam_data) %>%
  tibble::rownames_to_column(var = "ID") %>%
  dplyr::mutate(GP3 = as.character(GP3)) %>%
  dplyr::mutate(GP3 = ifelse(GP3 == "Nulliparous multigravida", 
                             paste("Nulliparous", 
                                   "multigravida", 
                                   sep = "\n"), GP3)) %>%
  dplyr::mutate(across(.cols = everything(), .fns = ~as.factor(.x)))

## Combine
plot_data <- ASV_table %>%
  dplyr::inner_join(meta_plot)



##### Define plot variables ##-----
# Define path to figures
path_to_fig <- here::here("Results", "Mscpt_plots")

## Check or create directory
if (!dir.exists(here::here(path_to_fig))) {
  dir.create(here::here(path_to_fig), recursive = T)
}



p1 <- plot_data %>%
  dplyr::select(ID, 
                taxa, 
                Abundance, 
                Ga_at_sample_weeks, 
                L.crispatus_prc, 
                Primipara, 
                GP3) %>%
  dplyr::mutate(across(starts_with(c("Ga", "L.")), 
                       .fns = ~as.numeric(as.character(.x)))) %>%
  dplyr::filter(taxa == "Lactobacillus.crispatus") %>%
  base::droplevels() %>%
  dplyr::arrange(Ga_at_sample_weeks, 
                 L.crispatus_prc) %>%
  dplyr::mutate(Pooled = "Pooled") %>%
  tidyr::pivot_longer(cols = c(Primipara, GP3, Pooled), 
                      names_to = "sub", 
                      values_to = "split") %>%
  dplyr::filter(if_else(sub == "GP3", !split == "Multiparous", TRUE)) %>%
  dplyr::mutate(split = factor(split, levels = c("Pooled",
                                                 "Nulliparous",
                                                 "Primigravida",
                                                 "Nulliparous\nmultigravida",
                                                 "Multiparous"))) 

gg_grobs <- c()

## Construct model and plot
p <- ggplot2::ggplot(p1, aes(x = Ga_at_sample_weeks, 
                             y = Abundance)) +
  
  ## add loess model
  ggplot2::stat_smooth(aes(colour = split, 
                           fill = split),
                       method = "loess", 
                       formula = y ~ x, 
                       size = 2.5,
                       fill = "grey90",
                       alpha = 0.50,
                       se = T,
                       fullrange = F,
                       show.legend = F) +
  
  ## add colours
  ggplot2::scale_fill_viridis_d(
    option = "F",
    end = 0.8,
    alpha = 0.6,
    name = "") +
  ggplot2::scale_colour_viridis_d(
    option = "F", 
    end = 0.8, 
    alpha = 0.9,
    name = "") +
  
  ## log10 transform y-axis
  ggplot2::scale_y_log10() +

  ## Adjust labels and theme
  ggplot2::ylab("Relative Abundance (log10)") +
  ggplot2::xlab("Gestational age (weeks)") +
  ggplot2::theme(axis.text.x = element_text(angle = 0, 
                                            colour = "grey50"))

## append object to list of grobs
gg_grobs <- append(gg_grobs, list(p))

box <- ggplot2::ggplot(p1, aes(x = split, 
                               y = Abundance,
                               colour = split,
                               fill = split)) +
  
  ## Add jitter plot
  ggplot2::geom_jitter(show.legend = F, 
                       size = 4.5,
                       alpha = 0.6) +
  
  ## half violin in the left
  gghalves::geom_half_violin(trim = F, 
                             side = "l",
                             draw_quantiles = c(0.25, 0.5, 0.75), 
                             scale = "area",
                             colour = "grey70",
                             size = .2,
                             alpha = 0.6,
                             show.legend = T) +
  
  ## half box on the right
  gghalves::geom_half_boxplot(side = "r", 
                              nudge = .05,
                              size = .2, 
                              colour = "grey70",
                              outlier.shape = NA, 
                              alpha = 0.6,
                              show.legend = F) +
  
  ## add discrete fill and colour
  ggplot2::scale_fill_viridis_d(
    option = "F", 
    end = 0.8, 
    alpha = 0.6,
    name = "") +
  ggplot2::scale_colour_viridis_d(
    option = "F", 
    end = 0.8, 
    alpha = 0.6,
    name = "") +
  
  ## add median
  ggplot2::stat_summary(fun = "median", 
                        geom = "point", 
                        colour = "skyblue", 
                        size = 3.5, 
                        show.legend = F) +
  
  ## log10 transform y-axis
  ggplot2::scale_y_log10() +
  
  ## adjust labels and theme
  ggplot2::theme(axis.text.x = element_blank(),
                 legend.position = "bottom",
                 legend.key.size = unit(1, "cm")) +
  ggplot2::xlab("") +
  ggplot2::ylab("Relative Abundance (log10)") +
  ggplot2::guides(fill = guide_legend(ncol = 5))

## append object to list of grobs
gg_grobs <- append(gg_grobs, list(box))


p1 <- plot_data %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::mutate(Pooled = "Pooled") %>%
  tidyr::pivot_longer(cols = c(Primipara, GP3, Pooled), 
                      names_to = "sub", 
                      values_to = "split") %>%
  dplyr::filter(if_else(sub == "GP3", !split == "Multiparous", TRUE)) %>%
  dplyr::mutate(split = factor(split, levels = c("Pooled",
                                                 "Nulliparous",
                                                 "Primigravida",
                                                 "Nulliparous\nmultigravida",
                                                 "Multiparous"))) %>%
  dplyr::group_by(taxa, split) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  dplyr::filter(taxa == "Lactobacillus crispatus")

## construct plot
bar <- ggplot2::ggplot(p1, aes(x = split, 
                               y = Abundance,
                               fill = split)) +
  
  ## add bar plot
  ggplot2::geom_bar(stat = "identity", 
                    position = position_stack(reverse = T),
                    show.legend = F) +
  
  ## add percentage within the bars
  ggplot2::geom_text(aes(label = paste(round(Abundance*100, digits = 2), "%"), 
                         fill = split),
                     size = 6, 
                     colour = "grey50",
                     alpha = 0.8,
                     position = position_stack(1.1, reverse = T), 
                     show.legend = F) +
  
  ## adjust labels and theme
  ggplot2::theme(axis.text.x = element_blank(),
                 axis.text.y = element_blank()) +
  ggplot2::xlab("") +
  ggplot2::ylab("") +
  ggplot2::scale_fill_viridis_d(option = "F", end = 0.8, name = "", alpha = 0.6)

## append object to list of grobs
gg_grobs <- append(gg_grobs, list(bar))


##### Arrange grobs and create illustration ##-----
## define layout design
layout_des <- rbind(c(1,1,3,3),
                    c(2,2,2,2))

## arrange grobs
g <- gridExtra::arrangeGrob(grobs = gg_grobs,
                            layout_matrix = layout_des,
                            padding = unit(0.5, "cm"),
                            top = "\n",
                            left = "\n",
                            right = "\n")

## add panel labels
p <- as_ggplot(g) + 
  cowplot::draw_plot_label(label = c("a", "b", "c"), 
                           x = c(0.02, 0.50, 0.02), 
                           y = c(0.98, 0.98, 0.47), 
                           size = 20,
                           fontface = "bold",
                           colour = "grey50") 

## save illustration
ggplot2::ggsave(path = here::here(path_to_fig),
                plot = p,
                filename = paste("Fig3_Lcrsp.pdf"), 
                device = "pdf", dpi = 320, height = 20, width = 35,
                units = "cm", limitsize = F)
