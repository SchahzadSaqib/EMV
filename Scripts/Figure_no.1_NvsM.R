##### Figure 1: Polar plot Nulliparous vs Multiparous


##### load packages and custom scripts ##-----
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
  dplyr::select("Lactobacillus.crispatus") %>%
  base::as.matrix()

## Add L.crispatus abundance to metadata
phyloseq::sample_data(EMMI_phy324.trf)$L.crispatus_prc <- L.crispatus_prc

## Extract L.iners abundance in each sample 
L.iners_prc <- data.frame(EMMI_phy324.trf@otu_table) %>%
  dplyr::select("Lactobacillus.iners") %>%
  base::as.matrix()

## Add L.iners abundance to metadata
phyloseq::sample_data(EMMI_phy324.trf)$L.iners_prc <- L.iners_prc

## Extract G.vaginalis abundance in each sample 
G.vag_prc <- data.frame(EMMI_phy324.trf@otu_table) %>%
  dplyr::select("Gardnerella.vaginalis") %>%
  base::as.matrix()

## Add G.vaginalis abundance to metadata
phyloseq::sample_data(EMMI_phy324.trf)$G.vag_prc <- G.vag_prc





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

## Compile data for plotting
plot_data <- ASV_table %>%
  dplyr::inner_join(meta_plot)



##### define plot variables ##-----
## Define path to figures
path_to_fig <- here::here("Results", "Mscpt_plots")

## Path to figures
if (!dir.exists(here::here(path_to_fig))) {
  dir.create(here::here(path_to_fig), recursive = T)
}

# ---------------------------------------------------------------------------- #
                 #### Polar plot - Nulliparous vs Multiparous ####
# ---------------------------------------------------------------------------- #
## gather plotting data 
p1 <- plot_data %>%
  dplyr::select(ID, 
                taxa, 
                Abundance, 
                Ga_at_sample_weeks, 
                Primipara, 
                GP3, 
                ends_with("_prc")) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::mutate(taxa = ifelse(Abundance < 0.05, "Others", taxa)) %>%
  dplyr::arrange(desc(Abundance)) %>%
  dplyr::mutate(across(starts_with(c("Ga")) | ends_with("_prc"), 
                       .fns = ~as.numeric(as.character(.x)))) %>%
  dplyr::group_by(taxa) %>%
  dplyr::add_tally() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(n = ifelse(taxa == "Others", 0, n)) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
  dplyr::mutate(value = 1.07) %>%
  dplyr::mutate(value2 = 1.14) %>%
  dplyr::mutate(Primipara = factor(Primipara, levels = unique(Primipara))) %>%
  dplyr::mutate(ID = factor(ID, levels = unique(ID)))


## define species colours
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
                                       begin = 0.65))

## colours for gradient outer band and segments
pal <- PNWColors::pnw_palette(name = "Moth", n = 100)
pal2 <- PNWColors::pnw_palette(name = "Sunset2", n = 7) %>%
  .[c(1, 6)]

## Set gaps between groups
empty_bar <- 10
nbreaks <- nlevels(p1$Primipara)
nObs <- nlevels(factor(plot_data$taxa))
to_add <- data.frame(matrix(NA, empty_bar*nbreaks*nObs, ncol(p1)) )
colnames(to_add) <- colnames(p1)
to_add$Primipara <- rep(levels(p1$Primipara), 
                        each = empty_bar*nObs)

## Define order of bars
p2 <- p1 %>%
  dplyr::bind_rows(to_add) %>%
  dplyr::mutate(taxa_abd = case_when((L.crispatus_prc >= 0.50) ~ 4,
                                     (L.iners_prc >= 0.50) ~ 3,
                                     (G.vag_prc >= 0.50) ~ 2,
                                     TRUE ~ 1)) %>%
  dplyr::arrange(Primipara, 
                 Ga_at_sample_weeks, 
                 taxa_abd, 
                 L.crispatus_prc, 
                 ID) %>%
  dplyr::mutate(labels = rep(seq(1, nlevels(ID) + empty_bar*nbreaks), 
                             each = nObs)) %>%
  dplyr::mutate(Primipara = as.character(Primipara)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ID = as.character(ID)) %>%
  dplyr::mutate(labels = factor(labels, 
                                levels = unique(labels))) %>%
  dplyr::mutate(Abundance = ifelse(is.na(Abundance), 
                                   0, 
                                   Abundance))

## Define label positions around polar plot
label_data <- p2 %>% 
  dplyr::group_by(labels, ID, Ga_at_sample_weeks) %>% 
  dplyr::summarise(tot = sum(Abundance)) %>%
  dplyr::mutate(labels = as.numeric(as.character(labels)))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$labels - 0.5) / number_of_bar
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle + 180, angle)

## prepare a data frame for base lines
base_data <- p2 %>% 
  dplyr::group_by(Primipara) %>% 
  dplyr::mutate(labels = as.numeric(as.character(labels))) %>%
  dplyr::summarise(start = min(labels), 
                   end = max(labels) - empty_bar) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(title = mean(c(start, end)))

## prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[c(nrow(grid_data), 
                                 1:nrow(grid_data) - 1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

## Define plot theme
plot_theme(text_size = 14)

## Make the plot
p <- ggplot2::ggplot(p2) + 
  
  ## control variables for polar plot
  ggplot2::coord_polar() +
  ggplot2::ylim(-.50,1.33) + ## control inner circle radius
  
  ## add outer band bar for gestational age
  ggplot2::geom_col(aes(x = labels, 
                        y = value, 
                        fill = Ga_at_sample_weeks),
                    width = 1,
                    position = "dodge", 
                    alpha = 0.7) +
  
  ## gradient fill manual
  ggplot2::scale_fill_gradientn(colours = pal, 
                                name = "Gestational age (weeks)") +
  
  ## start new scale
  ggnewscale::new_scale_fill() +
  
  ## add stacked bar for taxa
  ggplot2::geom_bar(aes(x = labels, 
                        y = Abundance, 
                        fill = taxa), 
                    stat = "identity",
                    width = 1,
                    position = position_stack(reverse = T)) +
  ggplot2::scale_fill_manual(values = mycols, 
                             name = "Species") +
  
  ## Add text for incremental y-axis tick labels
  ggplot2::annotate("text", x = max(as.numeric(as.character(p2$labels))) + 3, 
                    y = c(0, .25, .50, .75, 1.0), 
                    label = c("0", "25", "50", "75", "100") , 
                    color = "grey60", 
                    size = 7, 
                    angle = 0, 
                    fontface = "bold", 
                    hjust = 1) +

  ## Add samples names on top of bars
  ggplot2::geom_text(data = label_data, 
                     aes(x = labels, 
                         y = 1.10, 
                         label = ID, 
                         hjust = hjust), 
                     color = "grey60", 
                     size = 2.5, 
                     angle = label_data$angle, 
                     inherit.aes = FALSE) +
  
  ## start new colour scale
  ggnewscale::new_scale_colour() +
  
  ## Add segment at base of the plot
  ggplot2::geom_segment(data = base_data, 
                        aes(x = start, 
                            y = -0.03, 
                            xend = end, 
                            yend = -0.03,
                            colour = Primipara), 
                        alpha = 0.8, 
                        size = 0.8, 
                        inherit.aes = FALSE,
                        show.legend = F)  +
  
  ## label segments
  ggplot2::geom_text(data = base_data, 
                     aes(x = title, 
                         y = -0.45, 
                         label = Primipara,
                         colour = Primipara), 
                     hjust = c(-0.01, 1.01), 
                     alpha = 0.8, 
                     size = 6.5, 
                     fontface = "bold", 
                     inherit.aes = FALSE, 
                     show.legend = F) +
  
  ## Add colour
  ggplot2::scale_colour_manual(values = pal2) +
  
  ## adjust labels and theme
  ggplot2::theme(axis.text.x = element_blank(),
                 axis.text.y = element_blank(), 
                 legend.position = c(1.0, 0.5),
                 plot.margin = unit(c(-3, -3, -5, -12), "cm")) +
  ggplot2::xlab("") +
  ggplot2::ylab("") +
  ggplot2::guides(fill = guide_legend(ncol = 1, order = 1))
  

## Save as pdf
ggsave(p, 
       path = path_to_fig, 
       file = "Fig1_polar_NvsM.pdf", 
       device = "pdf", 
       dpi = 320, 
       width = 45, 
       height = 37, 
       units = "cm")