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

## Add update L.crispatus abundance to metadata
phyloseq::sample_data(EMMI_phy324.trf)$L.crispatus_prc <- L.crispatus_prc



## Extract L.iners abundance in each sample 
L.iners_prc <- data.frame(EMMI_phy324.trf@otu_table) %>%
  dplyr::select("Lactobacillus.iners") %>%
  base::as.matrix()

## Add update L.crispatus abundance to metadata
phyloseq::sample_data(EMMI_phy324.trf)$L.iners_prc <- L.iners_prc


## Extract L.crispatus abundance in each sample 
G.vag_prc <- data.frame(EMMI_phy324.trf@otu_table) %>%
  dplyr::select("Gardnerella.vaginalis") %>%
  base::as.matrix()

## Add update L.crispatus abundance to metadata
phyloseq::sample_data(EMMI_phy324.trf)$G.vag_prc <- G.vag_prc





##### Rearrange data for plotting ##-----
## Extract ASV table and elongate
ASV_table <- data.frame(EMMI_phy324.trf@otu_table) %>%
  tibble::rownames_to_column(var = "ID") %>%
  tidyr::pivot_longer(cols = !ID, names_to = "taxa", values_to = "Abundance")

## Extract meta table
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



##### define plot variables ##-----
## Define path to figures
path_to_fig <- here::here("Results", "Mscpt_plots")

## Check or create directory
if (!dir.exists(here::here(path_to_fig))) {
  dir.create(here::here(path_to_fig), recursive = T)
}

plot_theme(text_size = 14)



# ---------------------------------------------------------------------------- #
                    #### Supplementary Figure 1: GP3 ####
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
  dplyr::mutate(ID = factor(ID, levels = unique(ID))) %>%
  dplyr::mutate(value = 1.07) %>%
  dplyr::mutate(value2 = 1.14)


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
                                       begin = 0.65))


pal <- PNWColors::pnw_palette(name = "Moth", n = 100)
pal2 <- PNWColors::pnw_palette(name = "Sunset2", n = 7) %>%
  .[c(1, 4, 6)]  

## Set a number of 'empty bar' to add at the end of each group
empty_bar <- 10
nbreaks <- nlevels(p1$GP3)
nObs <- nlevels(factor(plot_data$taxa))
to_add <- data.frame(matrix(NA, empty_bar*nbreaks*nObs, ncol(p1)) )
colnames(to_add) <- colnames(p1)
to_add$GP3 <- rep(levels(p1$GP3), each = empty_bar*nObs)

## Order of bars 
p2 <- p1 %>%
  dplyr::bind_rows(to_add) %>%
  dplyr::mutate(taxa_abd = case_when((L.crispatus_prc >= 0.50) ~ 4,
                                     (L.iners_prc >= 0.50) ~ 3,
                                     (G.vag_prc >= 0.50) ~ 2,
                                     TRUE ~ 1)) %>%
  dplyr::arrange(GP3, Ga_at_sample_weeks, taxa_abd, L.crispatus_prc, ID) %>%
  dplyr::mutate(labels = rep(seq(1, nlevels(ID) + empty_bar*nbreaks), 
                             each = nObs)) %>%
  dplyr::mutate(GP3 = as.character(GP3)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ID = as.character(ID)) %>%
  dplyr::mutate(labels = factor(labels, levels = unique(labels))) %>%
  dplyr::mutate(Abundance = ifelse(is.na(Abundance), 0, Abundance))

## Get the name and the y position of each label
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
  dplyr::group_by(GP3) %>% 
  dplyr::mutate(labels = as.numeric(as.character(labels))) %>%
  dplyr::summarise(start = min(labels), end = max(labels) - empty_bar) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(title = mean(c(start, end))) %>%
  dplyr::mutate(GP3 = ifelse(GP3 == "Nulliparous multigravida", 
                             paste("Nulliparous", 
                                   "multigravida", 
                                   sep = "\n"), GP3))

## prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data) - 1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]


## Make the plot
p <- ggplot(p2) +      
  
  ## add bar for gestational age
  geom_col(aes(x = labels, 
               y = value, 
               fill = Ga_at_sample_weeks),
           position = "dodge",
           width = 1,
           alpha = 0.7) +
  
  
  scale_fill_gradientn(colours = pal, 
                       name = "Gestational age (weeks)") +
  
  ## start new scale
  ggnewscale::new_scale_fill() +
  
  ## add stacked bar for taxa
  geom_bar(aes(x = labels, 
               y = Abundance, 
               fill = taxa), 
           stat = "identity",
           width = 1) +
  scale_fill_manual(values = mycols, name = "Species") +
  
  
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(1.0, 0.5),
        plot.margin = unit(c(-3, -3, -5, -12), "cm")) +
  xlab("") +
  ylab("") +
  
  ## Add text showing the value of each 100/75/50/25 lines
  ggplot2::annotate("text", x = max(as.numeric(as.character(p2$labels))) + 3, 
                    y = c(0, .25, .50, .75, 1.0), 
                    label = c("0", "25", "50", "75", "100") , 
                    color = "grey60", 
                    size = 7, 
                    angle = 0, 
                    fontface = "bold", 
                    hjust = 1) +
  
  guides(fill = guide_legend(ncol = 1, order = 1)) +
  
  ylim(-.50,1.33) + # control inner circle radius
  coord_polar() +
  scale_x_discrete() +
  
  ## start new colour scale
  ggnewscale::new_scale_colour() +
  
  ## Add labels on top of each bar
  geom_text(data = label_data, 
            aes(x = labels, 
                y = 1.10, 
                label = ID, 
                hjust = hjust), 
            color = "grey60", 
            size = 2.5, 
            angle = label_data$angle, 
            inherit.aes = FALSE) +
  
  ## Add base line information
  geom_segment(data = base_data, 
               aes(x = start, 
                   y = -0.03, 
                   xend = end, 
                   yend = -0.03,
                   colour = GP3), 
               alpha = 0.8, 
               size = 0.6, 
               inherit.aes = FALSE,
               show.legend = F)  +
  geom_text(data = base_data, 
            aes(x = title, 
                y = -0.45, 
                label = GP3,
                colour = GP3), 
            hjust = c(-0.2,0.4,1.2),
            vjust = c(0,2.5,0.9),            
            alpha = 0.8, 
            size = 5, 
            fontface = "bold", 
            inherit.aes = FALSE,
            show.legend = F) +
  ggplot2::scale_colour_manual(values = pal2)

## Save as pdf
ggsave(p, 
       path = path_to_fig, 
       file = "Supplementary_Fig_S1.pdf", 
       device = "pdf", 
       dpi = 320, 
       width = 45, 
       height = 35, 
       units = "cm")





# ---------------------------------------------------------------------------- #
            #### Supplementary Figure 2: Previous delivery method ####
# ---------------------------------------------------------------------------- #
sample_tally <- plot_data %>%
  dplyr::select(ID, Previous_delivery_method, Parity_0or1) %>%
  dplyr::filter(Parity_0or1 == "Yes") %>%
  dplyr::group_by(Previous_delivery_method) %>%
  dplyr::distinct(ID, .keep_all = T) %>%
  dplyr::tally() %>%
  dplyr::ungroup()

p1 <- plot_data %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::filter(Parity_0or1 == "Yes") %>%
  dplyr::group_by(taxa, Previous_delivery_method) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  dplyr::mutate(taxa = if_else(Abundance < 0.002, "Others", taxa)) %>%
  dplyr::inner_join(sample_tally) %>%
  dplyr::arrange(desc(Abundance)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
  as.data.frame() 

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
                                       begin = 0.65)) %>%
  as.data.frame() %>%
  purrr::set_names("colours") %>%
  tibble::rownames_to_column(var = "taxa") %>%
  dplyr::mutate(lightness = shades::lightness(colours)) %>%
  dplyr::mutate(shades = if_else(lightness > 65, "light", "dark"))


p1 <- p1 %>%
  inner_join(mycols) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(Abundance)) %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
  dplyr::mutate(shades = factor(shades, levels = c("dark", "light"))) %>%
  dplyr::mutate(Previous_delivery_method = forcats::fct_relevel(
    Previous_delivery_method, "Nulliparous", "Vaginal")) %>%
  as.data.frame() 

## Define plot theme
text_size <- 10
plot_colour <- "grey50"
plot_theme(text_size = text_size)


bar_plot <- ggplot(p1, aes(x = Previous_delivery_method, 
                           y = Abundance, 
                           fill = taxa)) +
  geom_bar(stat = "identity", position = position_stack(reverse = T)) +
  
  ## Add label for sample tally  
  geom_text(data = sample_tally, 
            aes(x = Previous_delivery_method,
                y = 1.02,
                label = paste("n =", n, sep = " ")), 
            angle = 0, 
            hjust = "left",
            size = text_size*0.42,
            colour = plot_colour,
            inherit.aes = F) +
  
  
  ## flip coordinates
  coord_flip() +
  
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
  scale_fill_manual(values = mycols$colours, name = "Species") +
  
  ## scale colour for labels
  ggplot2::scale_colour_manual(values = c("grey85", "grey35"), 
                               name = "Species") +
  
  ## convert y-axis scale to percentage
  scale_y_continuous(labels = scales::percent, 
                     limits = c(0, 1.1), 
                     expand = c(0,0)) +
  ggplot2::xlab("") +
  ggplot2::ylab("") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 270)) +
  guides(fill = guide_legend(ncol = 8))

ggsave(bar_plot, 
       path = path_to_fig, 
       filename = paste("Supplementary_Fig_S2.pdf", 
                        sep = ""), 
       device = "pdf", 
       dpi = 320, 
       height = 15, 
       width = 47,
       units = "cm")




# ---------------------------------------------------------------------------- #
              #### Supplementary Figure 3: Parity (3 class) ####
# ---------------------------------------------------------------------------- #

## gather plotting data 
p1 <- plot_data %>%
  dplyr::select(ID, 
                taxa, 
                Abundance, 
                Parity_3class,
                ends_with("_prc")) %>%
  dplyr::mutate(Parity_3class = forcats::fct_recode(Parity_3class, 
                                                    "2+" = "2 or more")) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::mutate(taxa = ifelse(Abundance < 0.05, "Others", taxa)) %>%
  dplyr::arrange(desc(Abundance)) %>%
  dplyr::mutate(across(ends_with("_prc"), 
                       .fns = ~as.numeric(as.character(.x)))) %>%
  dplyr::group_by(taxa) %>%
  dplyr::add_tally() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(n = ifelse(taxa == "Others", 0, n)) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
  dplyr::mutate(Parity_3class = factor(Parity_3class, 
                                       levels = unique(Parity_3class))) %>%
  dplyr::mutate(ID = factor(ID, levels = unique(ID))) %>%
  dplyr::mutate(value = 1.07) %>%
  dplyr::mutate(value2 = 1.14)

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
                                       begin = 0.65))

pal <- PNWColors::pnw_palette(name = "Moth", n = 100)
pal2 <- PNWColors::pnw_palette(name = "Sunset2", n = 7) %>%
  .[c(1, 3, 6)]


## Set a number of 'empty bar' to add at the end of each group
empty_bar <- 10
nbreaks <- nlevels(p1$Parity_3class)
nObs <- nlevels(factor(plot_data$taxa))
to_add <- data.frame(matrix(NA, empty_bar*nbreaks*nObs, ncol(p1)) )
colnames(to_add) <- colnames(p1)
to_add$Parity_3class <- rep(levels(p1$Parity_3class), 
                            each = empty_bar*nObs)

## Order of bars 
p2 <- p1 %>%
  dplyr::bind_rows(to_add) %>%
  dplyr::mutate(taxa_abd = case_when((L.crispatus_prc >= 0.50) ~ 4,
                                     (L.iners_prc >= 0.50) ~ 3,
                                     (G.vag_prc >= 0.50) ~ 2,
                                     TRUE ~ 1)) %>%
  dplyr::arrange(Parity_3class, L.crispatus_prc, taxa_abd, ID) %>%
  dplyr::mutate(labels = rep(seq(1, nlevels(ID) + empty_bar*nbreaks), 
                             each = nObs)) %>%
  dplyr::mutate(Parity_3class = as.character(Parity_3class)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ID = as.character(ID)) %>%
  dplyr::mutate(labels = factor(labels, levels = unique(labels))) %>%
  dplyr::mutate(Abundance = ifelse(is.na(Abundance), 0, Abundance))

## Get the name and the y position of each label
label_data <- p2 %>% 
  dplyr::group_by(labels, ID) %>% 
  dplyr::summarise(tot = sum(Abundance)) %>%
  dplyr::mutate(labels = as.numeric(as.character(labels)))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$labels - 0.5) / number_of_bar
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle + 180, angle)

## prepare a data frame for base lines
base_data <- p2 %>% 
  dplyr::group_by(Parity_3class) %>% 
  dplyr::mutate(labels = as.numeric(as.character(labels))) %>%
  dplyr::summarise(start = min(labels), 
                   end = max(labels) - empty_bar) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(title = mean(c(start, end)))

## prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c(nrow(grid_data), 
                                  1:nrow(grid_data) - 1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

## Define plot theme
plot_theme(text_size = 14)

## Make the plot
p <- ggplot(p2) + 
  
  ## add stacked bar for taxa
  geom_bar(aes(x = labels, 
               y = Abundance, 
               fill = taxa), 
           stat = "identity",
           width = 1,
           position = position_stack(reverse = T)) +
  scale_fill_manual(values = mycols, name = "Species") +
  
  
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(1.0, 0.5),
        plot.margin = unit(c(-3, -3, -5, -12), "cm")) +
  xlab("") +
  ylab("") +
  
  ## Add text showing the value of each 100/75/50/25 lines
  ggplot2::annotate("text", x = max(as.numeric(as.character(p2$labels))) + 3, 
                    y = c(0, .25, .50, .75, 1.0), 
                    label = c("0", "25", "50", "75", "100") , 
                    color = "grey60", 
                    size = 7, 
                    angle = 0, 
                    fontface = "bold", 
                    hjust = 1) +
  
  guides(fill = guide_legend(ncol = 1, order = 1)) +
  
  ylim(-.50,1.33) + ## control inner circle radius
  coord_polar() +
  scale_x_discrete() +
  
  ## Add labels on top of each bar
  geom_text(data = label_data, 
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
  geom_segment(data = base_data, 
               aes(x = start, 
                   y = -0.03, 
                   xend = end, 
                   yend = -0.03,
                   colour = Parity_3class), 
               alpha = 0.8, 
               size = 0.8, 
               inherit.aes = FALSE,
               show.legend = F)  +
  
  ## label segments
  geom_text(data = base_data, 
            aes(x = title, 
                y = -0.45, 
                label = Parity_3class,
                colour = Parity_3class),
            position = position_stack(vjust = 0.7),
            alpha = 0.8, 
            size = 8, 
            fontface = "bold", 
            inherit.aes = FALSE, 
            show.legend = F) +
  ggplot2::scale_colour_manual(values = pal2)

## Save as pdf
ggsave(p, 
       path = path_to_fig, 
       file = "Supplementary_Fig_S3.pdf", 
       device = "pdf", 
       dpi = 320, 
       width = 45, 
       height = 37, 
       units = "cm")
