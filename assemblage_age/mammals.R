# Sat Jan 28 14:22:57 2023 ------------------------------
# new metric 'assemblage age' from herodotools package
# mammals

# libraries
library(daee)
library(BioGeoBEARS)
library(Herodotools)
library(tidyverse)
library(phytools)
library(picante)

# limits of map
coastline <- rnaturalearth::ne_coastline(returnclass = "sf")
map_limits <- list(
  x = c(-58, -35),
  y = c(-30, -5)
)

# load data
coords_id <- read.table("Shapefiles/id_coordinates.txt")
nrow(coords_id) # 432 sites
mammals_full_list <- read.table('02_metrics/full_data/list_mammals.txt', header=T)
mammals_full_phy <- read.tree("02_metrics/full_data/mammals.tre")


# adjustments in mammals phy
mammals.upham.drop <- drop.tip(mammals_full_phy, "_Anolis_carolinensis") #remove Anolis
sp_names_full <- mammals.upham.drop$tip.label
split_names <- sapply(strsplit(sp_names_full, split='_'), function(x) (c(x[1], x[2])))
split_names_t <- data.frame(t(split_names))
new_names <- paste(split_names_t$X1, split_names_t$X2, sep="_")

# check if it is correct
data.frame(original=mammals.upham.drop$tip.label[-1][1:10], fixed=new_names[1:10])
mammals.upham.drop$tip.label <- new_names
mammals.upham.drop #4175 tips

# adjustments in mammals list
#Trocando . por _
names(mammals_full_list) <- gsub(".", "_", fixed=T, names(mammals_full_list))
#Colocando em ordem alfabética
mammals_full_list <- mammals_full_list[ , order(names(mammals_full_list))]

#Mudando na composição
names(mammals_full_list) <- gsub("Lycalopex", "Pseudalopex", fixed=T, names(mammals_full_list)) 

# ordem alfabética
mammals_list <- mammals_full_list[ , order(names(mammals_full_list))]
head(mammals_list)

ncol(mammals_list) # 236 spp
mammals.upham.drop  # 196 spp

# prune with phylogeny
mammals_phy <- prune.sample(mammals_list, mammals.upham.drop) # 196 spp.

match.phylo.comm(mammals_phy, mammals_list)
ncol( mammals_list) # precisamos remover 40 espécies que estão a mais na composição

rem.col.phy <- c("Abrawayaomys_chebezi","Abrawayaomys_ruschii","Akodon_sanctipaulensis","Brucepattersonius_griserufescens","Brucepattersonius_guarani","Brucepattersonius_misionensis","Brucepattersonius_paradisus","Cabassous_tatouay","Callicebus_barbarabrownae","Callicebus_melanochir","Callithrix_flaviceps","Dasyprocta_iacki","Dasyprocta_prymnolopha","Dasypus_hybridus","Dasypus_septemcinctus","Graomys_chacoensis","Hylaeamys_oniscus","ID","Leontopithecus_caissara","Leontopithecus_chrysopygus","Leopardus_geoffroyi","Leopardus_guttulus","Marmosa_paraguayana","Monodelphis_unistriata","Nectomys_rattus","Oxymycterus_angularis","Oxymycterus_caparoae","Oxymycterus_hispidus","Oxymycterus_roberti","Phaenomys_ferrugineus","Phyllomys_kerri","Phyllomys_medius","Phyllomys_thomasi","Phyllomys_unicolor","Reithrodon_typicus","Sapajus_cay","Sylvilagus_tapetillus","Tolypeutes_tricinctus","Trinomys_mirapitanga","Wilfredomys_oenax")

#Retirando na lista
mammals_list_phy <- mammals_list[,!(names(mammals_list)%in% rem.col.phy)]
ncol(mammals_list_phy) #196
mammals_phy # 196 spp.


# squamata richness
richness_mammals <- rowSums(mammals_list_phy)

# just coordinates
coords <- coords_id[,-3]
names(coords)

# map squamata richness
mammals_richness <- 
  dplyr::bind_cols(coords, richness = richness_mammals) %>% 
  ggplot2::ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = richness_mammals)) + 
  rcartocolor::scale_fill_carto_c(name = "Richness", type = "quantitative", palette = "SunsetDark") +
  ggplot2::geom_sf(data = coastline) +
  ggplot2::coord_sf(xlim = map_limits$x, ylim = map_limits$y) +
  ggplot2::ggtitle("") +
  ggplot2::xlab("Longitude") +
  ggplot2::ylab("Latitude") +
  ggplot2::labs(fill = "Richness") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "mm"),
    legend.text = element_text(size = 6), 
    axis.text = element_text(size = 5),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)
  )
mammals_richness

# construindo evoregions
mammals_phy <- force.ultrametric(mammals_phy)

regions_mammals <- 
  Herodotools::calc_evoregions(
    comm = mammals_list_phy,
    phy = mammals_phy
  )

regions_mammals$cluster_evoregions
site_region_mammals <- regions_mammals$cluster_evoregions

# visualize evoregions
evoregion_mammals <- data.frame(
  coords, 
  site_region_mammals
)

mammals_evoregion <- terra::rast(evoregion_mammals)
mammals_evoregion

# Converting evoregion to a spatial polygon data frame, so it can be plotted
mammals_evoregion <- terra::as.polygons(mammals_evoregion) %>% 
  sf::st_as_sf()

# Assigning the same projection to both spatial objects
sf::st_crs(mammals_evoregion) <- sf::st_crs(coastline)

# Colours to plot evoregions
col_five_hues <- c(
           "#3d291a",
           "#a9344f",
           "#578a5b",
           "#83a6c4",
           "#504162",
           "#FF9B35")
           
map_evoregion_mammals <- 
  evoregion_mammals %>% 
  ggplot2::ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = site_region_mammals)) + 
  ggplot2::scale_fill_manual(
    name = "", 
    labels = LETTERS[1:6],
    values = rev(col_five_hues)
  ) +
  ggplot2::geom_sf(data = coastline) +
  ggplot2::geom_sf(
    data = mammals_evoregion, 
    color = "#040400",
    fill = NA, 
    size = 0.2) +
  ggplot2::coord_sf(xlim = map_limits$x, ylim = map_limits$y) +
  ggplot2::ggtitle("") + 
  ggplot2::theme_bw() +
  ggplot2::xlab("Longitude") +
  ggplot2::ylab("Latitude") +
  ggplot2::theme(
    legend.position = "bottom",
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.text = element_text(color = "black", size = 8)
  )
map_evoregion_mammals

# pertencimento a cada evoregion
# Selecting only axis with more than 5% of explained variance from evoregion output
axis_sel <- which(regions_mammals$PCPS$prop_explainded >= regions_mammals$PCPS$tresh_dist)
PCPS_thresh <- regions_mammals$PCPS$vectors[, axis_sel] 

# distance matrix using 4 significant PCPS axis accordingly to the 5% threshold 
dist_phylo_PCPS <- vegan::vegdist(PCPS_thresh, method = "euclidean")

# calculating affiliation values for each assemblage 
afi <- calc_affiliation_evoreg(phylo.comp.dist = dist_phylo_PCPS,
                               groups = site_region_mammals) 

# binding the information in a data frame
sites <- dplyr::bind_cols(coords, site_region =  site_region_mammals, afi)

map_joint_evoregion_mammals <- 
  evoregion_mammals %>% 
  ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = site_region_mammals), 
                       alpha = sites[, "afilliation"]) + 
  ggplot2::scale_fill_manual(
    name = "", 
    labels = LETTERS[1:6],
    values = rev(col_five_hues)
  ) +
  ggplot2::geom_sf(data = coastline, size = 0.4) +
  ggplot2::geom_sf(
    data = mammals_evoregion, 
    color = rev(col_five_hues),
    fill = NA, 
    size = 0.7) +
  ggplot2::coord_sf(xlim = map_limits$x, ylim = map_limits$y) +
  ggplot2::ggtitle("") + 
  guides(guide_legend(direction = "vertical")) +
  ggplot2::theme_bw() +
  ggplot2::xlab("Longitude") +
  ggplot2::ylab("Latitude") +
  ggplot2::theme(
    legend.position = "bottom",
    axis.title = element_text(size = 10)
  )

map_joint_evoregion_mammals

# temos que definir a ocorrência de cada espécie nas evoregiões. Para fazer isso, podemos usar a função get_region_occe obter um quadro de dados de espécies nas linhas e evoregiões nas colunas.
mammals_region <- get_region_occ_v2(comm = mammals_list_phy, site.region = site_region_mammals)
nrow(mammals_region) # 160
ncol(mammals_list_phy) # 196 - 36 spp a mais 

# O objeto criado na última etapa pode ser usado em uma função auxiliar no Herodotools para produzir facilmente o arquivo Phyllip necessário para executar a análise da reconstrução da área ancestral usando o BioGeoBEARS.
# save phyllip file
getwd()
Herodotools::get_tipranges_to_BioGeoBEARS(mammals_region,filename = "assemblage_age/geo_area_mammals_full.data",areanames = NULL) 

###------------------------------------------------------------------------
# depois da reconstrução de range ancestral 
# Aqui começamos a rodar assemblage_age
# tree file for 
mammals_phy_bio <- force.ultrametric(mammals_phy)
mammals_phy_bio <- ape::multi2di(mammals_phy_bio)

remove <- c("Akodon_montensis", "Bibimys_labiosus", "Callithrix_penicillata", "Cavia_aperea", "Cerdocyon_thous", "Chrysocyon_brachyurus", "Coendou_prehensilis", "Coendou_spinosus", "Cuniculus_paca", "Dasypus_novemcinctus", "Didelphis_albiventris", "Euphractus_sexcinctus", "Euryzygomatomys_spinosus", "Galictis_cuja", "Herpailurus_yagouaroundi", "Hydrochoerus_hydrochaeris", "Kannabateomys_amblyonyx", "Lontra_longicaudis", "Mazama_gouazoubira", "Monodelphis_dimidiata", "Monodelphis_domestica", "Myrmecophaga_tridactyla", "Nasua_nasua", "Necromys_lasiurus", "Nectomys_squamipes","Oligoryzomys_eliurus", "Oligoryzomys_nigripes", "Pecari_tajacu", "Procyon_cancrivorus", "Puma_concolor", "Sapajus_nigritus", "Sciurus_aestuans", "Speothos_venaticus", "Tamandua_tetradactyla", "Tapirus_terrestris", "Tayassu_pecari")

mammals_phy_bio <- drop.tip(mammals_phy_bio, remove)
ape::write.tree(mammals_phy_bio , 'assemblage_age/mammals_biogeo_full.new')

# idade da comunidade
# converting numbers to character
biogeo_area <- data.frame(biogeo = chartr("123456", "ABCDEF", evoregion_mammals$site_region_mammals)) 

# getting the ancestral range area for each node 
node_area <- 
  Herodotools::get_node_range_BioGeoBEARS(
    resDECj,
    phyllip.file = "assemblage_age/geo_area_mammals_full.data",
    mammals_phy_bio,
    max.range.size = 3 
  )

# remove species in list of species
ncol(mammals_list_phy) #196 spp
mammals_list_biogeo <- mammals_list_phy[,!(names(mammals_list_phy)%in% remove)]
ncol(mammals_list_biogeo) #160 spp

# calculating age arrival 
age_comm <- Herodotools::calc_age_arrival(W = mammals_list_biogeo, 
                                          tree = mammals_phy_bio, 
                                          ancestral.area = node_area, 
                                          biogeo = biogeo_area) 

sites <- dplyr::bind_cols(coords, site_region =  evoregion_mammals$site_region_mammals, age = age_comm$mean_age_per_assemblage)

min(age_comm$mean_age_per_assemblage)
max(age_comm$mean_age_per_assemblage)

map_age_mammals <- sites %>% 
  ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = mean_age_arrival)) + 
  rcartocolor::scale_fill_carto_c(type = "quantitative", 
                                  palette = "SunsetDark",
                                  direction = 1, 
                                  limits = c(0.0, 60),  ## max percent overall
                                  breaks = seq(0.0, 60, by = 20),
                                  labels = glue::glue("{seq(0.0, 60, by = 20)}")) +
  ggplot2::geom_sf(data = coastline, size = 0.4) +
  ggplot2::coord_sf(xlim = map_limits$x, ylim = map_limits$y) +
  ggplot2::ggtitle("") + 
  ggplot2::theme_bw() +
  ggplot2::labs(fill = "Mean age (Myr)") +
  ggplot2::guides(fill = guide_colorbar(barheight = unit(2.3, units = "mm"),  
                                        direction = "horizontal",
                                        ticks.colour = "grey20",
                                        title.position = "top",
                                        label.position = "bottom",
                                        title.hjust = 0.5)) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.title = element_blank(),
    legend.text = element_text(size = 10), 
    axis.text = element_text(size = 5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Thu Mar  9 21:34:19 2023 ------------------------------
# diversification rate in herodotools

# in situ
mammals_diversification <- 
  Herodotools::calc_insitu_diversification(W = mammals_list_biogeo, 
                                           tree = mammals_phy_bio, 
                                           ancestral.area = node_area, 
                                           biogeo = biogeo_area, 
                                           diversification = "jetz",
                                           type = "equal.splits")


# join dataset for plot
sites_new_mammals <- dplyr::bind_cols(coords,
                              site_region =  site_region,
                              age = age_comm$mean_age_per_assemblage,
                              diversification_model_based = mammals_diversification$model_based_Jetz_harmonic_mean_site,
                              diversification = mammals_diversification$Jetz_harmonic_mean_site)

# plot of harmonic-mean of diversification rate 
#map_diversification <- 
max(sites_new_mammals$diversification)
min(sites_new_mammals$diversification)
 
dr_mammals_map <- sites_new_mammals %>% 
  ggplot2::ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = diversification)) + 
  rcartocolor::scale_fill_carto_c(type = "quantitative", 
                                  palette = "SunsetDark",
                                  direction = 1, 
                                  limits = c(0.04, 0.11),  ## max percent overall
                                  breaks = c(0.04, 0.11, by = 0.02),
                                  labels = glue::glue("{c(0.04, 0.11, by = 0.02)}")) +
  ggplot2::geom_sf(data = coastline, size = 0.4) +
  ggplot2::coord_sf(xlim = map_limits$x, ylim = map_limits$y) +
  ggplot2::ggtitle("D") + 
  ggplot2::theme_bw() +
  ggplot2::labs(fill = "DR") +
  ggplot2::guides(fill = guide_colorbar(barheight = unit(3, units = "mm"),  
                                        direction = "horizontal",
                                        ticks.colour = "grey20",
                                        title.position = "top",
                                        label.position = "bottom",
                                        title.hjust = 0.5)) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.title = element_blank(),
    axis.text = element_text(size = 5)
  )

map_mammals_complete <- mammals_richness + map_joint_evoregion_mammals +
  map_age_mammals + dr_mammals_map

# Save dataset
dr_age_mammals <- cbind(coords_id,sites_new_mammals[,4:6])
write.table(dr_age_mammals,"assemblage_age/metrics_mammals_full.txt")
