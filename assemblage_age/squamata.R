# Sat Jan 28 10:21:01 2023 ------------------------------
# new metric 'assemblage age' from herodotools package
# squamata

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
squamata_full_list <- read.table('02_metrics/full_data/list_squamata.txt', header=T)
tree_squamata <- read.tree("02_metrics/full_data/squamata.tre")

coords_id <- read.table("Shapefiles/id_coordinates.txt")
nrow(coords_id) # 432 sites

# ordem alfabética
squamata_list <- squamata_full_list[ , order(names(squamata_full_list))] 
# trocando . por _
names(squamata_list) <- gsub(".", "_", fixed=T, names(squamata_list))

# prune with phylogeny
squamata_phy <- prune.sample(squamata_list, tree_squamata) # 487 spp.

match.phylo.comm(tree_squamata, squamata_list)
ncol(squamata_list) # precisamos remover 21 espécies que estão a mais na composição
rem.col.phy <- c("Acanthochelys_radiolata","Acanthochelys_spixii","Caiman_crocodilus","Caiman_latirostris","Caiman_yacare","Chelonoidis_carbonaria","Chelonoidis_denticulata","Hydromedusa_maximiliani","Hydromedusa_tectifera","ID","Kinosternon_scorpioides","Mesoclemmys_hogei","Mesoclemmys_tuberculata","Mesoclemmys_vanderhaegei","Paleosuchus_palpebrosus","Phrynops_geoffroanus","Phrynops_hilarii","Phrynops_tuberosus","Phrynops_williamsi","Trachemys_dorbigni","Tropidurus_imbituba") #apenas um lagarto, restante é testudine e crocodilia

#Retirando na lista
squamata_list_phy <- squamata_list[,!(names(squamata_list)%in% rem.col.phy)]
ncol(squamata_list_phy ) #382 spp
squamata_phy # 382 spp.

# squamata richness
richness_squamata <- rowSums(squamata_list_phy)

# just coordinates
coords <- coords_id[,-3]
names(coords)

# map squamata richness
squamata_richness <- 
  dplyr::bind_cols(coords, richness = richness_squamata) %>% 
  ggplot2::ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = richness_squamata)) + 
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
squamata_richness

# construindo evoregions
squamata_phy <- force.ultrametric(squamata_phy)

regions_squamata <- 
  Herodotools::calc_evoregions(
    comm = squamata_list_phy,
    phy = squamata_phy
  )

regions_squamata$cluster_evoregions
site_region_squamata <- regions_squamata$cluster_evoregions

# visualize evoregions
evoregion_squamata <- data.frame(
  coords, 
  site_region_squamata
)

site_region_squamata
evoregion_squamata$site_region_squamata

squamata_evoregion <- terra::rast(evoregion_squamata)
squamata_evoregion

# Converting evoregion to a spatial polygon data frame, so it can be plotted
squamata_evoregion <- terra::as.polygons(squamata_evoregion) %>% 
  sf::st_as_sf()

# Assigning the same projection to both spatial objects
sf::st_crs(squamata_evoregion) <- sf::st_crs(coastline)

# Colours to plot evoregions
col_five_hues <- c(
           "#3d291a",
           "#a9344f",
           "#578a5b",
           "#83a6c4")

map_evoregion_squamata <- 
  evoregion_squamata %>% 
  ggplot2::ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = site_region_squamata)) + 
  ggplot2::scale_fill_manual(
    name = "", 
    labels = LETTERS[1:4],
    values = rev(col_five_hues)
  ) +
  ggplot2::geom_sf(data = coastline) +
  ggplot2::geom_sf(
    data = squamata_evoregion, 
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
map_evoregion_squamata

# pertencimento a cada evoregion
# Selecting only axis with more than 5% of explained variance from evoregion output
axis_sel <- which(regions_squamata$PCPS$prop_explainded >= regions_squamata$PCPS$tresh_dist)
PCPS_thresh <- regions_squamata$PCPS$vectors[, axis_sel] 

# distance matrix using 4 significant PCPS axis accordingly to the 5% threshold 
dist_phylo_PCPS <- vegan::vegdist(PCPS_thresh, method = "euclidean")

# calculating affiliation values for each assemblage 
afi <- calc_affiliation_evoreg(phylo.comp.dist = dist_phylo_PCPS,
                               groups = site_region_squamata) 

# binding the information in a data frame
sites <- dplyr::bind_cols(coords, site_region =  site_region_squamata, afi)

map_joint_evoregion_squamata <- 
  evoregion_squamata %>% 
  ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = site_region_squamata), 
                       alpha = sites[, "afilliation"]) + 
  ggplot2::scale_fill_manual(
    name = "", 
    labels = LETTERS[1:6],
    values = rev(col_five_hues)
  ) +
  ggplot2::geom_sf(data = coastline, size = 0.4) +
  ggplot2::geom_sf(
    data = squamata_evoregion, 
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

map_joint_evoregion_squamata

# temos que definir a ocorrência de cada espécie nas evoregiões. Para fazer isso, podemos usar a função get_region_occe obter um quadro de dados de espécies nas linhas e evoregiões nas colunas.
squamata_region <- get_region_occ_v2(comm = squamata_list_phy, site.region = site_region_squamata)
nrow(squamata_region) # 382
ncol(squamata_list_phy) #382

# O objeto criado na última etapa pode ser usado em uma função auxiliar no Herodotools para produzir facilmente o arquivo Phyllip necessário para executar a análise da reconstrução da área ancestral usando o BioGeoBEARS.
# save phyllip file
getwd()
Herodotools::get_tipranges_to_BioGeoBEARS(squamata_region,filename = "assemblage_age/geo_area_squamata_full.data",areanames = NULL) 

###------------------------------------------------------------------------
# depois da reconstrução de range ancestral 
# Aqui começamos a rodar assemblage_age
# tree file for 
squamata_phy_bio <- force.ultrametric(squamata_phy)
squamata_phy_bio <- ape::multi2di(squamata_phy_bio)
ape::write.tree(squamata_phy_bio , 'assemblage_age/squamata_biogeo_full.new')

# idade da comunidade
# converting numbers to character
biogeo_area <- data.frame(biogeo = chartr("1234", "ABCD", evoregion_squamata$site_region_squamata)) 

# getting the ancestral range area for each node 
node_area <- 
  Herodotools::get_node_range_BioGeoBEARS(
    resDECj,
    phyllip.file = "assemblage_age/geo_area_squamata_full.data",
    squamata_phy_bio,
    max.range.size = 3 
  )


# calculating age arrival 
age_comm <- Herodotools::calc_age_arrival(W = squamata_list_phy, 
                                          tree = squamata_phy_bio, 
                                          ancestral.area = node_area, 
                                          biogeo = biogeo_area) 

sites <- dplyr::bind_cols(coords, site_region =  evoregion_squamata$site_region_squamata, age = age_comm$mean_age_per_assemblage)

max(sites$mean_age_arrival)

map_age <- sites %>% 
  ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = mean_age_arrival)) + 
  rcartocolor::scale_fill_carto_c(type = "quantitative", 
                                  palette = "SunsetDark",
                                  direction = 1, 
                                  limits = c(0, 120),  ## max percent overall
                                  breaks = seq(0, 120, by = 40),
                                  labels = glue::glue("{seq(0, 120, by = 40)}")) +
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
map_age

# Save dataset
#age_squamata <- cbind(coords_id,sites[,3:4])
#write.table(age_squamata,"assemblage_age/age_squamata_full.txt")
#read.table('assemblage_age/age_squamata.txt')

# Wed Mar  8 20:14:25 2023 ------------------------------
# diversification rate metrics

# in situ
squamata_diversification <- 
  Herodotools::calc_insitu_diversification(W = squamata_list_phy, 
                                           tree = squamata_phy_bio, 
                                           ancestral.area = node_area, 
                                           biogeo = biogeo_area, 
                                           diversification = "jetz",
                                           type = "equal.splits")


# join dataset for plot
sites_squamata <- dplyr::bind_cols(coords,
                          site_region =  site_region,
                          age = age_comm$mean_age_per_assemblage,
                          diversification_model_based = squamata_diversification$model_based_Jetz_harmonic_mean_site,
                          diversification = squamata_diversification$Jetz_harmonic_mean_site)

# Save dataset
dr_age_squamata <- cbind(coords_id,sites_squamata[,4:6])
write.table(dr_age_squamata,"assemblage_age/metrics_squamata_full.txt")

min(dr_age_squamata$diversification)
max(dr_age_squamata$diversification)


# plot of harmonic-mean of diversification rate 
map_diversification <- sites_squamata %>% 
  ggplot2::ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = diversification)) + 
  rcartocolor::scale_fill_carto_c(type = "quantitative", 
                                  palette = "SunsetDark",
                                  direction = 1, 
                                  limits = c(0.03, 0.07))+  ## max percent overall
  #breaks = seq(0, 0.030, by = 0.01))+
  #labels = glue::glue("{seq(0.010, 0.030, by = 0.010)}")) +
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

map_squamata_complete <- squamata_richness + map_joint_evoregion_squamata +
  map_age + map_diversification
