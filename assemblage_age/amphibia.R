# Thu Jan 26 20:49:12 2023 ------------------------------
# new metric 'assemblage age' from herodotools package
# amphibia 

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
  x = c(-70, -30),
  y = c(-40, 00)
)

# load data
getwd()
amphibia_full_list <- read.table('02_metrics/harmonized_data/list_amphibia.txt', header=T)
tree_amphibia <- read.tree("02_metrics/harmonized_data/amphibia.tre")
coords_id <- read.table("Shapefiles/id_coordinates.txt")
nrow(coords_id) # 432 sites

# adjustment in amphibia names
amphibia_list <- amphibia_full_list[ , order(names(amphibia_full_list))] #Ordem alfabética
#Arrumando nomes
names(amphibia_list) <- gsub("Scinax_x.signatus", "Scinax_x-signatus", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Scinax_v.signatus", "Scinax_v-signatus", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Brachycephalus_margariatus", "Brachycephalus_margaritatus", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Scinax_strigilata", "Scinax_strigilatus", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Pithecopus", "Phyllomedusa", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Ololygon", "Scinax", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Boana", "Hypsiboas", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Phyllomedusa_megacephalus", "Phyllomedusa_megacephala", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Phyllomedusa_nordestinus", "Phyllomedusa_nordestina", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Phyllomedusa_azureus", "Phyllomedusa_azurea", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Lithobates_palmipes", "Rana_palmipes", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Scinax_x.signatus", "Scinax_x-signatus", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Brachycephalus_margariatus", "Brachycephalus_margaritatus", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Aparasphenodon_ararapa", "Aparasphenodon_arapapa", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Agalychnis_aspera", "Hylomantis_aspera", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Physalaemus_nattereri", "Eupemphix_nattereri", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Agalychnis_granulosa", "Hylomantis_granulosa", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Scinax_v.signatus", "Scinax_v-signatus", fixed=T, names(amphibia_list))
names(amphibia_list) <- gsub("Scinax_strigilata", "Scinax_strigilatus", fixed=T, names(amphibia_list))

# prune with phylogeny
amphibia_phy <- prune.sample(amphibia_list, tree_amphibia) # 487 spp.
amphibia_phy # 421
ncol(amphibia_list) # 421

#match.phylo.comm(tree_amphibia, amphibia_list)
#ncol(amphibia_list) # precisamos remover 32 espécies que estão a mais na composição
#rem.col.phy <- c("Bokermannohyla_izecksohni","Brachycephalus_sulfuratus","Chiasmocleis_cordeiroi","Chiasmocleis_crucis","Crossodactylus_bokermanni","Crossodactylus_fransciscanus","Crossodactylus_timbuhy","Dendropsophus_baileyi","Dendropsophus_bromeliaceus", "Eleutherodactylus_bilineatus","Holoaden_pholeter","ID","Melanophryniscus_milanoi", "Melanophryniscus_xanthostomus","Scinax_strigilatus","Phyllodytes_megatympanum","Proceratophrys_fryi","Proceratophrys_mantiqueira","Proceratophrys_moratoi", "Proceratophrys_phyllostoma","Proceratophrys_pombali","Pseudopaludicola_atragula","Pseudopaludicola_pocoto","Scinax_caissara","Scinax_canastrensis","Scinax_centralis","Scinax_kautskyi","Scinax_melanodactylus","Scinax_rossaferesae","Scinax_skuki","Trachycephalus_typhonius","Vitreorana_baliomma")

#Retirando na lista
#amphibia_list_phy <- amphibia_list[,!(names(amphibia_list)%in% rem.col.phy)]
#ncol(amphibia_list_phy ) #487 spp, 432 sites
#amphibia_phy # 487 spp.

# amphibia richness
richness_amphibia <- rowSums(amphibia_list)

# just coordinates
coords <- coords_id[,-3]
names(coords)

# map amphibia richness
anura_richness <- 
  dplyr::bind_cols(coords, richness = richness_amphibia) %>% 
  ggplot2::ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = richness_amphibia)) + 
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
anura_richness

# construindo evoregions
amphibia_phy <- force.ultrametric(amphibia_phy)

regions <- 
  Herodotools::calc_evoregions(
    comm = amphibia_list,
    phy = amphibia_phy
  )

regions$cluster_evoregions
site_region <- regions$cluster_evoregions

# visualize evoregions
evoregion_df <- data.frame(
  coords, 
  site_region
)

r_evoregion <- terra::rast(evoregion_df)
r_evoregion

# Converting evoregion to a spatial polygon data frame, so it can be plotted
sf_evoregion <- terra::as.polygons(r_evoregion) %>% 
  sf::st_as_sf()

# Assigning the same projection to both spatial objects
sf::st_crs(sf_evoregion) <- sf::st_crs(coastline)

# Colours to plot evoregions
col_five_hues <- c(
           "#3d291a",
           "#a9344f",
           "#578a5b",
           "#83a6c4",
           "#fcc573",
           "#fc73ee"
)

map_evoregion <- 
  evoregion_df %>% 
  ggplot2::ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = site_region)) + 
  ggplot2::scale_fill_manual(
    name = "", 
    labels = LETTERS[1:6],
    values = rev(col_five_hues)
  ) +
  ggplot2::geom_sf(data = coastline) +
  ggplot2::geom_sf(
    data = sf_evoregion, 
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
map_evoregion

# pertencimento a cada evoregion
# Selecting only axis with more than 5% of explained variance from evoregion output
axis_sel <- which(regions$PCPS$prop_explainded >= regions$PCPS$tresh_dist)
PCPS_thresh <- regions$PCPS$vectors[, axis_sel] 

# distance matrix using 4 significant PCPS axis accordingly to the 5% threshold 
dist_phylo_PCPS <- vegan::vegdist(PCPS_thresh, method = "euclidean")

# calculating affiliation values for each assemblage 
afi <- calc_affiliation_evoreg(phylo.comp.dist = dist_phylo_PCPS,
                               groups = site_region) 

# binding the information in a data frame
sites <- dplyr::bind_cols(coords, site_region =  site_region, afi)

map_joint_evoregion_afilliation <- 
  evoregion_df %>% 
  ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = site_region), 
                       alpha = sites[, "afilliation"]) + 
  ggplot2::scale_fill_manual(
    name = "", 
    labels = LETTERS[1:6],
    values = rev(col_five_hues)
  ) +
  ggplot2::geom_sf(data = coastline, size = 0.4) +
  ggplot2::geom_sf(
    data = sf_evoregion, 
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

map_joint_evoregion_afilliation

# temos que definir a ocorrência de cada espécie nas evoregiões. Para fazer isso, podemos usar a função get_region_occe obter um quadro de dados de espécies nas linhas e evoregiões nas colunas.
a_region <- get_region_occ_v2(comm = amphibia_list, site.region = site_region)
nrow(a_region) # 409
ncol(amphibia_list) # 421
# por algum motivo a função esta excluindo 9 espécies, essas 9 espécies tem ampla distribuição, agora o porque disso estar acontecendo não faço a mínima ideia (?)

# O objeto criado na última etapa pode ser usado em uma função auxiliar no Herodotools para produzir facilmente o arquivo Phyllip necessário para executar a análise da reconstrução da área ancestral usando o BioGeoBEARS.
# save phyllip file
getwd()
Herodotools::get_tipranges_to_BioGeoBEARS(a_region,filename = "assemblage_age/geo_area_amphibia_harm.data",areanames = NULL)

# tree file for BioGeoBears
# We need need to remove some species and preparing the tree file 
amphibia_phy_bio <- force.ultrametric(tree_amphibia)
amphibia_phy_bio <- ape::multi2di(tree_amphibia)

remove <- c("Dendropsophus_minutus", "Hypsiboas_crepitans", "Hypsiboas_faber", "Hypsiboas_geographicus", "Leptodactylus_fuscus", "Leptodactylus_latrans", "Leptodactylus_mystacinus", "Physalaemus_cuvieri","Odontophrynus_americanus", "Rhinella_icterica", "Scinax_alter", "Scinax_squalirostris")

amphibia_phy_bio <- drop.tip(amphibia_phy_bio, remove)
ape::write.tree(amphibia_phy_bio , 'assemblage_age/amphibia_biogeo_harm.new')

# Assemblage age
# converting numbers to character
biogeo_area <- data.frame(biogeo = chartr("123456", "ABCDEF", evoregion_df$site_region)) 

# getting the ancestral range area for each node 
node_area <- 
  Herodotools::get_node_range_BioGeoBEARS(
    resDECj,
    phyllip.file = "assemblage_age/geo_area_amphibia_harm.data",
    amphibia_phy_bio,
    max.range.size = 3 
  )

# remove species in list of species
ncol(amphibia_list)
amphibia_list_biogeo <- amphibia_list[,!(names(amphibia_list)%in% remove)]
ncol(amphibia_list_biogeo)

# calculating age arrival 
age_comm <- Herodotools::calc_age_arrival(W = amphibia_list_biogeo, 
                                          tree = amphibia_phy_bio, 
                                          ancestral.area = node_area, 
                                          biogeo = biogeo_area) 
sites <- dplyr::bind_cols(coords, site_region =  site_region, age = age_comm$mean_age_per_assemblage)
max(sites$mean_age_arrival)
min(sites$mean_age_arrival)
#plot map
#map_age <- 
  sites %>% 
  ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = mean_age_arrival)) + 
  rcartocolor::scale_fill_carto_c(type = "quantitative", 
                                  palette = "SunsetDark",
                                  direction = 1, 
                                  limits = c(0.0, 35),  ## max percent overall
                                  breaks = seq(0.0, 35, by = 10),
                                  labels = glue::glue("{seq(0.0, 35, by = 10)}")) +
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
age_amphibia <- cbind(coords_id,sites[,3:4])
write.table(age_amphibia,"assemblage_age/age_amphibia_harm.txt")
read.table('assemblage_age/age_amphibia_harm.txt')
