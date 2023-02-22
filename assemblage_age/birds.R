# Sat Jan 28 10:21:01 2023 ------------------------------
# new metric 'assemblage age' from herodotools package
# birds

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
birds_full_list <- read.table('list_of_species/composition_birds.txt', header=T)

# tree data
tree_birds <- read.nexus("phylogeny/consensus_mcct.tre")
coords_id <- read.table("Shapefiles/id_coordinates.txt")
nrow(coords_id) # 432 sites

# ordem alfabética
birds_list <- birds_full_list[ , order(names(birds_full_list ))]
head(birds_list)
# trocando . por _
names(birds_list) <- gsub(".", "_", fixed=T, names(birds_list))

ncol(birds_list) # 816
tree_birds # 716 spp
# prune with phylogeny
birds_phy <- prune.sample(birds_list, mcct_birds)
# confering nomenclature names
names <- birds_mcct$tip.label
View(as.data.frame(names))
ncol(birds_list) # precisamos remover 121 espécies que estão a mais na composição

rem.col.phy <- c("Agelaioides_fringillarius","Amazilia_sapphirina","Anabacerthia_lichtensteini","Anodorhynchus_glaucus","Antrostomus_rufus","Antrostomus_sericocaudatus","Anumara_forbesi","Aramides_cajaneus","Asemospiza_fuliginosa","Asthenes_moreirae","Atticora_tibialis","Automolus_lammi","Basileuterus_auricapilla","Buteogallus_coronatus","Buteogallus_lacernulatus","Calidris_subruficollis","Cantorchilus_leucotis","Cantorchilus_longirostris","Castanozoster_thoracicus","Celeus_galeatus","Celeus_ochraceus","Celeus_tinnunculus","Ceratopipra_rubrocapilla","Cercomacroides_laeta","Chlorostilbon_notatus","Chordeiles_nacunda","Ciccaba_huhula","Ciccaba_virgata","Cichlocolaptes_holti","Cichlocolaptes_mazarbarnetti","Clibanornis_rectirostris","Colaptes_campestroides","Cyanocorax_coeruleus","Cyanoloxia_brissonii","Diopsittaca_cumanensis","Elaenia_sordida","Eupsittula_aurea","Eupsittula_cactorum","Formicivora_acutirostris","Formicivora_paludicola","Geranoaetus_albicaudatus","Griseotyrannus_aurantioatrocristatus","Guyramemua_affinis","Herpsilochmus_scapularis","Hirundinea_bellicosa","Hoploxypterus_cayanus","Hydropsalis_maculicaudus","Icterus_pyrrhopterus","ID","Islerothraupis_cristata","Leistes_superciliaris","Lipaugus_ater","Lipaugus_conditus","Mareca_sibilatrix","Microspingus_cabanisi","Microspingus_cinereus","Microspingus_lateralis","Myiodynastes_solitarius","Myiothlypis_flaveola","Myiothlypis_leucoblephara","Myiothlypis_leucophrys","Myiothlypis_rivularis","Myrmoderus_loricatus","Myrmoderus_ruficauda","Myrmoderus_squamosus","Nycticryphes_semicollaris","Nyctipolus_hirundinaceus","Ortalis_araucuan","Ortalis_squamata","Orthopsittaca_manilatus","Parabuteo_leucorrhous","Paraclaravis_geoffroyi","Phalcoboenus_chimango","Pheugopedius_genibarbis","Philohydor_lictor","Pionus_reichenowi","Pipraeidea_bonariensis","Pogonotriccus_eximius","Porphyriops_melanops","Pseudastur_polionotus","Pseudopipra_pipra","Psittacara_acuticaudatus","Psittacara_leucophthalmus","Pygochelidon_melanoleuca","Ramphastos_ariel","Ramphastos_culminatus","Rhopias_gularis","Rufirallus_viridis","Rupornis_magnirostris","Sarkidiornis_sylvicola", "Sclerurus_cearensis","Scytalopus_gonzagai","Scytalopus_petrophilus","Serpophaga_griseicapilla","Setopagis_parvula","Spatula_platalea","Spatula_versicolor","Spinus_magellanicus","Spinus_yarrellii","Sporophila_angolensis","Sporophila_beltoni","Sporophila_maximiliani","Sporophila_pileata","Stephanoxis_loddigesii","Sternula_superciliaris","Stigmatura_bahiae","Synallaxis_cinerea","Synallaxis_hellmayri","Systellura_longirostris","Tangara_brasiliensis","Tangara_cyanomelas","Tangara_flava","Tangara_ornata","Tangara_palmarum","Tangara_sayaca","Thlypopsis_pyrrhocoma","Tityra_braziliensis","Trogon_aurantius","Vireo_chivi","Xenops_rutilus","Xiphorhynchus_atlanticus") 

#Retirando na lista
birds_list_phy <- birds_list[,!(names(birds_list)%in% rem.col.phy)]
ncol(birds_list_phy ) #695
birds_phy # 695 spp.

View(birds_phy$edge.lengt == 0)

# squamata richness
richness_birds <- rowSums(birds_list_phy)

# just coordinates
coords <- coords_id[,-3]
names(coords)

# map squamata richness
birds_richness <- 
  dplyr::bind_cols(coords, richness = richness_birds) %>% 
  ggplot2::ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = richness_birds)) + 
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
birds_richness

# construindo evoregions
birds_phy <- force.ultrametric(birds_phy)

regions_birds <- 
  Herodotools::calc_evoregions(
    comm = birds_list_phy,
    phy = birds_phy
  )

regions_birds$cluster_evoregions
site_region_birds <- regions_birds$cluster_evoregions

# visualize evoregions
evoregion_birds <- data.frame(
  coords, 
  site_region_birds
)

birds_evoregion <- terra::rast(evoregion_birds)
birds_evoregion

# Converting evoregion to a spatial polygon data frame, so it can be plotted
birds_evoregion <- terra::as.polygons(birds_evoregion) %>% 
  sf::st_as_sf()

# Assigning the same projection to both spatial objects
sf::st_crs(birds_evoregion) <- sf::st_crs(coastline)

# Colours to plot evoregions
col_five_hues <- c(
  "#3d291a",
           "#a9344f",
           "#578a5b",
           "#83a6c4",
           "#f56cd9")
           
map_evoregion_birds <- 
  evoregion_birds %>% 
  ggplot2::ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = site_region_birds)) + 
  ggplot2::scale_fill_manual(
    name = "", 
    labels = LETTERS[1:5],
    values = rev(col_five_hues)
  ) +
  ggplot2::geom_sf(data = coastline) +
  ggplot2::geom_sf(
    data = birds_evoregion, 
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
map_evoregion_birds

# pertencimento a cada evoregion
# Selecting only axis with more than 5% of explained variance from evoregion output
axis_sel <- which(regions_birds$PCPS$prop_explainded >= regions_birds$PCPS$tresh_dist)
PCPS_thresh <- regions_birds$PCPS$vectors[, axis_sel] 

# distance matrix using 4 significant PCPS axis accordingly to the 5% threshold 
dist_phylo_PCPS <- vegan::vegdist(PCPS_thresh, method = "euclidean")

# calculating affiliation values for each assemblage 
afi <- calc_affiliation_evoreg(phylo.comp.dist = dist_phylo_PCPS,
                               groups = site_region_birds) 

# binding the information in a data frame
sites <- dplyr::bind_cols(coords, site_region =  site_region_birds, afi)

map_joint_evoregion_birds <- 
  evoregion_birds %>% 
  ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = site_region_birds), 
                       alpha = sites[, "afilliation"]) + 
  ggplot2::scale_fill_manual(
    name = "", 
    labels = LETTERS[1:5],
    values = rev(col_five_hues)
  ) +
  ggplot2::geom_sf(data = coastline, size = 0.4) +
  ggplot2::geom_sf(
    data = birds_evoregion, 
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

map_joint_evoregion_birds

# temos que definir a ocorrência de cada espécie nas evoregiões. Para fazer isso, podemos usar a função get_region_occe obter um quadro de dados de espécies nas linhas e evoregiões nas colunas.
birds_region <- get_region_occ_v2(comm = birds_list_phy, site.region = site_region_birds)
nrow(birds_region) # 695
ncol(birds_list_phy) #695
?get_region_occ
# O objeto criado na última etapa pode ser usado em uma função auxiliar no Herodotools para produzir facilmente o arquivo Phyllip necessário para executar a análise da reconstrução da área ancestral usando o BioGeoBEARS.
# save phyllip file
getwd()
Herodotools::get_tipranges_to_BioGeoBEARS(birds_region,filename = "assemblage_age/geo_area_birds.data",areanames = NULL) 

###------------------------------------------------------------------------
# depois da reconstrução de range ancestral 
# Aqui começamos a rodar assemblage_age
# tree file for 
birds_phy_bio <- force.ultrametric(birds_phy)
birds_phy_bio <- ape::multi2di(birds_phy_bio)

#remove <- c("Crotophaga_major", "Crypturellus_parvirostris", "Turdus_leucomelas")
#birds_phy_bio <- drop.tip(birds_phy_bio, remove)
ape::write.tree(birds_phy_bio , 'assemblage_age/birds_biogeo.new')

# idade da comunidade
# converting numbers to character
biogeo_area <- data.frame(biogeo = chartr("12345", "ABCDE", evoregion_birds$site_region_birds)) 

# getting the ancestral range area for each node 
node_area <- 
  Herodotools::get_node_range_BioGeoBEARS(
    resDECj,
    phyllip.file = "assemblage_age/geo_area_birds.data",
    birds_phy_bio,
    max.range.size = 3 
  )

# remove species in list of species
ncol(birds_list_phy) #695 spp
ncol(birds_list_phy) #695 spp

# calculating age arrival 
age_comm <- Herodotools::calc_age_arrival(W = birds_list_phy, 
                                          tree = birds_phy_bio, 
                                          ancestral.area = node_area, 
                                          biogeo = biogeo_area) 
min(age_comm$mean_age_per_assemblage)

sites <- dplyr::bind_cols(coords, site_region =  evoregion_birds$site_region_birds, age =age_comm$mean_age_per_assemblage)

#map_age <- 
map_age_birds <-  sites %>% 
  ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = V1, y = V2, fill = mean_age_arrival)) + 
  rcartocolor::scale_fill_carto_c(type = "quantitative", 
                                  palette = "SunsetDark",
                                  direction = 1, 
                                  limits = c(0.0, 25),  ## max percent overall
                                  breaks = seq(0.0, 25, by = 10),
                                  labels = glue::glue("{seq(0.0, 25, by = 10)}")) +
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


# Save dataset
age_birds <- cbind(coords_id,sites[,3:4])
write.table(age_birds,"assemblage_age/age_birds.txt")
read.table('assemblage_age/age_birds.txt')
