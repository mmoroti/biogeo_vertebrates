get_region_occ_v2 <- function (comm, site.region) 
{
  evoregion.data <- dplyr::bind_cols(comm, site.region = site.region) %>% 
    tidyr::pivot_longer(cols = 1:ncol(comm), names_to = "species", 
                        values_to = "presence") %>% dplyr::filter(.data$presence == 
                                                                    1) %>% dplyr::group_by(.data$species, site.region) %>% 
    dplyr::summarise(n = dplyr::n()) %>% dplyr::ungroup() %>% 
    dplyr::group_by(.data$species) %>% dplyr::mutate(species.total = sum(.data$n), 
                                                     prec.occupation = .data$n/.data$species.total) %>% dplyr::ungroup() %>% 
    dplyr::filter(.data$prec.occupation >= 0.25) %>% dplyr::mutate(area = LETTERS[site.region]) %>% 
    dplyr::select(.data$species, .data$area) %>% dplyr::mutate(value = 1) %>% 
    tidyr::pivot_wider(id_cols = .data$species, names_from = .data$area, 
                       names_sort = T, values_from = .data$value, values_fill = 0) %>% 
    as.data.frame()
  species.names <- evoregion.data[, 1]
  a.regions <- evoregion.data[, -1]
  rownames(a.regions) <- species.names
  return(a.regions)
}