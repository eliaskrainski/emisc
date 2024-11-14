#' Function to read Brazilian municipalities map for a year
#' @param y integer for year
#' @param rfolder the root folder
#' @return sf object with first and second column
#' as geocode and name
brmaps <- function(y, rfolder = "~/maps/br/") {
    stopifnot(length(y)==1)
    if(y %in% c(1872, 1900, 1911, 1920, 1933, 1940,
              1950, 1960, 1970, 1980, 1991)) {
        unzip(
            paste0(
                rfolder, "antigos_1872_1991/",
                "05_malha_municipal_", y, ".zip")
        )
        map <- sf::st_read(
            paste0(
                "05-malha\ municipal\ ", y, ".shp"),
            options = "ENCODING=WINDOWS-1252"
        )
        names(map)[1:2] <- c("geocode", "name")
        file.remove(
            paste0(
                "05-malha\ municipal\ ", y,
                c(".dbf", ".prj", ".sbn", ".sbx", ".shp", ".shp.xml", ".shx")
            )
        )
    }
    if(y == 2000) {
        map <- sf::st_read(
            paste0(rfolder, "censo2000/55MU500G.shp")
        )
        names(map)[c(5, 6)] <- c("geocode", "name")
    }
    if(y == 2010) {
        load(paste0(rfolder, "censo2010/br10m.RData"))
        map <- as(br10m, "sf")
        rm(br10m)
        names(map)[2:3] <- c("geocode", "name")
    }
    if(y == 2022) {
        map <- sf::st_read(
            paste0(rfolder, "censo2022/BR_Municipios_2022.shp"),
            options = "ENCODING=WINDOWS-1252"
        )

        names(map)[1:2] <- c("geocode", "name")
    }
    return(map)
}
