
setwd(here::here("examples", "cancer12br4"))
getwd()

if(!require(geobr)) {
    install.packages("geobr")
}

if(!file.exists("rds/sbrmap.rds")) {
    if(!file.exists("rds/brmap.rds")) {
        brmap <- read_municipality(year = 2022)
        saveRDS(object = brmap, 
                file = "rds/brmap.rds")
        sbrmap <- brmap[substr(brmap$code_muni, 1, 1) == "4", ]
        saveRDS(object = sbrmap,
                file = "rds/sbrmap.rds")
    }
}

sbrmap <- readRDS("rds/sbrmap.rds")
head(sbrmap,2)

library(spdep)

nbl <- poly2nb(sbrmap)
table(card(nbl))

nb2INLA(file = paste0("graph", nrow(sbrmap)),
        nb = nbl)
