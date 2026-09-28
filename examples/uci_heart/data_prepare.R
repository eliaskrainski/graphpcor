#################################################################
### 1) download
### 2) filter
### 3) prepare and save
#################################################################

## load the packages
library(ggplot2)
library(ggpubr)

## set the working directory
setwd(here::here("uci_heart"))
getwd()

###########################################################
## 1) get the UCI Heart Diseasex data from four hospitals
###########################################################
url <- "https://archive.ics.uci.edu/static/public/45/"
zipfl <- "heart+disease.zip"
if(!file.exists(zipfl))
    download.file(paste0(url, zipfl), zipfl)
unzip(zipfl)
dir()

### read data from 4 hospitals
anames <- c("age","sex","cp","trestbps","chol","fbs","restecg",
            "thalach","exang","oldpeak","slope","ca","thal","target")
names(anames) <- anames
locs <- c("cleveland", "hungarian", "switzerland", "va")
nlocs <- length(locs)
names(locs) <- c(locs[1:3], 'long_beach')
ldat0 <- lapply(locs, function(l)
    data.frame(local = l,
               read.csv(
                   file = paste0("processed.", l, ".data"),
                   header = FALSE,
                   na.strings = "?",
                   col.names = anames))
    ); names(ldat0) <- locs
sapply(ldat0, nrow)

### set chol==0 to NA
for(k in 1:nlocs) {
    ii <- which(ldat0[[k]]$chol==0)
    ldat0[[k]]$chol[ii] <- NA
}

### gplot
ggl0 <- lapply(anames, function(y) {
    dy <- do.call(
        "rbind", lapply(ldat0, function(d)
            d[c("local", y)]))
    names(dy)[2] <- 'y'
    ggplot(dy) + theme_minimal() +
        geom_histogram(aes(x=y)) + xlab(y) +
        facet_wrap(~local)
})

## ggarrange(plotlist = ggl0)

## selected variables (ordered)
snames <- c("chol", "age", "thalach", "oldpeak", "trestbps")
names(snames) <- snames
(p <- length(snames))

ggarrange(plotlist = ggl0[snames])

## remove lines with missing data on the selected variables
for(k in 1:nlocs) {
    ii <- which(complete.cases(ldat0[[k]][snames]))
    print(c(n0=nrow(ldat0[[k]]), n1=length(ii)))
    ldat0[[k]] <- ldat0[[k]][ii, ] 
}

## number of observations in each hospital
locals.n <- sapply(ldat0, nrow)
locals.n

## data.frame with 1st hospital and another with the other three
alldat <- do.call("rbind", ldat0[locals.n>0])
dat1 <- ldat0[[1]][snames]
dat2 <- do.call("rbind", ldat0[which(locals.n>0)[-1]])[snames]
c(nrow(dat1), nrow(dat2))

saveRDS(dat1, 'dat1.rds')
saveRDS(dat2, 'dat2.rds')

