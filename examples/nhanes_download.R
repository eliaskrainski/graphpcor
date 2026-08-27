###############################################################################
# Wang, Mukherjee & Park (2018)
#
# 18 urinary/blood metal biomarkers
# NHANES 2003-2004 through 2013-2014
#
# Outputs:
#   * harmonized 18-metal data for each NHANES cycle
#   * pooled 18-metal data
#   * cycle-specific Spearman correlation matrices
#   * pooled Spearman correlation matrix
#   * inverse-correlation (precision) matrices
#   * partial-correlation matrices
#
# IMPORTANT:
#   - NHANES already replaces observations below LLOD by LLOD/sqrt(2).
#   - Therefore we DO NOT perform a second LOD substitution.
#   - TMO / URXUTM is NOT one of the 18 Wang et al. biomarkers.
#
###############################################################################

## ============================================================================
## 0. Packages
## ============================================================================

packages <- c(
  "haven",
  "dplyr",
  "purrr",
  "tibble",
  "tidyr"
)

to_install <- packages[
  !packages %in% rownames(installed.packages())
]

if (length(to_install) > 0)
  install.packages(to_install)

library(haven)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)


## ============================================================================
## 1. Settings
## ============================================================================

options(timeout = 600)

data_dir <- "NHANES_Wang2018"

dir.create(
  data_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================================
## 2. The 18 biomarkers
## ============================================================================

metal_names <- c(
  "u_antimony",
  "u_arsenic",
  "u_arsenobetaine",
  "u_mma",
  "u_dma",
  "u_barium",
  "u_cadmium",
  "u_cobalt",
  "u_cesium",
  "u_molybdenum",
  "u_mercury",
  "u_lead",
  "u_thallium",
  "u_tungsten",
  "u_uranium",
  "b_cadmium",
  "b_lead",
  "b_mercury"
)


## ============================================================================
## 3. Exact NHANES files
##
## 2003-2004:
##   L06HM_C   urine metals
##   L06UAS_C  total + speciated arsenics
##   L06UHG_C  urine mercury
##   L06BMT_C  blood Cd/Pb/Hg
##
## 2005-2006:
##   UHM_D
##   UAS_D
##   UHG_D
##   PbCd_D
##
## 2007-2008:
##   UHM_E
##   UAS_E
##   UHG_E
##   PBCD_E
##
## 2009-2010:
##   UHM_F
##   UAS_F
##   UHG_F
##   PbCd_F
##
## 2011-2012:
##   UHM_G
##   UAS_G
##   UHG_G
##   PbCd_G
##
## 2013-2014:
##   UM_H     urine metals
##   UAS_H    speciated arsenics
##   UTAS_H   TOTAL urinary arsenic
##   UHG_H    urine mercury
##   PBCD_H   blood Cd/Pb/Hg
##
## ============================================================================

files <- tibble::tribble(
  ~cycle,       ~year, ~urine_metals, ~arsenic_file, ~total_arsenic_file,
  ~urine_mercury, ~blood_metals,

  "2003-2004", 2003,
  "L06HM_C", "L06UAS_C", "L06UAS_C",
  "L06UHG_C", "L06BMT_C",

  "2005-2006", 2005,
  "UHM_D", "UAS_D", "UAS_D",
  "UHG_D", "PBCD_D",

  "2007-2008", 2007,
  "UHM_E", "UAS_E", "UAS_E",
  "UHG_E", "PBCD_E",

  "2009-2010", 2009,
  "UHM_F", "UAS_F", "UAS_F",
  "UHG_F", "PBCD_F",

  "2011-2012", 2011,
  "UHM_G", "UAS_G", "UAS_G",
  "UHG_G", "PBCD_G",

  "2013-2014", 2013,
  "UM_H", "UAS_H", "UTAS_H",
  "UHG_H", "PBCD_H"
  )

## ============================================================================
## 4. CDC download URL
## ============================================================================

cdc_url <- function(year, file) {

  paste0(
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/",
    year,
    "/DataFiles/",
    file,
    ".XPT"
  )

}


## ============================================================================
## 5. Download/read XPT
## ============================================================================

get_nhanes_file <- function(year, file) {

  if (is.na(file))
    return(NULL)

  destination <- file.path(
    data_dir,
    paste0(file, ".XPT")
  )

  url <- cdc_url(year, file)

  if (!file.exists(destination)) {

    message("Downloading ", file, " ...")

    download.file(
      url = url,
      destfile = destination,
      mode = "wb",
      quiet = FALSE
    )

  } else {

    message("Using cached ", file)

  }

  dat <- read_xpt(destination)

  names(dat) <- toupper(names(dat))

  dat
}


## ============================================================================
## 6. Required raw CDC variables
## ============================================================================

required_raw <- list(

  urine_metals = c(
    "SEQN",
    "URXUSB",   # antimony
    "URXUBA",   # barium
    "URXUCD",   # cadmium
    "URXUCO",   # cobalt
    "URXUCS",   # cesium
    "URXUMO",   # molybdenum
    "URXUPB",   # lead
    "URXUTL",   # thallium
    "URXUTU",   # tungsten
    "URXUUR"    # uranium
  ),

  arsenic_species = c(
    "SEQN",
    "URXUAB",   # arsenobetaine
    "URXUDMA",  # DMA
    "URXUMMA"   # MMA
  ),

  total_arsenic = c(
    "SEQN",
    "URXUAS"
  ),

  urine_mercury = c(
    "SEQN",
    "URXUHG"
  ),

  blood_metals = c(
    "SEQN",
    "LBXBCD",
    "LBXBPB",
    "LBXTHG"
  )
)


## ============================================================================
## 7. Check raw variables
## ============================================================================

check_vars <- function(dat, required, file) {

  missing <- setdiff(
    required,
    names(dat)
  )

  if (length(missing) > 0) {

    stop(
      "\nMissing variables in ", file, ":\n",
      paste(missing, collapse = ", "),
      "\n"
    )

  }

  invisible(TRUE)
}


## ============================================================================
## 8. Extract one cycle
## ============================================================================

extract_cycle <- function(
    cycle,
    year,
    urine_metals,
    arsenic_file,
    total_arsenic_file,
    urine_mercury,
    blood_metals
) {

  message("\n============================================")
  message("Processing ", cycle)
  message("============================================")

  ## ------------------------------------------------------------
  ## Download/read files
  ## ------------------------------------------------------------

  um <- get_nhanes_file(
    year = year,
    file = urine_metals
  )

  ua <- get_nhanes_file(
    year = year,
    file = arsenic_file
  )

  uhg <- get_nhanes_file(
    year = year,
    file = urine_mercury
  )

  pbcd <- get_nhanes_file(
    year = year,
    file = blood_metals
  )

  ## Total arsenic:
  ## 2003-2012: same file as arsenic species
  ## 2013-2014: separate UTAS_H file

  if (identical(arsenic_file, total_arsenic_file)) {

    ua_total <- ua

  } else {

    ua_total <- get_nhanes_file(
      year = year,
      file = total_arsenic_file
    )

  }


  ## ------------------------------------------------------------
  ## Check required variables
  ## ------------------------------------------------------------

  check_vars(
    um,
    c(
      "SEQN",
      "URXUSB",
      "URXUBA",
      "URXUCD",
      "URXUCO",
      "URXUCS",
      "URXUMO",
      "URXUPB",
      "URXUTL",
      "URXUTU",
      "URXUUR"
    ),
    urine_metals
  )

  check_vars(
    ua,
    c(
      "SEQN",
      "URXUAB",
      "URXUDMA",
      "URXUMMA"
    ),
    arsenic_file
  )

  check_vars(
    ua_total,
    c(
      "SEQN",
      "URXUAS"
    ),
    total_arsenic_file
  )

  check_vars(
    uhg,
    c(
      "SEQN",
      "URXUHG"
    ),
    urine_mercury
  )

  check_vars(
    pbcd,
    c(
      "SEQN",
      "LBXBCD",
      "LBXBPB",
      "LBXTHG"
    ),
    blood_metals
  )


  ## ------------------------------------------------------------
  ## Urinary metals
  ## ------------------------------------------------------------

  um2 <- um |>
    dplyr::select(
      SEQN,
      URXUSB,
      URXUBA,
      URXUCD,
      URXUCO,
      URXUCS,
      URXUMO,
      URXUPB,
      URXUTL,
      URXUTU,
      URXUUR
    ) |>
    dplyr::rename(
      u_antimony   = URXUSB,
      u_barium     = URXUBA,
      u_cadmium    = URXUCD,
      u_cobalt     = URXUCO,
      u_cesium     = URXUCS,
      u_molybdenum = URXUMO,
      u_lead       = URXUPB,
      u_thallium   = URXUTL,
      u_tungsten   = URXUTU,
      u_uranium    = URXUUR
    )


  ## ------------------------------------------------------------
  ## Arsenic species
  ## ------------------------------------------------------------

  ua2 <- ua |>
    dplyr::select(
      SEQN,
      URXUAB,
      URXUDMA,
      URXUMMA
    ) |>
    dplyr::rename(
      u_arsenobetaine = URXUAB,
      u_dma           = URXUDMA,
      u_mma           = URXUMMA
    )


  ## ------------------------------------------------------------
  ## Total arsenic
  ## ------------------------------------------------------------

  ua_total2 <- ua_total |>
    dplyr::select(
      SEQN,
      URXUAS
    ) |>
    dplyr::rename(
      u_arsenic = URXUAS
    )


  ## ------------------------------------------------------------
  ## Urinary mercury
  ## ------------------------------------------------------------

  uhg2 <- uhg |>
    dplyr::select(
      SEQN,
      URXUHG
    ) |>
    dplyr::rename(
      u_mercury = URXUHG
    )


  ## ------------------------------------------------------------
  ## Blood metals
  ## ------------------------------------------------------------

  pbcd2 <- pbcd |>
    dplyr::select(
      SEQN,
      LBXBCD,
      LBXBPB,
      LBXTHG
    ) |>
    dplyr::rename(
      b_cadmium = LBXBCD,
      b_lead    = LBXBPB,
      b_mercury = LBXTHG
    )


  ## ------------------------------------------------------------
  ## Merge
  ## ------------------------------------------------------------

  out <- um2 |>
    dplyr::inner_join(
      ua2,
      by = "SEQN"
    ) |>
    dplyr::inner_join(
      ua_total2,
      by = "SEQN"
    ) |>
    dplyr::inner_join(
      uhg2,
      by = "SEQN"
    ) |>
    dplyr::inner_join(
      pbcd2,
      by = "SEQN"
    ) |>
    dplyr::mutate(
      cycle = cycle,
      year = year,
      .before = 1
    )


  ## ------------------------------------------------------------
  ## Ensure numeric
  ## ------------------------------------------------------------

  out <- out |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(metal_names),
        as.numeric
      )
    )


  ## ------------------------------------------------------------
  ## Return
  ## ------------------------------------------------------------

  out |>
    dplyr::select(
      cycle,
      year,
      SEQN,
      dplyr::all_of(metal_names)
    )
}

## ============================================================================
## 9. Download and process all six cycles
## ============================================================================

cycle_data <- purrr::pmap(
  files,
  extract_cycle
)

names(cycle_data) <- files$cycle


## ============================================================================
## 10. Confirm the 18-variable common set
## ============================================================================

available_by_cycle <- map(
  cycle_data,
  ~ names(.x)
)

common_variables <- Reduce(
  intersect,
  available_by_cycle
)

common_metal_variables <- intersect(
  common_variables,
  metal_names
)

cat("\n============================================\n")
cat("COMMON BIOMARKERS\n")
cat("============================================\n")

print(common_metal_variables)

cat(
  "\nNumber of common biomarkers:",
  length(common_metal_variables),
  "\n"
)

stopifnot(
  length(common_metal_variables) == 18
)


## ============================================================================
## 11. Missingness before complete-case filtering
## ============================================================================

missingness <- map_dfr(
  cycle_data,
  function(dat) {

    tibble(
      variable = metal_names,
      n = nrow(dat),
      n_missing = map_int(
        dat[metal_names],
        ~ sum(is.na(.x))
      )
    )

  },
  .id = "cycle"
)

missingness <- missingness |>
  mutate(
    percent_missing =
      100 * n_missing / n
  )

print(missingness)


## ============================================================================
## 12. LOD treatment
##
## IMPORTANT:
##
## CDC has already inserted LLOD/sqrt(2) into the measurement variable for
## observations below the detection limit.
##
## We therefore retain the values exactly as released.
##
## The LC variables are retained in the downloaded raw files but are not
## required for the correlation calculation.
## ============================================================================

cycle_data <- map(
  cycle_data,
  function(dat) {

    dat |>
      mutate(
        across(
          all_of(metal_names),
          as.numeric
        )
      )

  }
)


## ============================================================================
## 13. Complete cases for all 18 biomarkers
##
## Using complete cases gives a single sample for every element of each
## 18 x 18 correlation matrix.
## ============================================================================

cycle_complete <- map(
  cycle_data,
  function(dat) {

    dat |>
      filter(
        if_all(
          all_of(metal_names),
          ~ !is.na(.x)
        )
      )

  }
)


## ============================================================================
## 14. Sample sizes
## ============================================================================

sample_sizes <- tibble(
  cycle = names(cycle_data),

  n_merged = map_int(
    cycle_data,
    nrow
  ),

  n_complete_18 = map_int(
    cycle_complete,
    nrow
  )
)

print(sample_sizes)


## ============================================================================
## 15. Spearman correlation
## ============================================================================

spearman_cor <- function(dat) {

  X <- dat |>
    select(all_of(metal_names))

  R <- cor(
    X,
    method = "spearman",
    use = "complete.obs"
  )

  ## Remove tiny numerical asymmetry
  R <- (R + t(R)) / 2

  R
}


R_cycle <- map(
  cycle_complete,
  spearman_cor
)


## ============================================================================
## 16. Precision matrix
##
## Omega = R^{-1}
##
## This is the precision matrix associated with the rank-transformed
## variables, i.e. the inverse Spearman correlation matrix.
## ============================================================================

precision_from_cor <- function(R) {

  R <- (R + t(R)) / 2

  eigenvalues <- eigen(
    R,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  if (min(eigenvalues) <= 1e-10) {

    stop(
      "Correlation matrix is singular or nearly singular."
    )

  }

  solve(R)
}


Omega_cycle <- map(
  R_cycle,
  precision_from_cor
)


## ============================================================================
## 17. Partial correlations
##
## rho_ij | rest =
##
##       -Omega_ij /
##       sqrt(Omega_ii * Omega_jj)
##
## ============================================================================

partial_cor_from_precision <- function(Omega) {

  d <- sqrt(
    diag(Omega)
  )

  P <- -Omega / outer(
    d,
    d
  )

  diag(P) <- 1

  P
}


P_cycle <- map(
  Omega_cycle,
  partial_cor_from_precision
)


## ============================================================================
## 18. Pooled data
## ============================================================================

pooled <- bind_rows(
  cycle_data
)


pooled_complete <- pooled |>
  filter(
    if_all(
      all_of(metal_names),
      ~ !is.na(.x)
    )
  )


cat(
  "\nPooled observations after merging:",
  nrow(pooled),
  "\n"
)

cat(
  "Pooled complete 18-metal observations:",
  nrow(pooled_complete),
  "\n"
)


## ============================================================================
## 19. Pooled matrices
## ============================================================================

R_pooled <- spearman_cor(
  pooled_complete
)

Omega_pooled <- precision_from_cor(
  R_pooled
)

P_pooled <- partial_cor_from_precision(
  Omega_pooled
)


## ============================================================================
## 20. Add pooled results to lists
## ============================================================================

R_all <- c(
  R_cycle,
  list(
    pooled = R_pooled
  )
)

Omega_all <- c(
  Omega_cycle,
  list(
    pooled = Omega_pooled
  )
)

P_all <- c(
  P_cycle,
  list(
    pooled = P_pooled
  )
)


## ============================================================================
## 21. Display results
## ============================================================================

cat("\n============================================\n")
cat("POOLED SPEARMAN CORRELATION\n")
cat("============================================\n\n")

print(
  round(
    R_pooled,
    3
  )
)


cat("\n============================================\n")
cat("POOLED PRECISION MATRIX\n")
cat("============================================\n\n")

print(
  round(
    Omega_pooled,
    3
  )
)


cat("\n============================================\n")
cat("POOLED PARTIAL CORRELATION\n")
cat("============================================\n\n")

print(
  round(
    P_pooled,
    3
  )
)


## ============================================================================
## 22. Save complete data
## ============================================================================

saveRDS(
  cycle_data,
  file.path(
    data_dir,
    "NHANES_18metals_by_cycle.rds"
  )
)

saveRDS(
  cycle_complete,
  file.path(
    data_dir,
    "NHANES_18metals_complete_by_cycle.rds"
  )
)

## ============================================================================
## 23. Save matrices
## ============================================================================

saveRDS(
  R_all,
  file.path(
    data_dir,
    "Spearman_correlations.rds"
  )
)

saveRDS(
  Omega_all,
  file.path(
    data_dir,
    "Precision_matrices.rds"
  )
)

saveRDS(
  P_all,
  file.path(
    data_dir,
    "Partial_correlations.rds"
  )
)

## ============================================================================
## 24. CSV output
## ============================================================================

write_matrix_csv <- function(
    matrix_list,
    prefix) {

  walk(
    names(matrix_list),
    function(nm) {

      filename <- paste0(
        prefix,
        "_",
        gsub(
          "[^A-Za-z0-9]+",
          "_",
          nm
        ),
        ".csv"
      )

      write.csv(
        matrix_list[[nm]],
        file.path(
          data_dir,
          filename
        ),
        row.names = TRUE
      )

    }
  )
}


write_matrix_csv(
  R_all,
  "Spearman"
)

write_matrix_csv(
  Omega_all,
  "Precision"
)

write_matrix_csv(
  P_all,
  "PartialCorrelation"
)


## ============================================================================
## 25. Save sample-size and missingness diagnostics
## ============================================================================

write.csv(
  sample_sizes,
  file.path(
    data_dir,
    "sample_sizes.csv"
  ),
  row.names = FALSE
)

write.csv(
  missingness,
  file.path(
    data_dir,
    "missingness.csv"
  ),
  row.names = FALSE
)


## ============================================================================
## END
## ============================================================================

