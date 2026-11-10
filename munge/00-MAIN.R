# Project specific packages, functions and settings -----------------------

source(here::here("setup/setup.R"))

# Munge data --------------------------------------------------------------

source(here("munge/01-load-data.R"))
source(here("munge/02-fix-vars.R"))

# write.xlsx(esclab, "./data/meta-data/metavars_tmp.xlsx")
metaesc <- read.xlsx("./data/meta-data/metavars.xlsx")

# Cache/save data ---------------------------------------------------------

# Version number

version <- "004"

dir.create(paste0("./data/clean-data/v", version))

save(
  file = here(paste0("data/clean-data/v", version, "/escdata", version, ".RData")),
  list = c(
    "escdata",
    "metaesc"
  )
)
