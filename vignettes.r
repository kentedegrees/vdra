cat("\014")
devtools::install(pkg = ".", build_vignettes = TRUE)
Sys.time()
# library(vdra)
# vignette(package = "vdra")
browseVignettes("vdra")
# browseVignettes()
