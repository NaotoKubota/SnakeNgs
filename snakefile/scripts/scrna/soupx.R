suppressPackageStartupMessages({library(Matrix); library(SoupX); library(jsonlite)})
dir.create(dirname(snakemake@log[[1]]), recursive=TRUE, showWarnings=FALSE)
logfile <- file(snakemake@log[[1]], open="wt")
sink(logfile); sink(logfile, type="message")
tryCatch({
  src <- snakemake@input[[1]]
  cfg <- fromJSON(file.path(src, "settings.json"))
  set.seed(cfg$seed %% .Machine$integer.max)
  genes <- readLines(file.path(src, "genes.tsv"))
  barcodes <- readLines(file.path(src, "barcodes.tsv"))
  toc <- as(readMM(file.path(src, "toc.mtx")), "CsparseMatrix")
  rownames(toc) <- genes; colnames(toc) <- barcodes
  clusters <- read.delim(file.path(src, "clusters.tsv"), colClasses="character", check.names=FALSE)
  profile <- read.delim(file.path(src, "soup_profile.tsv"), row.names=1, check.names=FALSE)
  rho <- setNames(rep(0, length(barcodes)), barcodes)
  method <- "disabled"
  if (cfg$enabled) {
    # Only the aggregate kb empty-droplet profile is needed here. Supplying it
    # explicitly avoids transferring a millions-of-droplets matrix into R.
    sc <- SoupChannel(tod=toc, toc=toc, calcSoupProfile=FALSE)
    sc <- setSoupProfile(sc, profile[genes, c("est", "counts"), drop=FALSE])
    sc <- setClusters(sc, setNames(clusters$soupx_cluster, clusters$barcode))
    if (is.null(cfg$contamination_fraction)) {
      sc <- autoEstCont(sc, doPlot=FALSE)
      method <- "autoEstCont"
    } else {
      sc <- setContaminationFraction(sc, cfg$contamination_fraction)
      method <- "configured_fraction"
    }
    corrected <- adjustCounts(sc, roundToInt=TRUE)
    rho <- setNames(sc$metaData$rho, rownames(sc$metaData))[barcodes]
  } else {
    corrected <- toc
  }
  if (any(!is.finite(corrected@x)) || any(corrected@x < 0) || any(corrected@x != round(corrected@x)))
    stop("SoupX returned invalid integer counts")
  dir.create(dirname(snakemake@output[["matrix"]]), recursive=TRUE, showWarnings=FALSE)
  writeMM(corrected, snakemake@output[["matrix"]])
  writeLines(rownames(corrected), snakemake@output[["genes"]])
  writeLines(colnames(corrected), snakemake@output[["barcodes"]])
  write.table(data.frame(barcode=barcodes, rho=unname(rho)), snakemake@output[["rho"]],
              sep="\t", row.names=FALSE, quote=FALSE)
  write_json(list(sample=cfg$sample, enabled=cfg$enabled, method=method,
                  n_empty=cfg$n_empty, median_rho=median(rho),
                  removed_fraction=1-sum(corrected)/sum(toc),
                  soupx_version=as.character(packageVersion("SoupX")),
                  r_version=R.version.string), snakemake@output[["summary"]], auto_unbox=TRUE, pretty=TRUE)
  print(sessionInfo())
}, finally={sink(type="message"); sink(); close(logfile)})
