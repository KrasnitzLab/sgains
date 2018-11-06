library("SCclust")
# library("futile.logger")

options(echo=TRUE)

# flog.threshold(DEBUG)
# flog.info("running sgains scclust pipeline script")

args <- commandArgs(trailingOnly = TRUE)

print(args)

scgv_dir <- args[1]
case_name <- args[2]

varbin_dir <- args[3]
varbin_suffix <- args[4]

bin_boundaries <- args[5]

cytoband <- args[6]
nsim <- as.integer(args[7])
sharemin <- as.double(args[8])


# flog.debug("varbin_dir:", varbin_dir)

sgains_pipeline(
    scgv_dir, case_name,
    varbin_dir, varbin_suffix,
    bin_boundaries,
    cytoband,
    nsim=nsim,
    sharemin=sharemin)
