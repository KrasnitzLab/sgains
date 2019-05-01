# library("devtools")
# devtools::load_all("/home/lubo/Work/single-cell/SCclust")

library("SCclust")

options(echo=TRUE)

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
fdrthres <- as.double(args[9])
nshare <- as.integer(args[10])
climbtoshare <- as.integer(args[11])

# flog.debug("varbin_dir:", varbin_dir)

sgains_pipeline(
    scgv_dir, case_name,
    varbin_dir, varbin_suffix,
    bin_boundaries,
    cytoband,
    nsim=nsim,
    sharemin=sharemin,
    fdrthres=fdrthres,
    nshare=nshare,
    climbtoshare=climbtoshare)
