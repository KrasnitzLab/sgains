library("SCclust")

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

print(args)

bin_boundaries_filename <- args[1]
print(bin_boundaries_filename)

gc <- read.table(bin_boundaries_filename, header=TRUE)

file.names = args[2:length(args)]
segfile <- cbs.segment_files(file.names, gc, alpha = 0.05, nperm = 1000, undo.SD = 1.0, min.width = 5)

seg.quantal <- segfile$seg.quantal
ratio.quantal <- segfile$ratio.quantal

res1 <- preprocess_segfile(seg.quantal, gc, ploidies = TRUE)
breakpoint_table <- res1$breakpoint_table
ploidies_table <- res1$ploidies_table

smear_table <- findsmears(breakpoint_table, smear = 1, keepboundaries = FALSE, mask_XY = TRUE)

