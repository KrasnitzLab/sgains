library("SCclust")

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

print(args)

study_name <- args[1]
output_dirname <- args[2]

bin_boundaries_filename <- args[3]

file.names = args[4: length(args)]


gc <- read.table(bin_boundaries_filename, header=TRUE, as.is=T)

segfile <- cbs.segment_files(file.names, gc, alpha = 0.05, 
		nperm = 1000, undo.SD = 1.0, min.width = 5)

seg.quantal <- segfile$seg.quantal
ratio.quantal <- segfile$ratio.quantal

res1 <- preprocess_segfile(seg.quantal, gc, ploidies = TRUE)
breakpoint_table <- res1$breakpoint_table
ploidies_table <- res1$ploidies_table

smear_table <- findsmears(
		breakpoint_table, smear = 1, 
		keepboundaries = FALSE, mask_XY = TRUE)

res2 <- findpins(breakpoint_table, smear_table)
pins <- res2$pins
pinmat <- res2$pinmat
cell_names <- res2$cell_names
print(cell_names)

res3 <- simFisher_parallel(pins, pinmat, sim_round = 500)
true_fisherPV <- res3$true_fisherPV
sim_fisherPV <- res3$sim_fisherPV

res4 <- fdr_fisherPV(true_fisherPV, sim_fisherPV, 
		cell_names, lm_max = 0.001, graphic = F)
mat_fdr <- res4$mat_fdr
mat_dist <- res4$mat_dist

hc <- hclust_tree(pinmat, mat_fdr, mat_dist, hc_method = "average")
hc_clone <- find_clone(hc, fdr_thresh = -2, 
		share_min = 0.85, n_share = 3, 
		bymax = TRUE, 
		climb_from_size = 2, climb_to_share = 3, graphic = F)

sub_hc_clone <- find_subclone(
		hc_clone, pinmat, pins, min_node_size = 6, sim_round = 500, 
		lm_max = 0.001, hc_method = "average", base_share = 3, fdr_thresh = -2, 
		share_min = 0.90, 
		bymax = TRUE, climb_from_size = 2, climb_to_share = 3, graphic = F)

output_viewer(output_file_dir = output_dirname, 
		seg.quantal, ratio.quantal, pins, 
		pinmat, mat_dist, hc_clone, sub_hc_clone, 
		subcloneTooBig = 0.8, smear = 1, study=study_name,
		cell_names = cell_names)