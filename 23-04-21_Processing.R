
###############################################################################

# Building the processing pipeline from the beginning

###############################################################################

# Depending on the machine set the input and output directories

# Input directory (where are h5 files or hc4 directory located)
#input <- "./0_data/cov19_raw/" # Windows
input <- "~/data/cov19_raw/"    # Linux


# Output directory (where should the plots and the processed data go to?)
output <- paste("~/3_r_projects/Anders_Lab_Local/", format(Sys.time(), "%d-%m-%y"),
                "_processing/", sep="")
#output <- paste("~/pCloudDrive/3_r_projects/1-Anders_Lab/6_linked_charts/", 
#                format(Sys.time(), "%d-%m-%y"), "_processing/", sep="")
ifelse(!dir.exists(output), dir.create(output), FALSE)

output_data <- paste(output, "data/", sep="")
output_graphics <- paste(output, "graphics/", sep="")
ifelse(!dir.exists(output_data), dir.create(output_data), FALSE)
ifelse(!dir.exists(output_graphics), dir.create(output_graphics), FALSE)

###############################################################################

# Define global settings

# For the graphics
dpi = 300

# TODO
# What are global options one could vary, e.g. how the PCA space is constructed?

###############################################################################

# Load packages

library( tidyverse )
library( gridExtra )
library( Matrix )
library( sparseMatrixStats )

###############################################################################

# Get meta data

meta_dat <- read_tsv( url( "https://raw.githubusercontent.com/zhangzlab/covid_balf/master/meta.txt" ) )

# Select threshold used by the authors
thresholds <- unique( meta_dat [, 7:10] )

# Select sample info (sample ID, disease status -> HC = healthy control, M = moderate, S = severe)
sample_info <- select( meta_dat, sample = sample_new, sample_c_id = sample, 
                       group, disease )

# Assuming that everything was downloaded accordingly from GEO database and 
# already extracted and stored in the data folder.
# HC1 - S6 (without HC4): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145926
# HC4 (from another study): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3660650
files <- list.files(input)
sample_info %>% rowwise %>% mutate(file_name = paste(input, str_trim(files[grep(
  paste("[[:punct:][:alnum:]]*", sample_c_id, "[[:punct:][:alnum:]]*", "matrix", sep=""),
  files)]), sep="")) -> sample_info

###############################################################################

# Read in the count matrices

nCov <- list()
map(sample_info$file_name, function(file){
  print(file)
  if (str_ends(file, ".h5")){
    Seurat::Read10X_h5(file)
  }
  else {
    Seurat::Read10X(data.dir="./0_data/cov19_raw/hc4/")
  }
}) -> raw_counts

names(raw_counts) <- sample_info$sample

# Add nCoV to all samples which do not have the entry yet

for( i in 1:length( raw_counts ) ) {
  if( nrow( raw_counts[[i]] ) == 33538 )
    raw_counts[[i]] <- rbind( raw_counts[[i]], nCoV=0 )
}

map_lgl( raw_counts[-1], ~ identical( rownames(.), rownames( raw_counts[[1]] ) ) ) %>% 
  all %>% stopifnot

###############################################################################

# Make the Kneeplot

raw_counts %>% map(function(x){colSums(x != 0)}) %>% unlist -> unique_feat

png(file=paste(output_graphics, "1_kneeplot.png", sep=""), 
    width = dpi*8,height = dpi*6, units = "px", res = dpi, type='cairo')

ggplot() +
  geom_point( aes( x=1:length( unique_feat ), y=sort( unique_feat,decreasing=TRUE ), 
                   stroke=0.05 ) ) +
  scale_x_log10() + scale_y_log10() +
  geom_hline( yintercept = 1000 ) +
  labs( x = "Ranked non-empty droplet", y = "Number of detected feature" )

dev.off()

###############################################################################

# QC: Check features
# Exclude all features which are not present in at least 3 cells

# First visualize the distribution of the features 

x_min <- 0
x_max <- max(log10(unlist(map(raw_counts, rowSums))))

make_plots_1 = function(sample) {
  ggplot(data.frame(feat=log10(rowSums(raw_counts[[sample]])))) +
    geom_histogram(aes(x=feat), bins=100, color="black", fill="blue") + 
    labs(title=sample, x="Log10 Feature Sum") + lims(x=c(x_min, x_max))
}

plots <- comprehenr::to_list(for (smp in sample_info$sample) make_plots_1(smp))

png(file=paste(output_graphics, "2_sum_feat_before.png", sep=""), 
    width = dpi*12,height = dpi*10, units = "px", res = dpi, type='cairo')
do.call("grid.arrange", c(plots, ncol=4))
dev.off()

bef <- dim(raw_counts$HC1)[1]

raw_counts %>% map( function( x ) {rowSums( x != 0 ) >= 3 } ) -> above
which( Reduce( "&", above ) ) -> above_in_all

# Save gene_names for later and initialize vector to store reduced count matrices
gene_names <- rownames( raw_counts$HC1[above_in_all,] )
red_counts <- vector( mode="list", length = length( raw_counts ) )

# Use for loop to subset all the rows form the raw counts that fulfill the requirement
for(i in 1:length(raw_counts)) red_counts[[i]] <- raw_counts[[i]][above_in_all, ]
names(red_counts) <- names(raw_counts)
raw_counts <- red_counts
rm(red_counts)

# Save so we can compare later.
aft <- dim(raw_counts$HC1)[1]

# Visualize the features again

x_min <- min(log10(unlist(map(raw_counts, rowSums))))
x_max <- max(log10(unlist(map(raw_counts, rowSums))))

plots <- comprehenr::to_list(for (smp in sample_info$sample) make_plots_1(smp))

png(file=paste(output_graphics, "3_sum_feat_after.png", sep=""), 
    width = dpi*12,height = dpi*10, units = "px", res = dpi, type='cairo')
do.call("grid.arrange", c(plots, ncol=4))
dev.off()

###############################################################################

# QC: Check droplets (cells)

# 1. Number of detected features per cell per sample

png(file=paste(output_graphics, "4_features_before.png", sep=""), 
    width = dpi*8,height = dpi*6, units = "px", res = dpi, type='cairo')

raw_counts %>% map(function(x){colSums(x != 0)}) -> unique_feat

unique_feat %>% unlist()%>% 
  data.frame() %>% rownames_to_column() %>% rowwise() %>%
  mutate(sample = str_split(rowname, "\\.")[[1]][1]) %>%
  rename(count = ".") %>%
  ggplot(., aes(x=sample, y=count)) +
  geom_boxplot() +
  labs(x = "Sample", y = "Feature count", 
       title = "Number of detected features per cell before QC")

dev.off()

# 2. Number of counted UMIs per cell per sample

png(file=paste(output_graphics, "5_UMIs_before.png", sep=""), 
    width = dpi*8,height = dpi*6, units = "px", res = dpi, type='cairo')

raw_counts %>% map(function(x){colSums(x)}) -> umi_counts

umi_counts %>% unlist()%>% 
  data.frame() %>% rownames_to_column() %>% rowwise() %>%
  mutate(sample = str_split(rowname, "\\.")[[1]][1]) %>%
  rename(count = ".") %>%
  ggplot(., aes(x=sample, y=count)) +
  geom_boxplot() +
  labs(x = "Sample", y = "UMI count", 
       title = "Number of counted UMIs per cell before QC")

dev.off()

# 3. Percentage of mitochondrial RNA

png(file=paste(output_graphics, "6_MT_before.png", sep=""), 
    width = dpi*8,height = dpi*6, units = "px", res = dpi, type='cairo')

# Assuming order of features is the same in every cell, which is indeed true
raw_counts %>% map(function(x){rownames(x) %>% grep("^MT-", ., ignore.case = FALSE)}) -> mito_index

raw_counts %>% map(function(x){colSums(x[mito_index$HC1, , drop=FALSE]) / colSums(x)*100}) -> mito_frac

mito_frac %>% unlist()%>% 
  data.frame() %>% rownames_to_column() %>% rowwise() %>%
  mutate(sample = str_split(rowname, "\\.")[[1]][1]) %>%
  rename(count = ".") %>%
  ggplot(., aes(x=sample, y=count)) +
  geom_boxplot() +
  labs(x = "Sample", y = "Mito percentage", 
       title = "Mitochondrial percentage before QC")

dev.off()

# Combine counts and add sample same to barcode, so each cell label is unique

counts <- raw_counts
rm(raw_counts) # save memory

for( name in names(counts) )
  colnames(counts[[name]]) %>% 
  str_remove( "-1" ) %>% 
  str_c( name, "-", . ) -> colnames(counts[[name]])

counts %>% do.call( cbind, . ) -> counts  

# Subset cells according to the thresholds used by the autors

# Checking the distribution of the cells over the samples before subsetting
table(unlist(str_split(colnames(counts), "-"))[seq(1, dim(counts)[2]*2, 2)]) %>%
  data.frame() %>% t() -> before

# Q: Use above or equal to ro above?
counts <- counts[, colSums(counts) > thresholds$nCount_RNA_threshold & 
                   colSums(counts!=0) > thresholds$nFeature_RNA_low &
                   colSums(counts!=0) < thresholds$nFeature_RNA_high & 
                   unlist(mito_frac) < thresholds$percent.mito]

table(unlist(str_split(colnames(counts), "-"))[seq(1, dim(counts)[2]*2, 2)])  %>%
  data.frame() %>% t() -> after

summary <- rbind(before, after[2,])
colnames(summary) <- summary[1,]
summary <- summary[-1,]
rownames(summary) <- c("Before", "After")
sample_info <- cbind(sample_info, t(summary))

# Check the QC again

# 1. Number of detected features per cell per sample

png(file=paste(output_graphics, "7_features_after.png", sep=""), 
    width = dpi*8,height = dpi*6, units = "px", res = dpi, type='cairo')

data.frame(barcode = colnames(counts), feature = colSums(counts!=0)) %>%
  rowwise() %>% mutate(sample = str_split(barcode, "-")[[1]][1]) %>%
  ggplot(., aes(x = sample, y = feature)) + geom_boxplot() +
  labs(x = "Sample", y = "Feature count", 
       title = "Number of detected features per cell after QC")

dev.off()

# 2. Number of counted UMIs per cell per sample

png(file=paste(output_graphics, "8_UMIs_after.png", sep=""), 
    width = dpi*8,height = dpi*6, units = "px", res = dpi, type='cairo')

data.frame(barcode = colnames(counts), feature = colSums(counts)) %>%
  rowwise() %>% mutate(sample = str_split(barcode, "-")[[1]][1]) %>%
  ggplot(., aes(x = sample, y = feature)) + geom_boxplot() +
  labs(x = "Sample", y = "UMI count", 
       title = "Number of counted UMIs per cell after QC")

dev.off()

# 3. Percentage of mitochondrial RNA

png(file=paste(output_graphics, "9_MT_after.png", sep=""), 
    width = dpi*8,height = dpi*6, units = "px", res = dpi, type='cairo')

rownames(counts) %>% grep("^MT-", ., ignore.case = FALSE) -> mito_index

data.frame(barcode = colnames(counts), 
           feature = (colSums(counts[mito_index,])/colSums(counts))*100) %>%
  rowwise() %>% mutate(sample = str_split(barcode, "-")[[1]][1]) %>%
  ggplot(., aes(x = sample, y = feature)) + geom_boxplot() +
  labs(x = "Sample", y = "Mito percentage", 
       title = "Mitochondrial percentage after QC")

dev.off()

###############################################################################

# Calculate the fraction -> Scale Normalization (accounting for the library size)

t( t(counts) / colSums(counts) ) -> fracs

###############################################################################

# Create dataframe to store metadata for each cell

tibble( colname = colnames( fracs ) ) %>%
  mutate( sample = str_extract( colname, "^[^-]*") ) %>%
  mutate( disease = ifelse( str_starts( colname, "HC"), FALSE, TRUE ) )-> cells
cells$total <- colSums(counts)

###############################################################################

# Save the data

save(counts, file = paste(output_data,"counts.RData", sep=""), compress = FALSE)
save(fracs, file = paste(output_data,"fracs.RData", sep=""), compress = FALSE)
save(cells, file = paste(output_data,"cells.RData", sep=""), compress = FALSE)

###############################################################################

# Calculate gene means and gene variance for each gene in each sample

sample_info$sample %>% 
  map(function(x)rowMeans( fracs[ , str_starts( colnames(fracs), x ) ] ) ) %>%
  do.call( cbind, . ) -> gene_means
colnames(gene_means) <- sample_info$sample

sample_info$sample %>% 
  map(function(x)rowVars( fracs[ , str_starts( colnames(fracs), x ) ] ) ) %>%
  do.call( cbind, . ) -> gene_vars
rownames(gene_vars) <- rownames(fracs)
colnames(gene_vars) <- sample_info$sample

###############################################################################

# Select highly variable gene

# Variance stabilizing transformation: Dividing variance by the mean 
# assuming Poisson noise (were variance = expect value)

# Set the threshold
# Why this threshold -> Returning about 1500 highly variable genes

thr <- 0.0006

hvg_per_smp <- comprehenr::to_vec(for (i in 1:dim(gene_means)[2]) 
  sum( (gene_means[,i] > 0) & gene_vars[,i]/gene_means[,i] > thr ))

sample_info <- cbind(sample_info, hvg_per_smp)

png(file=paste(output_graphics, "10_hvg.png", sep=""), 
    width = dpi*16,height = dpi*12, units = "px", res = dpi, type='cairo')

# Join the means and the variance, both extented (pivot_longer) according ot the sample
inner_join(
  gene_means %>% 
    as_tibble( rownames="gene" ) %>% 
    pivot_longer( -gene, names_to="sample", values_to="mean" ), 
  gene_vars %>% 
    as_tibble( rownames="gene" ) %>% 
    pivot_longer( -gene, names_to="sample", values_to="var" ) ) %>%
  ggplot(.) + 
  geom_point( aes( x=mean, y=var/mean ), size=.1, alpha=.2 ) + facet_wrap( ~ sample ) + 
  scale_x_log10() + scale_y_log10() + geom_hline( yintercept = thr, col="blue" ) +
  labs(x="Mean", y="Variance/Mean")

dev.off()

# Select genes, which are above the threshold in at least two cells
table( rowSums( gene_means > 0 & gene_vars / gene_means > thr ) )

( rowSums( gene_means > 0 & gene_vars / gene_means > thr ) > 1 ) %>% 
  which %>% names -> hvg

###############################################################################

# Linear dimensional reduction -> PCA

# Rationale of doing PCA on principal components form sample-sepcific PCA
# -> Do not want to capture batch effects in the first PCs

# For each sample do a PCA and retain the first 40 principal components 
# (we center the data, but do not scale them)

# We stabilize the variance of the hv genes by taking the hyperbolic sine of the
# square root of the fraction scaled a scaling factor of 10,000
sc_factor <- 1e4

# Using the irlba package so set a seed for reproducible results
set.seed(100)

map( sample_info$sample, function(sm) {
  cat( sm, " ")
  asinh( sqrt( sc_factor * fracs[ hvg, cells$sample == sm ] ) ) %>% 
    t %>% irlba::prcomp_irlba( n=40, retx=TRUE, center=TRUE, scale.=FALSE ) 
} ) -> pcs

names(pcs) <- sample_info$sample

# Scale the principal components from each sample by their variance (square root)
# Bind all the principal components and use as input for another PCA
# All PCs are defined in the feature space spanned by the hv genes

map(pcs, function(x) t( x$rotation %*% diag( x$sdev ) ) ) %>% 
  do.call( rbind, . ) %>%
  prcomp() -> pc

# Visualize scree plot

png(file=paste(output_graphics, "11_scree_plot.png", sep=""), 
    width = dpi*8,height = dpi*6, units = "px", res = dpi, type='cairo')

ggplot() + geom_point(aes(x=1:40, y=pc$sdev[1:40]^2))

dev.off()

# Map all the cells into the mutual pc space using the first 40 principal components
# (which are again defined in the feature space, basically projecting the data onto the components)

map( sample_info$sample, function(sm) {
  list( scale=FALSE, center=pcs[[sm]]$center, sdev=pc$sdev[1:40], 
        rotation=pc$rotation[,1:40] ) -> a
  class(a) <- "prcomp"
  predict( a , t( asinh( sqrt( sc_factor * fracs[ hvg, cells$sample == sm ] ) ) ) ) } ) -> a
x <- do.call( rbind, a )

# Regress out the correlation between the PC coordinates and the library size
# Why does this correlation exist? -> Because especially for low expressed
# genes there is a high correlation of the library size and the fraction

# Correlation before (for first 9 principal components)

plots <- vector("list", length = 9)
for (i in 1:9) {
  cf_i <- coefficients(lm(x[, i] ~ log(colSums(counts))))
  p1 <- ggplot() +
    geom_point(aes(x = log(colSums(counts)), y = x[, i]), size=.2) +
    geom_line(aes(x = log(colSums(counts)), y = log(colSums(counts))*cf_i[2]+cf_i[1]), 
              color="red") +
    labs(title=paste("Slope: ", cf_i[2]))
  plots[[i]] <- p1
}

png(file=paste(output_graphics, "12_PC_correlation_before.png", sep=""), 
    width = dpi*10,height = dpi*10, units = "px", res = dpi, type='cairo')

print(do.call("grid.arrange", c(plots, ncol=3)))

dev.off()

# Regress out

x %>% apply( 2, function(x) residuals( lm( x ~ log(colSums( counts ) ) ) ) ) -> pcx

# Correlation after

plots <- vector("list", length = 9)
for (i in 1:9) {
  cf_i <- coefficients(lm(pcx[, i] ~ log(colSums(counts))))
  p1 <- ggplot() +
    geom_point(aes(x = log(colSums(counts)), y = pcx[, i]), size=.2) +
    geom_line(aes(x = log(colSums(counts)), y = log(colSums(counts))*cf_i[2]+cf_i[1]), 
              color="red") +
    labs(title=paste("Slope: ", cf_i[2]))
  plots[[i]] <- p1
}

png(file=paste(output_graphics, "13_PC_correlation_after.png", sep=""), 
    width = dpi*10,height = dpi*10, units = "px", res = dpi, type='cairo')

do.call("grid.arrange", c(plots, ncol=3))

dev.off()

# Do the same for the sample specific pcs space (not visualizing the results)
# We already have the coordinates of the cells mapped into the sample specific PC space

map(pcs, function(pc) pc$x) -> coordinates
x <- do.call( rbind, coordinates )
x %>% apply( 2, function(x) residuals( lm( x ~ log(colSums( counts ) ) ) ) ) -> pcs

###############################################################################

# Validate that the sample-specific pcs space is similar to the mutual pcx space

# Compare the distances of 10000 randomly choosen pairs of cells measured first 
# in the pcs and in the pcx space, and look how well they are correlated

make_plots_2 = function(sample_num, limit=40) {
  
  sample_str <- sample_info$sample[sample_num]
  
  a1 <- sample.int(sum(cells$sample==sample_str), 10000, replace=TRUE)
  a2 <- sample.int(sum(cells$sample==sample_str), 10000, replace=TRUE)
  
  x <- sqrt(rowSums( 
    ( pcs[which(cells$sample==sample_str)[a1],] - 
        pcs[which(cells$sample==sample_str)[a2],] )^2 ))
  
  y <- sqrt(rowSums( 
    ( pcx[which(cells$sample==sample_str)[a1],] - 
        pcx[which(cells$sample==sample_str)[a2],] )^2 ))
  
  coeff <- round(lm(y ~ x)[[1]][2],3)
  
  plt <- ggplot(data.frame(x, y), aes(x=x, y=y)) + geom_point(cex = .2) +
    labs(x="PCS (sample-specific PC space)", y="PCX (mututal PC space)",
         title=paste(sample_str, "|", coeff)) +
    lims(x=c(0,limit), y=c(0,limit)) + geom_abline(intercept=0, slope=1, col="red")
  
  return(plt)
}

plots <- comprehenr::to_list(for (i in 1:dim(sample_info)[1]) make_plots_2(i))

png(file=paste(output_graphics, "14_PCS_PCX_comparison.png", sep=""), 
    width = dpi*12,height = dpi*12, units = "px", res = dpi, type='cairo')

do.call("grid.arrange", c(plots, ncol=4))

dev.off()

# Visualize again just looking at the interesting region (distance up to 12)

plots <- comprehenr::to_list(for (i in 1:dim(sample_info)[1]) make_plots_2(i, limit=12))

png(file=paste(output_graphics, "15_PCS_PCX_comparison_zoomed.png", sep=""), 
    width = dpi*12,height = dpi*12, units = "px", res = dpi, type='cairo')

do.call("grid.arrange", c(plots, ncol=4))

dev.off()

# Q: Is it a problem that there is not perfect correlation in all samples?

png(file=paste(output_graphics, "16_PCX_samples.png", sep=""), 
    width = dpi*12,height = dpi*12, units = "px", res = dpi, type='cairo')

ggplot(data.frame(x=pcx[,1], y=pcx[,2]), aes(x=x, y=y)) +
  geom_point(aes(col=factor(cells$sample)))

dev.off()

###############################################################################

# Calculate 2D embeddings

# Importantly, both tSNE and UMAP must be initiualized by PCA

# tSNE PCX
set.seed(100)

map( sample_info$sample, function(sm) {
  cat( sm, " " )
  Rtsne::Rtsne( pcx[ cells$sample==sm, ], pca=FALSE, dims=2, 
                pca_center=TRUE, pca_scale=FALSE, normalize=TRUE)
} ) -> tsnes
names(tsnes) <- sample_info$sample

# Combine all coordinates of the tsnes into one data frame
tsnes %>% map( ~ .$Y ) %>% do.call( rbind, . ) -> tsne_pcx
rownames(tsne_pcx) <- cells$colname
colnames(tsne_pcx) <- c( "tsne1", "tsne2")

png(file=paste(output_graphics, "17_tSNE_pcx.png", sep=""), 
    width = dpi*10,height = dpi*10, units = "px", res = dpi, type='cairo')

cells %>% cbind(tsne_pcx) %>% select(sample, tsne1, tsne2) %>%
  ggplot(.) + geom_point(aes(x=tsne1, y=tsne2), size=.2) + facet_wrap(~ sample)

dev.off()

# tSNE PCS
set.seed(100)

map( sample_info$sample, function(sm) {
  cat( sm, " " )
  Rtsne::Rtsne( pcs[ cells$sample==sm, ], pca=FALSE, dims=2, 
                pca_center=TRUE, pca_scale=FALSE, normalize=TRUE)
} ) -> tsnes
names(tsnes) <- sample_info$sample

# Combine all coordinates of the tsnes into one data frame
tsnes %>% map( ~ .$Y ) %>% do.call( rbind, . ) -> tsne_pcs
rownames(tsne_pcs) <- cells$colname
colnames(tsne_pcs) <- c( "tsne1", "tsne2")

png(file=paste(output_graphics, "18_tSNE_pcs.png", sep=""), 
    width = dpi*10,height = dpi*10, units = "px", res = dpi, type='cairo')

cells %>% cbind(tsne_pcs) %>% select(sample, tsne1, tsne2) %>%
  ggplot(.) + geom_point(aes(x=tsne1, y=tsne2), size=.2) + facet_wrap(~ sample)

dev.off()

# UMAP PCX
set.seed(100)

map( sample_info$sample, function(sm) {
  cat( sm, " " )
  # Setting pca to 40, should ensure that no PCA is performed (because we have 40 PCs)
  uwot::umap( pcx[ cells$sample==sm, ], n_neighbors=15, n_components=2,
              nn_method="annoy", scale=FALSE, pca=40 ) } ) -> umaps
names(umaps) <- sample_info$sample

# Combine all coordinates of the umaps into one data frame
umaps %>% map(function(x)x) %>% do.call( rbind, . ) -> umap_pcx
rownames(umap_pcx) <- cells$colname
colnames(umap_pcx) <- c( "umap1", "umap2")

png(file=paste(output_graphics, "19_UMAP_pcx.png", sep=""), 
    width = dpi*10,height = dpi*10, units = "px", res = dpi, type='cairo')

cells %>% cbind(umap_pcx) %>% select(sample, umap1, umap2) %>%
  ggplot(.) + geom_point(aes(x=umap1, y=umap2), size=.2) + facet_wrap(~ sample)

dev.off()

# UMAP PCS
set.seed(100)

map( sample_info$sample, function(sm) {
  cat( sm, " " )
  # Setting pca to 40, should ensure that no PCA is performed (because we have 40 PCs)
  uwot::umap( pcs[ cells$sample==sm, ], n_neighbors=15, n_components=2,
              nn_method="annoy", scale=FALSE, pca=40 ) } ) -> umaps
names(umaps) <- sample_info$sample

# Combine all coordinates of the umaps into one data frame
umaps %>% map(function(x)x) %>% do.call( rbind, . ) -> umap_pcs
rownames(umap_pcs) <- cells$colname
colnames(umap_pcs) <- c( "umap1", "umap2")

png(file=paste(output_graphics, "20_UMAP_pcs.png", sep=""), 
    width = dpi*10,height = dpi*10, units = "px", res = dpi, type='cairo')

cells %>% cbind(umap_pcs) %>% select(sample, umap1, umap2) %>%
  ggplot(.) + geom_point(aes(x=umap1, y=umap2), size=.2) + facet_wrap(~ sample)

dev.off()

###############################################################################

# Calculate 100 nearest neighbors within each sample using the pcx space

sample_info$sample %>% map( function(sm) {
  cat( sm, " " )
  RANN::nn2( pcx[ cells$sample==sm, ], k=100 ) } ) -> a

# Adjust the indices

for( i in 2:length(a) ) {
  cat(names(a)[i], " ")
  a[[i]]$nn.idx <- a[[i]]$nn.idx + sum( map_int( a, ~nrow(.$nn.idx))[1:(i-1)] )
}

# Merge everything into one list by rbind.
list( 
  idx = ( map( a, ~.$nn.idx ) %>% do.call( rbind, . ) ),
  dist = ( map( a, ~.$nn.dist ) %>% do.call( rbind, . ) ) ) -> all_NN_pcx

# Stop if not first column of nn$idx is the same as 1:number of rows of nn$idx
stopifnot( identical( all_NN_pcx$idx[,1], 1:nrow(all_NN_pcx$idx) ) )


# Calculate 100 nearest neighbors within each sample using the pcs space

sample_info$sample %>% map( function(sm) {
  cat( sm, " " )
  RANN::nn2( pcs[ cells$sample==sm, ], k=100 ) } ) -> a

# Adjust the indices

for( i in 2:length(a) ) {
  cat(names(a)[i], " ")
  a[[i]]$nn.idx <- a[[i]]$nn.idx + sum( map_int( a, ~nrow(.$nn.idx))[1:(i-1)] )
}

# Merge everything into one list by rbind.
list( 
  idx = ( map( a, ~.$nn.idx ) %>% do.call( rbind, . ) ),
  dist = ( map( a, ~.$nn.dist ) %>% do.call( rbind, . ) ) ) -> all_NN_pcs

# Stop if not first column of nn$idx is the same as 1:number of rows of nn$idx
stopifnot( identical( all_NN_pcs$idx[,1], 1:nrow(all_NN_pcs$idx) ) )

###############################################################################

# Cumulative correcting factor (useful for later)

cumulative <- vector(mode = "numeric", length=dim(sample_info)[1])
for (i in 1:dim(sample_info)[1]) {
  print(i-1)
  cumulative[i] <- sum(as.numeric(sample_info$After)[0:(i-1)])
}
sample_info <- cbind(sample_info, cumulative)

###############################################################################

# Calculate the Jaccard Indices for all samples (could take very long this way)

# How to speed up the calculation
# Multithreading? -> https://furrr.futureverse.org/

# Option 1
calc_jaccard = function(dims, all_NN_subset) {
  
  # Initialize matrix
  mt <- matrix(nrow=dims, ncol=dims)
  
  # Speed it up by just doing half the matrix (which also contains all the information)
  for (i in 1:dims) {
    for (j in i:dims) {
      mt[i, j] <- sum(all_NN_subset$idx[i,] %in% all_NN_subset$idx[j,])
      # mt[i, j] <- length(union(all_NN_subset$idx[i,], all_NN_subset$idx[j,]))
    }
  }
  return(mt/100)
}

# Compile functions (does that even help -> TODO benchmark different methods)
calc_jaccard_cmd <- compiler::cmpfun(calc_jaccard)

# Wrapper for the function
wrapper_jaccard_pcx = function(sm) {
  
  # Progress
  cat( sm, " " )
  
  # Get the subset information
  subset <- cells$sample==sm
  dims <- sum(subset)
  correction <- sample_info$cumulative[grep(sm, sample_info$sample)]
  all_NN_subset <- list(idx = (all_NN_pcx$idx[subset,] - correction), 
                        dist = all_NN_pcx$dist[subset,])
  
  jac_i <- calc_jaccard_cmd(dims=dims, all_NN_subset=all_NN_subset)
  return(jac_i)
}

# Wrapper for the function
wrapper_jaccard_pcs = function(sm) {
  
  # Progress
  cat( sm, " " )
  
  # Get the subset information
  subset <- cells$sample==sm
  dims <- sum(subset)
  correction <- sample_info$cumulative[grep(sm, sample_info$sample)]
  all_NN_subset <- list(idx = (all_NN_pcs$idx[subset,] - correction), 
                        dist = all_NN_pcs$dist[subset,])
  
  jac_i <- calc_jaccard_cmd(dims=dims, all_NN_subset=all_NN_subset)
  return(jac_i)
}

# Using furrr
#library(furrr)
#plan(multisession, workers = 13)

# PCX space

sample_info$sample %>% map(wrapper_jaccard_pcx) -> jaccard_pcx

names(jaccard_pcx) <- sample_info$sample

# PCS space

sample_info$sample %>% map(wrapper_jaccard_pcs) -> jaccard_pcs

names(jaccard_pcs) <- sample_info$sample

###############################################################################

# Process the annotation data

ref <- read.table("./3_covid_balf/all.cell.annotation.meta.txt", header=TRUE) %>% 
  mutate(colname = str_c(sample_new, "-", str_replace(ID, "_\\d*", ""))) %>%
  rename(sample_old = sample)

# Join the ref with our cell data frame
cells %>% left_join(ref, by="colname") -> cells_ref

# Build dictionary mapping genes to celltypes
celltype_marker <- c("Epithelial", "Epithelial", "Macrophages", "Neutrophil", 
                     "mDC", "pDC", "Mast", "T", "NK", "B", "Plasma")
names(celltype_marker) <- c("TPPP3", "KRT18", "CD68", "FCGR3B", "CD1C",
                            "LILRA4", "TPSB2", "CD3D", "KLRD1", "MS4A1", "IGHG4")

###############################################################################

# Save the data

save(hvg, file = paste(output_data,"hvg.RData", sep=""), compress = FALSE)
save(pcs, file = paste(output_data,"pcx.RData", sep=""), compress = FALSE)
save(pcx, file = paste(output_data,"pcs.RData", sep=""), compress = FALSE)
save(tsne_pcx, file = paste(output_data,"tsne_pcx.RData", sep=""), compress = FALSE)
save(tsne_pcs, file = paste(output_data,"tsne_pcs.RData", sep=""), compress = FALSE)
save(umap_pcx, file = paste(output_data,"umap_pcx.RData", sep=""), compress = FALSE)
save(umap_pcs, file = paste(output_data,"umap_pcs.RData", sep=""), compress = FALSE)
save(all_NN_pcx, file = paste(output_data,"all_NN_pcx.RData", sep=""), compress = FALSE)
save(all_NN_pcs, file = paste(output_data,"all_NN_pcs.RData", sep=""), compress = FALSE)
save(jaccard_pcx, file = paste(output_data,"jaccard_pcx.RData", sep=""), compress = FALSE)
save(jaccard_pcs, file = paste(output_data,"jaccard_pcs.RData", sep=""), compress = FALSE)
save(sample_info, file = paste(output_data, "sample_info.RData", sep=""), compress=FALSE)
save(cells_ref, file = paste(output_data, "cells_ref.RData", sep=""), compress=FALSE)
save(celltype_marker, file = paste(output_data, "celltype_marker.RData", sep=""), compress=FALSE)

###############################################################################

# Saving meta data

cat(" Sample meta data: \n",file=paste(output, "info.txt", sep=""))
write.table(sample_info, paste(output, "info.txt", sep=""), quote=FALSE, 
            append=TRUE, row.names=FALSE)
cat("\n QC Thresholds: \n",file=paste(output, "info.txt", sep=""), append=TRUE)
write.table(thresholds, paste(output, "info.txt", sep=""), quote=FALSE, 
            append=TRUE, row.names=FALSE)
cat(paste("\n QC Features (in at least 3 cells): \n", 
          bef, "->", aft), 
    file=paste(output, "info.txt", sep=""), append=TRUE)

