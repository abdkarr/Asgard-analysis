library(here)
library(Seurat)
library(SeuratObject)

library(devtools)

load_all(here("..", "Asgard"))

OUTPUT_DIR <- here("leukemia", "outputs")

# Set FALSE if you don't want to regress cell cycle
do_cell_cycle_regression = TRUE

# Load Data
sobjects <- load_luekemia()

# Pre-processing
for (i in seq_along(sobjects)) {
    sobjects[[i]] <- NormalizeData(sobjects[[i]], verbose=FALSE)
    sobjects[[i]] <- FindVariableFeatures(sobjects[[i]], nfeatures=2000, 
                                          verbose=FALSE)
}

# Integration
anchors <- FindIntegrationAnchors(sobjects, anchor.features = 2000, dims = 1:15)
integrated <- IntegrateData(anchors, dims = 1:15)

if (do_cell_cycle_regression) {
    integrated <- CellCycleScoring(integrated, s.features = cc.genes$s.genes, 
                                   g2m.features = cc.genes$g2m.genes, 
                                   set.ident = TRUE)
    integrated <- ScaleData(integrated, vars.to.regress = c("S.Score", "G2M.Score"), 
                            features = rownames(integrated))
    integrated <- RunPCA(integrated, npcs = 15, verbose = FALSE)
} else { 
    integrated <- ScaleData(integrated, verbose = FALSE)
}

# PCA
integrated <- Seurat::RunPCA(integrated, npcs = 15, verbose = FALSE)

# Save PCA: In the downstream analysis, we only need PCA of integrated data
# So we are going to save PCA results
save_dir <- file.path(OUTPUT_DIR, "integration")
dir.create(save_dir, recursive=TRUE)

if (do_cell_cycle_regression) {
    pca_file <- file.path(save_dir, "seurat_pca_cell_cycle_regressed.csv")
} else {
    pca_file <- file.path(save_dir, "seurat_pca.csv")
}

X_pca = Embeddings(integrated, reduction="pca")
write.table(X_pca, pca_file, quote=FALSE, sep=",", col.names=FALSE)