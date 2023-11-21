library(here)
library(Seurat)
library(SeuratObject)
library(SingleR)
library(celldex)
library(Asgard)

OUTPUT_DIR <- here("leukemia", "outputs")

# Load Data
sobjects <- load_luekemia()

hpca <- celldex::HumanPrimaryCellAtlasData()

for (i in seq_along(sobjects)) {
    # Log-normalization
    sobjects[[i]] <- Seurat::NormalizeData(sobjects[[i]], verbose=FALSE)

    # Convert to SCE
    sce <- Seurat::as.SingleCellExperiment(sobjects[[i]])

    # Annotate
    cell_types <- SingleR::SingleR(test=sce, ref=hpca, labels = hpca$label.main)
    cell_types <- data.frame(
        row.names = row.names(cell_types), celltype = cell_types$labels
    )

    # Save annotations
    sample_name <- names(sobjects)[i]
    save_dir <- file.path(OUTPUT_DIR, "annotations")
    dir.create(save_dir, recursive=TRUE)
    fname <- file.path(save_dir,  sprintf("%s_cell_types.csv", sample_name))
    write.table(cell_types, fname, quote=FALSE, sep=",")
}