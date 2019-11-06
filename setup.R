full <- read.csv("scrna_ather_g4_count_tmm.csv")
counts <- t(full[,2:ncol(full)])
rownames(counts) <- sapply(strsplit(rownames(counts), split = '_'), '[', 1)
colnames(counts) <- as.character(full[,1])
counts <- counts - rowMeans(counts)
counts <- round(counts,3)

bcode_info <- readxl::read_xlsx("42255_2019_102_MOESM5_ESM.xlsx",skip=3) %>% 
  mutate(VAE_ClusterLabel = factor(VAE_ClusterLabel))

library(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl")
anno <- getBM(attributes = c("mgi_symbol","description",
                             "ensembl_gene_id"),
              filters = "ensembl_gene_id",
              values = rownames(counts),mart=mart) %>% 
  separate(description, into=c("desc","meta"),sep="\\[") %>% 
  dplyr::select(-meta)
saveRDS(anno, file="anno.rds")
saveRDS(counts, file="counts.rds")
saveRDS(bcode_info, file="barcodes.rds")
