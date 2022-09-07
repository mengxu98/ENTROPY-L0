

rm(list = ls())
gc()

library("readxl")
library("SCEnt")
library("Matrix")
library("tidyr")
library("Seurat")
library("L0Learn")
library("Metrics")
source("as_matrix_cpp.R")
source("entropy_SCEnt_function.R")
#
if (Sys.info()[1] == "Windows") {
  savepath <- "../results/"
} else if (Sys.info()[1] == "Darwin") {
  savepath <- "../results/"
} else {
  savepath <- "/data/mengxu/data/L0/"
}

if (TRUE) {
  if (FALSE) {
    # Use 200000 cells with harmony

    stage <- c("normal", "1", "2", "3", "4")
    seu_list <- list()
    for (i in 1:length(stage)) {
      s <- stage[i]
      load(paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_SCT_harmony_PCA.Rdata"))
      seu_list[[i]] <- seu_obj_data
    }
  } else {
    # Use 61271 cells with harmony

    load("../data/lung_L0_data_harmony_new.Rdata")
    celltype <- scRNA_harmony$new_celltype %>% levels()
    scRNA_harmony@assays$RNA <- c() #先去掉counts矩阵，保证内存
    #rm(scRNA_harmony)
    gc()
  }
} else {
  # Use 61271 cells but no harmony
  load(paste0("/data/mengxu/data/L0/lung_L0_data.Rdata"))
  seu_list <- Seurat::SplitObject(seu_obj, split.by = "stage")
  rm(seu_obj)
  gc()
}

#counts_L0_seu_list <- list()
L0_information <- c()

for (c in 1:length(celltype)) {
  
  #load("../data/lung_L0_data_harmony_new.Rdata")
  scRNA_harmony_sel <- subset(x = scRNA_harmony, subset = celltype == celltype[c])
  seu_list <- SplitObject(scRNA_harmony_sel, split.by = "stage")
  rm(scRNA_harmony_sel)
  #rm(scRNA_harmony)
  gc()

  samples_combine <- c()
  counts_combine <- c()

  for (i in 1:length(seu_list)) {
    
    message("\n----- Now process the data of ", "stage-", i - 1, " Cell: ", celltype[c], "! -----\n")
    
    matrix <- seu_list[[i]]@assays$SCT@data
    
    # if (length(seu_list[[i]]@assays) == 1) {
    #   # 修改，如果没有进行harmony，则进行harmony，同时使用SCT
    # 
    #   matrix <- seu_list[[i]]@assays$RNA@counts
    #   # matrix <- seu_list[[i]]@assays$RNA@data %>% as.matrix() %>% as.data.frame()
    # } else {
    # 
    #   # matrix <- seu_list[[i]]@assays$SCT@counts %>% as.matrix() %>% as.data.frame()
    #   matrix <- seu_list[[i]]@assays$SCT@data
    # }
    
    if (TRUE) {
      counts <- as.matrix(matrix)
      # counts <- matrix
      # counts <- as.matrix(matrix[,1:200])
      rm(matrix)
      gc()

      if (FALSE) {
        max_counts <- max(counts)

        # if(max_counts > 10000){
        #   counts=log(counts+1)/log(10)
        # }else{
        #   counts=log(counts+1)/log(2)
        # }

        counts <- count_filter(counts)
        # counts <- scale(counts) #bad
        
      }

      if (TRUE) {

        # remove ERCC genes and Ribosomal genes
        ercc.genes <- grep("^ERCC", rownames(counts), value = TRUE)
        rb.genes <- grep("^RPL|^RPS|^MRPL|^MRPS", rownames(counts), value = TRUE)
        mt.genes <- grep("^MT", rownames(counts), value = TRUE)
        hb.genes <- grep("^HBA|^HBB", rownames(counts), value = TRUE)

        counts <- counts[which(!(rownames(counts) %in%
                                   c(ercc.genes, rb.genes, mt.genes, hb.genes))), ]
        
      }

      samples <- colnames(counts) %>% as.data.frame()
      samples$label <- i - 1
      # Combine -----------------------------------------------------------------
      samples_combine <- rbind(
        samples_combine,
        samples
      )

      if (is.null(counts_combine)) {
        counts_combine <- rbind(
          counts_combine,
          t(counts)
        )
      } else {
        gene_inter <- intersect(rownames(counts), colnames(counts_combine))
        counts <- counts[gene_inter, ]
        counts_combine <- counts_combine[, gene_inter]
        counts_combine <- rbind(
          counts_combine,
          t(counts)
        )
      }
      
    }

    entropy <- compute_entropy(counts,
                               unit = "log2",
                               method_entropy = "KL.plugin", 
                               # Now method can select one of KL.plugin(Default), ChaoShen, Dirichlet or entropy
                               compute_var = TRUE,
                               normalise = TRUE,
                               entropy_type = "Heterogeneity", # or Homogeneity
                               trans = FALSE,
                               select_num = 2000
    )
    # entropy <- Entropy(counts)
    write.csv(
      entropy,
      paste(savepath, "Entropy_stage-", i - 1, ".csv", sep = "")
    )
    
  }
  
  # Process L0 data ---------------------------------------------------------
  if (TRUE) {
    rm(entropy)
    rm(samples)
    rm(seu_list)
    rm(counts)
    gc()

    counts_raw <- as.matrix(counts_combine)

    row.names(counts_raw) <- row.names(counts_combine)

    Y_label <- samples_combine$label %>%
      as.numeric() %>%
      as.vector()
    names(Y_label) <- samples_combine$.

    entropy_normal <- read.csv(paste0(savepath, "Entropy_stage-0.csv"),
                               row.names = 1
    )
    genes_normal <- entropy_normal$gene %>% as.data.frame()
    entropy_stage1 <- read.csv(paste0(savepath, "Entropy_stage-1.csv"),
                               row.names = 1
    )
    genes_stage1 <- entropy_stage1$gene %>% as.data.frame()
    entropy_stage2 <- read.csv(paste0(savepath, "Entropy_stage-2.csv"),
                               row.names = 1
    )
    genes_stage2 <- entropy_stage2$gene %>% as.data.frame()
    entropy_stage3 <- read.csv(paste0(savepath, "Entropy_stage-3.csv"),
                               row.names = 1
    )
    genes_stage3 <- entropy_stage3$gene %>% as.data.frame()
    entropy_stage4 <- read.csv(paste0(savepath, "Entropy_stage-4.csv"),
                               row.names = 1
    )
    genes_stage4 <- entropy_stage4$gene %>% as.data.frame()
  }

  if (FALSE) {
    genes_union <- union(genes_normal$., genes_stage1$.)
    genes_union <- union(genes_union, genes_stage2$.)
    genes_union <- union(genes_union, genes_stage3$.)
    genes_union <- union(genes_union, genes_stage4$.)
    genes_union <- intersect(genes_union, colnames(counts_raw))

    Y_label <- Y_label[-which(Y_label == "4")]
    row.names(Y_label) <- row.names(Y_label)
  } else {
    genes_union <- union(genes_normal$., genes_stage1$.)
    genes_union <- union(genes_union, genes_stage2$.)
    genes_union <- union(genes_union, genes_stage3$.)
    genes_union <- union(genes_union, genes_stage4$.)
    genes_union <- intersect(genes_union, colnames(counts_raw))
  }
  
  rm(entropy_normal, entropy_stage1, entropy_stage2, entropy_stage3, entropy_stage4,
     genes_normal, genes_stage1, genes_stage2, genes_stage3, genes_stage4)
  
  count_L0 <- counts_raw[names(Y_label), genes_union]
  row.names(count_L0) <- row.names(count_L0)
  rm(counts_raw)
  gc()

  # L0 test-----------------------------------------------------------------------
  if (TRUE) {
    maxSNVSize <- 0.05 * min(ncol(count_L0), nrow(count_L0))
    penalty <- c("L0", "L0L1", "L0L2")
    algorithm <- c("CD", "CDPSI")

    p <- 1
    a <- 1
    cat(
      "\n----- Now run:",
      paste(penalty[p], algorithm[a], sep = "_"),
      "-----\n"
    )

    pre_time <- proc.time()
    if (TRUE) {
      fit_L0 <- L0Learn.fit(count_L0, Y_label,
                            penalty = penalty[p],
                            algorithm = algorithm[a],
                            maxSuppSize = maxSNVSize
      )
      # png(
      #   filename = paste0(sever_path, "fit_L0.png"),
      #   width = 1000,
      #   height = 600,
      #   units = "px"
      # )
      # plot(fit_L0, showLines = T)
      # dev.off()
      print(fit_L0)
    }
    if (FALSE) {
      cvfit_L0 <- L0Learn.cvfit(count_L0, Y_label,
                                penalty = penalty[p],
                                algorithm = algorithm[a],
                                nFolds = 10,
                                maxSuppSize = maxSNVSize
      )
      # png(
      #   filename = paste0(sever_path, "cvfit_L0.png"),
      #   width = 500,
      #   height = 400,
      #   units = "px"
      # )
      # plot(cvfit_L0, showLines = T)
      # dev.off()
      plot(cvfit_L0, showLines = T)
      print(cvfit_L0)
    }
    time <- proc.time() - pre_time

    # Extract coefficient at middle lambda
    fit_L0_information <- as.data.frame(print(fit_L0))
    fit_L0_information <- fit_L0_information[order(fit_L0_information$suppSize, decreasing = TRUE), ]
    lambda_L0 <- fit_L0_information$lambda[1]
    gamma_L0 <- fit_L0_information$gamma[1]

    pred_L0 <- predict(fit_L0,
                       newx = count_L0,
                       lambda = lambda_L0,
                       gamma = gamma_L0
    ) %>% as.vector()

    coef_L0 <- coef(fit_L0,
                    lambda = lambda_L0,
                    gamma = gamma_L0
    ) %>% as.vector()

    coef_L0 <- coef_L0[-1]
    coef_L0 <- which(coef_L0 != 0)
    coef_L0 <- colnames(count_L0)[coef_L0]

    if (length(coef_L0) == 1) {
      X_Y <- cbind(count_L0[, coef_L0], Y_label)
      colnames(X_Y)[1] <- coef_L0
      X_Y <- as.data.frame(X_Y)
    } else {
      X_Y <- cbind(count_L0[, coef_L0], Y_label) %>% as.data.frame()
    }
    lmfit <- lm(Y_label ~ ., data = X_Y)
    fit_coef_L0 <- summary(lmfit)
    ###
    write.csv(
      fit_coef_L0$coefficients,
      file <- paste0(
        paste0(savepath, "feature_selected_byL0_"),
        celltype[c],
        "_",
        penalty[p],
        algorithm[a],
        #"_",
        #format(Sys.time(), "%b%d_%H_%M_%S"),
        ".csv"
      )
    )

    ###
    R2_L0 <- MuMIn::r.squaredGLMM(lmfit)[1]
    res_rmse <- rmse(Y_label, pred_L0)
    res_rse <- rse(Y_label, pred_L0)
    Rsquare_L0 <- 1 - res_rse

    L0_information_mid <- c(
      celltype[c],
      paste(penalty[p], algorithm[a], sep = "_"),
      # R2_L0,
      Rsquare_L0,
      time[3],
      length(coef_L0)
    )
    L0_information <- rbind.data.frame(L0_information, L0_information_mid)
    names(L0_information) <- c(
      "CellType",
      "Method",
      # "R2",
      "Rsquare",
      "Time(s)",
      "NO._features"
    )
    cat("\n-----", L0_information_mid, "-----\n")
    cat("\n-----", paste(penalty[p], algorithm[a], sep = "_"), "final -----\n")

  }

}

#
load("../data/lung_L0_data_harmony_new.Rdata")

celltype <- scRNA_harmony$new_celltype %>% levels()

seu_list <- SplitObject(scRNA_harmony, split.by = "new_celltype")
seu_list2 <- list()

penalty <- c("L0", "L0L1", "L0L2")
algorithm <- c("CD", "CDPSI")

p <- 1
a <- 1

for (i in 1:length(seu_list)) {

  #scRNA_harmony_sel <- subset(x = scRNA_harmony, subset = celltype == celltype[c])

  genes <- read.csv(file = paste0(savepath, 
                                  "feature_selected_byL0_",
                                  celltype[i],
                                  "_",
                                  penalty[p],
                                  algorithm[a],
                                  ".csv"
                                  )
                    )
  
  genes <- genes$X
  seu <- seu_list[[i]]

  seu <- seu[genes,]
  
  seu_list2[[i]] <- seu
  
  # counts_L0_seu <- CreateSeuratObject(
  #   counts = mat,
  #   project = "Entropy-L0",
  #   min.features = 200,
  #   min.cells = 3
  # )
  # counts_L0_seu$celltype <- celltype[c]
  # counts_L0_seu$stage <- paste0("stage-", i - 1)
  # counts_L0_seu_list[[c]] <- counts_L0_seu

}

seu_list_s1 <- SplitObject(seu_list2[[1]], split.by = "stage")
seu_list_s2 <- SplitObject(seu_list2[[2]], split.by = "stage")
seu_list_s3 <- SplitObject(seu_list2[[3]], split.by = "stage")
seu_list_s4 <- SplitObject(seu_list2[[4]], split.by = "stage")
seu_list_s5 <- SplitObject(seu_list2[[5]], split.by = "stage")
seu_list_s6 <- SplitObject(seu_list2[[6]], split.by = "stage")
seu_list_s7 <- SplitObject(seu_list2[[7]], split.by = "stage")
seu_list_s8 <- SplitObject(seu_list2[[8]], split.by = "stage")
seu_list_s9 <- SplitObject(seu_list2[[9]], split.by = "stage")
seu_list_s10 <- SplitObject(seu_list2[[10]], split.by = "stage")

rm(seu_list, seu_list2, seu, scRNA_harmony)
gc()

sss <- merge(seu_list_s1[[1]],
             y = c(seu_list_s2[[1]], 
                   seu_list_s3[[1]], 
                   seu_list_s4[[1]], 
                   seu_list_s5[[1]], 
                   seu_list_s6[[1]], 
                   seu_list_s7[[1]], 
                   seu_list_s8[[1]], 
                   seu_list_s9[[1]], 
                   seu_list_s10[[1]]))

library("CellChat")
cellchat <- createCellChat(object = sss, meta = sss@meta.data, group.by = "new_celltype")
levels(cellchat@idents)


###
CellChatDB <- CellChatDB.human # use  CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") # Cell-Cell Contact, ECM-Receptor,  Secreted Signaling
cellchat@DB <- CellChatDB.use
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

#future::plan("multiprocess", workers = 8)
future::plan("multisession", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(2, 3), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pathways.show <- c("CXCL")
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells
vertex.receiver <- seq(1, 4) # a numeric vector.
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow = c(1, 1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow = c(1, 1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.
