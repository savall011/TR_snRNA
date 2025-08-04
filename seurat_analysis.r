# 1. 加载包
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(monocle3)

# 2. 定义文件与时间点对应关系
h5_files <- c(
  "filtered_feature_bc_matrix_1.h5",  # 0h
  "filtered_feature_bc_matrix_2.h5",  # 2h
  "filtered_feature_bc_matrix_3.h5",  # 8h
  "filtered_feature_bc_matrix_4.h5",  # 1d
  "filtered_feature_bc_matrix_5.h5",  # 3d
  "filtered_feature_bc_matrix_6.h5"   # 7d
)
timepoints <- c("0h","2h","8h","1d","3d","7d")

# 3. 读取并创建 Seurat 对象列表
seurat_list <- lapply(seq_along(h5_files), function(i) {
  message("Loading ", h5_files[i], " as timepoint ", timepoints[i])
  mat <- Read10X_h5(h5_files[i])
  so  <- CreateSeuratObject(counts = mat,
                             project = "TimeCourse",
                             min.cells = 3,
                             min.features = 200)
  so$timepoint <- timepoints[i]
  so
})
names(seurat_list) <- timepoints

# 4. 质控统计 & 过滤（根据需要调整阈值）
for (tp in timepoints) {
  so <- seurat_list[[tp]]
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
  VlnPlot(so, c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3) + ggtitle(tp)
  seurat_list[[tp]] <- subset(so,
    subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10
  )
}

# 5. SCTransform 归一化并鉴定变量基因
seurat_list <- lapply(seurat_list, SCTransform, verbose = FALSE)

# 6. 对齐与整合
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list,
                                  anchor.features = features,
                                  verbose = FALSE)
anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                  normalization.method = "SCT",
                                  anchor.features = features,
                                  verbose = FALSE)
seurat_integrated <- IntegrateData(anchorset = anchors,
                                   normalization.method = "SCT",
                                   verbose = FALSE)

# 7. 降维 & 聚类
DefaultAssay(seurat_integrated) <- "integrated"
seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE)
ElbowPlot(seurat_integrated, ndims = 30)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:20, verbose = FALSE)
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:20, verbose = FALSE)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.5, verbose = FALSE)

# 可视化
p1 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "timepoint", label = TRUE)
p2 <- DimPlot(seurat_integrated, reduction = "umap", label = TRUE)
print(p1 + p2)

# 8. 差异基因分析
markers <- FindAllMarkers(seurat_integrated,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat_integrated, features = top10$gene) + NoLegend()

# 9. 保存 Seurat 对象与标记基因
saveRDS(seurat_integrated, file = "seurat_timecourse_integrated.rds")
write.csv(markers, file = "seurat_timecourse_markers.csv", row.names = FALSE)

# ---------------------------------------
# 10. Monocle3 轨迹分析
# 将 Seurat 对象转换为 CellDataSet
cds <- as.cell_data_set(seurat_integrated)
# 保留 Seurat 上的 UMAP 坐标
cds <- cluster_cells(cds, reduction_method = "UMAP")  # 可选：重新聚类
cds <- learn_graph(cds, use_partition = TRUE)
# 指定起始细胞群（可根据时间或 cluster 设定）
# 例如：选择 0h 中的细胞为根节点
root_cells <- colnames(cds)[cds@colData$timepoint == "0h"]
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds, root_cells))

# 绘制伪时序
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE) +
  ggtitle("Pseudotime Trajectory")

# 根据 pseudotime 做差异基因
deg_pseudotime <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
top_deg <- deg_pseudotime %>% arrange(q_value) %>% head(100)
# 绘制热图
plot_genes_in_pseudotime(cds, genes = top_deg$gene_id) +
  theme(axis.text.y = element_text(size = 6))

# 保存 Monocle 结果
saveRDS(cds, file = "monocle3_cds.rds")
write.csv(as.data.frame(deg_pseudotime), file = "monocle3_deg_pseudotime.csv", row.names = FALSE)

# =============================================================================
# End of script
# =============================================================================
