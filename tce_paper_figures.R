#########################################################
## "T cell engagers control solid tumors through ########  
## IL2-driven effector differentiation of CD8 T-cells" ##
## Manuscript Companion Code - Chris McGinnis, PhD ######
## Ansuman Satpathy Lab, Stanford University, 2025 ######
#########################################################

library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ComplexHeatmap)

##############################################################################################################################
## Figure 4: TCE-IL2REH combination therapy promotes anti-tumor TIL responses and differentially affects CD4+ and CD8+ TILs ##
##############################################################################################################################
load('seu_cd4_clean.Robj')
load('seu_cd8_clean.Robj')
load('cd4_sub_freq.Robj')
load('cd8_sub_freq.Robj')
load('cd8_tex_zscores.Robj')
load('cd8_act_zscores.Robj')

## Figure 4B: CD4 and CD8 subtype UMAPs
DimPlot(seu_cd4_clean, cols=c('black','tan4','cadetblue4','dodgerblue','cadetblue2')) + NoLegend() + NoAxes()
DimPlot(seu_cd8_clean, cols=c('goldenrod','lightcoral','red','tan4','darkred','black')) + NoLegend() + NoAxes()


## Figure 4C: CD4 subtype proportion barchart and sample UMAPs
# cd4_sub_freq <- table(seu_cd4_clean@meta.data$sample, seu_cd4_clean@meta.data$subtype)
# cd4_sub_freq <- round(cd4_sub_freq/rowSums(cd4_sub_freq), digits=4)
# cd4_sub_freq <- melt(cd4_sub_freq)
# colnames(cd4_sub_freq) <- c('sample','subtype','freq')
# cd4_sub_freq$sample <- factor(cd4_sub_freq$sample, levels=c('hbs','il2','tce','tce_il2'))
# cd4_sub_freq$subtype <- factor(cd4_sub_freq$subtype, paste0('cd4_',c('naive','mem','act','treg','prolif')))
# 3x4.5 pdf
ggplot(cd4_sub_freq, aes(x=sample, y=freq, fill=subtype)) + geom_col(color='black') + theme_classic() + scale_fill_manual(values=c('cadetblue2','dodgerblue','cadetblue4','tan4','black')) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
DimPlot(seu_cd4_clean, cells.highlight = colnames(seu_cd4_clean)[which(seu_cd4_clean@meta.data$sample == 'hbs')], cols.highlight = 'black', sizes.highlight = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_cd4_clean, cells.highlight = colnames(seu_cd4_clean)[which(seu_cd4_clean@meta.data$sample == 'il2')], cols.highlight = 'dodgerblue3', sizes.highlight = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_cd4_clean, cells.highlight = colnames(seu_cd4_clean)[which(seu_cd4_clean@meta.data$sample == 'tce')], cols.highlight = 'orange2', sizes.highlight = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_cd4_clean, cells.highlight = colnames(seu_cd4_clean)[which(seu_cd4_clean@meta.data$sample == 'tce_il2')], cols.highlight = 'orchid4', sizes.highlight = 0.1) + NoLegend() + NoAxes()


## Figure 4D: CD8 subtype proportion barchart
# cd8_sub_freq <- table(seu_cd8_clean@meta.data$sample, seu_cd8_clean@meta.data$subtype)
# cd8_sub_freq <- round(cd8_sub_freq/rowSums(cd8_sub_freq), digits=4)
# cd8_sub_freq <- melt(cd8_sub_freq)
# colnames(cd8_sub_freq) <- c('sample','subtype','freq')
# cd8_sub_freq$sample <- factor(cd8_sub_freq$sample, levels=c('hbs','il2','tce','tce_il2'))
# cd8_sub_freq$subtype <- factor(cd8_sub_freq$subtype, paste0('cd8_',c('naive','act','tpex','tex','tex_prolif','ifn')))
# 3x4.5 pdf
ggplot(cd8_sub_freq, aes(x=sample, y=freq, fill=subtype)) + geom_col(color='black') + theme_classic() + scale_fill_manual(values=c('lightcoral','red','tan4','darkred','black','goldenrod')) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))


## Figure 4E: Tex z-score heatmap
# cd8_tex_hbs <- colnames(seu_cd8_clean)[which(seu_cd8_clean@meta.data$subtype == 'cd8_tex' & seu_cd8_clean@meta.data$sample == 'hbs')]
# cd8_tex_il2 <- colnames(seu_cd8_clean)[which(seu_cd8_clean@meta.data$subtype == 'cd8_tex' & seu_cd8_clean@meta.data$sample == 'il2')]
# cd8_tex_tce <- colnames(seu_cd8_clean)[which(seu_cd8_clean@meta.data$subtype == 'cd8_tex' & seu_cd8_clean@meta.data$sample == 'tce')]
# cd8_tex_tce_il2 <- colnames(seu_cd8_clean)[which(seu_cd8_clean@meta.data$subtype == 'cd8_tex' & seu_cd8_clean@meta.data$sample == 'tce_il2')]
# cd8_tex_degs_il2 <- FindMarkers(seu_cd8_clean, ident.1 = cd8_tex_hbs, ident.2 = cd8_tex_il2, logfc.threshold = log(2), min.pct = 0.8)
# cd8_tex_degs_tce <- FindMarkers(seu_cd8_clean, ident.1 = cd8_tex_hbs, ident.2 = cd8_tex_tce, logfc.threshold = log(2), min.pct = 0.8)
# cd8_tex_degs_tce_il2 <- FindMarkers(seu_cd8_clean, ident.1 = cd8_tex_hbs, ident.2 = cd8_tex_tce_il2, logfc.threshold = log(2), min.pct = 0.8)
# cd8_tex_degs <- unique(c(rownames(cd8_tex_degs_il2), rownames(cd8_tex_degs_tce), rownames(cd8_tex_degs_tce_il2)))
# cd8_tex_zscores <- as.data.frame(matrix(0L, nrow=4, ncol=length(cd8_tex_degs)))
# rownames(cd8_tex_zscores) <- c('cd8_tex_hbs','cd8_tex_il2','cd8_tex_tce','cd8_tex_tce_il2')
# colnames(cd8_tex_zscores) <- cd8_tex_degs
# cd8_tex_cells <- c(cd8_tex_hbs,cd8_tex_il2,cd8_tex_tce,cd8_tex_tce_il2)
# cd8_tex_genes_mean <- rowMeans(seu_cd8_clean@assays$RNA$data[cd8_tex_degs, cd8_tex_cells])
# cd8_tex_genes_sd <- apply(seu_cd8_clean@assays$RNA$data[cd8_tex_degs, cd8_tex_cells], 1, FUN = function(x) { sd(x) } )
# for (i in rownames(cd8_tex_zscores)) {
#   cells <- get(i, envir = .GlobalEnv)
#   temp <-  rowMeans(seu_cd8_clean@assays$RNA$data[cd8_tex_degs, cells])
#   cd8_tex_zscores[i, ] <- (temp-cd8_tex_genes_mean)/cd8_tex_genes_sd
# }
Heatmap(cd8_tex_zscores, col = rev(brewer.pal(name='RdBu',n=9)), cluster_rows = F, column_names_gp = gpar(fontsize=8), row_names_gp = gpar(fontsize=8))


## Fig. 4F: Tex DEG violin plots
# seu_cd8_clean <- SetIdent(seu_cd8_clean, value=seu_cd8_clean@meta.data$tex_sample)
# 3x3 pdf
VlnPlot(seu_cd8_clean, 'Tox', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_clean, 'Tigit', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_clean, 'Gzmb', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))


## Figure 4G: Activated CD8+ TIL z-score heatmap
# cd8_act_hbs <- colnames(seu_cd8_clean)[which(seu_cd8_clean@meta.data$subtype == 'cd8_act' & seu_cd8_clean@meta.data$sample == 'hbs')]
# cd8_act_il2 <- colnames(seu_cd8_clean)[which(seu_cd8_clean@meta.data$subtype == 'cd8_act' & seu_cd8_clean@meta.data$sample == 'il2')]
# cd8_act_tce <- colnames(seu_cd8_clean)[which(seu_cd8_clean@meta.data$subtype == 'cd8_act' & seu_cd8_clean@meta.data$sample == 'tce')]
# cd8_act_tce_il2 <- colnames(seu_cd8_clean)[which(seu_cd8_clean@meta.data$subtype == 'cd8_act' & seu_cd8_clean@meta.data$sample == 'tce_il2')]
# cd8_act_degs_il2 <- FindMarkers(seu_cd8_clean, ident.1 = cd8_act_hbs, ident.2 = cd8_act_il2, logfc.threshold = log(2.5), min.pct = 0.8)
# cd8_act_degs_tce <- FindMarkers(seu_cd8_clean, ident.1 = cd8_act_hbs, ident.2 = cd8_act_tce, logfc.threshold = log(2.5), min.pct = 0.8)
# cd8_act_degs_tce_il2 <- FindMarkers(seu_cd8_clean, ident.1 = cd8_act_hbs, ident.2 = cd8_act_tce_il2, logfc.threshold = log(2.5), min.pct = 0.8)
# cd8_act_degs <- unique(c(rownames(cd8_act_degs_il2), rownames(cd8_act_degs_tce), rownames(cd8_act_degs_tce_il2)))
# cd8_act_zscores <- as.data.frame(matrix(0L, nrow=4, ncol=length(cd8_act_degs)))
# rownames(cd8_act_zscores) <- c('cd8_act_hbs','cd8_act_il2','cd8_act_tce','cd8_act_tce_il2')
# colnames(cd8_act_zscores) <- cd8_act_degs
# cd8_act_cells <- c(cd8_act_hbs,cd8_act_il2,cd8_act_tce,cd8_act_tce_il2)
# cd8_act_genes_mean <- rowMeans(seu_cd8_clean@assays$RNA$data[cd8_act_degs, cd8_act_cells])
# cd8_act_genes_sd <- apply(seu_cd8_clean@assays$RNA$data[cd8_act_degs, cd8_act_cells], 1, FUN = function(x) { sd(x) } )
# for (i in rownames(cd8_act_zscores)) {
#   cells <- get(i, envir = .GlobalEnv)
#   temp <-  rowMeans(seu_cd8_clean@assays$RNA$data[cd8_act_degs, cells])
#   cd8_act_zscores[i, ] <- (temp-cd8_act_genes_mean)/cd8_act_genes_sd
# }
Heatmap(cd8_act_zscores, col = rev(brewer.pal(name='RdBu',n=9)), cluster_rows = F, column_names_gp = gpar(fontsize=8), row_names_gp = gpar(fontsize=8))


## Fig. 4H: Activated CD8+ TIL DEG violin plots
# seu_cd8_clean <- SetIdent(seu_cd8_clean, value=seu_cd8_clean@meta.data$act_sample)
# 3x3 pdf
VlnPlot(seu_cd8_clean, 'Gzmb', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_clean, 'Tnfrsf18', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_clean, 'Tnfrsf9', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))





###########################################################################################################################
## Figure 5: TCE-mediated tumor growth control is associated with clonal replacement followed by activation of CD8+ TILs ##
###########################################################################################################################
load('sample_tcr_freq.Robj')
load('cd8_clone_mat.Robj')
load('clone_counts_df.Robj')
load('seu_cd8_tcr.Robj')
load('seu_cd8_tcr_sub.Robj')
load('subtype_cd8_tcr_freq.Robj')
load('tex_prolif_zscores')
load('gly_zscores.Robj')
load('costim_eff_zscores.Robj')
load('tol_zscores.Robj')

## Figure 5A: Sample-level TCR diversity and clone size bar charts
# sample_tcr_freq <- table(tcr_meta$sample, tcr_meta$cloneSize)
# sample_tcr_freq <- round(sample_tcr_freq/rowSums(sample_tcr_freq),digits=4)
# sample_tcr_freq <- melt(sample_tcr_freq)
# colnames(sample_tcr_freq) <- c('sample','cloneSize','freq')
# sample_tcr_freq$sample <- factor(sample_tcr_freq$sample, levels=c('hbs','il2','tce','tce_il2'))
# sample_tcr_freq$cloneSize <- factor(sample_tcr_freq$cloneSize, levels=rev(c('rare','small','medium','large')))
# 4.5x3 pdf
ggplot(sample_tcr_freq, aes(x=sample, y=freq, fill=cloneSize)) + geom_col(color='black') + theme_classic() + scale_fill_manual(values=c('darkred','maroon','red','lightcoral')) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))

# shannon_mat <- as.data.frame(matrix(0L, nrow=100, ncol=4))
# colnames(shannon_mat) <- c('hbs','il2','tce','tce_il2')
# for (i in 1:100) {
#   print(i)
#   ind <- c(sample(which(tcr_meta$sample == 'hbs'), 2000),
#            sample(which(tcr_meta$sample == 'il2'), 2000),
#            sample(which(tcr_meta$sample == 'tce'), 2000),
#            sample(which(tcr_meta$sample == 'tce_il2'), 2000))
#   g <- plot_diversity(input = tcr_meta[ind,], data_col = "clonotype_id", cluster_col = "sample", method = abdiv::shannon)
#   shannon_mat[i,] <- g$data$diversity
# }
# shannon_summary <- as.data.frame(matrix(0L, nrow=4, ncol=3))
# colnames(shannon_summary) <- c('sample','mean','sd')
# shannon_summary$sample <- c('hbs','il2','tce','tce_il2')
# shannon_summary$mean <- c(mean(shannon_mat[,1]), mean(shannon_mat[,2]), mean(shannon_mat[,3]), mean(shannon_mat[,4]))
# shannon_summary$sd <- c(sd(shannon_mat[,1]), sd(shannon_mat[,2]), sd(shannon_mat[,3]), sd(shannon_mat[,4]))
# ggplot(shannon_summary, aes(x=sample, y=mean, fill=sample)) + geom_col(color="black") + theme_classic() + 
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) + 
#   theme(legend.position = 'none', axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm")) + 
#   scale_fill_manual(values=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8))


## Figure 5B: Subtype clone frequency heatmap and count histogram
# clone_counts <- table(cd8_tcr_meta$clonotype_id)
# clone_counts <- clone_counts[which(clone_counts >= 2)]
# cd8_clone_mat <- table(cd8_tcr_meta$clonotype_id, cd8_tcr_meta$subtype)
# cd8_clone_mat <- cd8_clone_mat[names(clone_counts), ]
# cd8_clone_mat <- round(cd8_clone_mat/rowSums(cd8_clone_mat), digits=4)
Heatmap(cd8_clone_mat, col = brewer.pal(name='Reds',n=9), show_row_names = F)
# g <- Heatmap(cd8_clone_mat, col = brewer.pal(name='Reds',n=9), show_row_names = F)
# g <- prepare(g)
# clone_counts_sub <- clone_counts[g@row_order]
# clone_counts_df <- as.data.frame(matrix(0L, nrow=length(clone_counts_sub), ncol=2))
# colnames(clone_counts_df) <- c('clone','ncell')
# clone_counts_df$clone <- names(clone_counts_sub)
# clone_counts_df$ncell <- clone_counts_sub
# clone_counts_df$clone <- factor(clone_counts_df$clone, levels=rev(clone_counts_df$clone))
ggplot(clone_counts_df, aes(x=clone, y=log2(ncell))) + geom_col(fill='black',color='black') + theme_classic() + theme(axis.text.x = element_blank())

# HBS: nCell = 3143; nClone = 1014
# TCE: nCell = 6682; nClone = 3889
# IL2: nCell = 10605; nClone = 4180
# Combo: nCell = 10001; nClone = 6543

## Figure 5C: CD8+ TIL UMAP w/ subtypes and clone sizes
DimPlot(seu_cd8_tcr, cols=c('goldenrod','black','darkred','tan4','red','lightcoral'), group.by = 'subtype') + NoLegend() + NoAxes() + theme(plot.title = element_blank())
DimPlot(seu_cd8_tcr, group.by = 'cloneSize', cols=c('darkred','maroon','lightcoral','red')) + NoLegend() + NoAxes() + theme(plot.title = element_blank())


## Figure 5D: Subsetted CD8+ TIL UMAP w/ subtypes, clone sizes, and treatment group projections
DimPlot(seu_cd8_tcr_sub, group.by = 'cloneSize', cols=c('darkred','maroon','lightcoral','red'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
DimPlot(seu_cd8_tcr_sub, group.by = 'subtype', cols=c('black','darkred','red'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
DimPlot(seu_cd8_tcr_sub, cells.highlight = colnames(seu_cd8_tcr_sub)[which(seu_cd8_tcr_sub@meta.data$sample == 'hbs')], cols.highlight = 'black', sizes.highlight = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_cd8_tcr_sub, cells.highlight = colnames(seu_cd8_tcr_sub)[which(seu_cd8_tcr_sub@meta.data$sample == 'il2')], cols.highlight = 'dodgerblue3', sizes.highlight = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_cd8_tcr_sub, cells.highlight = colnames(seu_cd8_tcr_sub)[which(seu_cd8_tcr_sub@meta.data$sample == 'tce')], cols.highlight = 'orange2', sizes.highlight = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_cd8_tcr_sub, cells.highlight = colnames(seu_cd8_tcr_sub)[which(seu_cd8_tcr_sub@meta.data$sample == 'tce_il2')], cols.highlight = 'orchid4', sizes.highlight = 0.1) + NoLegend() + NoAxes()


## Figure 5E: CD8 subtype clone size barchart
# subtype_cd8_tcr_freq <- table(seu_cd8_tcr_sub@meta.data$subtype_sample, seu_cd8_tcr_sub@meta.data$cloneSize)
# subtype_cd8_tcr_freq <- round(subtype_cd8_tcr_freq/rowSums(subtype_cd8_tcr_freq),digits=4)
# subtype_cd8_tcr_freq <- melt(subtype_cd8_tcr_freq)
# colnames(subtype_cd8_tcr_freq) <- c('subtype','cloneSize','freq')
# subtype_cd8_tcr_freq$cloneSize <- factor(subtype_cd8_tcr_freq$cloneSize, levels=rev(c('rare','small','medium','large')))
# subtype_cd8_tcr_freq$cloneSize <- factor(subtype_cd8_tcr_freq$cloneSize, levels=rev(c('rare','small','medium','large')))
ggplot(subtype_cd8_tcr_freq[grep('cd8_act', subtype_cd8_tcr_freq$subtype), ], aes(x=subtype, y=freq, fill=cloneSize)) + geom_col(color='black') + theme_classic() + scale_fill_manual(values=c('darkred','maroon','red','lightcoral')) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
ggplot(subtype_cd8_tcr_freq[grep('prolif', subtype_cd8_tcr_freq$subtype), ], aes(x=subtype, y=freq, fill=cloneSize)) + geom_col(color='black') + theme_classic() + scale_fill_manual(values=c('darkred','maroon','red','lightcoral')) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
ggplot(subtype_cd8_tcr_freq[grep('tex_[h|i|t|]', subtype_cd8_tcr_freq$subtype), ], aes(x=subtype, y=freq, fill=cloneSize)) + geom_col(color='black') + theme_classic() + scale_fill_manual(values=c('darkred','maroon','red','lightcoral')) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))


## Figure 5F: Tex-prolif expanded vs diverse marker gene z-score heatmap
# tex_prolif_degs <- FindMarkers(seu_cd8_tcr_sub, ident.1 = colnames(seu_cd8_tcr_sub)[which(seu_cd8_tcr_sub@meta.data$tex_prolif_sub == 'diverse')], ident.2 = colnames(seu_cd8_tcr_sub)[which(seu_cd8_tcr_sub@meta.data$tex_prolif_sub == 'biased')], logfc.threshold=log(2), min.pct=0.7)
# tex_prolif_degs <- rownames(tex_prolif_degs)
# tex_prolif_zscores <- as.data.frame(matrix(0L, nrow=2, ncol=length(tex_prolif_degs)))
# rownames(tex_prolif_zscores) <- c('biased','diverse')
# colnames(tex_prolif_zscores) <- tex_prolif_degs
# cells <- colnames(seu_cd8_tcr_sub)[which(seu_cd8_tcr_sub@meta.data$tex_prolif_sub != 0)]
# tex_prolif_genes_mean <- rowMeans(seu_cd8_tcr_sub@assays$RNA$data[tex_prolif_degs, cells])
# tex_prolif_genes_sd <- apply(seu_cd8_tcr_sub@assays$RNA$data[tex_prolif_degs, cells], 1, FUN = function(x) { sd(x) } )
# for (i in rownames(tex_prolif_zscores)) {
#   cells <- colnames(seu_cd8_tcr_sub)[which(seu_cd8_tcr_sub@meta.data$tex_prolif_sub == i)]
#   temp <-  rowMeans(seu_cd8_tcr_sub@assays$RNA$data[tex_prolif_degs, cells])
#   tex_prolif_zscores[i, ] <- (temp-tex_prolif_genes_mean)/tex_prolif_genes_sd
# }
Heatmap(tex_prolif_zscores, col = rev(brewer.pal(name='RdBu',n=9)), cluster_rows = F, column_names_gp = gpar(fontsize=8), row_names_gp = gpar(fontsize=8))


## Figure 5G: Tex-prolif expanded vs diverse marker gene violin plots
# seu_cd8_tcr_sub <- SetIdent(seu_cd8_tcr_sub, value=seu_cd8_tcr_sub@meta.data$tex_prolif_sample)
VlnPlot(seu_cd8_tcr_sub, 'Tox', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_tcr_sub, 'Id2', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_tcr_sub, 'Pdcd1', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_tcr_sub, 'Pim1', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_tcr_sub, 'Cd69', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_tcr_sub, 'Nfkb1', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))


## Figure 5I: Clonally-diverse CD8 TIL glycolysis, effector/co-stim, and tolerance marker z-score heatmaps
# WIKIPATHWAYS: WP_GLYCOLYSIS_AND_GLUCONEOGENESIS
# gly_genes <- c('Gck','Ldhal6b','Aldoa','Aldoc','Dld','Eno1','Eno2','Eno3','Fbp2','Fbp1','G6pc1','Gapdh','Gapdhs','Got1','Got2','Gpi1','Hk1','Hk2','Ldha','Ldhb','Ldhc','Mdh2','Mdh1','Pck1','Pcx','Pdha1','Pdha2','Pfkl','Pfkm','Pgam1','Pgk1','Pgk2','Pkm','Pklr','Slc2a1','Slc2a2','Slc2a3','Slc2a4','Hk3','Tpi1','Aldob','Dlat','Pdhx','Mpc1','Pgam2','Pfkp','Slc2a5','Pdhb','Mpc2')
# gly_genes <- gly_genes[which(gly_genes %in% VariableFeatures(seu_cd8_tcr))]
# gly_zscores <- as.data.frame(matrix(0L, nrow=2, ncol=length(gly_genes)))
# rownames(gly_zscores) <- c('tce_rare','tce_il2_rare')
# colnames(gly_zscores) <- gly_genes
# cells <- colnames(seu_cd8_tcr)[which(seu_cd8_tcr@meta.data$sample_cloneSize %in% c('tce_rare','tce_il2_rare'))]
# gly_genes_mean <- rowMeans(seu_cd8_tcr@assays$RNA$data[gly_genes, cells])
# gly_genes_sd <- apply(seu_cd8_tcr@assays$RNA$data[gly_genes, cells], 1, FUN = function(x) { sd(x) } )
# for (i in rownames(gly_zscores)) {
#   cells <- colnames(seu_cd8_tcr)[which(seu_cd8_tcr@meta.data$sample_cloneSize == i)]
#   temp <-  rowMeans(seu_cd8_tcr@assays$RNA$data[gly_genes, cells])
#   gly_zscores[i, ] <- (temp-gly_genes_mean)/gly_genes_sd
# }
Heatmap(gly_zscores, col = rev(brewer.pal(name='RdBu',n=9)), cluster_rows = F, column_names_gp = gpar(fontsize=8), row_names_gp = gpar(fontsize=8))

# Effector + co-stim signature
# costim_eff_genes <- c('Cd27','Cd28','Tnfrsf4','Tnfrsf9','Tnfrsf18','Gzma','Gzmb','Gzmk','Casp3', 'Zbtb32','Ptma','Gadd45b','Serpinb9','Prf1','Ifng')
# costim_eff_genes <- costim_eff_genes[which(costim_eff_genes %in% VariableFeatures(seu_cd8_tcr))]
# costim_eff_zscores <- as.data.frame(matrix(0L, nrow=2, ncol=length(costim_eff_genes)))
# rownames(costim_eff_zscores) <- c('tce_rare','tce_il2_rare')
# colnames(costim_eff_zscores) <- costim_eff_genes
# cells <- colnames(seu_cd8_tcr)[which(seu_cd8_tcr@meta.data$sample_cloneSize %in% c('tce_rare','tce_il2_rare'))]
# costim_eff_genes_mean <- rowMeans(seu_cd8_tcr@assays$RNA$data[costim_eff_genes, cells])
# costim_eff_genes_sd <- apply(seu_cd8_tcr@assays$RNA$data[costim_eff_genes, cells], 1, FUN = function(x) { sd(x) } )
# for (i in rownames(costim_eff_zscores)) {
#   cells <- colnames(seu_cd8_tcr)[which(seu_cd8_tcr@meta.data$sample_cloneSize == i)]
#   temp <-  rowMeans(seu_cd8_tcr@assays$RNA$data[costim_eff_genes, cells])
#   costim_eff_zscores[i, ] <- (temp-costim_eff_genes_mean)/costim_eff_genes_sd
# }
Heatmap(costim_eff_zscores, col = rev(brewer.pal(name='RdBu',n=9)), cluster_rows = F, column_names_gp = gpar(fontsize=8), row_names_gp = gpar(fontsize=8))

# CD8+ T-cell tolerance signature 
# tol_genes <- c('Ass1','Btla','Capg','Capn3','Cd200','Coro1b','Dapl1','Dgkz','Endod1','Gm19220','Gm23935','Ift80','Ikzf2','Izumo1r','Mxd4','Nrgn','Nrn1','Rgcc','Rpl31-ps8','Rxra','Scd2','Smc4','Snora57','Tns1','Tspan32','Wnt10a')
# tol_genes <- tol_genes[which(tol_genes %in% VariableFeatures(seu_cd8_tcr))]
# tol_zscores <- as.data.frame(matrix(0L, nrow=2, ncol=length(tol_genes)))
# rownames(tol_zscores) <- c('tce_rare','tce_il2_rare')
# colnames(tol_zscores) <- tol_genes
# cells <- colnames(seu_cd8_tcr)[which(seu_cd8_tcr@meta.data$sample_cloneSize %in% c('tce_rare','tce_il2_rare'))]
# tol_genes_mean <- rowMeans(seu_cd8_tcr@assays$RNA$data[tol_genes, cells])
# tol_genes_sd <- apply(seu_cd8_tcr@assays$RNA$data[tol_genes, cells], 1, FUN = function(x) { sd(x) } )
# for (i in rownames(tol_zscores)) {
#   cells <- colnames(seu_cd8_tcr)[which(seu_cd8_tcr@meta.data$sample_cloneSize == i)]
#   temp <-  rowMeans(seu_cd8_tcr@assays$RNA$data[tol_genes, cells])
#   tol_zscores[i, ] <- (temp-tol_genes_mean)/tol_genes_sd
# }
Heatmap(tol_zscores, col = rev(brewer.pal(name='RdBu',n=9)), cluster_rows = F, column_names_gp = gpar(fontsize=8), row_names_gp = gpar(fontsize=8))





#############################
## Supplemental Figure XXX ##
#############################
load("seu_cd4_clean.Robj")
load("seu_cd8_clean.Robj")
load("cd4_cd8_freq.Robj")
load("cd8_sub_freq.Robj")
load("subtype_cd4_tcr_freq.Robj")
load('cd4_genes.Robj')
load('cd8_genes.Robj')

## Figure SXA: CD4 and CD8 subtype marker dot plots
# cd4_markers <- FindAllMarkers(seu_cd4_clean, only.pos = T, logfc.threshold = log(2), min.pct = 0.5)
# cd4_genes <- NULL
# for (i in paste0('cd4_',c('naive','mem','act','treg','prolif'))) { cd4_genes <- c(cd4_genes, cd4_markers$gene[which(cd4_markers$cluster == i)[1:3]]) }
# cd4_genes <- unique(cd4_genes)
# 5x3 pdf
DotPlot(seu_cd4_clean, features=cd4_genes, cols='RdBu', dot.scale = 3) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))

# cd8_markers <- FindAllMarkers(seu_cd8_clean, only.pos = T, logfc.threshold = log(2), min.pct = 0.5)
# cd8_genes <- NULL
# for (i in paste0('cd8_',c('naive','act','tpex','tex','tex_prolif','ifn'))) { cd8_genes <- c(cd8_genes, cd8_markers$gene[which(cd8_markers$cluster == i)[1:3]]) }
# cd8_genes <- unique(cd8_genes)
# 5.25x3 pdf
DotPlot(seu_cd8_clean, features=cd8_genes, cols='RdBu', dot.scale = 3) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))

## Fig. SXB: CD4/CD8 proportion bar chart
# cd4_cd8_freq <- table(seu_clean@meta.data$sample, seu_clean@meta.data$cd4cd8)
# cd4_cd8_freq <- round(cd4_cd8_freq/rowSums(cd4_cd8_freq), digits=4)
# cd4_cd8_freq <- melt(cd4_cd8_freq)
# colnames(cd4_cd8_freq) <- c('sample','cd4cd8','freq')
# cd4_cd8_freq$sample <- factor(cd4_cd8_freq$sample, levels=c('hbs','il2','tce','tce_il2'))
ggplot(cd4_cd8_freq, aes(x=sample, y=freq, fill=cd4cd8)) + geom_col(color='black') + theme_classic() + scale_fill_manual(values=c('black','grey40')) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))


## Fig. SXC: Activated CD8+ TIL and Tex proportion staggered barcharts
ggplot(cd8_sub_freq[which(cd8_sub_freq$subtype %in% c('cd8_tex','cd8_act') & cd8_sub_freq$sample != 'tce'), ], aes(x=subtype, y=freq, fill=sample)) + geom_col(color='black', position=position_dodge()) + theme_classic() + scale_fill_manual(values=alpha(c('black','dodgerblue3','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))


## Fig. SXD: Glycolysis gene score violin plot
# WIKIPATHWAYS: WP_GLYCOLYSIS_AND_GLUCONEOGENESIS
# gly_genes <- c('Gck','Ldhal6b','Aldoa','Aldoc','Dld','Eno1','Eno2','Eno3','Fbp2','Fbp1','G6pc1','Gapdh','Gapdhs','Got1','Got2','Gpi1','Hk1','Hk2','Ldha','Ldhb','Ldhc','Mdh2','Mdh1','Pck1','Pcx','Pdha1','Pdha2','Pfkl','Pfkm','Pgam1','Pgk1','Pgk2','Pkm','Pklr','Slc2a1','Slc2a2','Slc2a3','Slc2a4','Hk3','Tpi1','Aldob','Dlat','Pdhx','Mpc1','Pgam2','Pfkp','Slc2a5','Pdhb','Mpc2')
# gly_genes <- gly_genes[which(gly_genes %in% VariableFeatures(seu_cd8_clean))]
# seu_cd8_clean <- AddModuleScore(seu_cd8_clean, features = list(gly_genes), ctrl=5, name='glycolysis')
# seu_cd8_clean <- SetIdent(seu_cd8_clean, value=seu_cd8_clean@meta.data$act_sample)
VlnPlot(seu_cd8_clean, 'glycolysis1', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))


## Fig. SXE: CD4 subtype clonal diversity bar charts
# subtype_cd4_tcr_freq <- table(cd4_tcr_meta$subtype, cd4_tcr_meta$cloneSize)
# subtype_cd4_tcr_freq <- round(subtype_cd4_tcr_freq/rowSums(subtype_cd4_tcr_freq),digits=4)
# subtype_cd4_tcr_freq <- melt(subtype_cd4_tcr_freq)
# colnames(subtype_cd4_tcr_freq) <- c('subtype','cloneSize','freq')
# subtype_cd4_tcr_freq$subtype <- factor(subtype_cd4_tcr_freq$subtype, levels=paste0('cd4_',c('naive','mem','act','treg','prolif')))
# subtype_cd4_tcr_freq$cloneSize <- factor(subtype_cd4_tcr_freq$cloneSize, levels=rev(c('rare','small','medium','large')))
ggplot(subtype_cd4_tcr_freq, aes(x=subtype, y=freq, fill=cloneSize)) + geom_col(color='black') + theme_classic() + scale_fill_manual(values=c('darkred','maroon','red','lightcoral')) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))

