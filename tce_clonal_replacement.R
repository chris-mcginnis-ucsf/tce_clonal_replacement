#########################################################
## "T cell engagers control solid tumors through ########  
## IL2-driven effector differentiation of CD8 T-cells" ##
## Manuscript Companion Code - Chris McGinnis, PhD ######
#########################################################

##############
## Figure 1 ##
##############
load('FIG1_trp2_abund.Robj')
load('FIG1_trp2_tact_lys.Robj')
load('FIG1_dll3_lys.Robj')

## Fig. 1A: Trp2 abundance bar chart
ggplot(trp2_abund[which(trp2_abund$stain == '+100nM_SPM14'),], aes(x=cell_line, y=digg, fill=cell_line)) + 
  geom_bar(stat = 'summary', fun = 'mean', color='black') + 
  geom_jitter(width = 0.2, size=0.8) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.5) +
  scale_y_continuous(trans = "pseudo_log",  breaks = c(100, 1000, 10000)) +
  coord_cartesian(ylim = c(100, NA)) + 
  theme_classic() + 
  scale_fill_manual(values=alpha(c('lightcoral','red','cadetblue2','maroon','dodgerblue','navy','darkred'),0.8)) + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

## Fig. 1B: Trp2 TCE design -- Exported image
## Fig. 1C: Trp2 gel -- Exported images
## Fig. 1D-E: In vitro co-culture analysis ~ T-cell activation, Trp2 Db/IgG TCE scatter-line plot 
ggplot(trp2_tact_lys[which(trp2_tact_lys$ab == 'db' & trp2_tact_lys$target %in% c('none',paste0('MC38_',c('low','mid','high','GFP')))), ], aes(x=conc, y=pCD137, color=target)) +
  stat_summary(geom="point") + 
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.1) +
  geom_smooth(method = "drm", method.args=list(fct=L.4()), se = FALSE, linewidth=1) +
  scale_x_log10(labels=label_number(drop0trailing = TRUE)) +
  theme_classic() +
  scale_color_manual(values=c('seagreen','navy','cadetblue2','dodgerblue','grey')) +
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

ggplot(trp2_tact_lys[which(trp2_tact_lys$ab == 'db' & trp2_tact_lys$target %in% c('none','MC38_GFP',paste0('B16',c('','_IFN','_IFN_TRP2')))), ], aes(x=conc, y=pCD137, color=target)) +
  stat_summary(geom="point") + 
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.1) +
  geom_smooth(method = "drm", method.args=list(fct=L.4()), se = FALSE, linewidth=1) +
  scale_x_log10(labels=label_number(drop0trailing = TRUE)) +
  theme_classic() +
  scale_color_manual(values=c('lightcoral','maroon','darkred','seagreen','grey')) +
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

ggplot(trp2_tact_lys[which(trp2_tact_lys$ab == 'igg' & trp2_tact_lys$target %in% c('none',paste0('MC38_',c('low','mid','high','GFP')))), ], aes(x=conc, y=pCD137, color=target)) +
  stat_summary(geom="point") + 
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.1) +
  geom_smooth(method = "drm", method.args=list(fct=L.4()), se = FALSE, linewidth=1) +
  scale_x_log10(labels=label_number(drop0trailing = TRUE)) +
  theme_classic() +
  scale_color_manual(values=c('seagreen','navy','cadetblue2','dodgerblue','grey')) +
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm")) + ylim(c(0,100))

ggplot(trp2_tact_lys[which(trp2_tact_lys$ab == 'igg' & trp2_tact_lys$target %in% c('none','MC38_GFP',paste0('B16',c('','_IFN','_IFN_TRP2')))), ], aes(x=conc, y=pCD137, color=target)) +
  stat_summary(geom="point") + 
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.1) +
  geom_smooth(method = "drm", method.args=list(fct=L.4()), se = FALSE, linewidth=1) +
  scale_x_log10(labels=label_number(drop0trailing = TRUE)) +
  theme_classic() +
  scale_color_manual(values=c('lightcoral','maroon','darkred','seagreen','grey')) +
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm")) + ylim(c(0,100))

## Fig. 1F: In vitro co-culture analysis ~ Cell lysis, Trp2 Db TCE scatter-line plot
ggplot(trp2_tact_lys[which(trp2_tact_lys$ab == 'db' & trp2_tact_lys$target %in% paste0('MC38_',c('low','mid','high','GFP'))), ], aes(x=conc, y=specific_lysis, color=target))+
  stat_summary(geom="point") + 
  stat_summary(geom="errorbar", width=0.1) +
  geom_smooth(method = "drm", method.args=list(fct=L.4()), se = FALSE, linewidth=1) +
  scale_x_log10(labels=label_number(drop0trailing = TRUE)) +
  theme_classic() +
  scale_color_manual(values=c('seagreen','navy','cadetblue2','dodgerblue','grey')) +
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm")) + ylim(c(0,100))

ggplot(trp2_tact_lys[which(trp2_tact_lys$ab == 'db' & trp2_tact_lys$target %in% c('MC38_GFP',paste0('B16',c('','_IFN','_IFN_TRP2')))), ], aes(x=conc, y=specific_lysis, color=target))+
  stat_summary(geom="point") + 
  stat_summary(geom="errorbar", width=0.1) +
  geom_smooth(method = "drm", method.args=list(fct=L.4()), se = FALSE, linewidth=1) +
  scale_x_log10(labels=label_number(drop0trailing = TRUE)) +
  theme_classic() +
  scale_color_manual(values=c('lightcoral','maroon','darkred','seagreen','grey')) +
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm")) + ylim(c(0,100))

## Fig. 1G: In vitro co-culture analysis ~ Cell lysis, DLL3 Db TCE scatter-line plot
ggplot(dll3_lys, aes(x=conc, y=specific_lysis, color=target))+
  stat_summary(geom="point") + 
  stat_summary(geom="errorbar", width=0.1) +
  geom_smooth(method = "drm", method.args=list(fct=L.4()), se = FALSE, linewidth=1) +
  scale_x_log10(labels=label_number(drop0trailing = TRUE)) +
  theme_classic() +
  scale_color_manual(values=c('darkorchid4','plum2','darkorchid2')) +
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm")) + ylim(c(0,100))

##############
## Figure 2 ##
##############
load('FIG2_til_b16_trp2.Robj')
load('FIG2_til_kp1_dll3.Robj')
load('FIG2_h2kb_kp1_dll3.Robj')

## Fig. 2A: Trp2 TCE treatment schematic -- Exported image
## Fig. 2B: Tumor growth dynamics (size) ~ High antigen (MC38) Trp2 TCE vs HBS -- PRISM analysis
## Fig. 2C: Tumor growth dynamics (proportion) ~ High antigen (MC38) Trp2 TCE vs HBS -- PRISM analysis
## Fig. 2D: Tumor growth dynamics (size) ~ Low antigen (B16) 5mg/kg Trp2 TCE vs 10mg/kg Trp2 TCE vs HBS -- PRISM analysis
## Fig. 2E: TIL numbers, low antigen (B16) 5mg/kg Trp2 TCE vs 10mg/kg Trp2 TCE vs HBS
ggplot(til_b16_trp2, aes(x=sample, y=tils_per_mg_tumor, fill=sample)) +
  geom_boxplot(outliers=F) + geom_jitter(width=0.1) +
  scale_y_continuous(trans = "pseudo_log",  breaks = c(30, 100, 300, 1000)) +
  theme_classic() + scale_fill_manual(values=c('grey30','goldenrod','tan4')) +
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

## Fig. 2F: DLL3 TCE treatment schematic -- Exported image
## Fig. 2G: Tumor growth dynamics (size) ~ High antigen (16T) DLL3 TCE vs HBS -- PRISM analysis
## Fig. 2H: Tumor growth dynamics (proportion) ~ High antigen (16T) DLL3 TCE vs HBS -- PRISM analysis
## Fig. 2I: Tumor growth dynamics (size) ~ Low antigen (KP1) 2mg/kg DLL3 TCE vs 5mg/kg DLL3 TCE vs 10mg/kg DLL3 TCE vs HBS -- PRISM analysis
## Fig. 2J: TIL numbers, low antigen (KP1) 2mg/kg DLL3 TCE vs 5mg/kg DLL3 TCE vs 10mg/kg DLL3 TCE vs HBS 
ggplot(til_kp1_dll3, aes(x=sample, y=tils_per_mg_tumor, fill=sample)) + 
  geom_boxplot(outliers=F) + geom_jitter(width=0.1) +
  scale_y_continuous(trans = "pseudo_log",  breaks = c(30, 100, 300, 1000, 3000)) +
  theme_classic() + scale_fill_manual(values=c('grey30','khaki3','goldenrod','tan4')) + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

## Fig. 2K: Tumor H2-Kb, low antigen (KP1) 2mg/kg DLL3 TCE vs 5mg/kg DLL3 TCE vs 10mg/kg DLL3 TCE vs HBS 
ggplot(h2kb_kp1_dll3, aes(x=sample, y=h2kb, fill=sample)) + 
  geom_boxplot(outliers=F) + geom_jitter(width=0.1) +
  scale_y_continuous(trans = "pseudo_log",  breaks = c(1e2, 1e3, 1e4, 1e5)) +
  theme_classic() + scale_fill_manual(values=c('grey30','khaki3','goldenrod','tan4')) + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

##############
## Figure 3 ##
##############
load('FIG3_ot1_tce_cyto.Robj')
load('FIG3_cd25_kp1_dll3.Robj')
load('FIG3_til_kp1_dll3_il2.Robj')
load('FIG3_gzmb_kp1_dll3_il2.Robj')

## Fig. 3A: OT-I Ova experimental schematic -- exported image
## Fig. 3B: OT-I Ova +/- Trp2 TCE IFNg secretion
ggplot(ot1_tce_cyto, aes(x=sample, y=ifng, fill=sample)) + 
  geom_bar(stat = 'summary', fun = 'mean', color='black') + 
  geom_jitter(width = 0.2, size=0.8) + ylim(c(0,45)) + 
  theme_classic() + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm")) + 
  scale_fill_manual(values=c('goldenrod','grey30','tan4'))

## Fig. 3C: OT-I Ova +/- Trp2 TCE IL2 secretion
ggplot(ot1_tce_cyto, aes(x=sample, y=il2, fill=sample)) + 
  geom_bar(stat = 'summary', fun = 'mean', color='black') + 
  geom_jitter(width = 0.2, size=0.8) + ylim(c(0,45)) + 
  theme_classic() + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm")) + 
  scale_fill_manual(values=c('goldenrod','grey30','tan4'))

## Fig. 3D: Tumor growth dynamics (size) ~ B16 Trp2 TCE vs wtIL2 vs Combo vs HBS -- PRISM analysis
## Fig. 3E-F: CD4+ and CD8+ T cell CD25 expression across organs +/- DLL3 TCE -- Flowjo analysis
## Fig. 3G: CD25+ CD8+ proportions in LN, spleen, and tumor ~ KP1 2mg/kg DLL3 TCE vs 5mg/kg DLL3 TCE vs 10mg/kg DLL3 TCE vs HBS 
ggplot(cd25_kp1_dll3[which(cd25_kp1_dll3$tissue == 'LN'), ], aes(x=sample, y=pCD25pos_cd8, fill=sample)) + 
  geom_boxplot(outliers=F) + geom_jitter(width=0.1) + ylim(c(0,100)) +
  theme_classic() + scale_fill_manual(values=c('grey30','khaki3','goldenrod','tan4')) + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

ggplot(cd25_kp1_dll3[which(cd25_kp1_dll3$tissue == 'Spleen'), ], aes(x=sample, y=pCD25pos_cd8, fill=sample)) + 
  geom_boxplot(outliers=F) + geom_jitter(width=0.1) + ylim(c(0,100)) +
  theme_classic() + scale_fill_manual(values=c('grey30','khaki3','goldenrod','tan4')) + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

ggplot(cd25_kp1_dll3[which(cd25_kp1_dll3$tissue == 'Tumor'), ], aes(x=sample, y=pCD25pos_cd8, fill=sample)) + 
  geom_boxplot(outliers=F) + geom_jitter(width=0.1) + ylim(c(0,100)) +
  theme_classic() + scale_fill_manual(values=c('grey30','khaki3','goldenrod','tan4')) + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

## Fig. 3H: CD25+ CD4+ proportions, KP1 tumors, 2mg/kg DLL3 TCE vs 5mg/kg DLL3 TCE vs 10mg/kg DLL3 TCE vs HBS 
ggplot(cd25_kp1_dll3[which(cd25_kp1_dll3$tissue == 'LN'), ], aes(x=sample, y=pCD25pos_cd4, fill=sample)) + 
  geom_boxplot(outliers=F) + geom_jitter(width=0.1) + ylim(c(0,100)) +
  theme_classic() + scale_fill_manual(values=c('grey30','khaki3','goldenrod','tan4')) + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

ggplot(cd25_kp1_dll3[which(cd25_kp1_dll3$tissue == 'Spleen'), ], aes(x=sample, y=pCD25pos_cd4, fill=sample)) + 
  geom_boxplot(outliers=F) + geom_jitter(width=0.1) + ylim(c(0,100)) +
  theme_classic() + scale_fill_manual(values=c('grey30','khaki3','goldenrod','tan4')) + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

ggplot(cd25_kp1_dll3[which(cd25_kp1_dll3$tissue == 'Tumor'), ], aes(x=sample, y=pCD25pos_cd4, fill=sample)) + 
  geom_boxplot(outliers=F) + geom_jitter(width=0.1) + ylim(c(0,100)) +
  theme_classic() + scale_fill_manual(values=c('grey30','khaki3','goldenrod','tan4')) + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

## Fig. 3I: B16/KP1 Trp2 TCE vs IL2reh vs Combo vs HBS experimental schematic -- exported image
## Fig. 3J: Tumor growth dynamics (size) ~ B16 Trp2 TCE vs IL2reh vs Combo vs HBS -- PRISM analysis
## Fig. 3L: TIL proportions, KP1 DLL3 TCE vs IL2reh vs Combo vs HBS
ggplot(til_kp1_dll3_il2, aes(x=sample, y=pTcell, fill=sample)) + 
  geom_boxplot(outliers=F) + geom_jitter(width=0.1) + ylim(c(0,100)) +
  theme_classic() + scale_fill_manual(values=c('grey30','dodgerblue3','goldenrod','orchid4')) + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

## Fig. 3M: CD8+ TIL proportions, KP1 DLL3 TCE vs IL2reh vs Combo vs HBS
ggplot(til_kp1_dll3_il2, aes(x=sample, y=pCD8, fill=sample)) + 
  geom_boxplot(outliers=F) + geom_jitter(width=0.1) + ylim(c(0,100)) +
  theme_classic() + scale_fill_manual(values=c('grey30','dodgerblue3','goldenrod','orchid4')) + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

## Fig. 3N: GZMB MFI, KP1 DLL3 TCE vs IL2reh vs Combo vs HBS
ggplot(gzmb_kp1_dll3_il2, aes(x=sample, y=gzmb_mfi, fill=sample)) + 
  geom_boxplot(outliers=F) + geom_jitter(width=0.1) + ylim(c(0,6e5)) +
  theme_classic() + scale_fill_manual(values=c('grey30','dodgerblue3','goldenrod','orchid4')) + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

##############
## Figure 4 ##
##############
load('FIG4_seu_cd8_clean.Robj')
load('FIG4_seu_cd4_clean.Robj')
load('FIG4_cd4_sub_freq.Robj')
load('FIG4_cd8_sub_freq.Robj')
load('FIG4_cd8_tex_zscores.Robj')
load('FIG4_cd8_act_zscores.Robj')

## Fig. 4B: CD4 and CD8 subtype UMAPs
DimPlot(seu_cd4_clean, cols=c('black','tan4','cadetblue4','dodgerblue','cadetblue2')) + NoLegend() + NoAxes()
DimPlot(seu_cd8_clean, cols=c('goldenrod','lightcoral','red','tan4','darkred','black')) + NoLegend() + NoAxes()

## Fig. 4C: CD4 subtype proportion barchart and sample UMAPs
ggplot(cd4_sub_freq, aes(x=sample, y=freq, fill=subtype)) + geom_col(color='black') + theme_classic() + scale_fill_manual(values=c('cadetblue2','dodgerblue','cadetblue4','tan4','black')) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
DimPlot(seu_cd4_clean, cells.highlight = colnames(seu_cd4_clean)[which(seu_cd4_clean@meta.data$sample == 'hbs')], cols.highlight = 'black', sizes.highlight = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_cd4_clean, cells.highlight = colnames(seu_cd4_clean)[which(seu_cd4_clean@meta.data$sample == 'il2')], cols.highlight = 'dodgerblue3', sizes.highlight = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_cd4_clean, cells.highlight = colnames(seu_cd4_clean)[which(seu_cd4_clean@meta.data$sample == 'tce')], cols.highlight = 'orange2', sizes.highlight = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_cd4_clean, cells.highlight = colnames(seu_cd4_clean)[which(seu_cd4_clean@meta.data$sample == 'tce_il2')], cols.highlight = 'orchid4', sizes.highlight = 0.1) + NoLegend() + NoAxes()

## Fig. 4D: CD8 subtype proportion barchart
ggplot(cd8_sub_freq, aes(x=sample, y=freq, fill=subtype)) + geom_col(color='black') + theme_classic() + scale_fill_manual(values=c('lightcoral','red','tan4','darkred','black','goldenrod')) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))

## Fig. 4E: Tex z-score heatmap
Heatmap(cd8_tex_zscores, col = rev(brewer.pal(name='RdBu',n=9)), cluster_rows = F, column_names_gp = gpar(fontsize=8), row_names_gp = gpar(fontsize=8))

## Fig. 4F: Tex DEG violin plots
VlnPlot(seu_cd8_clean, 'Tox', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_clean, 'Tigit', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_clean, 'Gzmb', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))

## Fig. 4G: Activated CD8+ TIL z-score heatmap
Heatmap(cd8_act_zscores, col = rev(brewer.pal(name='RdBu',n=9)), cluster_rows = F, column_names_gp = gpar(fontsize=8), row_names_gp = gpar(fontsize=8))

## Fig. 4H: Activated CD8+ TIL DEG violin plots
VlnPlot(seu_cd8_clean, 'Gzmb', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_clean, 'Tnfrsf18', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_clean, 'Tnfrsf9', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))

##############
## Figure 5 ##
##############
load('FIG5_sample_tcr_freq.Robj')
load('FIG5_cd8_clone_mat.Robj')
load('FIG5_clone_counts_df.Robj')
load('FIG5_seu_cd8_tcr.Robj')
load('FIG5_seu_cd8_tcr_sub.Robj')
load('FIG5_subtype_cd8_tcr_freq.Robj')
load('FIG5_tex_prolif_zscores')
load('FIG5_gly_zscores.Robj')
load('FIG5_costim_eff_zscores.Robj')
load('FIG5_tol_zscores.Robj')

## Fig. 5A: Sample-level TCR diversity and clone size bar charts
ggplot(sample_tcr_freq, aes(x=sample, y=freq, fill=cloneSize)) + geom_col(color='black') + theme_classic() + scale_fill_manual(values=c('darkred','maroon','red','lightcoral')) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))

## Fig. 5B: Subtype clone frequency heatmap and count histogram
Heatmap(cd8_clone_mat, col = brewer.pal(name='Reds',n=9), show_row_names = F)
ggplot(clone_counts_df, aes(x=clone, y=log2(ncell))) + geom_col(fill='black',color='black') + theme_classic() + theme(axis.text.x = element_blank())

## Fig. 5C: CD8+ TIL UMAP w/ subtypes and clone sizes
DimPlot(seu_cd8_tcr, cols=c('goldenrod','black','darkred','tan4','red','lightcoral'), group.by = 'subtype') + NoLegend() + NoAxes() + theme(plot.title = element_blank())
DimPlot(seu_cd8_tcr, group.by = 'cloneSize', cols=c('darkred','maroon','lightcoral','red')) + NoLegend() + NoAxes() + theme(plot.title = element_blank())

## Fig. 5D: Subsetted CD8+ TIL UMAP w/ subtypes, clone sizes, and treatment group projections
DimPlot(seu_cd8_tcr_sub, group.by = 'cloneSize', cols=c('darkred','maroon','lightcoral','red'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
DimPlot(seu_cd8_tcr_sub, group.by = 'subtype', cols=c('black','darkred','red'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
DimPlot(seu_cd8_tcr_sub, cells.highlight = colnames(seu_cd8_tcr_sub)[which(seu_cd8_tcr_sub@meta.data$sample == 'hbs')], cols.highlight = 'black', sizes.highlight = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_cd8_tcr_sub, cells.highlight = colnames(seu_cd8_tcr_sub)[which(seu_cd8_tcr_sub@meta.data$sample == 'il2')], cols.highlight = 'dodgerblue3', sizes.highlight = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_cd8_tcr_sub, cells.highlight = colnames(seu_cd8_tcr_sub)[which(seu_cd8_tcr_sub@meta.data$sample == 'tce')], cols.highlight = 'orange2', sizes.highlight = 0.1) + NoLegend() + NoAxes()
DimPlot(seu_cd8_tcr_sub, cells.highlight = colnames(seu_cd8_tcr_sub)[which(seu_cd8_tcr_sub@meta.data$sample == 'tce_il2')], cols.highlight = 'orchid4', sizes.highlight = 0.1) + NoLegend() + NoAxes()

## Fig. 5E: CD8 subtype clone size barchart
ggplot(subtype_cd8_tcr_freq[grep('cd8_act', subtype_cd8_tcr_freq$subtype), ], aes(x=subtype, y=freq, fill=cloneSize)) + geom_col(color='black') + theme_classic() + scale_fill_manual(values=c('darkred','maroon','red','lightcoral')) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
ggplot(subtype_cd8_tcr_freq[grep('prolif', subtype_cd8_tcr_freq$subtype), ], aes(x=subtype, y=freq, fill=cloneSize)) + geom_col(color='black') + theme_classic() + scale_fill_manual(values=c('darkred','maroon','red','lightcoral')) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
ggplot(subtype_cd8_tcr_freq[grep('tex_[h|i|t|]', subtype_cd8_tcr_freq$subtype), ], aes(x=subtype, y=freq, fill=cloneSize)) + geom_col(color='black') + theme_classic() + scale_fill_manual(values=c('darkred','maroon','red','lightcoral')) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))

## Figure 5F: Tex-prolif expanded vs diverse marker gene z-score heatmap
Heatmap(tex_prolif_zscores, col = rev(brewer.pal(name='RdBu',n=9)), cluster_rows = F, column_names_gp = gpar(fontsize=8), row_names_gp = gpar(fontsize=8))

## Figure 5G: Tex-prolif expanded vs diverse marker gene violin plots
VlnPlot(seu_cd8_tcr_sub, 'Tox', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_tcr_sub, 'Id2', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_tcr_sub, 'Pdcd1', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_tcr_sub, 'Pim1', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_tcr_sub, 'Cd69', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
VlnPlot(seu_cd8_tcr_sub, 'Nfkb1', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))

## Figure 5I: Clonally-diverse CD8 TIL glycolysis, effector/co-stim, and tolerance marker z-score heatmaps
Heatmap(gly_zscores, col = rev(brewer.pal(name='RdBu',n=9)), cluster_rows = F, column_names_gp = gpar(fontsize=8), row_names_gp = gpar(fontsize=8))
Heatmap(costim_eff_zscores, col = rev(brewer.pal(name='RdBu',n=9)), cluster_rows = F, column_names_gp = gpar(fontsize=8), row_names_gp = gpar(fontsize=8))
Heatmap(tol_zscores, col = rev(brewer.pal(name='RdBu',n=9)), cluster_rows = F, column_names_gp = gpar(fontsize=8), row_names_gp = gpar(fontsize=8))

##########################
## Supplemental Figures ##
##########################
###############
## Figure S1 ##
###############
load('FIG1_trp2_tact_lys.Robj')
load('FIGS1_dll3_abund.Robj')

## Fig. S1A: TRP2 TCE chromatography traces -- Exported image 
## Fig. S1B: In vitro co-culture gating strategy -- Made in Flowjo
## Fig. S1C: In vitro co-culture analysis ~ Cell lysis @ 32nM, Trp2 Db vs IgG TCE bar charts
ggplot(trp2_tact_lys[which(trp2_tact_lys$conc == 32 & trp2_tact_lys$target != 'none'), ], aes(x=target, y=specific_lysis, fill=ab)) + 
  geom_bar(stat = 'summary', fun = 'mean', color='black', position=position_dodge(width = 1)) + 
  geom_point(size=0.8, position=position_dodge(width = 1)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.5, position=position_dodge(width=1)) +
  theme_classic() + 
  scale_fill_manual(values=alpha(c('grey30','cadetblue4'),0.8)) + 
  theme(legend.position = 'none') + ylim(c(0,100))

## Fig. S1D: In vitro co-culture analysis ~ T-cell activation @ 32nM, Trp2 Db vs IgG TCE bar charts
ggplot(trp2_tact_lys[which(trp2_tact_lys$conc == 32), ], aes(x=target, y=pCD137, fill=ab)) + 
  geom_bar(stat = 'summary', fun = 'mean', color='black', position=position_dodge(width = 1)) + 
  geom_point(size=0.8, position=position_dodge(width = 1)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.5, position=position_dodge(width=1)) +
  theme_classic() + 
  scale_fill_manual(values=alpha(c('grey30','cadetblue4'),0.8)) + 
  theme(legend.position = 'none') + ylim(c(0,100))

## Fig. S1E: DLL3 TCE chromatography trace -- Exported image
## Fig. S1F: DLL3 expression on 16T and KP1 -- Made in flowjo
## Fig. S1G: DLL3 abundance bar chart
ggplot(dll3_abund, aes(x=cell_line, y=abc, fill=cell_line)) + 
  geom_bar(stat = 'summary', fun = 'mean', color='black') + 
  geom_jitter(width = 0.2, size=0.8) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.5) +
  scale_fill_manual(values=alpha(c('plum2','darkorchid2','darkorchid4'),0.8)) +
  theme_classic() + 
  theme(legend.position = 'none')

###############
## Figure S2 ##
###############

## Fig. S2A: Tumor growth dynamics (size) ~ High antigen (MC38) Trp2 TCE vs HBS, individual mice -- PRISM analysis
## Fig. S2B: Tumor growth dynamics (size) ~ Low antigen (B16) 5mg/kg Trp2 TCE vs 10mg/kg Trp2 TCE vs HBS, individual mice -- PRISM analysis
## Fig. S2C: Tumor growth dynamics (size) ~ High antigen (16T) DLL3 TCE vs HBS, individual mice -- PRISM analysis
## Fig. S2D: Tumor growth dynamics (size) ~ Low antigen (KP1) 2mg/kg DLL3 TCE vs 5mg/kg DLL3 TCE vs 10mg/kg DLL3 TCE vs HBS, individual mice -- PRISM analysis
## Fig. S2E: KP1 MHC-I expression +/- IFNg -- PRISM analysis

###############
## Figure S3 ##
###############
# load('FIG3_tumor_growth_b16_trp2_il2.Robj')
# load('FIG3_tumor_growth_kp1_dll3_il2.Robj')
# load('FIG3_tumor_growth_b16_trp2_wtil2.Robj')
load('FIGS3_survival_summary_b16_tce_wtil2.Robj')
load('FIGS3_exh_kp1_dll3_il2.Robj')

## Fig. S3A: Survival curves ~ B16 Trp2 TCE with varying wtIL2 doses
ggplot(survival_summary_b16_tce_wtil2[which(survival_summary_b16_tce_wtil2$treat == 'tce_12ug_wtIL2'),], aes(x=day, y=pAlive)) + 
  geom_step(direction='vh', color='tan4', linewidth=1) + ylim(c(0,1)) +
  theme_classic() + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

ggplot(survival_summary_b16_tce_wtil2[which(survival_summary_b16_tce_wtil2$treat == 'tce_6ug_wtIL2'),], aes(x=day, y=pAlive)) + 
  geom_step(direction='vh', color='goldenrod', linewidth=1) + ylim(c(0,1)) +
  theme_classic() + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

ggplot(survival_summary_b16_tce_wtil2[which(survival_summary_b16_tce_wtil2$treat == 'tce_3ug_wtIL2'),], aes(x=day, y=pAlive)) + 
  geom_step(direction='vh', color='khaki3', linewidth=1) + ylim(c(0,1)) +
  theme_classic() + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

ggplot(survival_summary_b16_tce_wtil2[which(survival_summary_b16_tce_wtil2$treat == 'tce_0_wtIL2'),], aes(x=day, y=pAlive)) + 
  geom_step(direction='vh', color='grey30', linewidth=1) + ylim(c(0,1)) +
  theme_classic() + 
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

## Fig. S3B: Tumor growth dynamics (size) ~ B16 Trp2 TCE vs wtIL2 vs Combo vs HBS, individual mice -- PRISM analysis
## Fig. S3C: Tumor growth dynamics (size) ~ B16 Trp2 TCE vs IL2reh vs Combo vs HBS, individual mice -- PRISM analysis
## Fig. S3D: Tumor growth dynamics (size) ~ KP1 DLL3 TCE vs IL2reh vs Combo vs HBS, individual mice -- PRISM analysis
## Fig. S3E: Exhaustion marker prpoportions, KP1 DLL3 TCE vs IL2reh vs Combo vs HBS
ggplot(exh_kp1_dll3_il2, aes(x=sample, y=pd1_mfi, fill=sample)) +
  geom_boxplot(outliers=F) + geom_jitter(width=0.1) + 
  theme_classic() + scale_fill_manual(values=c('grey30','dodgerblue3','goldenrod','orchid4')) +
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

ggplot(exh_kp1_dll3_il2, aes(x=sample, y=tim3_mfi, fill=sample)) +
  geom_boxplot(outliers=F) + geom_jitter(width=0.1) + 
  theme_classic() + scale_fill_manual(values=c('grey30','dodgerblue3','goldenrod','orchid4')) +
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

ggplot(exh_kp1_dll3_il2, aes(x=sample, y=lag3_mfi, fill=sample)) +
  geom_boxplot(outliers=F) + geom_jitter(width=0.1) + 
  theme_classic() + scale_fill_manual(values=c('grey30','dodgerblue3','goldenrod','orchid4')) +
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

ggplot(exh_kp1_dll3_il2, aes(x=sample, y=cd39_mfi, fill=sample)) +
  geom_boxplot(outliers=F) + geom_jitter(width=0.1) + 
  theme_classic() + scale_fill_manual(values=c('grey30','dodgerblue3','goldenrod','orchid4')) +
  theme(legend.position = 'none', axis.line = element_line(linewidth = 1), axis.ticks = element_line(linewidth = 1.1), axis.ticks.length = unit(0.35, "cm"))

###############
## Figure S4 ##
###############
load('FIGS4_body_weight_b16_tce_wtil2.Robj')
load('FIGS4_body_weight_b16_tce_il2reh.Robj')
load('FIGS4_body_weight_kp1_tce_il2reh.Robj')
load('FIGS4_vitiligo_freq.Robj')

## Fig. S4A: Body weight curves ~ B16 Trp2 TCE vs wtIL2 vs Combo vs HBS 
ggplot(body_weight_b16_tce_wtil2[which(body_weight_b16_tce_wtil2$Group == 'hbs'),], aes(x=Day, y=WeightPercent, color=MouseID)) + 
  geom_line() + theme_classic() + ylim(c(50,130)) + geom_hline(yintercept = 100) + scale_color_manual(values=rep('grey30',3)) +
  stat_summary(aes(group=Group), fun.y=mean, geom='line', linewidth=1.5, color='grey30') +
  stat_summary(aes(group=Group), fun.y=mean, geom='point', size=2, color='grey30') +
  stat_summary(aes(group=Group), geom='errorbar', fun.data='mean_se', width=0.5, color='grey30') + NoLegend()
  
ggplot(body_weight_b16_tce_wtil2[which(body_weight_b16_tce_wtil2$Group == 'wtil2'),], aes(x=Day, y=WeightPercent, color=MouseID)) + 
  geom_line() + theme_classic() + ylim(c(50,130)) + geom_hline(yintercept = 100) + scale_color_manual(values=rep('dodgerblue3',3)) +
  stat_summary(aes(group=Group), fun.y=mean, geom='line', linewidth=1.5, color='dodgerblue3') +
  stat_summary(aes(group=Group), fun.y=mean, geom='point', size=2, color='dodgerblue3') +
  stat_summary(aes(group=Group), geom='errorbar', fun.data='mean_se', width=0.5, color='dodgerblue3') + NoLegend()

ggplot(body_weight_b16_tce_wtil2[which(body_weight_b16_tce_wtil2$Group == 'tce'),], aes(x=Day, y=WeightPercent, color=MouseID)) + 
  geom_line() + theme_classic() + ylim(c(50,130)) + geom_hline(yintercept = 100) + scale_color_manual(values=rep('goldenrod',6)) +
  stat_summary(aes(group=Group), fun.y=mean, geom='line', linewidth=1.5, color='goldenrod') +
  stat_summary(aes(group=Group), fun.y=mean, geom='point', size=2, color='goldenrod') +
  stat_summary(aes(group=Group), geom='errorbar', fun.data='mean_se', width=0.5, color='goldenrod') + NoLegend()

ggplot(body_weight_b16_tce_wtil2[which(body_weight_b16_tce_wtil2$Group == 'tce_wtil2'),], aes(x=Day, y=WeightPercent, color=MouseID)) + 
  geom_line() + theme_classic() + ylim(c(50,130)) + geom_hline(yintercept = 100) + scale_color_manual(values=rep('orchid4',6)) +
  stat_summary(aes(group=Group), fun.y=mean, geom='line', linewidth=1.5, color='orchid4') +
  stat_summary(aes(group=Group), fun.y=mean, geom='point', size=2, color='orchid4') +
  stat_summary(aes(group=Group), geom='errorbar', fun.data='mean_se', width=0.5, color='orchid4') + NoLegend()

## Fig. S4B: Body weight curves ~ B16 Trp2 TCE vs IL2reh vs Combo vs HBS 
ggplot(body_weight_b16_tce_il2reh[which(body_weight_b16_tce_il2reh$Group == 'hbs'),], aes(x=Day, y=WeightPercent, color=MouseID)) + 
  geom_line() + theme_classic() + ylim(c(50,130)) + geom_hline(yintercept = 100) + scale_color_manual(values=rep('grey30',7)) +
  stat_summary(aes(group=Group), fun.y=mean, geom='line', linewidth=1.5, color='grey30') +
  stat_summary(aes(group=Group), fun.y=mean, geom='point', size=2, color='grey30') +
  stat_summary(aes(group=Group), geom='errorbar', fun.data='mean_se', width=0.5, color='grey30') + NoLegend()

ggplot(body_weight_b16_tce_il2reh[which(body_weight_b16_tce_il2reh$Group == 'il2reh'),], aes(x=Day, y=WeightPercent, color=MouseID)) + 
  geom_line() + theme_classic() + ylim(c(50,130)) + geom_hline(yintercept = 100) + scale_color_manual(values=rep('dodgerblue3',7)) +
  stat_summary(aes(group=Group), fun.y=mean, geom='line', linewidth=1.5, color='dodgerblue3') +
  stat_summary(aes(group=Group), fun.y=mean, geom='point', size=2, color='dodgerblue3') +
  stat_summary(aes(group=Group), geom='errorbar', fun.data='mean_se', width=0.5, color='dodgerblue3') + NoLegend()

ggplot(body_weight_b16_tce_il2reh[which(body_weight_b16_tce_il2reh$Group == 'tce'),], aes(x=Day, y=WeightPercent, color=MouseID)) + 
  geom_line() + theme_classic() + ylim(c(50,130)) + geom_hline(yintercept = 100) + scale_color_manual(values=rep('goldenrod',7)) +
  stat_summary(aes(group=Group), fun.y=mean, geom='line', linewidth=1.5, color='goldenrod') +
  stat_summary(aes(group=Group), fun.y=mean, geom='point', size=2, color='goldenrod') +
  stat_summary(aes(group=Group), geom='errorbar', fun.data='mean_se', width=0.5, color='goldenrod') + NoLegend()

ggplot(body_weight_b16_tce_il2reh[which(body_weight_b16_tce_il2reh$Group == 'tce_il2reh'),], aes(x=Day, y=WeightPercent, color=MouseID)) + 
  geom_line() + theme_classic() + ylim(c(50,130)) + geom_hline(yintercept = 100) + scale_color_manual(values=rep('orchid4',7)) +
  stat_summary(aes(group=Group), fun.y=mean, geom='line', linewidth=1.5, color='orchid4') +
  stat_summary(aes(group=Group), fun.y=mean, geom='point', size=2, color='orchid4') +
  stat_summary(aes(group=Group), geom='errorbar', fun.data='mean_se', width=0.5, color='orchid4') + NoLegend()

## Fig. S4C: Body weight curves ~ KP1 DLL3 TCE vs IL2reh vs Combo vs HBS 
ggplot(body_weight_kp1_tce_il2reh[which(body_weight_kp1_tce_il2reh$Group == 'hbs'),], aes(x=Day, y=WeightPercent, color=MouseID)) + 
  geom_line() + theme_classic() + ylim(c(50,130)) + geom_hline(yintercept = 100) + scale_color_manual(values=rep('grey30',10)) +
  stat_summary(aes(group=Group), fun.y=mean, geom='line', linewidth=1.5, color='grey30') +
  stat_summary(aes(group=Group), fun.y=mean, geom='point', size=2, color='grey30') +
  stat_summary(aes(group=Group), geom='errorbar', fun.data='mean_se', width=0.5, color='grey30') + NoLegend()

ggplot(body_weight_kp1_tce_il2reh[which(body_weight_kp1_tce_il2reh$Group == 'il2reh'),], aes(x=Day, y=WeightPercent, color=MouseID)) + 
  geom_line() + theme_classic() + ylim(c(50,130)) + geom_hline(yintercept = 100) + scale_color_manual(values=rep('dodgerblue3',10)) +
  stat_summary(aes(group=Group), fun.y=mean, geom='line', linewidth=1.5, color='dodgerblue3') +
  stat_summary(aes(group=Group), fun.y=mean, geom='point', size=2, color='dodgerblue3') +
  stat_summary(aes(group=Group), geom='errorbar', fun.data='mean_se', width=0.5, color='dodgerblue3') + NoLegend()

ggplot(body_weight_kp1_tce_il2reh[which(body_weight_kp1_tce_il2reh$Group == 'tce'),], aes(x=Day, y=WeightPercent, color=MouseID)) + 
  geom_line() + theme_classic() + ylim(c(50,130)) + geom_hline(yintercept = 100) + scale_color_manual(values=rep('goldenrod',10)) +
  stat_summary(aes(group=Group), fun.y=mean, geom='line', linewidth=1.5, color='goldenrod') +
  stat_summary(aes(group=Group), fun.y=mean, geom='point', size=2, color='goldenrod') +
  stat_summary(aes(group=Group), geom='errorbar', fun.data='mean_se', width=0.5, color='goldenrod') + NoLegend()

ggplot(body_weight_kp1_tce_il2reh[which(body_weight_kp1_tce_il2reh$Group == 'tce_il2reh'),], aes(x=Day, y=WeightPercent, color=MouseID)) + 
  geom_line() + theme_classic() + ylim(c(50,130)) + geom_hline(yintercept = 100) + scale_color_manual(values=rep('orchid4',10)) +
  stat_summary(aes(group=Group), fun.y=mean, geom='line', linewidth=1.5, color='orchid4') +
  stat_summary(aes(group=Group), fun.y=mean, geom='point', size=2, color='orchid4') +
  stat_summary(aes(group=Group), geom='errorbar', fun.data='mean_se', width=0.5, color='orchid4') + NoLegend()

## Fig. S4D: Mouse vitiligo pictures -- Exported image
## Fig. S4E: Vitiligo frequency bar charts -- B16 Trp2 TCE vs wtIL2/IL2reh vs Combo vs HBS
ggplot(vitiligo_freq, aes(x=Group, y=VitFrequency)) + 
  geom_bar(stat = 'summary', fun = 'mean', color='black', fill=alpha('black',0.75)) + 
  geom_point() + theme_classic()
  
###############
## Figure S5 ##
###############
load('FIG4_seu_cd4_clean.Robj')
load('FIG4_seu_cd8_clean.Robj')
load('FIGS5_cd4_genes.Robj')
load('FIGS5_cd8_genes.Robj')
load('FIGS5_cd4_cd8_freq.Robj')
load('FIG4_cd8_sub_freq.Robj')
load('FIG5_seu_cd8_tcr_sub.Robj')

## Fig. S5A: CD4 and CD8 subtype marker dot plots
DotPlot(seu_cd4_clean, features=cd4_genes, cols='RdBu', dot.scale = 3) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))
DotPlot(seu_cd8_clean, features=cd8_genes, cols='RdBu', dot.scale = 3) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))

## Fig. S5B: CD4/CD8 proportion bar chart
ggplot(cd4_cd8_freq, aes(x=sample, y=freq, fill=cd4cd8)) + geom_col(color='black') + theme_classic() + scale_fill_manual(values=c('black','grey40')) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))

## Fig. S5C: Activated CD8+ TIL and Tex proportion staggered barcharts
ggplot(cd8_sub_freq[which(cd8_sub_freq$subtype %in% c('cd8_tex','cd8_act') & cd8_sub_freq$sample != 'tce'), ], aes(x=subtype, y=freq, fill=sample)) + geom_col(color='black', position=position_dodge()) + theme_classic() + scale_fill_manual(values=alpha(c('black','dodgerblue3','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))

## Fig. S5D: Glycolysis gene score violin plot
VlnPlot(seu_cd8_clean, 'glycolysis1', idents = c('hbs','il2','tce','tce_il2'), pt.size = 0, cols=alpha(c('black','dodgerblue3','orange2','orchid4'),0.8)) + NoLegend() + theme(axis.line = element_line(size=0.5), axis.ticks = element_line(size=0.75), axis.ticks.length=unit(.2, "cm"))

## Fig. S5E: Cxcr6 expression across CD8 TIL subtypes
VlnPlot(seu_cd8_clean, 'Cxcr6', pt.size = 0, cols=c('lightcoral','red','maroon','darkred')) + NoLegend() + theme(plot.title = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.text.x = element_blank()) 

## Fig. S5F: Cxcr6 expression across Tex prolif clonality bins
VlnPlot(seu_cd8_tcr_sub, idents=c('cd8_tex_prolif_rare','cd8_tex_prolif_small','cd8_tex_prolif_medium','cd8_tex_prolif_large'), 'Cxcr6', pt.size = 0, cols=c('lightcoral','red','maroon','darkred')) + NoLegend() + theme(plot.title = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.text.x = element_blank()) 