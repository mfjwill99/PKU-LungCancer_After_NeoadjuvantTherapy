suppressMessages(library(Seurat))
suppressMessages(library(monocle))
library(DDRTree)
library(tidyverse)
library(ggthemes)
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
library(ggplot2)
library(paletteer)
library(nichenetr)
library(ggsignif)
library(circlize)



p <- DimPlot(data_new, reduction = 'umap') +
  scale_color_manual(values = palette) 
p
ggsave(paste(output.dir, paste('Umap_sub', ".jpeg", sep=""), sep = '/'), p, height = 4.2, width=8)

p <- DimPlot(data_new, reduction = 'umap', split.by = 'Group') +
  scale_color_manual(values = palette) 
p
ggsave(paste(output.dir, paste('Umap_sub_group', ".jpeg", sep=""), sep = '/'), p, height = 4, width=15)



markers <- c('CD3E', 'CD8A','CD4','FOXP3','TIGIT','SELL','IL7R','CCR7','CD69','TRDC',
             'CD68','FOLR2', 'ITGAX','ITGAM','PLAC8','MKI67','CTSK','CD79A','CD27','GPR183',
             'IGHD','CRIP1','IGHM','HSPA1A','TEX9','FCRL5','RGS13','GNLY','EPCAM','COL1A1','S100A8','PECAM1','JCHAIN',
             'LAMP3','TPSB2','ACTA2','CLEC4C')

pdf(paste(output.dir, paste('heatmap', ".pdf", sep=""), sep = '/'), height = 6.5, width=7.5)
AverageHeatmap(object = data_new,
               markerGene =  markers,
               column_names_rot = 60,
               annoCol = F,
               #  myanCol = palette,
               row_title = "",
               cluster_rows = F) 
dev.off()




## select MF
MF <- obj.subset
features <- c('FABP4','MCEMP1',
              'C1QA','C1QB','APOE','SIRPA','SPP1','FOLR2','FOLR3',
              'HLA-DRA', 'HLA-DQA1', 'HLA-DPA1','CD74','PLAC8',
              'TNF','CXCL2', 'CXCL5','CXCL9','CXCL10','CXCL11','IL10','IDO1','CTSK','ITGAX','ITGAM','MKI67','TOP2A')
row_split <- c(rep('AlveolarMps',2), rep('Immunosuppressive',7), 
               rep('AntigenPresenting',5), rep('Inflammatory',11),rep('Proliferating',2))


row_split = factor(row_split, levels = c('AlveolarMps','Immunosuppressive','AntigenPresenting','Inflammatory','Proliferating') )
source("Heat_Dot_data.R")
### set colnames order
plot_ord <- levels(obj.subset)

MF <- obj.subset
levels(MF)
data.plot <- Heat_Dot_data(object=MF, features=features)


# 
exp.mat <- data.plot %>% select(features.plot,id,avg.exp.scaled) %>% spread(id,avg.exp.scaled)
rownames(exp.mat) <- exp.mat$features.plot
exp.mat$features.plot <- NULL
exp.mat <- exp.mat[,plot_ord]
per.mat <- data.plot %>% select(features.plot,id,pct.exp) %>% spread(id,pct.exp)
rownames(per.mat) <- per.mat$features.plot
per.mat$features.plot <- NULL
per.mat <- per.mat[,plot_ord]/100

### plot heatmap
library(ComplexHeatmap)
library(circlize) ## color 
col_fun <- colorRamp2(c(-1.5, 0, 2.5), c("#118ab2", "#fdffb6", "#e63946"))
# split heatmap
col_split = c(rep("Monocyte",2),rep("Macrophage",6))
col_split =factor(col_split,levels = c("Monocyte","Macrophage"))
cols.ha <- palette[41:45]
names(cols.ha) <- levels(row_split)
ha = HeatmapAnnotation(df = data.frame(Marker=row_split),which = "row",
                       col = list(Marker = cols.ha))

pdf(paste0(dir,"/MP_maker_heat_CCI.pdf"),width = 5.5,height = 11)
Heatmap(exp.mat, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill){
          grid.rect(x = x, y = y, width = width, height = height,gp = gpar(col = "grey", fill = NA))
          grid.circle(x = x, y = y,r=per.mat[i,j]/2 * max(unit.c(width, height)),
                      gp = gpar(fill = col_fun(exp.mat[i, j]), col = NA))},
        column_names_side = "top",row_names_side = "right",
        row_split = row_split, # column_split = col_split,
        row_gap = unit(3, "mm"),column_gap =unit(3, "mm"), 
        left_annotation = ha,
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()






##### ----------- pseudotime ------------##########

sample.markers <- FindAllMarkers(obj.subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

i = 200

data <- as(as.matrix(obj.subset@assays$integrated@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = obj.subset@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

my_cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size()
)

my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds)
my_cds <- detectGenes(my_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(my_cds),
                                    num_cells_expressed >= 10))
print(head(pData(my_cds)))

# #### step5:  ####
disp_table <- dispersionTable(my_cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id

my_cds2 <- setOrderingFilter(my_cds, ordering_genes = disp.genes)
my_cds2 <- estimateSizeFactors(my_cds2)
my_cds2 <- reduceDimension(my_cds2, pseudo_expr = 1, verbose=T)
my_cds2 <- orderCells(my_cds2, reverse = TRUE)


#### step6: visualize the trajectory 
cols <- ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value
pdf(paste(output.dir, paste('4.3_pseudotime_MP_', i, 'Seurat.pdf', sep = ''),  sep = '/'), width = 6.3, height = 5)
print(
  plot_cell_trajectory(my_cds2[, ], 
                       color_by = 'subcluster_new', 
                       show_branch_points = TRUE,
                       alpha = 0.4,
                       cell_size = 0.1,
                       cell_link_size = 0.5) + 
    scale_color_manual(values = alpha(cols, 0.4)) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 15),
          legend.position = 'right') +
    guides(colour = guide_legend(override.aes = list(size = 3)))
)
dev.off()

pdf(paste(output.dir, paste('4.3_pseudotime_MP_', i, '2Seurat3.pdf', sep = ''),  sep = '/'), width = 4.3, height = 3)
print(
  plot_cell_trajectory(my_cds2[, ], 
                       color_by = 'subcluster_new', 
                       show_branch_points = TRUE,
                       alpha = 0.4,
                       cell_size = 0.1,
                       cell_link_size = 0.5) + 
    scale_color_manual(values = alpha(cols, 0.4)) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 15),
          legend.position = 'right') +
    guides(colour = guide_legend(override.aes = list(size = 3)))
)
dev.off()







### gene heatmap
cells_subset <- colnames(my_cds2)

Binner <- function(cells_subset){
  df <- data.frame(pData(my_cds2[,cells_subset]))
  df <- df[,c("Pseudotime", "State",'subcluster_new','Group')]
  df <- df[order(df$Pseudotime, decreasing = F),]
  len <- length(df$Pseudotime)
  bin<-round(len/100)
  subcluster_new <- c()
  Group <- c()
  State <- c()
  value1 <- c()
  value2 <- c()
  value3 <- c()
  for(i in 0:99){
    if(i < 99){
      start <- 1+(bin*i)
      stop <- bin+(bin*i)
      value1 <- median(as.numeric(as.vector(df$State[c(start:stop)])))
      State <- c(State, value1)
      
      value2 <- median(as.vector(df$subcluster_new[c(start:stop)]))
      subcluster_new <- c(subcluster_new, value2)
      
      value3 <- median(as.vector(df$Group[c(start:stop)]))
      Group <- c(Group, value3)
    }
    else{
      State <- c(State, value1)
      subcluster_new <- c(subcluster_new, value2)
      Group <- c(Group, value3)
    }
  }
  return(data.frame(State,subcluster_new,Group))
}
bin <- Binner(colnames(my_cds2))
bin$State <- as.character(bin$State)

col_group <- ggthemes_data$tableau$`color-palettes`$regular$`Color Blind`$value[1:4]
names(col_group) <- c('Normal','T-','N-','T+')
col_state <- ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value
names(col_state) <- 1:10

ann_colors <- list(
  Group = col_group,
  subcluster_new = palette,
  State = col_state
)




ct <- 'ResMac'

my_cds2_sub <- my_cds2[, my_cds2$subcluster_new %in% ct]
diff_test_res <- differentialGeneTest(my_cds2_sub,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01))

p <- plot_pseudotime_heatmap(my_cds2_sub[sig_gene_names,],
                             num_clusters = 2,
                             cores = 1, #add_annotation_col = bin,
                             show_rownames = T, return_heatmap = T)
p
ggsave(paste(output.dir, paste('Pseudotime_HP_',ct, "final.pdf", sep=""), sep = '/'), p, height = 2.3, width=4)


colnames(plotdata)
plotdata <- plotdata %>% filter(Group %in% c('T-','T+'))

col_group <- ggthemes::ggthemes_data$tableau$`color-palettes`$regular$`Color Blind`$value
p <- ggplot(plotdata,aes(x = Pseudotime,color =Group,  fill = Group)) +
  geom_density(alpha = 0.1) +
  scale_fill_manual(values = col_group)+
  scale_color_manual(values = col_group) +
  theme_classic() 
p
ggsave(paste(output.dir, paste('Pseudotime_density_',ct, "final.jpg", sep=""), sep = '/'), p, height = 1.2, width=4)







## nichenet ##############
ident = 'nichenet'
output.dir <- paste(dir, ident, sep = '/')
if(! dir.exists(output.dir)) {dir.create(output.dir)}


### nichenet ####
load('ligand_target_matrix.Rdata')
load('lr_network.Rdata')
load('weighted_networks.Rdata')



sender_celltypes =  c('ResMac') 
#### GOI receiver: find all markers for each cluster
receiver = c('CD4Tem') 
receiver = c('CD8Tem')

receiver = c('CD4Tem','CD8Tem')
n = 200

Idents(data_new) <- data_new$celltype.all
temp <- subset(data_new, idents = receiver)
Idents(temp) <- temp$Group
markers.t <- FindMarkers(temp, ident.1 = 'T-', ident.2 = 'T+')


geneset_oi <- markers.t %>% filter(avg_log2FC > 0) %>% row.names()
geneset_oi <- intersect(geneset_oi, id.gene)



## receiver
list_expressed_genes_receiver = receiver %>% unique() %>% lapply(get_expressed_genes, data_new, 0.10, assay_oi = 'RNA') # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_receiver = list_expressed_genes_receiver %>% unlist() %>% unique()

background_expressed_genes = expressed_genes_receiver %>% 
  .[. %in% rownames(ligand_target_matrix)]

## sender
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, data_new, 0.10, assay_oi = 'RNA') # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()



#### STEP6: Define a set of potential ligands ####
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% 
  filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
  pull(from) %>%  
  unique()

#### STRP7: Perform NicheNet ligand activity analysis ####
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank( dplyr::desc(pearson)))

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
best_upstream_ligands <- union(best_upstream_ligands, ligand_activities$test_ligand[grep('^CC|^CX', ligand_activities$test_ligand)])


### Heatmap ####
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network
ggsave(paste(output.dir, paste('2.2_heatmap_', paste0(receiver, collapse = '.'), ".jpeg", sep=""), sep = '/'), p_ligand_target_network, height = 6, width=8)
write.csv(vis_ligand_target, file = paste(output.dir, paste0(paste0(receiver, collapse = '.'), '.csv'), sep = '/'))


  
  





