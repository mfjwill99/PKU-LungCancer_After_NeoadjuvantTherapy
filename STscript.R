
library(Seurat)
library(SPATA2)
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(SingleCellExperiment)
library(tidyverse)
library("CellTrek")
library("dplyr")
library(magrittr)




## SpatialDimPlot ####
p <- SpatialDimPlot(data_celltrek) +
  scale_fill_manual(values = palette)
ggsave(paste('./celltrek', paste('Celltrek1', ".jpeg", sep=""), sep = '/'), p, height = 4.5, width=8)




##  Colocalization Network ####

levels(data_celltrek)

# PT
# 1 PT+
# 4 PT-

glut_cell <- c('CD8Tem','CD4Tem','Treg','ExhaustedCD4Trm','CD4Tcm','CD8Teff',
               'gdT','NaiveCD8T','ExhaustedCD8Trm',
               
               'ResMac','CD11c_Mac','CD11b_Mac','PLAC8_Mac',
               'MKI67_Mac','CTSK_Mac',
               
               'IGHD_Bcell','CRIP1_Bcell','IGHM_Bcell','HSPA1A_Bcell','TEX9_Bcell',
               'FCRL5_Bcell','RGS13_Bcell')


# 2 mLN-
# 3 mLN-

glut_cell <- c('CD8Tem','CD4Tem','Treg','ExhaustedCD4Trm','CD4Tcm','CD8Teff',
               'gdT','NaiveCD8T','ExhaustedCD8Trm',
               
               'IGHD_Bcell','CRIP1_Bcell','IGHM_Bcell','HSPA1A_Bcell','TEX9_Bcell',
               'FCRL5_Bcell','RGS13_Bcell')




brain_celltrek_glut <- subset(data_celltrek, idents = glut_cell)
brain_celltrek_glut$celltype.all <- Idents(brain_celltrek_glut)

brain_sgraph_KL <- CellTrek::scoloc(brain_celltrek_glut, col_cell='celltype.all', use_method='KL', eps=1e-50)


## We extract the minimum spanning tree (MST) result from the graph
brain_sgraph_KL_mst_cons <- brain_sgraph_KL$mst_cons
brain_sgraph_KL_mst_cons[match(glut_cell, rownames(brain_sgraph_KL_mst_cons)),
                         match(glut_cell, colnames(brain_sgraph_KL_mst_cons))]  %>% write.csv('ST1_net1.csv')


mat <- read.csv('ST1_net1.csv', row.names = 1)

library(igraph)
net <- graph_from_adjacency_matrix(as.matrix(mat), 
                                   mode = "undirected", 
                                   weighted = TRUE, 
                                   diag = FALSE)
isolated_nodes <- which(degree(net, mode = "all") == 0)
if (length(isolated_nodes) > 0) {
  net <- delete_vertices(net, isolated_nodes)
}

colors <- palette[names(palette) %in% V(net)$name]
V(net)$color <- colors

E(net)$width <- E(net)$weight * 10
set.seed(123) # 保证布局可重复

pdf('ST1_Net1.pdf')
plot(net, 
     vertex.label.color = "black",  # 标签颜色
     vertex.label.dist = 1.5,   
     vertex.label.cex = 1.2,    # 标签大小
     vertex.size = 8,          # 节点大小
     layout = layout_with_fr,   # Fruchterman-Reingold布局
     main = "Cell-Cell Colocalization Network")
dev.off()



