library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)
library(ggthemes)
library(ggpubr)
library(plyr)
library(ggbeeswarm)
library(scales)
library(cowplot)
library(viridis)
library(scater)
library(BiocParallel)
library(harmony)
library(SingleCellExperiment)
library(SpatialExperiment)
library(tidyverse)
library(imcRtools)
library(diffcyt)
library(dittoSeq)
library(bluster)
library(argparser)
library(FastPG)
library(glue)
set.seed(100)




platform<-'IMC'
csv.dir<-'../CSV/'
info.file<-'../info.csv'
marker.file<-'../marker.csv'
cofactor<-1
q.nor.n=101


# output dir
output.dir<-'./output2501'
if(! dir.exists(output.dir)) {dir.create(output.dir)}


# # demo stat
# group.stat.list <- list(group1 = list(stat = 't.test', stat.label = list(  c('Neg','Pos'))))
# group.label.list <- list(group1 = c('Neg','Pos'))

## 0. QC ##########################################################
print('start section 0: QC')
output.dir0 <-  paste(output.dir, '0_QC', sep = '/')
if(! dir.exists(output.dir0)) {dir.create(output.dir0)}

## read marker 
panel <- read.csv(marker.file)
colnames(panel) <- c('channel', 'marker','lineager')


# read csv
files <- list.files(csv.dir ,pattern='.csv$', full=TRUE)
#x <- files[1]
csv.list <- lapply(files, function(x) {
  temp <- read.csv(x) # %>% select(!contains('DAPI'))
  col.temp <- colnames(temp)[grep('X[0-9]+', colnames(temp))] 
  col.channel <- str_split( col.temp , pattern = '_', simplify = T)[,1] %>% str_sub(2, nchar(.))
  colnames(temp)[grep('X[0-9]+', colnames(temp))]  <- panel$marker[match(col.channel, panel$channel)]
  temp <- temp[, !is.na(colnames(temp))]
  temp
})
names(csv.list) <- files


cur_features <- do.call('rbind.fill', csv.list)


## read sample info ####
info <- read.csv(info.file)
colnames(info)[1] <- 'sample_id'


if (!'roi_location' %in% (info %>% colnames())){
  filess <- list.files(csv.dir ,pattern='.csv$')
  
  TName <- str_extract_all(filess,pattern = 'T[0-9]+', simplify = T)
  ROIName <- str_extract_all(filess,pattern = 'ROI_[A-Z0-9]+', simplify = T)
  ROIID <- str_extract_all(filess,pattern = 'ROI[0-9]+', simplify = T)
  
  matchdata<-data.frame(roi_id=paste0(TName,'_',str_replace_all(ROIID[,1],'_','')),
                        roi_location=paste0(TName,'_',str_replace_all(ROIName[,1],'_','')))
  info<-merge(info,matchdata)
}



### 0.1. creat spe object ####
marker <- panel %>% filter(! lineager == 0) %>% pull(marker)
counts <- cur_features %>% select(marker) 

meta <- cur_features %>% select(roi, AreaID,  CellID) %>%
  mutate(roi_id = roi, .keep = 'unused') %>%
  left_join(info) ### add sample info
coords <- cur_features %>% select(contains('position'))
colnames(coords) <- c("Pos_X", "Pos_Y")



spe <- SpatialExperiment(assays = list(counts = t(counts)),
                         colData = meta, 
                         sample_id = as.character(meta$sample_id), 
                         image_id = as.character(meta$roi_id),
                         spatialCoords = as.matrix(coords))
colnames(spe) <- paste0(spe$roi_id, ".", spe$CellID)


# define channel for cluster
rowData(spe)$use_channel <- rownames(spe) %in% (panel %>% filter(lineager == '1') %>% pull(marker))

### 0.2 counts transformation and normalization ####

# transformation
counts.tsf <- asinh(counts(spe)/cofactor)
# normalization
counts.nor <- counts.tsf
for (i in 1:nrow(counts.nor)) {
  print(i)
  counts.nor[i,] <- q.nor(x = counts.nor[i,], n = q.nor.n)
  counts.nor[i,] <- minMax(counts.nor[i,])
}
assay(spe, "exprs") <- counts.nor



### 0.3 batch correction ####
set.seed(220225)


# before
spe <- runUMAP(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs")


# svd分解
mat <- t(assay(spe, "exprs")) [,rowData(spe)$use_channel]
svd_res <- svd(mat)
variance_explained <- (svd_res$d^2)/sum(svd_res$d^2)
npc <- which(cumsum(variance_explained) >= 0.95)[1]

if(length(unique(spe$batch_id)) < 2) {
  # PCA: after without correction
  harmony_emb <- prcomp(t(mat)) ### batch_id
  reducedDim(spe, "harmony") <- harmony_emb$rotation[, 1:npc]
  spe <- runUMAP(spe, dimred = "harmony", name = "UMAP_harmony")
  
} else {
  # HARMONY: after with correction
  harmony_emb <- HarmonyMatrix(mat, as.factor(spe$batch_id), do_pca = T, npcs = npc) ### batch_id
  reducedDim(spe, "harmony") <- harmony_emb
  spe <- runUMAP(spe, dimred = "harmony", name = "UMAP_harmony")
}




## 1. reduction and cluster #############################################################
output.dir1 <-  paste(output.dir, '1_Cluster', sep = '/')
if(! dir.exists(output.dir1)) {dir.create(output.dir1)}
dir.create(output.dir1)


# 
### 1.1 cluster ####
expr.mat <- reducedDim(spe, "harmony")
cluster_PhenoGraph <- FastPG::fastCluster(as.matrix(expr.mat),k = 30, num_threads = 100)
spe$pg_cluster <- as.factor(cluster_PhenoGraph$communities)


cols_clst <- ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value
cols_clst <- colorRampPalette(cols_clst)(length(levels(spe$pg_cluster)))
names(cols_clst) <- levels(spe$pg_cluster)


colnames(spe) <- paste(colnames(spe), 1:ncol(spe))

temp <- data.frame(spe$pg_cluster, reducedDim(spe, "UMAP_harmony"))
colnames(temp) <- c('pg_cluster','UMAP_1','UMAP_2')
lc.cent = temp %>% group_by(pg_cluster) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarise_all(median)
p <- dittoDimPlot(spe, var = "pg_cluster", 
                  reduction.use = "UMAP_harmony", size = 0.2) + 
  # ggtitle("Phenograph clusters expression on UMAP") +
  guides(color=guide_legend(ncol=2, override.aes = list(size = 3))) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.position = 'none') +
  geom_text(data = lc.cent, aes(x = UMAP_1, y = UMAP_2, label = as.character(pg_cluster))) +
  scale_color_manual(values = cols_clst) +
  labs(x = 'UMAP1', y = 'UMAP2', title = '')
p
ggsave(paste(output.dir1, '/1.1_UMAP/Umap.num.pdf', sep = '/'), p, height = 7, width=7)
ggsave(paste(output.dir1, '/1.1_UMAP/Umap.num.png', sep = '/'), p, height = 7, width=7,dpi=300)




###### 3.1 Annotation #############################################################
output.dir1 <-  paste(output.dir, '2_CellType', sep = '/')
if(! dir.exists(output.dir1)) {dir.create(output.dir1)}



spe$celltype<-mapvalues(spe$pg_cluster,from = celltype$new_pg_cluster,to = celltype$celltype)
spe <- spe[, ! spe$celltype == 'rm']
spe$celltype<-factor(spe$celltype,levels = unique(celltype$celltype)[which(unique(celltype$celltype)!='rm')])

spe$new_pg_cluster<-mapvalues(spe$pg_cluster,from = celltype$pg_cluster,to = celltype$new_pg_cluster)
celltype_number<-unique(celltype$celltype)[which(unique(celltype$celltype)!='rm')] %>% length
print(celltype_number)

celltype_umap<-paste0(c(1:celltype_number),':',unique(celltype$celltype)[which(unique(celltype$celltype)!='rm')])
spe$celltype_umap<-mapvalues(spe$celltype,from = unique(celltype$celltype)[which(unique(celltype$celltype)!='rm')],to = celltype_umap)
spe$celltype_umap<-factor(spe$celltype_umap,levels = celltype_umap[!grepl('rm',celltype_umap)])


cols_clst <- ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value
cols_clst <- colorRampPalette(cols_clst)(length(levels(spe$celltype)))
names(cols_clst) <- levels(spe$celltype)

cols_clst_umap <- ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value
cols_clst_umap <- colorRampPalette(cols_clst)(length(levels(spe$celltype_umap)))
names(cols_clst_umap) <- levels(spe$celltype_umap)
spe$pg_cluster_umap<-str_split_fixed(spe$celltype_umap,pattern = ':',n = 2)[,1]
spe$pg_cluster_umap<-factor(spe$pg_cluster_umap %>% as.numeric())


temp <- data.frame(spe$pg_cluster_umap, reducedDim(spe, "UMAP_harmony"))
colnames(temp) <- c('new.id','UMAP_1','UMAP_2')
lc.cent = temp %>% group_by(new.id) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarise_all(median)
p <- dittoDimPlot(spe, var = "celltype_umap", 
                  reduction.use = "UMAP_harmony", size = 0.2) + 
  ggtitle("") +
  guides(color=guide_legend(ncol=2, override.aes = list(size = 2))) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_blank()) +
  scale_color_manual(values = cols_clst_umap) + coord_fixed()+
  geom_text(data = lc.cent, aes(x = UMAP_1, y = UMAP_2, label = as.character(new.id))) +
  labs(x = 'UMAP1', y = 'UMAP2')
dir.create(paste0(output.dir1,'/2.1_UMAP'))

ggsave(paste(output.dir1, '/2.1_UMAP/UmapAnno.num.pdf', sep = '/'), p, height = 4, width=6)
ggsave(paste(output.dir1, '/2.1_UMAP/UmapAnno.num.png', sep = '/'), p, height = 4, width=6,dpi=300)




### 2.3 marker median heatmap ####
hmdata0 <- cbind(t(assay(spe, "exprs"))[,marker], Anno = spe$celltype) %>% 
  as_tibble() %>% group_by(Anno) %>% 
  summarise_all(median) %>% column_to_rownames(var = 'Anno')

cols <- colorRampPalette(rev(brewer.pal(7, 'RdYlBu')))(101)
# Heatmap
temp <- spe$celltype %>% table() %>% as.data.frame() 
row_ha = rowAnnotation(counts = anno_barplot(temp$Freq, gp = gpar(fill = cols_clst)))
rownames(hmdata0) <- spe$celltype %>% levels()

png(paste(output.dir1, '2.2_Heatmap/heatmap.median.anno.png', sep = '/'), width = 7000, height =3500,res=72*15)
Heatmap(hmdata0, 
        col =  cols,
        name = 'Expression value',
        cluster_columns = F,
        cluster_rows = F,
        #  heatmap_width = unit(19, 'cm'),
        #  heatmap_height = unit(12, 'cm'),
        
        right_annotation = row_ha,
        heatmap_legend_param = list(legend_height = unit(4, "cm"),title_position = "lefttop-rot")
)
dev.off()






#### marker ######

stat.dir.sub <- paste(output.dir,"/marker.dot",sep = "")
if(! dir.exists(stat.dir.sub)) {dir.create(stat.dir.sub)}

marker.fun <- c("B7H4", "CD137", "LAG3", "PD1", "PD.L1" ,"TIGIT",'TIM3','VISTA')

tempp <- colData(spe)
exprtempp <- t(assay(spe, "exprs")) %>% as.data.frame()
exprtempp <- cbind(exprtempp,celltype = tempp$celltype)


## all marker
data <- exprtempp
col <- colnames(data)[1:(ncol(data)-1)]

plot_data <- data %>% as_tibble %>% group_by(celltype)  %>% summarise_all(median) %>% 
  pivot_longer(cols = col,  names_to = c("Marker"), values_to = "Expression")
plot_data <- plot_data %>% filter( celltype %in% c('ProlT', 'NK','ProlNK'))
head(plot_data)


plot_data$Marker <- factor(plot_data$Marker, levels = c(marker.fun, setdiff(rownames(spe), marker.fun)))

p=ggplot(plot_data,aes(x = Marker, 
                       y = celltype, # 按照富集度大小排序
                       size = Expression, 
                       color = celltype)) +
  geom_point(shape = 16)+ 
  labs(x = '', y = '')+           # 设置x，y轴的名称
  scale_color_calc()+
  # scale_colour_manual(values = cols_clst) +
  theme_bw() + 
  theme(title=element_text(size=12, face = 'bold'))+ 
  theme(panel.grid = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_text(angle=90,hjust = 1, size = 12,color = "black"),
        axis.text.y=element_text(size = 12,color = "black")) 
p
ggsave(paste(stat.dir.sub, paste('marker.dot.NKT', ".pdf", sep=""), sep = '/'), p, height = 2.5, width=8.5)
ggsave(paste(stat.dir.sub, paste('marker.dot.NKT.legend', ".pdf", sep=""), sep = '/'), p, height = 8.5, width=12.5)





## CN ####
spe <- buildSpatialGraph(spe, img_id = "roi_id", type = "knn", k = knn_k)
spe <- aggregateNeighbors(spe, colPairName = "knn_interaction_graph", 
                          aggregate_by = "metadata", count_by = "celltype")


roi.num <- spe$roi_id %>% unique() %>% length

CN.num <- 15
## CN number #####
for(CN.num in c(14)) {
  print(CN.num)
  # CN.num = as.numeric(CN.num)
  
  dir.create(paste0(output.dir3,'/CN_',CN.num))
  dir.create(paste0(output.dir3,'/CN_',CN.num,'/4.3_CN_Heatmap','_',CN.num))
  
  cn_1 <- kmeans(spe$aggregatedNeighbors, centers = CN.num)
  
  spe$cn_celltypes <- as.factor(cn_1$cluster)
  
  saveRDS(spe, file = paste0(output.dir3,'/CN_',CN.num,'.rds'))  # *
  
  if (is.null(spe$roi_location)){
    spe$roi_location <- spe$roi_id
  }
  
  p <- plotSpatial(spe,
                   node_color_by = "cn_celltypes",
                   img_id = "roi_id",
                   ncols = ceiling(sqrt(roi.num)),
                   nrows = ceiling(sqrt(roi.num)),
                   node_size_fix = 0.5) +
    scale_color_tableau('Classic 20')
  ggsave(paste(output.dir3, paste0('4_CN_',CN.num,'.pdf'), sep = '/'), p, height = max(ceiling(sqrt(roi.num))*2, 10), width=max(ceiling(sqrt(roi.num))*2, 10))
  
  
  ## heatmap
  for_plot <- colData(spe) %>% as_tibble() %>%
    group_by(cn_celltypes, celltype) %>%
    dplyr::summarize(count = n()) %>%
    mutate(freq = count / sum(count)) %>%
    pivot_wider(id_cols = cn_celltypes, names_from = celltype, 
                values_from = freq, values_fill = 0) %>%
    ungroup() %>%
    select(-cn_celltypes)
  rownames(for_plot) <- paste0('CN', rownames(for_plot))
  for_plot %>% write.csv(paste(output.dir3, paste0('/CN_',CN.num,'/4.3_CN_Heatmap','_',CN.num,'/','CN_heatmap_',CN.num,'_rowscale.csv'), sep = '/'))
  
  p <- pheatmap::pheatmap(for_plot, color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
                          scale = "row")
  ggsave(paste(output.dir3, paste0('/CN_',CN.num,'/4.3_CN_Heatmap','_',CN.num,'/','CN_heatmap_',CN.num,'_rowscale.pdf'), sep = '/'), p, height = 5, width = 2.5 + ncol(for_plot)/10)
  ggsave(paste(output.dir3, paste0('/CN_',CN.num,'/4.3_CN_Heatmap','_',CN.num,'/','CN_heatmap_',CN.num,'_rowscale.png'), sep = '/'), p, height = 5, width = 2.5 + ncol(for_plot)/10,dpi=300)
  
  #paste0(output.dir3,'/CN_',CN.num,'/4.3_CN_Heatmap','_',CN.num)
  
  for_plot <- colData(spe) %>% as_tibble() %>%
    group_by( celltype,cn_celltypes) %>%
    dplyr::summarize(count = n()) %>%
    mutate(freq = count / sum(count)) %>%
    pivot_wider(id_cols = cn_celltypes, names_from = celltype, 
                values_from = freq, values_fill = 0) %>%
    ungroup() %>%
    select(-cn_celltypes)
  rownames(for_plot) <- paste0('CN', rownames(for_plot))
  for_plot %>% write.csv(paste(output.dir3, paste0('/CN_',CN.num,'/4.3_CN_Heatmap','_',CN.num,'/','CN_heatmap_',CN.num,'_colscale.csv'), sep = '/'))
  
  p <- pheatmap::pheatmap(for_plot, color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
                          scale = "column")
  ggsave(paste(output.dir3, paste0('/CN_',CN.num,'/4.3_CN_Heatmap','_',CN.num,'/','CN_heatmap_',CN.num,'_colscale.pdf'), sep = '/'), p, height = 5, width = 2.5 + ncol(for_plot)/10)
  ggsave(paste(output.dir3, paste0('/CN_',CN.num,'/4.3_CN_Heatmap','_',CN.num,'/','CN_heatmap_',CN.num,'_colscale.png'), sep = '/'), p, height = 5, width = 2.5 + ncol(for_plot)/10,dpi=300)
  
}







## 3. cell ratio stat ######################################################################

ttest.data <- colData(spe) %>% as_tibble() %>% group_by(roi_id) %>% dplyr::select(celltype)%>% 
  table() %>% normalize()%>% as.data.frame() %>% left_join(clinical)
table(ttest.data$celltype, ttest.data$group4.2)



### 2.2 t test #####
g = group.num
for (g in group.num) {
  
  stat.dir <- paste0(output.dir2, '/3.3_celltype.ratio2.', g)
  if(! dir.exists(stat.dir)) {dir.create(stat.dir)}
  
  ttest.data.sub <- ttest.data %>% select(roi_id:Freq,g)
  colnames(ttest.data.sub)[4] <- 'group'
  ttest.data.sub <- ttest.data.sub %>% filter(! group == '')
  
  ## t test
  
  for (i in unique(ttest.data.sub$celltype)) {
    print(i)
    ttest.data_i <- ttest.data.sub %>% filter(celltype == i)
    
    y_position <- max(ttest.data_i$Freq) * .8  * seq(2,0,by=-0.15)[1:length(group.stat.list[g][[1]][['stat.label']])]
    
    ttest.data_i$group <- factor(ttest.data_i$group, levels = group.label.list[[1]])
    
    p <- ggplot(ttest.data_i, aes(x = group, y = Freq)) +
      geom_boxplot(aes(color = group), outlier.alpha = 0, alpha = 0,notch = F) +
      geom_signif(comparisons = group.stat.list[g][[1]][['stat.label']], test = 't.test',#vjust = 1.5,
                  textsize = 3, y_position = y_position) +
      theme_classic() +
      labs(y = 'Ratio', x = '', title = i) +
      guides(color=guide_legend(ncol=1,byrow=F)) +
      scale_color_tableau('Tableau 20') +
      scale_fill_tableau('Tableau 20') +
      theme(strip.background = element_blank(),
            strip.text = element_text(angle = 0, size = 10),
            legend.position = 'none',
            axis.text.x = element_text(angle = 60, hjust = 1),
            plot.title = element_text(size = 12),
            axis.text = element_text(size = 10))
    temp.i <- str_replace_all(i,':','_')
    temp.i <- str_replace_all(temp.i,'/','_')
    ggsave(paste(stat.dir,  paste0('ratio_cluster_',temp.i,'.pdf'), sep = '/'), p, height = 4, width=3,limitsize = FALSE)
    ggsave(paste(stat.dir,  paste0('ratio_cluster_',temp.i,'.png'), sep = '/'), p, height = 4, width=3,dpi=300,limitsize = FALSE)
    
  }
  
  
  
}




