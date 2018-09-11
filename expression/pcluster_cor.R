#2016-7-28 for pcluster
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggdendro))
suppressMessages(library(ggplot2))
suppressMessages(library(argparser))
options(stringsAsFactors = F)
#read data

args<-commandArgs(T)
file_path <- args[1]
name <- args[2]
method <- args[3]
sperm <- args[4]
out_dir <- args[5]
#color_theme <- args[6]
#samples <- args[6]
#genes <- args[7]

# for test
#file_path <- './test.matrix.txt'
#name <- 'test'
#method <- 'pearson'
#sperm <- 'new_name_map_germ_sample'
#out_dir <- './'


data <- read.table(file_path,header = T,row.names = 1,com = "", check.names = F)
name_if <- read.table(sperm, header = F, row.names = 2)
#sample_df <- read.delim(samples, header = F)
#gene_df <- read.delim(genes, header = F) 
#data <- data[gene_df$V1,]
#data <- data[rowSums(data) > 0, ]

#test 50 samples
#data.small <- data %>% sample_n(48)
data.mat <- as.matrix(data)
data.mat.log <- log2(data.mat+1)
sample.name <- names(data)

names <- name_if[sample.name, 'V1']
#names2 <- name_if[sample.name, 'sample_c']
#names3 <- name_if[sample.name, 'sperm']

col_num = length(unique(names))

annotation_col = data.frame(
  names = sample.name,
  tissues = factor(names)
#  samples = factor(names2),
#  sperms = factor(names3)
)


#base cluster
pcluster <- function(data,types){
  
  pcluster_theme <- theme_bw()+ theme(axis.line.y=element_blank(),
                                      axis.ticks.y=element_blank(),
                                      axis.text.y=element_blank(),
                                      axis.ticks.x=element_blank(),
                                      axis.text.x=element_blank(),                                      
                                      panel.background=element_rect(fill="white"),
                                      panel.grid=element_blank())

    if (method == 'pearson') {
        sample_cor = cor(data, method=method, use='pairwise.complete.obs')
        ##data.mat.log.scale <- scale(data)
        sample_dist = as.dist(1-sample_cor)
        sample_dist = sample_dist*100
    } else {
        data = t(scale(t(data), scale=F))
        sample_dist = dist(t(data), method=method)
    }
  fit_hc <- hclust(sample_dist)
  dendr <- dendro_data(fit_hc,type="rectangle") 
  text.df <- merge(label(dendr),annotation_col,by.x="label",by.y="names")
  
  max_y <- max(segment(dendr)$yend)
  add_y <- ceiling(max_y / 15)
  max_y <- ceiling(max_y + max_y/15)
#  col1 <- brewer.pal(9,"Set1")
#  col2 <- brewer.pal(8,"Set2")
#  col3 <- brewer.pal(9,"Set3")
#  ann_colors <- c(col1, col3)[1:col_num]
                                        #  names(ann_colors) <- unique(names)
  ## color_cfg <- read.table(color_theme, header = F,
    ##                         comment.char = "", stringsAsFactors = F)
  group_colors <- colorRampPalette(brewer.pal( 9 , "Set1" ) )(col_num)
  #group_colors <- color_cfg$V2
  names(group_colors) <- unique(names)
  p <- ggplot() 
  p <- p + geom_segment(data=segment(dendr), aes(x=x, y=y+add_y, xend=xend, yend=yend+add_y),size = 0.2)
#  p <- p + geom_bar(data = text.df,stat = "identity",aes(x = x,y = add_y-2,fill = tissues),width = 1)
#  p <- p + geom_bar(data = text.df,stat = "identity",aes(x = x,y = -(add_y-2)-2,fill = samples),width = 1)
  p <- p + geom_bar(data = text.df,stat = "identity",aes(x = x, y = add_y, fill = tissues),width = 1) + scale_fill_manual(values = group_colors)
  p <- p + geom_text(data = text.df,aes(x=x,y = -add_y/5,label=label),angle = 90,vjust = 0,hjust = 1,size = 1)
  p <- p + labs(fill = "")+ylab("") + xlab('')
  p <- p + ylim(c(-add_y-10,max_y + 10)) + pcluster_theme
  out_name = paste(types,method,'pcluster',sep = ".")
  ggsave(paste(out_dir,'/',out_name,'.pdf',sep = ""),plot = p,width = 12,height = 6)
  ggsave(paste(out_dir,'/',out_name,'.png',sep = ""),type = 'cairo-png',plot = p,width = 12,height = 6)
}
#---
pcluster(data.mat.log,name)

