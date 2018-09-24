#2016-7-28 for pcluster
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggdendro))
suppressMessages(library(ggplot2))
suppressMessages(library(argparser))
options(stringsAsFactors = F)


p <- arg_parser("perform cluster analysis")
p <- add_argument(p, '--file_path', help = 'expression matrix file.')
p <- add_argument(p, '--name', help = 'output file name.')
p <- add_argument(p, '--method', help = 'method for dist/cor.')
p <- add_argument(p, '--sample_inf', help = 'sample information with sample names and group names')
p <- add_argument(p, '--out_dir', help = 'output directory')
p <- add_argument(p, '--genes', help = 'user provided gene list', default = NULL)
argv <- parse_args(p)


file_path <- argv$file_path
name <- argv$name
method <- argv$method
sperm <- argv$sample_inf
out_dir <- argv$out_dir
genes <- argv$genes
#color_theme <- args[6]
#samples <- args[6]
#genes <- args[7]

data <- read.table(file_path,header = T,row.names = 1,com = "", check.names = F)

if (! is.na(genes)) {
    gene_df <- read.table(genes, header=F)
    data <- data[gene_df$V1,]
    print(dim(data))
}

name_if <- read.table(sperm, header = F, row.names = 2)
data.mat <- as.matrix(data)
data.mat.log <- log2(data.mat+1)
sample.name <- names(data)

names <- name_if[sample.name, 'V1']
col_num = length(unique(names))

annotation_col = data.frame(
  names = sample.name,
  tissues = factor(names)
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
  group_colors <- colorRampPalette(brewer.pal( 9 , "Set1" ) )(col_num)
  names(group_colors) <- unique(names)
  p <- ggplot() 
  p <- p + geom_segment(data=segment(dendr), aes(x=x, y=y+add_y, xend=xend, yend=yend+add_y),size = 0.2)
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

