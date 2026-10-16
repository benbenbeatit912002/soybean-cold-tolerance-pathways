#先祖何在設計

install.packages('igraph')
library(igraph)

datdir='C:/Users/chris/Desktop/程式語言/研究法/大豆基因/network/crosstalk/0221版本'
datfile = paste(datdir, 'jc_point_short.txt', sep='/')
jc_point_short= read.table(datfile,header = F,sep='\t')

datdir='C:/Users/chris/Desktop/程式語言/研究法/大豆基因/network/crosstalk/0221版本'
datfile = paste(datdir, 'jc_point_middle.txt', sep='/')
jc_point_middle= read.table(datfile,header = F,sep='\t')

datfile = paste(datdir, 'both.txt', sep='/')
both=read.table(datfile,header = F,sep='\t')
both=as.character(as.matrix(both))

datfile = paste(datdir, 'only_short.txt', sep='/')
only_short= read.table(datfile,header = F,sep='\t')
datfile = paste(datdir, 'only_middle.txt', sep='/')
only_middle= read.table(datfile,header = F,sep='\t')
#################
links=as.matrix(jc_point_short)
ct_gene=only_short[,4]
breaks=seq(0, 15, by = 3) ##將0~18 用3分割
ct_gene_size=(cut(ct_gene, breaks, labels = FALSE))*7 ##數據分割


node_names = as.character(links[1, -1])
links = links[-1, -1]
loc=which(links<0.004456328)   #(0.0044563280,0.5)
links[loc]=0

net =graph.adjacency(links, mode = "undirected", weighted = TRUE)
n=nrow(only_short)
net = simplify(net)
V(net)$name = node_names


{
  n=nrow(only_short)
  k=0
  score_width=numeric()  
  for (i in 1:n) {
    for (j in 1:n){
      
      score_width[k]=links[i,j]
      k=k+1
    }
  }
  
  
  
  score_width=unique(score_width)
  loc=which(score_width==1 | score_width==0 | is.na(score_width) )
  score_width=score_width[-loc]
  
  
  score_width = as.numeric(score_width)
  breaks_width=seq(0, 1, by = 0.1) ##將0~1 用0.1分割
  score_width=cut(score_width, breaks_width, labels = FALSE) ##數據分割
} #width

#####
links=as.matrix(jc_point_middle)
ct_gene=only_middle[,4]
breaks=seq(0, 120, by = 24) ##將0~120 用20分割
ct_gene_size2=(cut(ct_gene, breaks, labels = FALSE))*7 ##數據分割


node_names = as.character(links[1, -1])
links = links[-1, -1]
loc=which(links<0.14160401)  #(0.14160401,0.5)
links[loc]=0


net2 =graph.adjacency(links, mode = "undirected", weighted = TRUE)
n=nrow(only_middle)

net2 = simplify(net2)
V(net2)$name = node_names

{
  k=0
  score_width2=numeric()  
  for (i in 1:n) {
    for (j in 1:n){
      
      score_width2[k]=links[i,j]
      k=k+1
    }
  } 
  score_width2=unique(score_width2)
  loc=which(score_width2==1 | score_width2==0 | is.na(score_width2) )
  score_width2=score_width2[-loc]
  
  score_width2 = as.numeric(score_width2)
  
  breaks_width=seq(0, 1, by = 0.1) ##將0~1 用0.1分割
  score_width2=cut(score_width2, breaks_width, labels = FALSE) ##數據分割
}

######################################
netcb=net+net2
#####  # 恢復原始結構
{
  
  loc_netcb=as.data.frame(get.edgelist(netcb))
  loc_netcb$combined <- paste(loc_netcb$V1, loc_netcb$V2, sep = " ")
  loc_net=as.data.frame(get.edgelist(net))
  loc_net$combined <- paste(loc_net$V1, loc_net$V2, sep = " ")
  loc_net2=as.data.frame(get.edgelist(net2))
  loc_net2$combined <- paste(loc_net2$V1, loc_net2$V2, sep = " ")
  
  
  
  
  k = integer(nrow(loc_netcb))
  for (i in 1:nrow(loc_netcb)) {
    row_match = (loc_net[, 1] == loc_netcb[i, 1] & loc_net[, 2] == loc_netcb[i, 2]) |
      (loc_net[, 1] == loc_netcb[i, 2] & loc_net[, 2] == loc_netcb[i, 1])
    match_index = which(row_match)
    if (length(match_index) > 0) {
      k[i] = match_index[1]  # 只取第一个匹配的位置
    }
  }
  
  k2 = integer(nrow(loc_netcb))
  for (i in 1:nrow(loc_netcb)) {
    row_match = (loc_net2[, 1] == loc_netcb[i, 1] & loc_net2[, 2] == loc_netcb[i, 2]) |
      (loc_net2[, 1] == loc_netcb[i, 2] & loc_net2[, 2] == loc_netcb[i, 1])
    match_index = which(row_match)
    if (length(match_index) > 0) {
      k2[i] = match_index[1]  # 只取第一个匹配的位置
    }
  }
  loc1=which(k!=0) 
  loc2=which(k2!=0) 
  loc3=47                                   ######自己找兩個金色交###### get.edgelist(netcb)
  length(loc2)
  lo=which(k==0)
  lo2=which(k2==0)
  k=k[-lo]  #
  k2=k2[-lo2]#
} #處理  
E(netcb)[loc1]$width =score_width[k]*3
E(netcb)[loc2]$width =score_width2[k2]*3

V(netcb)$color=ifelse(V(netcb)$name %in% both,'gold2',
                      ifelse(V(netcb)$name%in%only_short[c(13,9,5,11),2],'lightskyblue' ,#自己找p值小的，從短期p-value的excel中找 淺色代表p值小的
                             ifelse(V(netcb)$name%in% only_middle[,2],'#FBAED2','lightskyblue3')))
#F64A8A
V(netcb)$label.cex= 2
E(netcb)$color[loc1] = ifelse(E(netcb)$width[loc1] <=6, 'gray80',
                              ifelse(E(netcb)$width[loc1] <= 12, 'gray70',
                                     ifelse(E(netcb)$width[loc1] <= 30, 'gray65')))
E(netcb)$color[loc2] = ifelse(E(netcb)$width[loc2] <=6, 'wheat1',
                              ifelse(E(netcb)$width[loc2] <= 12, 'bisque2',
                                     ifelse(E(netcb)$width[loc2] <= 30, 'burlywood')))
E(netcb)$color[loc3]='darkseagreen3'
ct_gene_size2=ct_gene_size2[-c(9,10)] #刪除重複的點 從中期的點中找中間兩點V(net2)
V(netcb)$size = c(ct_gene_size,ct_gene_size2)
plot_tk = tkplot(netcb, width = 1000, height = 1000, tkplot.bg = "transparent",vertex.frame.color = NA)