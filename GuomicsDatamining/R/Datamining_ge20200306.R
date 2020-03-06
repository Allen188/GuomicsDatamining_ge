#' ge.na.ratio
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
ge.na.ratio <- function(x){
  sum(is.na(x))/dim(x)[1]/dim(x)[2]
}

#' ge.split
#'
#' @param data
#' @param split
#' @param which
#'
#' @return
#' @export
#'
#' @examples
ge.split <- function(data,split,which=1){
  sapply(data,function(v){strsplit(v,split)[[1]][which]})
}

#' ge.na.ratio
#'
#' @param data
#' @param sep
#' @param header
#'
#' @return
#' @export
#'
#' @examples
ge.readtable <- function(data,sep = "\t",header = T){
  read.table(data,sep = sep,header = header,stringsAsFactors = F)
}

#' ge.writetable
#'
#' @param data
#' @param filename
#' @param sep
#' @param col.names
#' @param row.names
#' @param quote
#'
#' @return
#' @export
#'
#' @examples
ge.writetable <- function(data,filename ,sep = "\t",col.names = T,row.names = T,quote = F){
  write.table(data,filename,sep=sep,col.names = col.names,row.names = row.names,quote = quote)
}

#' ge.remove.techrep
#'
#' @param data
#' @param pattern
#' @param method
#'
#' @return
#' @export
#'
#' @examples
ge.remove.techrep <- function(data,pattern="_repB",method="mean"){
  repB <- names(data)[grepl(pattern, names(data))]
  for (i in repB) {
    repa <- str_split(i,pattern)[[1]][1]
    df1 <- data[grepl(repa, names(data)),]
    data <- data[!grepl(repa, names(data)),]
    new_mean <- apply(df1, 1, function(x){ifelse(sum(is.na(x))==2,NA, mean(as.numeric(x),na.rm=T))} )
    data <- cbind(data,new_mean)
    names(data)[ncol(data)] <- repa
  }
  return(data)
}

#' ge.plot.techrep.correlation
#'
#' @param cor1
#' @param cor2
#' @param name
#'
#' @return
#' @export
#'
#' @examples
ge.plot.techrep.correlation <- function(cor1,cor2,name="pearson_correlation"){
  pdf(paste0(name,".pdf"))
  r <- cor(cor1, cor2, use = "pairwise.complete.obs")
  smoothScatter(cor1, cor2, nrpoints = 100,cex = 2,
                colramp = colorRampPalette(c(blues9,"orange", "red")),
                main = name, xlab = "repA", ylab = "repB")
  abline(lm(cor1 ~ cor2), col="red", lwd=2, lty=2)
  text(min(cor1,na.rm = T)*1.2,max(cor2,na.rm = T)*0.9,labels =paste0( "r =", as.character(round(r,4))),cex = 1.2)
  dev.off()
}

#' ge.plot.pool.correlation
#'
#' @param data
#' @param name
#' @param method
#'
#' @return
#' @export
#'
#' @examples
ge.plot.pool.correlation <- function(data,name="bio_cor",method="circle"){
  df_cor <- data.frame(data)
  pdf(paste0(name,".pdf"))
  mycor=cor(df_cor, use = "pairwise.complete.obs")
  corrplot(mycor, method=method,type = "upper",tl.col = "black",tl.srt = 45, tl.cex = 0.5)
  dev.off()
}

#' ge.plot.boxplot
#'
#' @param data
#' @param x
#' @param y
#' @param type
#' @param filename
#'
#' @return
#' @export
#'
#' @examples
ge.plot.boxplot <- function(data,x,y,type,filename){
  a <- ggplot(data=data, aes(x =x, y =y ,color=type,group=type)) +
    geom_jitter(alpha = 0.3,size=3) +
    geom_boxplot(alpha = .5,size=1)+
    labs(x="sample",y="value",fill= "type")+
    theme_bw() +
    theme(panel.border = element_blank())+
    theme(axis.line = element_line(size=1, colour = "black")) +
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,colour = "black"),text = element_text(size = 15,colour = "black"))+
    theme(axis.text.x = element_text( hjust = 1,angle = 45))
  ggsave(paste0(filename, ".pdf"),plot=a,width=8,height=8)
}

#' ge.plot.valcano
#'
#' @param data
#' @param title
#' @param fd
#' @param pvalue
#'
#' @return
#' @export
#'
#' @examples
ge.plot.valcano <- function(data, title,fd=1,pvalue=0.05){
  df8 <- data
  pdf(paste0(title, "_volcano.pdf"))
  plot(df8$fd, -log10(df8$P_value_adjust), col="#00000033", pch=19,
       xlab=paste("log2 (fold change)"),
       ylab="-log10 (P_value_adjust)",
       main= title)

  up <- subset(df8, df8$P_value_adjust < pvalue & df8$fd > fd)
  down <- subset(df8, df8$P_value_adjust< pvalue & df8$fd< -1*fd)
  write.csv(up,file = paste0(title, "_up.csv"))
  write.csv(down,file = paste0(title, "_down.csv"))
  points(up$fd, -log10(up$P_value_adjust), col=1, bg = brewer.pal(9, "YlOrRd")[6], pch=21, cex=1.5)
  points(down$fd, -log10(down$P_value_adjust), col = 1, bg = brewer.pal(11,"RdBu")[9], pch = 21,cex=1.5)
  abline(h=-log10(pvalue),v=c(-1*fd,fd),lty=2,lwd=1)

  dev.off()
}

ge.plot.pca <- function(data,type,title=""){
  library(ggbiplot)
  df10 <- data
  df10[is.na(df10)] <- 0
  names <-type
  df10 <- t(apply(df10, 1, scale))
  colnames(df10) <- names
  df.pr <- prcomp(t(df10))
  a<- ggbiplot(df.pr, obs.scale = 1, var.scale = 10, groups =names,alpha = 0,varname.size= 1, ellipse =F, circle = F,var.axes = F)+
    geom_point(aes(colour=names),size = 3,alpha=0.8)+
    # geom_point(aes(shape=df1$column),size = 3,alpha=1/2)+
    scale_color_manual(name="type",values=c("#537e35","#e17832","#f5b819","#5992c6","#282f89"))+
    theme(legend.direction = 'horizontal',legend.position = 'top',legend.text = element_text(size = 15,color = "black"), legend.title = element_text(size=15,color="black") ,panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+ theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(plot.subtitle=element_text(size=30, hjust=0, color="black"))+
    theme(axis.title.x=element_text(size=17, hjust=0.5, color="black"))+
    theme(axis.title.y=element_text(size=17, hjust=0.5, color="black"))
  ggsave(paste0(title,"_pca.pdf"),plot =a ,device = NULL)
}


ge.plot.tsne <- function(data,type,title=""){
  library(Rtsne)
  df10 <- data
  df10[is.na(df10)] <- 0
  names <-type
  df10 <- t(apply(df10, 1, scale))
  colnames(df10) <- names

  color <- factor(names) #,levels = c("red","#74A9CF")
  pdf(paste0(title,"_TNSE.pdf"))
  df11.tsne <- Rtsne(t(df10), dims = 2, perplexity = (ncol(data)-1)/3-1, verbose = T , check_duplicates = FALSE)
  plot(df11.tsne$Y,col=color, main = "tsne", pch = 20,cex=2,cex.axis=2,cex.lab=2)

  plot(df11.tsne$Y, type = "n", main = "tsne", pch = 20)
  text(df11.tsne$Y, labels = names(data), col= "DimGrey")
  dev.off()
}

ge.plot.umap<- function(data,type,title=""){
  library(umap)
  df10 <- data
  df10[is.na(df10)] <- 0
  names <-type
  df10 <- t(apply(df10, 1, scale))
  colnames(df10) <- names

  color <- factor(names) #,levels = c("red","#74A9CF")
  pdf(paste0(title,"_UMAP.pdf"))
  df.umap <- umap(t(df10))
  plot(df.umap$layout,col = color, main = "umap", pch = 20,cex=2,cex.axis=2,cex.lab=2)
  dev.off()
}


ge.mfuzz.cselection <- function(data,range=seq(5,50,5),repeats = 5){
  df3a<-as.matrix(data)
  df3Ex<- ExpressionSet(assayData = df3a)
  if(interactive()){
    df3F <- filter.NA(df3Ex)
    df3F <- fill.NA(df3F)
    df3F <- standardise(df3F)
  }

  df3F <- filter.NA(df3F)
  m<-mestimate(df3F)
  cselection(df3F,m=m,crange = range,repeats = repeats,visu = T)
  return(df3F)
}

ge.mfuzz.getresult <- function(data, pic,filename){
  cl <- mfuzz(data,c=pic,m=1.25)
  dir.create(path=filename,recursive = TRUE)
  pdf(paste0(filename,".pdf"))
  mfuzz.plot2(data, cl=cl,mfrow=c(4,4),centre=TRUE,x11=F,centre.lwd=0.2)#min.mem=0.99
  dev.off()

  for(i in 1:pic){
    potname<-names(cl$cluster[unname(cl$cluster)==i])
    write.csv(cl[[4]][potname,i],paste0(filename,"/mfuzz_",i,".csv"))
  }
}




