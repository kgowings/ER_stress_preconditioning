#Open GSEAplot.txt in R and run the following code: 
#Example
  #Term	pvalue	ES	No_genes
  #male meiosis I	0.015876766	0.518760295	8
  #negative regulation of growth	0.01509331	0.521337929	6



require(ggplot2)
GSEAplot <- read.delim("GSEAplot.txt")
data <- GSEAplot
data$Term<-as.character(data$Term)
data$Term <- factor(data$Term, levels=unique(data$Term))
ggplot(data=data,aes(x=ES,y=Term, size=No_genes))+
  geom_point(aes(color=pvalue, size=No_genes))+
  scale_size_continuous(range = c(2.5,7.25))+
  scale_color_gradient(low="red", high="blue")+
  xlab("Enrichment Score") + 
  ylab("Enriched Ontological Terms") 
  #+ ggtitle("GSEA for Stress Prconditioning GWAS")


RStudio 2022.02.2+485 "Prairie Trillium" Release (8acbd38b0d4ca3c86c570cf4112a8180c48cc6fb, 2022-04-19) for Windows
Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) QtWebEngine/5.12.8 Chrome/69.0.3497.128 Safari/537.36

ggplot2 package version 3.3.6
