####load packages
library(readxl)
library(pheatmap)
library(ggplot2)
library(fpc)
library(ggfortify)
library(reshape2)
library(ggrepel)
library(VennDiagram)
library(tidyverse)
library(ComplexHeatmap)
library(pROC)

###fig2b
fig2b<-read_excel("Source Data Fig. 2.xlsx", sheet = "fig2b")
fig2b<-data.frame(fig2b)
rownames(fig2b)<-fig2b$Uniprot
fig2b<-fig2b[,-1]
fig2b<-t(fig2b)

pheatmap(fig2b,show_rownames=FALSE,show_colnames=FALSE,fontsize=12,scale="none",clustering_method="median",
         cluster_cols=F,cluster_rows=FALSE,
         color = colorRampPalette(colors=c("#67001f","#b5172a","#d7604d","#f5a581","#fddbc6","#cce6f1","#90c5dd","#4793c2","#1f66ac","#073060"))(100)
)

###fig2c
fig2c<-read_excel("Source Data Fig. 2.xlsx", sheet = "fig2c")
data <- fig2c[,2:4]
nk=2:10
set.seed(22)
wss<-sapply(nk,function(k){kmeans(data,centers = k)$tot.withinss})
wss
plot(nk,wss,type = "l",xlab="Number of k",ylab="Within sum of squares")
sw=sapply(nk,function(k){cluster.stats(dist(data),kmeans(data,centers = k)$cluster)$avg.silwidth})
sw
plot(nk,sw,type = "l",xlab="Number of clusters",ylab = "Average silhouette width")
fit<-kmeans(data,4)
fit
cluster_output <- data.frame(mydata,fit$cluster)

set.seed(22)
autoplot(kmeans(data, 4),data=data,size=3,label=FALSE,color="k",frame=TRUE)+theme_bw()+
  theme (panel.grid.major=element_line(colour=NA),
         panel.background = element_rect(fill = "transparent",colour = NA),
         panel.grid.minor = element_blank())+
  scale_color_manual(values=c("#7c4caf","#1576ae","#fd9200","#d33524"))+
  scale_fill_manual(values=c("white","white","white","white"))

###fig2d
fig2d<-read_excel("Source Data Fig. 2.xlsx", sheet = "fig2d")
lineplotdata<-melt(fig2d,id.vars=c("Uniprot","fit.cluster"),measures.vars=c("x1","x2","x3"),variable.name="X",value.name="value" )
lineplotdata$timepoint<-unclass(lineplotdata$X)

lineplotdata1<-subset(lineplotdata,fit.cluster==1)
lineplotdata2<-subset(lineplotdata,fit.cluster==2)
lineplotdata3<-subset(lineplotdata,fit.cluster==3)
lineplotdata4<-subset(lineplotdata,fit.cluster==4)

ggplot(data = lineplotdata1,aes(x=timepoint,y=value,group = Uniprot))+
  geom_line(color="#d33524")+
  xlab("Timepoint")+
  ylab("Protein levels (z-scores)")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(.075,.915),
        axis.text.x = element_text(size = 12,color="black"),
        axis.text.y = element_text(size = 12,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        legend.box.background = element_rect(color="black"))+
  scale_x_continuous(limits=c(1,3),breaks = seq(1,3,1))+
  scale_y_continuous(limits = c(-0.51,0.51))

ggplot(data = lineplotdata2,aes(x=timepoint,y=value,group = Uniprot))+
  geom_line(color="#fd9200")+
  xlab("Timepoint")+
  ylab("Protein levels (z-scores)")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(.075,.915),
        axis.text.x = element_text(size = 12,color="black"),
        axis.text.y = element_text(size = 12,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        legend.box.background = element_rect(color="black"))+
  scale_x_continuous(limits=c(1,3),breaks = seq(1,3,1))+
  scale_y_continuous(limits = c(-0.51,0.51))

ggplot(data = lineplotdata3,aes(x=timepoint,y=value,group = Uniprot))+
  geom_line(color="#7c4caf")+
  xlab("Timepoint")+
  ylab("Protein levels (z-scores)")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(.075,.915),
        axis.text.x = element_text(size = 12,color="black"),
        axis.text.y = element_text(size = 12,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        legend.box.background = element_rect(color="black"))+
  scale_x_continuous(limits=c(1,3),breaks = seq(1,3,1))+
  scale_y_continuous(limits = c(-0.51,0.51))


ggplot(data = lineplotdata4,aes(x=timepoint,y=value,group = Uniprot))+
  geom_line(color="#1576ae")+
  xlab("Timepoint")+
  ylab("Protein levels (z-scores)")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(.075,.915),
        axis.text.x = element_text(size = 12,color="black"),
        axis.text.y = element_text(size = 12,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        legend.box.background = element_rect(color="black"))+
  scale_x_continuous(limits=c(1,3),breaks = seq(1,3,1))+
  scale_y_continuous(limits = c(-0.51,0.51))

###fig2e
fig2e<-read_excel("Source Data Fig. 2.xlsx", sheet = "fig2e")
ggplot(fig2e,aes(x=coeftime,y=logqtime))+geom_point(pch=21,color="black",aes(size=abs(coeftime),bg=group),show.legend = FALSE)+
  scale_fill_manual(values=c("#d33524","#fd9200","#7c4caf","#1576ae","darkgrey"))+ 
  geom_text_repel(aes(label=label),size=3,direction="both",force)+
  scale_x_continuous(limits=c(-0.3,0.3),breaks=c(-0.30,-0.20,-0.10,0,0.10,0.20,0.30),labels=c(expression("-0.30"),expression("-0.20"),expression("-0.10"),expression("0"),expression("0.10"),expression("0.20"),expression("0.30")))+
  scale_y_continuous(limits=c(0,31),breaks=c(0,1.301029995664,10,20,30),labels=c(expression("0"),expression("0.05"),expression("1x10"^"-10"),expression("1x10"^"-20"),expression("1x10"^"-30")))+
  theme_bw() +
  theme(panel.grid=element_blank(),legend.box.background = element_rect(color="black"),panel.background = element_blank())+#axis.line = element_line(size=0.1))+
  theme(axis.text=element_text(size=6))+
  geom_hline(aes(yintercept=1.30103),color="darkgrey",linetype="dashed",size=0.5)+
  xlab("Effect size (time)")+ylab("Q-value")+theme(axis.title=element_text(size=6)) 

###fig2f
fig2f<-read_excel("Source Data Fig. 2.xlsx", sheet = "fig2f")

ggplot(fig2f,aes(x=logq,y=order,fill=cluster_factor))+
  geom_bar(stat="identity",width=0.7)+
  scale_fill_manual(values=c("#d33524","#fd9200","#7c4caf","#1576ae"),labels = c("1","2","3","4"))+
  theme_bw()+
  theme(panel.grid=element_blank(),panel.background = element_blank(),axis.line = element_line(color="black"),axis.text = element_text(color="black"))+
  scale_x_continuous(limits=c(0,50),breaks=seq(0,50,10))+
  scale_y_discrete(labels=fig2f$name_label)+
  ylab("")+
  xlab("-log10 (Q-value)")

###fig3a
fig3a<-read_excel("Source Data Fig. 3.xlsx", sheet = "fig3a")
ggplot(fig3a,aes(x=coefage_discovery,y=-log10(q_discovery)))+
  geom_point(data=subset(fig3a,group=="null"),aes(size=abs(coefage_discovery)),pch=21,color="black",fill="darkgrey",show.legend = FALSE)+
  geom_point(data=subset(fig3a,group=="positive"),aes(size=abs(coefage_discovery)),pch=21,color="black",fill="#d33524",show.legend = FALSE)+
  geom_point(data=subset(fig3a,group=="negative"),aes(size=abs(coefage_discovery)),pch=21,color="black",fill="#1576ae",show.legend = FALSE)+
  geom_text_repel(aes(label=label),size=4,direction="both",force)+
  scale_x_continuous(limits=c(-0.034,0.033),breaks=c(-0.03,-0.02,-0.01,0,0.01,0.02,0.03),labels=c(expression("-0.03"),expression("-0.02"),expression("-0.01"),expression("0"),expression("0.01"),expression("0.02"),expression("0.03")))+
  scale_y_continuous(limits=c(0,45),breaks=c(0,10,20,30,40),labels=c(expression("1"),expression("1x10"^"-10"),expression("1x10"^"-20"),expression("1x10"^"-30"),expression("1x10"^"-40")))+
  theme_bw() +
  theme(panel.grid=element_blank(),legend.box.background = element_rect(color="black"),panel.background = element_blank())+#axis.line = element_line(size=0.1))+   ##去除网格线##去除外层边框##去除背景色##设置刻度线格式
  geom_hline(aes(yintercept=1.30103),color="darkgrey",linetype="dashed",size=0.5)+
  xlab("Effect size (age)")+ylab("Q-value (age)")  ##设置坐标轴名称与字体

###fig3b
fig3b<-read_excel("Source Data Fig. 3.xlsx", sheet = "fig3b")
ggplot(fig3b,aes(x=coefage_validation,y=-log10(q_validation)))+
  geom_point(data=subset(fig3b,group=="null"),aes(size=abs(coefage_validation)),pch=21,color="black",fill="darkgrey",show.legend = FALSE)+
  geom_point(data=subset(fig3b,group=="positive"),aes(size=abs(coefage_validation)),pch=21,color="black",fill="#d33524",show.legend = FALSE)+
  geom_point(data=subset(fig3b,group=="negative"),aes(size=abs(coefage_validation)),pch=21,color="black",fill="#1576ae",show.legend = FALSE)+
  geom_text_repel(aes(label=label),direction="both",force)+
  scale_x_continuous(limits=c(-0.034,0.033),breaks=c(-0.03,-0.02,-0.01,0,0.01,0.02,0.03),labels=c(expression("-0.03"),expression("-0.02"),expression("-0.01"),expression("0"),expression("0.01"),expression("0.02"),expression("0.03")))+
  scale_y_continuous(limits=c(0,40),breaks=c(0,10,20,30,40),labels=c(expression("1"),expression("1x10"^"-10"),expression("1x10"^"-20"),expression("1x10"^"-30"),expression("1x10"^"-40")))+
  theme_bw() +
  theme(panel.grid=element_blank(),legend.box.background = element_rect(color="black"),panel.background = element_blank())+#axis.line = element_line(size=0.1))+   ##去除网格线##去除外层边框##去除背景色##设置刻度线格式
  geom_hline(aes(yintercept=1.30103),color="darkgrey",linetype="dashed",size=0.5)+
  xlab("Effect size (age)")+ylab("Q-value (age)")  ##设置坐标轴名称与字体

###code for fig 3c###
fig3c<-read_excel("Source Data Fig. 3.xlsx", sheet = "fig3c")
ggplot(fig3c,aes(x=coefage_discovery,y=coefage_validation))+
  geom_point(data=subset(fig3c,group=="null"),size=3,pch=21,color="black",fill="darkgrey")+
  geom_point(data=subset(fig3c,group=="positive"),size=3,pch=21,color="black",fill="#d33524")+
  geom_point(data=subset(fig3c,group=="negative"),size=3,pch=21,color="black",fill="#1576ae")+
  geom_text_repel(aes(label=label),size=4,direction="both",force)+
  scale_x_continuous(limits=c(-0.035,0.035),breaks=seq(-0.03,0.03,0.01))+
  scale_y_continuous(limits=c(-0.035,0.035),breaks=seq(-0.03,0.03,0.01))+
  theme_bw() +
  theme(panel.grid=element_blank(),legend.box.background = element_rect(color="black"),panel.background = element_blank())+#axis.line = element_line(size=0.1))+   ##去除网格线##去除外层边框##去除背景色##设置刻度线格式
  geom_hline(aes(yintercept=0),color="darkgrey",linetype="dashed",size=0.5)+geom_vline(aes(xintercept=0),color="darkgrey",linetype="dashed",size=0.5)+ ##设置参考线
  xlab("Effect size (Primary analysis)")+ylab("Effect size (Validation)")  ##设置坐标轴名称与字体

###code for fig 3d###
fig3d<-read_excel("Source Data Fig. 3.xlsx", sheet = "fig3d")
ggplot(fig3d,aes(x=coefsex_discovery,y=-log10(qsex_discovery)))+
  geom_point(data=subset(fig3d,group=="null"),aes(size=abs(coefsex_discovery)),pch=21,color="black",fill="darkgrey",show.legend = FALSE)+
  geom_point(data=subset(fig3d,group=="male enriched"),aes(size=abs(coefsex_discovery)),pch=21,color="black",fill="#c47644",show.legend = FALSE)+
  geom_point(data=subset(fig3d,group=="female enriched"),aes(size=abs(coefsex_discovery)),pch=21,color="black",fill="#9b4ae0",show.legend = FALSE)+
  geom_text_repel(aes(label=label),size=4,direction="both",force)+
  scale_x_continuous(limits=c(-1.05,1.05),breaks=seq(-1,1,0.2))+
  scale_y_continuous(limits=c(0,180),breaks=c(0,30,60,90,120,150,180),labels=c(expression("1"),expression("1x10"^"-30"),expression("1x10"^"-60"),expression("1x10"^"-90"),expression("1x10"^"-120"),expression("1x10"^"-150"),expression("1x10"^"-180")))+
  theme_bw() +
  theme(panel.grid=element_blank(),legend.box.background = element_rect(color="black"),panel.background = element_blank())+#axis.line = element_line(size=0.1))+   ##去除网格线##去除外层边框##去除背景色##设置刻度线格式
  geom_hline(aes(yintercept=1.30103),color="darkgrey",linetype="dashed",size=0.5)+
  xlab("Effect size (sex)")+ylab("Q-value (sex)") ##设置坐标轴名称与字体

###code for fig 3e###
fig3e<-read_excel("Source Data Fig. 3.xlsx", sheet = "fig3e")
ggplot(fig3e,aes(x=coefsex_discovery,y=coefsex_validation))+
  geom_point(data=subset(fig3e,group=="null"),size=3,pch=21,color="black",fill="darkgrey")+
  geom_point(data=subset(fig3e,group=="male enriched"),size=3,pch=21,color="black",fill="#c47644")+
  geom_point(data=subset(fig3e,group=="female enriched"),size=3,pch=21,color="black",fill="#9b4ae0")+
  geom_text_repel(aes(label=label),size=4,direction="both",force)+
  scale_x_continuous(limits=c(-1.05,1.05),breaks=seq(-1,1,0.2))+
  scale_y_continuous(limits=c(-1.05,1.05),breaks=seq(-1,1,0.2))+
  theme_bw() +
  theme(panel.grid=element_blank(),legend.box.background = element_rect(color="black"),panel.background = element_blank())+#axis.line = element_line(size=0.1))+   ##去除网格线##去除外层边框##去除背景色##设置刻度线格式
  geom_hline(aes(yintercept=0),color="darkgrey",linetype="dashed",size=0.5)+geom_vline(aes(xintercept=0),color="darkgrey",linetype="dashed",size=0.5)+ ##设置参考线
  xlab("Effect size (Primary analysis)")+ylab("Effect size (Validation)")  ##设置坐标轴名

###fig3f
fig3f<-read_excel("Source Data Fig. 3.xlsx", sheet = "fig3f")
fig3f<-data.frame(fig3f)
Ageproteins<-fig3f[which(fig3f$agerelated=="yes"),1]
Sexproteins<-fig3f[which(fig3f$sexrelated=="yes"),1]
p=venn.diagram(list("Age-related proteins"=Ageproteins,"Sex-related proteins"=Sexproteins),cat.cex=1,cex=1,main=NULL,col=c("#1576ae","#30BD00"),filename=NULL)
pdf("venn.pdf")
grid.draw(p)
dev.off()

###fig3g
fig3g<-read_excel("Source Data Fig. 3.xlsx", sheet = "fig3g")
fig3g %>% mutate(protein = fct_reorder(protein,order)) %>%
  ggplot(aes(x=Participants,y=protein))+
  geom_point(pch=21,color="black",aes(size=-log10(q_age),fill=coef_age))+
  scale_fill_gradient2(low="#1576ae",mid="white",high="#d33524")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  theme(axis.text.x = element_text(angle = 90)) 

###fig4b
fig4b<-read_excel("Source Data Fig. 4.xlsx", sheet = "fig4b")

ggplot(fig4b,aes(x=logq,y=order,fill=network_factor))+
  geom_bar(stat="identity",width=0.7)+
  scale_fill_manual(values=c("#e2aa2d","#0aee82","#754feb","#f428dc"),labels = c("1","2","3","4"))+ 
  scale_x_continuous(limits=c(0,20),breaks=seq(0,20,5))+
  scale_y_discrete(labels=fig4b$label)+
  theme_bw()+
  theme(panel.grid=element_blank(),panel.background = element_blank(),axis.line = element_line(color="black"),axis.text = element_text(color="black"))+
  ylab("")+
  xlab("-log10 (fdr)")


###fig4c left panel
fig4c_left<-read_excel("Source Data Fig. 4.xlsx", sheet = "fig4c_left")
data<-data.frame(fig4c_left)
rownames(data)<-data$protein
data<-data[,3:34]
data<-t(data)

bk <- c(seq(-0.3,-0.01,by=0.001),seq(0.01,0.3,by=0.001))
anno_row <- data.frame(Trait = c(rep("Anthropometric", 4), rep("Lipid metabololism", 4), rep("Glucose metabololism",3), 
                                 rep("Immunity", 3), rep("Liver and kideney index", 6), rep("MMSE", 12)))
rownames(anno_row)<-row.names(data)
anno_row1 <- data.frame(Trait = c(rep("Anthropometric", 4), rep("Lipid metabololism", 4), rep("Glucose metabololism",3), 
                                  rep("Immunity", 3), rep("Live", 3), rep("Renal", 3),rep("MMSE", 12)))
rownames(anno_row1)<-row.names(data)
anno_col<-data.frame(Protein = c(rep("Network1", 18), rep("Network2", 17), rep("Network3",12),
                                 rep("Network4", 9), rep("Others", 30)))
rownames(anno_col)<-colnames(data)

ann_col = list(
  Protein = c(network1="#d55e00",network2="#e69f00",network3="#009e73",network4="#56b4e9",network5="darkgrey")
)

pheatmap(data,clustering_method="complete",
         cluster_rows=T, 
         annotation_row = anno_row1,annotation_col = anno_col,
         row_split=anno_row1$Trait,##column_split=anno_col$Protein,
         color = c(colorRampPalette(colors = c("#1576ae","#4793c2","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","#d7604d","#d33524"))(length(bk)/2)),
         breaks=bk
)

###fig4c right panel
line_bar<-read_excel("Source Data Fig. 4.xlsx", sheet = "fig4c_right")

p <- ggplot(line_bar)  +  
  geom_bar(aes(x=as.factor(order), y=sig),stat="identity", fill="cyan",colour="#006000")+ 
  labs(x="Traits",y="r correlation")+
  ###scale_y_continuous(limits=c(-0.26,1),breaks=seq(-0.26,1,0.1))+
  theme_bw()+#去掉背景灰色
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),#以上theme中代码用于去除网格线且保留坐标轴边框
        legend.box.background = element_rect(color="black"))#为图例田间边框线
p

p1<-p+geom_line(aes(x=order, y=r*40),color="grey")+
  geom_point(aes(x=order, y=r*40),shape=21, color="black", fill="#69b3a2", size=6)+
  scale_y_continuous(name="Number of significant proteins",sec.axis = sec_axis(~./40,name = 'Pearsons r'))

p1

###fig5
fig5<-read_excel("Source Data Fig. 5.xlsx")
dys<-subset(fig5,disease=="dyslipidemia")
hyper<-subset(fig5,disease=="hypertension")
t2d<-subset(fig5,disease=="T2D")
fattyliver<-subset(fig5,disease=="fatty liver")
hepa<-subset(fig5,disease=="hepatitis")
renal<-subset(fig5,disease=="renal disease")

ggplot(dys,aes(y=reorder(Protein,loghr), x=loghr))+
  geom_vline(aes(xintercept=0),color="darkgrey",linetype="dashed",size=0.5)+
  scale_color_manual(values=c("#d33524","#d33524","#1576ae","#1576ae"),labels=c("fdr<0.05","p<0.05","p<0.05","fdr<0.05"))+ 
  scale_shape_manual(values=c(16,1,1,16),labels=c("fdr<0.05","p<0.05","p<0.05","fdr<0.05"))+
  geom_errorbar(aes(xmin=logll, xmax=logul, color=group), width = 0.2,size=0.5)+
  geom_point(size=4,aes(shape=group,color=group))+
  scale_x_continuous(limits=c(-0.30102999566398114,0.30102999566398114),breaks=c(-0.30102999566398114,0,0.30102999566398114),labels=c("0.5","1.0","2.0"))+
  xlab("HR (95% CI)")+
  ylab(NULL)+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.5),
        axis.text=element_text(size=10),
        axis.text.x = element_text(hjust=1,vjust=0.5,color="black"),
        legend.title=element_blank(),
        legend.position = "right") 

ggplot(hyper,aes(y=reorder(Protein,loghr), x=loghr))+
  geom_vline(aes(xintercept=0),color="darkgrey",linetype="dashed",size=0.5)+
  scale_color_manual(values=c("#d33524","#1576ae","#1576ae"),labels=c("p<0.05","p<0.05","fdr<0.05"))+ 
  scale_shape_manual(values=c(1,1,16),labels=c("p<0.05","p<0.05","fdr<0.05"))+
  geom_errorbar(aes(xmin=logll, xmax=logul, color=group), width = 0.2,size=0.5)+
  geom_point(size=4,aes(shape=group,color=group))+
  scale_x_continuous(limits=c(-0.30102999566398114,0.30102999566398114),breaks=c(-0.30102999566398114,0,0.30102999566398114),labels=c("0.5","1.0","2.0"))+
  xlab("HR (95% CI)")+
  ylab(NULL)+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.5),
        axis.text=element_text(size=10),
        axis.text.x = element_text(hjust=1,vjust=0.5,color="black"),
        legend.title=element_blank(),
        legend.position = "none") 


ggplot(t2d,aes(y=reorder(Protein,loghr), x=loghr))+
  geom_vline(aes(xintercept=0),color="darkgrey",linetype="dashed",size=0.5)+ ##设置参考线
  scale_color_manual(values=c("#d33524","#d33524","#1576ae","#1576ae"),labels=c("fdr<0.05","p<0.05","p<0.05","fdr<0.05"))+ 
  scale_shape_manual(values=c(16,1,1,16),labels=c("fdr<0.05","p<0.05","p<0.05","fdr<0.05"))+
  geom_errorbar(aes(xmin=logll, xmax=logul, color=group), width = 0.2,size=0.5)+
  geom_point(size=4,aes(shape=group,color=group))+
  scale_x_continuous(limits=c(-0.335,0.30102999566398114),breaks=c(-0.30102999566398114,0,0.30102999566398114),labels=c("0.5","1.0","2.0"))+
  xlab("HR (95% CI)")+
  ylab(NULL)+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.5),
        axis.text=element_text(size=10),
        axis.text.x = element_text(hjust=1,vjust=0.5,color="black"),
        legend.position = "none") 


ggplot(fattyliver,aes(y=reorder(Protein,loghr), x=loghr))+
  geom_vline(aes(xintercept=0),color="darkgrey",linetype="dashed",size=0.5)+ ##设置参考线
  scale_color_manual(values=c("#d33524","#d33524","#1576ae","#1576ae"),labels=c("fdr<0.05","p<0.05","p<0.05","fdr<0.05"))+ 
  scale_shape_manual(values=c(16,1,1,16),labels=c("fdr<0.05","p<0.05","p<0.05","fdr<0.05"))+
  geom_errorbar(aes(xmin=logll, xmax=logul, color=group), width = 0.2,size=0.5)+
  geom_point(size=4,aes(shape=group,color=group))+
  scale_x_continuous(limits=c(-0.30102999566398114,0.30102999566398114),breaks=c(-0.30102999566398114,0,0.30102999566398114),labels=c("0.5","1.0","2.0"))+
  xlab("HR (95% CI)")+
  ylab(NULL)+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.5),
        axis.text=element_text(size=10),
        axis.text.x = element_text(hjust=1,vjust=0.5,color="black"),
        legend.position = "none") ##去除网格线##去除外层边框##去除背景色##设置刻度线格式

ggplot(hepa,aes(y=reorder(Protein,loghr), x=loghr))+
  geom_vline(aes(xintercept=0),color="darkgrey",linetype="dashed",size=0.5)+ ##设置参考线
  scale_color_manual(values=c("#d33524","#d33524","#1576ae","#1576ae"),labels=c("fdr<0.05","p<0.05","p<0.05","fdr<0.05"))+ 
  scale_shape_manual(values=c(16,1,1,16),labels=c("fdr<0.05","p<0.05","p<0.05","fdr<0.05"))+
  geom_errorbar(aes(xmin=logll, xmax=logul, color=group), width = 0.2,size=0.5)+
  geom_point(size=4,aes(shape=group,color=group))+
  scale_x_continuous(limits=c(-1,1.05),breaks=c(-1,0,1),labels=c("0.1","1.0","10"))+
  xlab("HR (95% CI)")+
  ylab(NULL)+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.5),
        axis.text=element_text(size=10),
        axis.text.x = element_text(hjust=1,vjust=0.5,color="black"),
        legend.position = "none")

ggplot(renal,aes(y=reorder(Protein,loghr), x=loghr))+
  geom_vline(aes(xintercept=0),color="darkgrey",linetype="dashed",size=0.5)+
  scale_color_manual(values=c("#d33524","#1576ae","#1576ae"),labels=c("p<0.05","p<0.05","fdr<0.05"))+ 
  scale_shape_manual(values=c(1,1,16),labels=c("p<0.05","p<0.05","fdr<0.05"))+
  geom_errorbar(aes(xmin=logll, xmax=logul, color=group), width = 0.2,size=0.5)+
  geom_point(size=4,aes(shape=group,color=group))+
  scale_x_continuous(limits=c(-1,1.05),breaks=c(-1,0,1),labels=c("0.1","1.0","10"))+
  xlab("HR (95% CI)")+
  ylab(NULL)+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.5),
        axis.text=element_text(size=10),
        axis.text.x = element_text(hjust=1,vjust=0.5,color="black"),
        legend.position = "none") 

###fig6a
data<-read_excel("Source Data Fig. 6.xlsx", sheet = "fig6a")
roc_allproteins<-roc(data$healthy,as.numeric(data$predict_408))
roc_ageingproteins<-roc(data$healthy,as.numeric(data$predict_86))
roc_ageingproteins_top<-roc(data$healthy,as.numeric(data$predict_22))
plot(roc_allproteins,col="#fd9200",grid=c(0.2, 0.2))
plot(roc_ageingproteins,add=TRUE,col="#d33524")
plot(roc_ageingproteins_top,add=TRUE,col="#1576ae")
legend("bottomright", legend=c("AUC=0.72(408 proteins)","AUC=0.70 (86 ageing proteins)","AUC=0.70 (22 ageing proteins)"),col=c("#fd9200","#d33524","#1576ae"),lty=1)

###fig6b
data<-read_excel("Source Data Fig. 6.xlsx", sheet = "fig6b")
data %>% mutate(protein = fct_reorder(protein,MeanDecreaseAccuracy)) %>%
  ggplot () +
  geom_bar(aes(x=protein,y=MeanDecreaseAccuracy,fill=MeanDecreaseAccuracy),stat="identity") +
  coord_flip() +
  scale_y_continuous(limits=c(0,22),breaks=c(0,10,20))+
  scale_fill_gradient(low="#1576ae",high="#d33524")+
  xlab("")+
  theme_bw()+
  theme(panel.grid=element_blank(),panel.background = element_blank(),legend.position="none")

###fig6c
heatcoef<-read_excel("Source Data Fig. 6.xlsx", sheet = "fig6c")
heatcoef<-data.frame(heatcoef)
row.names(heatcoef)<-heatcoef$protein
heatcoef<-heatcoef[,c(2,5)]

bk <- c(seq(-0.5,-0.0001,by=0.01),seq(0.0001,0.5,by=0.01))
pheatmap(as.matrix(heatcoef),show_rownames=TRUE,show_colnames=FALSE,fontsize=12,scale="none",
         cluster_cols=F,cluster_rows=FALSE,
         color = c(colorRampPalette(colors = c("#4793c2","#90c5dd","#cce6f1","#f2f2f2"))(length(bk)/2),colorRampPalette(colors = c("#f2f2f2","#fddbc6","#f5a581","#d7604d"))(length(bk)/2)),
         breaks=bk)

###fig6d
data<-read_excel("Source Data Fig. 6.xlsx", sheet = "fig6d")
roc_full<-roc(data$healthy,as.numeric(data$predict_full))
roc_ageingproteins_top<-roc(data$healthy,as.numeric(data$predict_22))
roc_agesexbmi<-roc(data$healthy,as.numeric(data$predict_agesexbmi))

plot(roc_full,col="#d33524",grid=c(0.2, 0.2))
plot(roc_ageingproteins_top,add=TRUE,col="#1576ae")
plot(roc_agesexbmi,add=TRUE,col="#fd9200")
legend("bottomright",lty=1,legend=c("AUC=0.72 (Age,sex,BMI+22 ageing proteins)","AUC=0.70 (PHAS/22 ageing proteins)","AUC=0.63 (Age,sex,BMI)"),col=c("#d33524","#1576ae","#fd9200"))

###fig6e
data<-read_excel("Source Data Fig. 6.xlsx", sheet = "fig6e")
data$logq_discovery=-log10(data$q_discovery)

ggplot(data,aes(x=coef_discovery,y=logq_discovery))+geom_point(pch=21,color="black",aes(size=abs(coef_discovery),bg=group),show.legend = FALSE)+
  geom_text_repel(aes(label=trait),size=3,direction="both",force)+
  scale_x_continuous(limits=c(-1.5,1),breaks=c(-1.5,-1.0,-0.5,0,0.5,1.0),labels=c(expression("-1.5"),expression("-1.0"),expression("-0.5"),expression("0"),expression("0.5"),expression("1.0")))+
  theme_bw() +
  theme(panel.grid=element_blank(),panel.background = element_blank(),legend.box.background = element_rect(color="black"),legend.position = "right")+ 
  geom_hline(aes(yintercept=1.30103),color="darkgrey",linetype="dashed",size=0.5)+
  xlab("Effect size (PHAS)")+
  ylab("Q-value")

###fig6f
data<-read_excel("Source Data Fig. 6.xlsx", sheet = "fig6f")
ggplot(data,aes(x=coef_discovery,y=coef_validation))+
  geom_point(size=4,pch=21,color="black",aes(fill=group))+
  geom_text_repel(aes(label=label),size=4,direction="both",force)+
  theme_bw() +
  theme(legend.position = "none")+
  theme(panel.grid=element_blank(),legend.box.background = element_rect(color="black"),panel.background = element_blank())+#axis.line = element_line(size=0.1))+   ##去除网格线##去除外层边框##去除背景色##设置刻度线格式
  geom_hline(aes(yintercept=0),color="darkgrey",linetype="dashed",size=0.5)+geom_vline(aes(xintercept=0),color="darkgrey",linetype="dashed",size=0.5)+ ##设置参考线
  xlab("Effect size (Primary analysis)")+ylab("Effect size (Validation)")  ##设置坐标轴名称与字体

###fig6g
data<-read_excel("Source Data Fig. 6.xlsx", sheet = "fig6g")
ggplot(data,aes(y=reorder(Diseases,order), x=loghr))+
  geom_vline(aes(xintercept=0),color="darkgrey",linetype="dashed",size=0.5)+
  scale_color_manual(values=c("#1576ae","black"))+ 
  scale_shape_manual(values=c(16,1))+
  geom_errorbar(aes(xmin=logll, xmax=logul, color=group), width = 0.2,size=0.5)+
  geom_point(size=4,aes(shape=group,color=group))+
  scale_x_continuous(limits=c(-1,0.6989700043360187),breaks=c(-1,0,0.6989700043360187),labels=c("0.1","1","5"))+
  xlab("HR (95% CI)")+
  ylab(NULL)+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.5),
        axis.text=element_text(size=10),
        axis.text.x = element_text(hjust=1,vjust=0.5,color="black"),
        legend.position = "None") 

ggplot(data,aes(y=reorder(Diseases,order), x=loghr2))+
  geom_vline(aes(xintercept=0),color="darkgrey",linetype="dashed",size=0.5)+ ##设置参考线
  scale_color_manual(values=c("#1576ae","black"))+ 
  scale_shape_manual(values=c(16,1))+
  geom_errorbar(aes(xmin=logll2, xmax=logul2, color=group2), width = 0.2,size=0.5)+
  geom_point(size=4,aes(shape=group2,color=group2))+
  scale_x_continuous(limits=c(-1.222,0.6989700043360187),breaks=c(-1,0,0.6989700043360187),labels=c("0.1","1","5"))+
  xlab("HR (95% CI)")+
  ylab(NULL)+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.5),
        axis.text=element_text(size=10),
        axis.text.x = element_text(hjust=1,vjust=0.5,color="black"),
        legend.position = "None")

###fig7a
data<-read_excel("Source Data Fig. 7.xlsx", sheet = "fig7a")

data %>% mutate(var = fct_reorder(var,desc(n))) %>%
  ggplot () +
  geom_bar(aes(x=var,y=permanova,fill=as.factor(n)),width=0.6,stat="identity") +
  coord_flip() +
  scale_y_continuous(limits=c(0,11),breaks=c(0,5,10))+
  scale_fill_manual(values=c("#1576ae","#95c5f0","#cccccc","#ff7f00","#d33524"))+
  xlab("")+
  theme_bw()+
  theme(panel.grid=element_blank(),panel.background = element_blank(),legend.position="none")

data %>% mutate(var = fct_reorder(var,desc(n))) %>%
  ggplot () +
  geom_bar(aes(x=var,y=lm,fill=as.factor(n)),width=0.6,stat="identity") +
  coord_flip() +
  scale_y_continuous(limits=c(0,10),breaks=c(0,5,10))+
  scale_fill_manual(values=c("#1576ae","#95c5f0","#cccccc","#ff7f00","#d33524"))+
  xlab("")+
  theme_bw()+
  theme(panel.grid=element_blank(),panel.background = element_blank(),legend.position="none")


###fig7b
variance_protein<-read_excel("Source Data Fig. 7.xlsx", sheet = "fig7b")

variance_protein_age<-subset(variance_protein,group==1)
variance_protein_age<-variance_protein_age[,c(7,3:6,2)]
variance_protein_age<-melt(variance_protein_age,id.vars='Protein')

variance_protein_diet<-subset(variance_protein,group==3)
variance_protein_diet<-variance_protein_diet[,c(7,2,3,5,6,4)]
variance_protein_diet<-melt(variance_protein_diet,id.vars='Protein')

variance_protein_meta<-subset(variance_protein,group==4)
variance_protein_meta<-variance_protein_meta[,c(7,2,3,4,6,5)]
variance_protein_meta<-melt(variance_protein_meta,id.vars='Protein')

variance_protein_genome<-subset(variance_protein,group==5)
variance_protein_genome<-variance_protein_genome[,c(7,2:6)]
variance_protein_genome<-melt(variance_protein_genome,id.vars='Protein')

variance_protein_age%>%mutate(n=row_number())%>%mutate(Protein = fct_reorder(Protein,n))%>%
  ggplot(aes(x=Protein,value,fill=variable))+
  geom_bar(stat="identity",position="stack", color="black", width=0.6,size=0.25)+
  scale_fill_manual(values=c("#95c5f0","#cccccc","#ff7f00","#d33524","#1576ae"))+
  labs(y = "Adjusted R2")+
  scale_y_continuous(breaks=seq(0,40,10),limits=c(0,40),expand=c(0,0))+
  xlab("")+
  theme_bw()+
  theme(panel.grid=element_blank(),panel.background = element_blank(),legend.position="none")

variance_protein_diet%>%mutate(n=row_number())%>%mutate(Protein = fct_reorder(Protein,n))%>%
  ggplot(aes(x=Protein,value,fill=variable))+
  geom_bar(stat="identity",position="stack", color="black", width=0.6,size=0.25)+
  scale_fill_manual(values=c("#1576ae","#95c5f0","#ff7f00","#d33524","#cccccc"))+
  labs(y = "Adjusted R2")+
  scale_y_continuous(breaks=seq(0,40,10),limits=c(0,40),expand=c(0,0))+
  xlab("")+
  theme_bw()+
  theme(panel.grid=element_blank(),panel.background = element_blank(),legend.position="none")

variance_protein_meta%>%mutate(n=row_number())%>%mutate(Protein = fct_reorder(Protein,n))%>%
  ggplot(aes(x=Protein,value,fill=variable))+
  geom_bar(stat="identity",position="stack", color="black", width=0.6,size=0.25)+
  scale_fill_manual(values=c("#1576ae","#95c5f0","#cccccc","#d33524","#ff7f00"))+
  labs(y = "Adjusted R2")+
  scale_y_continuous(breaks=seq(0,40,10),limits=c(0,40),expand=c(0,0))+
  xlab("")+
  theme_bw()+
  theme(panel.grid=element_blank(),panel.background = element_blank(),legend.position="none")

variance_protein_genome%>%mutate(n=row_number())%>%mutate(Protein = fct_reorder(Protein,n))%>%
  ggplot(aes(x=Protein,value,fill=variable))+
  geom_bar(stat="identity",position="stack", color="black", width=0.6,size=0.25)+
  scale_fill_manual(values=c("#1576ae","#95c5f0","#cccccc","#ff7f00","#d33524"))+
  labs(y = "Adjusted R2")+
  scale_y_continuous(breaks=seq(0,40,10),limits=c(0,40),expand=c(0,0))+
  xlab("")+
  theme_bw()+
  theme(panel.grid=element_blank(),panel.background = element_blank(),legend.position="none")

###fig7c
data<-read_excel("Source Data Fig. 7.xlsx", sheet = "fig7c")

ggplot(data,aes(x=coef,y=logp))+
  geom_point(data=subset(data,group=="c"),aes(size=logq),pch=21,color="black",fill="darkgrey",show.legend = FALSE)+
  geom_point(data=subset(data,group=="a"),aes(size=logq),pch=21,color="black",fill="#d33524",show.legend = FALSE)+
  geom_point(data=subset(data,group=="b"),aes(size=logq),pch=21,color="black",fill="#1576ae",show.legend = FALSE)+
  geom_text_repel(aes(label=label),direction="both",force)+
  scale_x_continuous(limits=c(-0.012,0.012),breaks=c(-0.01,0,0.01))+
  scale_y_continuous(limits=c(0.823909,3),breaks=c(0,1,2,3),labels=c(expression("1"),expression("1x10"^"-1"),expression("1x10"^"-2"),expression("1x10"^"-3")))+
  theme_bw() +
  theme(panel.grid=element_blank(),legend.box.background = element_rect(color="black"),panel.background = element_blank())+
  geom_hline(aes(yintercept=1.30103),color="darkgrey",linetype="dashed",size=0.5)+
  xlab("Effect size (Microbioal species)")+ylab("Q-value")


###fig7d
data<-read_excel("Source Data Fig. 7.xlsx", sheet = "fig7d")
ggplot(data)+
  geom_point(aes(y=phas,x=metascore),shape=21,color="black",fill="#0d7091",size=4)+
  geom_smooth(aes(y=phas,x=metascore),method="lm",color="red",fill="grey",se=TRUE)+
  labs(x="Microbial score")+
  labs(y="Proteomic healthy ageing score")+
  scale_y_continuous(limits=c(0,0.8),breaks = c(0,0.4,0.8))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.position = "bottom")


###fig7e
data<-read_excel("Source Data Fig. 7.xlsx", sheet = "fig7e")
ggplot(data)+
  geom_point(aes(y=phas,x=metascore),shape=21,color="black",fill="#0d7091",size=4)+
  geom_smooth(aes(y=phas,x=metascore),method="lm",color="red",fill="grey",se=TRUE)+
  labs(x="Microbial score")+
  labs(y="Proteomic healthy ageing score")+
  scale_y_continuous(limits=c(0,0.8),breaks = c(0,0.4,0.8))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.position = "bottom")
