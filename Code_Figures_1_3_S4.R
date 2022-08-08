.libPaths("C:/project/libraryOld/")
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(vegan)
library(caret)
library(caret)

# functions ---------------------------------------------------------------

getwilcoxTest <- function(caseVec, contVec, imputeVal=10) {
  wtest = wilcox.test(caseVec, contVec, alternative = "two.sided")
  return(wtest$p.value)
}
getwilcoxTestMat <- function(caseMat, contMat, imputeVal=10) {
  caseIds = colnames(caseMat)
  contIds = colnames(contMat)
  caseContMat = cbind(caseMat, contMat)
  caseInd = which(colnames(caseContMat) %in% caseIds)
  contInd = which(colnames(caseContMat) %in% contIds)
  wilcoxTest = apply(caseContMat, 1, function(currVec) {
    caseVec = currVec[caseInd]
    contVec = currVec[contInd]
    ES = getwilcoxTest(caseVec, contVec)
    return(ES)
  })
  return(wilcoxTest)
}
################################ Input data ################################################################
inputData=readRDS("Fig_1.rds")
inputData1=readRDS("Fig_3_S4.rds")
################################  Figure 1 #################################################################
inp=inputData$Reactobiome
outp <- inputData$metadata
inp=inp[rownames(outp),]
all(rownames(inp)==rownames(outp))
iris.scal <- cmdscale(dist(log2(inp+1)),k=2)
PcoA2D_AT <- as.data.frame(iris.scal)

colnames(PcoA2D_AT)=c("PCoA1","PCoA2")
rownames(PcoA2D_AT)=rownames(inp)
####fig 1b
ggplot(outp, aes(x = MAGMA_richness, y = Reaction_richness))+geom_point(shape=21, color="#4682B4", size=3,alpha = 0.6)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

####fig 1c
inp=as.matrix(inp)
class(inp)="numeric"

data1=as.data.frame(inp)
a=nearZeroVar(data1)
inp=inp[,-a]
rem1=c()
for (i in 1:nrow(inp)) {
  a=inp[i,]
  g=i+1
  h=c()
  cnt=1
  if (g <=nrow(inp)) {
  for (j in c(g:nrow(inp))) {
    b=inp[j,]
    c=a==b
    h[cnt]=all(c==TRUE)
    cnt=cnt+1
  }
  }
  rem1[i]=any(h==TRUE)
}
rem2=which(rem1)
inp1=inp[-rem2,]

example_NMDS=metaMDS(log2(inp1+1), # Our community-by-species matrix
                     k=2) # The number of reduced dimensions
stressplot(example_NMDS,pch = 16,p.col = "grey")

###fig 1d
ggplot(PcoA2D_AT, aes(x = PCoA1, y = PCoA2))+geom_point(shape=21, color="#4682B4", size=3,alpha = 0.6)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

PCoA1=tapply(PcoA2D_AT$PCoA1,outp$Geography,mean)
PCoA2=tapply(PcoA2D_AT$PCoA2,outp$Geography,mean)
var1=tapply(PcoA2D_AT$PCoA1,outp$Geography,var)
var2=tapply(PcoA2D_AT$PCoA2,outp$Geography,var)
geo=data.frame(cbind(PCoA1,PCoA2,var1,var2))
geo[,5]=sqrt(geo$var1+geo$var2)

ggplot(geo, aes(x = PCoA1, y = PCoA2,label=rownames(geo)))+geom_point(shape=21, color="Black",fill="navy", size=geo$V5/3,alpha = 0.3)+
  scale_fill_manual(values=c("#4682B4","#1E90FF","#0000CD",	"#AFEEEE","#5F9EA0"))+
  geom_text(aes(label=rownames(geo)),hjust=0, vjust=0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

####fig 1e

PcoA2D_AT2=cbind(PcoA2D_AT,outp$Reaction_richness)
colnames(PcoA2D_AT2)[3]="Reaction_richness"

ggplot(PcoA2D_AT2, aes(x = PCoA1, y = Reaction_richness))+geom_point(shape=21, color="grey",size=2,alpha = 0.9)+
  geom_smooth(size=2)+
  ylim(3000,4200)+
  # scale_fill_manual(values=c("#4682B4","#1E90FF","#0000CD",	"#AFEEEE","#5F9EA0"))+
  #scale_fill_continuous(type = "viridis")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

###Figure 1f 
ggplot(PcoA2D_AT, aes(x = PCoA1, y = PCoA2,colour = outp$Reactotype,fill= outp$Reactotype))+geom_point(shape=21, color="black", size=3,alpha = 0.6,group=outp$Reactotype)+
  scale_fill_manual(values=c('#F0F9E866','#BAE4BC66','#7BCCC466','#43A2CA66','#0868AC66'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

mosaicplot(Reactotype ~ enteroType, data = outp,las=1,cex.axis=1,col=colorRampPalette(brewer.pal(9,"Blues"))(3))

###Figure 1g
ggplot(outp, aes(x=Reaction_richness,fill=Reactotype,colour=Reactotype))+geom_density(size=1)+theme_classic()+ xlim(3200,4136) +
  scale_color_manual(values=c('#F0F9E8','#BAE4BC','#7BCCC4','#43A2CA','#0868AC'))+
  scale_fill_manual(values=c('#F0F9E866','#BAE4BC66','#7BCCC466','#43A2CA66','#0868AC66'))

###Figure 1h, 
mosaicplot(Reactotype ~ Continent, data = outp,las=1,cex.axis=1,col=colorRampPalette(brewer.pal(9,"Blues"))(5))


######################################### Fig 3    ##############################################################################

###Fig 3b
fluxpc <- inputData1$FluxCP1_reto1reto4
flux <- inputData1$Flux1_reto1reto4
rownames(fluxpc)==rownames(flux)

flu=as.data.frame(t(flux[1:2,]))

ggplot(data=flu,aes(x=group,y=Biomass_Bacteria,group=group,colour = group))+ geom_boxplot(notch = T,outlier.fill =NA,outlier.stroke = F,
                                                                                          outlier.size = 0.3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 45,hjust = 1))

fluxpc1=as.matrix(fluxpc[-1,])
flux1=as.matrix(flux[-1,])
flux1p=flux1;flux1c=flux1
fluxpc1p=fluxpc1;fluxpc1c=fluxpc1
flux1p[flux1p<0]=0;flux1c[flux1c>0]=0
fluxpc1p[fluxpc1p<0]=0;fluxpc1c[fluxpc1c>0]=0
class(fluxpc1p)="numeric"
class(fluxpc1c)="numeric"
class(flux1p)="numeric"
class(flux1c)="numeric"

fluxpc1p=fluxpc1p[rowSums(fluxpc1p)>0,]
fluxpc1c=fluxpc1c[rowSums(abs(fluxpc1c))>0,]
flux1p=flux1p[rowSums(flux1p)>0,]
flux1c=flux1c[rowSums(abs(flux1c))>0,]

rownames(flux1p)=paste("prod_",rownames(flux1p),sep="")
rownames(flux1c)=paste("cons_",rownames(flux1c),sep="")

rownames(fluxpc1p)=paste("prod_",rownames(fluxpc1p),sep="")
rownames(fluxpc1c)=paste("cons_",rownames(fluxpc1c),sep="")

fpSum=log2(abs(apply(flux1p[,flux[1,]=="increased"],1,sum)))-log2(abs(apply(flux1p[,flux[1,]=="decreased"],1,sum)))
fpmean=log2(abs(apply(flux1p[,flux[1,]=="increased"],1,mean)))-log2(abs(apply(flux1p[,flux[1,]=="decreased"],1,mean)))
fcSum=log2(abs(apply(flux1c[,flux[1,]=="increased"],1,sum)))-log2(abs(apply(flux1c[,flux[1,]=="decreased"],1,sum)))
fcmean=log2(abs(apply(flux1c[,flux[1,]=="increased"],1,mean)))-log2(abs(apply(flux1c[,flux[1,]=="decreased"],1,mean)))

fpSum=abs(fpSum)
fpmean=abs(fpmean)

fcSum=-abs(fcSum)
fcmean=-abs(fcmean)

f1=rbind(cbind(fpSum,fpmean),cbind(fcSum,fcmean))


p=as.data.frame(log2(abs(apply(fluxpc1p[,fluxpc[1,]=="increased"],1,sum))/47)-log2(abs(apply(fluxpc1p[,fluxpc[1,]=="decreased"],1,sum))/48))
c=as.data.frame(log2(abs(apply(fluxpc1c[,fluxpc[1,]=="increased"],1,sum))/47)-log2(abs(apply(fluxpc1c[,fluxpc[1,]=="decreased"],1,sum))/48))

colnames(p)="LogFC"
colnames(c)="LogFC"
pc=rbind(p,c)
f1=f1[rownames(pc),]
rownames(f1)==rownames(pc)
final=cbind(f1,pc)
final[final==Inf]=1
final[final==-Inf]=-1
ggplot(final, aes(x = LogFC, y = fpmean))+geom_point(size =2,alpha=.5,shape=16,colour="blue")+ylim(-3,3)+scale_x_continuous()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


####figure 3c
fluxpc <- read.table("Biomass.txt", header=T, sep = "\t", row.names = 1, as.is=TRUE)
library(reshape2)
reza=melt(fluxpc)

ggplot(data=reza,aes(x=variable,y=value,colour = group))+ geom_boxplot(notch = F,outlier.fill =NA,outlier.stroke = F,
                                                                       outlier.size = 0.3)+ylim(0,3.2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 45,hjust = 1))


###Fig 3d
yy <-inputData1$FluxPersonalized

yy$significance <- -log10(yy$FDR)
yy$differentially_expressed <- yy$FDR < 0.0001
ggplot(yy, aes(x = V4,y=significance,col=differentially_expressed)) + geom_point()+xlab("Flux ratio")+
  scale_y_continuous(trans = 'log2')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("-log10(Padj)")

################################ supplementary fig 4####################################
####supplementary fig 4a
PSE_output_coverage=inputData1$PSE_output_coverage_0.01
metadata <- inputData1$metadata
metadata=metadata[colnames(PSE_output_coverage),]

PSE_output_coverage=as.matrix(t(PSE_output_coverage))


PSE_output_coverage_up=PSE_output_coverage[metadata$Reactotype=="Reto-1",]
PSE_output_coverage_down=PSE_output_coverage[metadata$Reactotype=="Reto-4",]

caseMat=as.matrix(t(PSE_output_coverage_up))
class(caseMat)="numeric"
contMat=as.matrix(t(PSE_output_coverage_down))
class(contMat)="numeric"
wt1=as.data.frame(getwilcoxTestMat(caseMat,contMat))

s=apply(caseMat,1,mean)-apply(contMat,1,mean)
s1=apply(caseMat,1,mean)/apply(contMat,1,mean)
s2=log2(apply(caseMat,1,mean))-log2(apply(contMat,1,mean))

wt1[,2]=apply(caseMat,1,mean)
wt1[,3]=apply(contMat,1,mean)
wt1[,4]=s
wt1[,5]=s1
wt1[,6]=s2

#calculating FDR
g=as.vector(wt1[,1])
a=p.adjust(g, method = c("BH"), n = length(g))
#append FDR to Sif
sif=cbind(wt1,FDR=a)

yy <- sif

yy$significance <- -log10(yy$FDR)
yy$differentially_expressed <- yy$FDR < 0.01
ggplot(yy, aes(x = V6,y=significance,col=differentially_expressed)) + geom_point()+xlab("Coverage of enrichment (logFoldChange)")+scale_y_continuous(trans = 'log2')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("-log10(Padj)")

########supplementary fig 4d

PSE_output_coverage=inputData1$PSE_output_coverage_0.001
metadata <- inputData1$metadata
metadata=na.omit(metadata)
tra1=intersect(rownames(metadata),colnames(PSE_output_coverage))
metadata=metadata[tra1,];PSE_output_coverage=PSE_output_coverage[,tra1]
PSE_output_coverage=as.matrix(t(PSE_output_coverage))
PSE_output_coverage_up=PSE_output_coverage[metadata$Continent=="Tr",]
PSE_output_coverage_down=PSE_output_coverage[metadata$Reactotype=="Reto-4",]

caseMat=as.matrix(t(PSE_output_coverage_up))
class(caseMat)="numeric"
contMat=as.matrix(t(PSE_output_coverage_down))
class(contMat)="numeric"
wt1=as.data.frame(getwilcoxTestMat(caseMat,contMat))

s=apply(caseMat,1,mean)-apply(contMat,1,mean)
s1=apply(caseMat,1,mean)/apply(contMat,1,mean)
s2=log2(apply(caseMat,1,mean))-log2(apply(contMat,1,mean))
s3=log(apply(caseMat,1,mean))-log(apply(contMat,1,mean))
wt1[,2]=apply(caseMat,1,mean)
wt1[,3]=apply(contMat,1,mean)
wt1[,4]=s
wt1[,5]=s1
wt1[,6]=s2
wt1[,7]=s3
#calculating FDR
g=as.vector(wt1[,1])
a=p.adjust(g, method = c("BH"), n = length(g))
#append FDR to Sif
sif=cbind(wt1,FDR=a)
yy <- sif

yy$significance <- -log10(yy$FDR)
yy$differentially_expressed <- yy$FDR < 0.01
ggplot(yy, aes(x = V6,y=significance,col=differentially_expressed)) + geom_point()+xlab("Coverage of enrichment (logFoldChange)")+#scale_y_continuous(trans = 'log2')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("-log10(Padj)")







