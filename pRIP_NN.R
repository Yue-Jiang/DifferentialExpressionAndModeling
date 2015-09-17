#load stuff
library(mixtools)
library(seqinr)
library(ape)
library(edgeR)
library(DESeq)
library(lattice)
library(MASS)
library(gplots)

MyGray <- rgb(t(col2rgb("black")), alpha=60, maxColorValue=255)
MyTransparent <- rgb(t(col2rgb("black")), alpha=0, maxColorValue=255)
MyGreen <- rgb(t(col2rgb("green")), alpha=150, maxColorValue=255)
MyRed <- rgb(t(col2rgb("red")), alpha=150, maxColorValue=255)
MyBlue <- rgb(t(col2rgb("blue")), alpha=150, maxColorValue=255)
MyPurple <- rgb(t(col2rgb("purple")), alpha=150, maxColorValue=255)

#differential expression analysis
DiffExp <- function (targets, countsTable) {
  Treat <- factor(targets$Treatment);Subject <- factor(targets$Subject);design <- model.matrix(~Subject+Treat)
  cds <- newCountDataSet(countsTable,Treat);cds <- estimateSizeFactors(cds);cds <- estimateDispersions(cds);d <- nbinomTest(cds,"0","1")
  e.litter <- DGEList(counts=countsTable)
  e.litter <- estimateGLMCommonDisp(e.litter,design)
  e.litter <- estimateGLMTrendedDisp(e.litter,design)
  e.litter <- estimateGLMTagwiseDisp(e.litter,design)
  fit <- glmFit(e.litter, design);lrt <- glmLRT(fit);diff <- topTags(lrt,n=dim(lrt)[1])$table
  result <- merge(merge(diff,countsTable,by=0,sort=F),d, by.x="Row.names", by.y="id",sort=F)
  colnames(result)[1] <- "id"
  return(result)
}
result.No_Ac <- DiffExp(targets.Ac,countsTable.No_Ac)
result.No_1pAc <- DiffExp(targets.Ac,countsTable.No_1pAc)
result.No_TMT <- DiffExp(targets.TMT,countsTable.No_TMT)
result.No_1pTMT <- DiffExp(targets.TMT,countsTable.No_1pTMT)

result.No_Ac$id[result.No_Ac$id=="Olfr175-ps1"] <- "Olfr174"
result.No_1pAc$id[result.No_1pAc$id=="Olfr175-ps1"] <- "Olfr174"
result.No_TMT$id[result.No_TMT$id=="Olfr175-ps1"] <- "Olfr174"
result.No_1pTMT$id[result.No_1pTMT$id=="Olfr175-ps1"] <- "Olfr174"

#fig2
filter.OR <- function(result,cpm=0,add=0.1) { #filter ORs by CPM, correct FDR 
  p1=0.05;p2=0.001;p3=0.001
  result.OR <- subset(result,grepl('Olfr',id))#|grepl('Taar',id)|grepl('Vmn',id))
  result.nonOR <- subset(result,!grepl('Olfr',id))
  A.cut <- cpm*sum(result$baseMeanA)/1000000; B.cut <- cpm*sum(result$baseMeanB)/1000000
  cut <- cpm*sum(result$baseMean)/1000000
  #  A.cut <- B.cut <- cpm
  sub.result.OR <- subset(result.OR, baseMeanA>=A.cut&baseMeanB>=B.cut)
  sub.result.nonOR <- subset(result.nonOR, baseMeanA>=A.cut&baseMeanB>=B.cut)
  EGR1 <- subset(result.nonOR, id=="Egr1")
  CFOS <- subset(result.nonOR, id=="Fos")
  M72 <- subset(result.OR,id=="Olfr160")
  plot(log10(sub.result.nonOR$baseMeanA+add),log10(sub.result.nonOR$baseMeanB+add),type="n",frame.plot=F,xlab=NA,ylab=NA)
  points(log10(sub.result.nonOR$baseMeanA+add),log10(sub.result.nonOR$baseMeanB+add),col=MyGray,pch=16)
  title(xlab=expression('Log'[10]*'(mean read counts, no-odor control)'), ylab=expression('Log'[10]*'(mean read counts, stimulated)') ,cex.lab=1.5)
  #title(xlab="Mean read counts, no-odor control, log10", ylab="Mean read counts, stimulated, log10",cex.lab=1.5)
  points(log10(sub.result.OR$baseMeanA+add),log10(sub.result.OR$baseMeanB+add),bg=MyRed,col=MyGray2,pch=21)
  points(log10(EGR1$baseMeanA+add),log10(EGR1$baseMeanB+add),bg="cyan",pch=21);points(log10(CFOS$baseMeanA+add),log10(CFOS$baseMeanB+add),bg="green",pch=21)
  abline(0,1,lwd=3,lty=2,col="gray")
  legend("bottomright",c("ORs","c-Fos","Egr1","Others"),pt.bg=c(MyRed,"green","cyan",MyTransparent),col=c("black","black","black",MyGray),pch=c(21,21,21,16),bty="n",cex=1.2)
  sub.result.OR$FDR <- p.adjust(sub.result.OR$PValue,method="fdr")
  p1sub <- subset(sub.result.OR,(FDR>p2&FDR<p1)&logFC>0); p2sub <- subset(sub.result.OR,(FDR>p3&FDR<p2)&logFC>0); p3sub <- subset(sub.result.OR,FDR<p3&logFC>0)
  
  MyGray <- "gray"
  plot(log10(sub.result.OR$baseMeanA+add),log10(sub.result.OR$baseMeanB+add),type="n",frame.plot=F,xlab=NA,ylab=NA)
  points(log10(sub.result.OR$baseMeanA+add),log10(sub.result.OR$baseMeanB+add),
         bg=ifelse(sub.result.OR$FDR<p1&sub.result.OR$logFC>0, MyTransparent, MyGray),
         col=ifelse(sub.result.OR$FDR<p1&sub.result.OR$logFC>0, MyTransparent, MyGray),pch=21,cex=1.2)
  title(xlab=expression('Log'[10]*'(mean read counts, no-odor control)'), ylab=expression('Log'[10]*'(mean read counts, stimulated)') ,cex.lab=1.5)
  #title(xlab="Mean read counts, no-odor control, log10", ylab="Mean read counts, stimulated, log10",cex.lab=1.5)
  points(log10(p1sub$baseMeanA+add),log10(p1sub$baseMeanB+add),col=MyGray2,bg=MyBlue,pch=21,cex=1.2)
  points(log10(p2sub$baseMeanA+add),log10(p2sub$baseMeanB+add),col=MyGray2,bg=MyPurple,pch=21,cex=1.2)
  points(log10(p3sub$baseMeanA+add),log10(p3sub$baseMeanB+add),col=MyGray2,bg=MyRed,pch=21,cex=1.2)
  abline(0,1,lwd=3,lty=2)
  #  points(log10(M72$baseMeanA),log10(M72$baseMeanB),col=MyGray2,bg="cyan",pch=21)
  #legend("bottomright",c(paste("FDR adjusted p<",p1),paste("FDR adjusted p<",p2),paste("FDR adjusted p<",p3),paste("#ORs enriched:",dim(subset(sub.result.OR,FDR<vivo.cut&logFC>0))[1])),bg=c(MyBlue,MyPurple,MyRed,MyTransparent),pch=21,bty="n",cex=0.9)
  #  legend("bottomright",c(paste("FDR adjusted p<",p1),paste("FDR adjusted p<",p3),paste("#ORs enriched:",dim(subset(sub.result.OR,FDR<vivo.cut&logFC>0))[1])),col=c(MyGray2,MyGray2,MyTransparent),pt.bg=c(MyBlue,MyRed,MyTransparent),pch=21,bty="n",cex=0.9)
  #  legend("bottomright",c(paste("p<",p1),paste("p<",p3),paste("#ORs enriched:",dim(subset(sub.result.OR,FDR<p1&logFC>0))[1])),col=c(MyGray2,MyGray2,MyTransparent),pt.bg=c(MyBlue,MyRed,MyTransparent),pch=21,cex=1.2,bty="n")
  #legend("bottomright",c(paste("p<",p3," (n = ",dim(subset(sub.result.OR,FDR<p3&logFC>0))[1],")",sep=""),paste("p< ",p1," (n=",dim(subset(sub.result.OR,FDR<p1&logFC>0))[1],")",sep="")),col=c(MyGray2,MyGray2),pt.bg=c(MyRed,MyBlue),pch=21,cex=1.2,bty="n")
  lgnd <- c(bquote(italic(p) < .(p3)*' '*(n == .(dim(subset(sub.result.OR,FDR<p3&logFC>0))[1]))),
            bquote(italic(p) < .(p1)*' '*(n == .(dim(subset(sub.result.OR,FDR<p1&logFC>0))[1]))))
  legend("bottomright", as.expression(lgnd),,col=c(MyGray2,MyGray2),pt.bg=c(MyRed,MyBlue),pch=21,cex=1.2,bty="n")
  MyGray <- rgb(t(col2rgb("black")), alpha=60, maxColorValue=255)
  plot(sub.result.OR$logFC,-log10(sub.result.OR$FDR),type="n",frame.plot=F,xlab=NA,ylab=NA)
  points(sub.result.OR$logFC,-log10(sub.result.OR$FDR),
         bg=ifelse(sub.result.OR$FDR<p1&sub.result.OR$logFC>0, MyTransparent, MyGray),
         col=ifelse(sub.result.OR$FDR<p1&sub.result.OR$logFC>0, MyTransparent, MyGray),pch=21,cex=1.2)
  title(xlab=expression('Log'[2]*'(fold change)'), ylab=expression('-Log'[10]*'(p-value)') ,cex.lab=1.5)
  #title(xlab="Mean read counts, no-odor control, log10", ylab="Mean read counts, stimulated, log10",cex.lab=1.5)
  points(p1sub$logFC,-log10(p1sub$FDR),col=MyGray2,bg=MyBlue,pch=21,cex=1.2)
  points(p2sub$logFC,-log10(p2sub$FDR),col=MyGray2,bg=MyPurple,pch=21,cex=1.2)
  points(p3sub$logFC,-log10(p3sub$FDR),col=MyGray2,bg=MyRed,pch=21,cex=1.2)
  legend("topleft", as.expression(lgnd),,col=c(MyGray2,MyGray2),pt.bg=c(MyRed,MyBlue),pch=21,cex=1.2,bty="n")
  #legend("topleft",c(paste("p<",p3," (n=",dim(subset(sub.result.OR,FDR<p3&logFC>0))[1],")",sep=""),paste("p<",p1," (n=",dim(subset(sub.result.OR,FDR<p1&logFC>0))[1],")",sep="")),col=c(MyGray2,MyGray2),pt.bg=c(MyRed,MyBlue),pch=21,cex=1.2,bty="n")
  
  
  return (sub.result.OR)
}
###
pdf("figure2.pdf",5.2,5)
par(mar=c(4,4.4,2,2)+0.1)
#for (i in c(0.001,0.01,0.1,1)) {
result.No_Ac.OR <- filter.OR (result.No_Ac,add=1,cpm=0)
result.No_1pAc.OR <- filter.OR (result.No_1pAc,add=1,cpm=0)

result.No_TMT.OR <- filter.OR (result.No_TMT,add=1)
result.No_1pTMT.OR <- filter.OR (result.No_1pTMT,add=1)
dev.off()

pdf("figure_HP.pdf",5.2,5)
par(mar=c(4,4.4,2,2)+0.1)
#for (i in c(0.001,0.01,0.1,1)) {
result.No_HP.OR <- filter.OR (result.No_HP,add=1,cpm=0)
result.No_1pHP.OR <- filter.OR (result.No_1pHP,add=1,cpm=0)
dev.off()

################################################## shared ORs, 1% vs 100%
shared.OR <- function(result.OR,result1p.OR) {
  merged.OR <- merge(result1p.OR,result.OR,by="id")
  p1=0.05;p2=0.01;p3=0.001
  plot(log10(merged.OR$FDR.x),log10(merged.OR$FDR.y),col=MyGray2,bg="gray",pch=21,frame.plot=F,
       xlab=NA,ylab=NA,
       #         xlim=c(0,-9),ylim=c(0,-9))
       xlim=c(0,log10(min(result1p.OR$FDR))-3),ylim=c(0,log10(min(result.OR$FDR))-3))
  title(xlab=expression('Log'[10]*'(p-value,1% stimulation)'), ylab=expression('Log'[10]*'(p-value,100% stimulation)') ,cex.lab=1.5)
  #    title(xlab="p-value,1% stimulation,log10",ylab="p-value,100% stimulation,log10",cex.lab=1.5)
  abline(v=log10(p1),lty=2,col="blue",lwd=2);abline(h=log10(p1),lty=2,col="blue",lwd=2)
  #    abline(v=log10(p2),lty=2,col="purple",lwd=2);abline(h=log10(p2),lty=2,col="purple",lwd=2)
  abline(v=log10(p3),lty=2,col="red",lwd=2);abline(h=log10(p3),lty=2,col="red",lwd=2)
  #legend("topright",lty=2,lwd=2,col=c("red","blue"),c("p=0.001","p=0.05"),bty="n",cex=1.2)
  legend("topright",lty=2,lwd=2,col=c("red","blue"),c(expression(italic(p)==0.001),expression(italic(p)==0.05)),bty="n",cex=1.2)
  #    sub1 <- subset(merged.OR,(FDR.x<0.05|FDR.y<0.05)&(logFC.x>0|logFC.y>0))
  #    sub2 <- subset(merged.OR,(FDR.x<0.001|FDR.y<0.001)&(logFC.x>0|logFC.y>0))
  #    sub3 <- subset(merged.OR,(FDR.x<0.05&(FDR.y>0.05|logFC.y<=0))&(logFC.x>0))
  #    plot(sub1$logFC.x,sub1$logFC.y,col=MyGray2,bg="blue",pch=21,frame.plot=F,
  #    xlab=NA,ylab=NA)
  #    points(sub2$logFC.x,sub2$logFC.y,col=MyGray2,bg="red",pch=21)
  #    points(sub3$logFC.x,sub3$logFC.y,col=MyGray2,bg="green",pch=21)
  #     #xlim=c(0,log10(min(result1p.OR$FDR))),ylim=c(0,log10(min(result.OR$FDR))))
  #    title(xlab="logFC,1% stimulation",ylab="logFC,100% stimulation",cex.lab=1.5)
  #    abline(v=log10(p1),lty=2,col="blue",lwd=2);abline(h=log10(p1),lty=2,col="blue",lwd=2)
  #    abline(v=log10(p2),lty=2,col="purple",lwd=2);abline(h=log10(p2),lty=2,col="purple",lwd=2)
  #    abline(v=log10(p3),lty=2,col="red",lwd=2);abline(h=log10(p3),lty=2,col="red",lwd=2)
  #    legend("topright",lty=2,col=c("red","blue"),c("p=0.001","p=0.05"),bty="n",cex=1.5)
}
###
pdf("figure_sharedOR.pdf",5.2,5)
par(mar=c(4,4.4,2,2)+0.1)
shared.OR(result.No_Ac.OR,result.No_1pAc.OR)
shared.OR(result.No_TMT.OR,result.No_1pTMT.OR)
dev.off()

#fig3
Box <- function(result.OR,luc,vivo.cut,add=1) {
  MyGray="gray"
  p1=0.05;p2=0.001;p3=0.001
  p1sub <- subset(result.OR,(FDR>p2&FDR<p1)&logFC>0); p2sub <- subset(result.OR,(FDR>p3&FDR<p2)&logFC>0); p3sub <- subset(result.OR,FDR<p3&logFC>0)
  plot(log10(result.OR$baseMeanA+add),log10(result.OR$baseMeanB+add),
       bg=ifelse(result.OR$FDR<p1&result.OR$logFC>0, MyTransparent, MyGray),
       col=ifelse(result.OR$FDR<p1&result.OR$logFC>0, MyTransparent, MyGray),pch=21,frame.plot=F,xlab=NA,ylab=NA)
  title(xlab=expression('Log'[10]*'(mean read counts, no-odor control)'), ylab=expression('Log'[10]*'(mean read counts, stimulated)') ,cex.lab=1.5)
  #title(xlab="Mean read counts, unstimulated, log10", ylab="Mean read counts, stimulated, log10",cex.lab=1.5)
  points(log10(p1sub$baseMeanA+add),log10(p1sub$baseMeanB+add),col=MyTransparent,bg=MyBlue,pch=21)
  points(log10(p2sub$baseMeanA+add),log10(p2sub$baseMeanB+add),col=MyTransparent,bg=MyPurple,pch=21)
  points(log10(p3sub$baseMeanA+add),log10(p3sub$baseMeanB+add),col=MyTransparent,bg=MyRed,pch=21)
  tested <- subset(result.OR,id%in%luc$OLFR)
  points(log10(subset(tested,id=="Olfr160")$baseMeanA+add),log10(subset(tested,id=="Olfr160")$baseMeanB+add),col=MyTransparent,bg="cyan",pch=21)
  points(log10(tested$baseMeanA+add),log10(tested$baseMeanB+add),pch=0,lwd=1)  
  legend("bottomright",c("Tested in vitro", "Olfr160 (M72)",as.expression(bquote('FDR adjusted '* italic(p) < .(p1))),as.expression(bquote('FDR adjusted '* italic(p) < .(p3)))),col=c("black",MyGray2,MyGray2,MyGray2),
         pt.bg=c(MyTransparent,"cyan",MyBlue,MyRed),pch=c(0,21,21,21),pt.lwd=c(1,1,1,1),bty="n",cex=0.9)
  
  list <- list(a=subset(luc,OLFR%in%subset(tested,FDR<vivo.cut&logFC>0)$id)$f1,b=subset(luc,OLFR%in%subset(tested,FDR>=vivo.cut|logFC<=0)$id)$f1,
               c=subset(luc,OLFR%in%subset(tested,FDR<vivo.cut&logFC>0)$id)$f2,d=subset(luc,OLFR%in%subset(tested,FDR>=vivo.cut|logFC<=0)$id)$f2,
               e=subset(luc,OLFR%in%subset(tested,FDR<vivo.cut&logFC>0)$id)$f3,f=subset(luc,OLFR%in%subset(tested,FDR>=vivo.cut|logFC<=0)$id)$f3)
  q <- data.frame(lapply(list,quantile,probs=c(0.05,0.25,0.5,0.75,0.95)))
  boxplot(list,col=c("red","gray","red","gray","red","gray"),border=MyTransparent,ylab=NA,names=NA,ylim=c(-5,120),frame.plot=F,xaxt="n")
  title(ylab="Normalized in vitro activation",cex.lab=1.5)
  axis(1,at=c(1.5,3.5,5.5),labels=c(expression(3*mu*"M"),expression(30*mu*"M"),expression(300*mu*"M")),tick=F,cex.axis=1.5)
  for (j in 1:6) {
    arrows(j,q[2,j],j,q[1,j],angle=90,col=ifelse(j%%2==0,"gray","red"),lwd=3);arrows(j,q[4,j],j,q[5,j],angle=90,col=ifelse(j%%2==0,"gray","red"),lwd=3)
    segments(j-0.4,q[3,j],j+0.4,q[3,j],lwd=3)
  }
  #  legend("topleft",c(paste("FDR adjusted p<",vivo.cut," (n=",dim(subset(luc,OLFR%in%subset(tested,FDR<vivo.cut&logFC>0)$id))[1],")", sep = ""),paste("FDR adjusted p>=",vivo.cut," (n=",dim(subset(luc,OLFR%in%subset(tested,FDR>=vivo.cut|logFC<=0)$id))[1],")", sep = "")),fill=c("red","gray"),border=MyTransparent,bty="n")
  legend("topleft",c(paste("Enriched (n = ",dim(subset(luc,OLFR%in%subset(tested,FDR<vivo.cut&logFC>0)$id))[1],")", sep = ""),paste("Not enriched (n = ",dim(subset(luc,OLFR%in%subset(tested,FDR>=vivo.cut|logFC<=0)$id))[1],")", sep = "")),fill=c("red","gray"),border=MyTransparent,bty="n",cex=1.2)
  #  text(1.2,84,paste("n:",dim(subset(luc,OLFR%in%subset(tested,FDR<vivo.cut&logFC>0)$id))[1],"vs",dim(subset(luc,OLFR%in%subset(tested,FDR>=vivo.cut|logFC<=0)$id))[1]))
  
  
  signif.code <- function (p) {if (p<=0.01) {return ("***")}; if (p>0.01&p<=0.05) {return ("*")}; if (p>0.05) {return ("N.S.")}}
  segments(1,max(q[5,1],q[5,2])+3,2,max(q[5,1],q[5,2])+3,lwd=2)
  text(1.5,max(q[5,1],q[5,2])+7,signif.code(wilcox.test(list$a,list$b)$p.value),cex=ifelse(signif.code(wilcox.test(list$a,list$b)$p.value)=="N.S.",1,2))
  segments(3,max(q[5,3],q[5,4])+3,4,max(q[5,3],q[5,4])+3,lwd=2)
  text(3.5,max(q[5,3],q[5,4])+7,signif.code(wilcox.test(list$c,list$d)$p.value),cex=ifelse(signif.code(wilcox.test(list$c,list$d)$p.value)=="N.S.",1,2))
  segments(5,max(q[5,5],q[5,6])+3,6,max(q[5,5],q[5,6])+3,lwd=2)
  text(5.5,max(q[5,5],q[5,6])+7,signif.code(wilcox.test(list$e,list$f)$p.value),cex=ifelse(signif.code(wilcox.test(list$e,list$f)$p.value)=="N.S.",1,2))
  
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  text(5,9,paste("FDR adjusted p:",vivo.cut,"; #ORs enriched:",dim(subset(result.OR,FDR<vivo.cut&logFC>0))[1]))
  text(5,6,paste("#ORs tested in vitro:",dim(tested)[1]))
  text(5,4,paste("wilcox n: red:",dim(subset(luc,OLFR%in%subset(tested,FDR<vivo.cut&logFC>0)$id))[1],"; gray:",dim(subset(luc,OLFR%in%subset(tested,FDR>=vivo.cut|logFC<=0)$id))[1]))
  text(5,3,paste("3uM wilcox p-value:",wilcox.test(list$a,list$b)$p.value));text(5,2,paste("30uM wilcox p-value:",wilcox.test(list$c,list$d)$p.value));text(5,1,paste("300uM wilcox p-value:",wilcox.test(list$e,list$f)$p.value))
}
###
pdf("figure3.pdf",5.2,5)
par(mar=c(4,4.4,2,2)+0.1)
Box(result.No_Ac.OR,agg.all.luc.dedup.clean,0.05)
Box(result.No_1pAc.OR,agg.all.luc.dedup.clean,0.05)
dev.off()
################################################## ROC curves, returns receptors confirmed in vitro
ROC <- function(result.OR,luc,vivo.cut,vitro.p) {
  set.seed(1);mixmodl.luc <- normalmixEM(luc$f3,mu=c(0,50));plot(mixmodl.luc,which=2,breaks=100)
  vitro.cut <- qnorm(p=(1-vitro.p), mixmodl.luc$mu[1],mixmodl.luc$sigma[1])
  hist(luc$f3,breaks=100,freq=F,xlab="Normalized fold of induction",ylab="density",main="Fold of induction, 300uM Acetophenone")
  abline(v=vitro.cut,col="blue")
  text(100,0.2,paste("in vitro cut-off:",round(vitro.cut,3)))
  acr <- subset(luc,f3>vitro.cut)
  no.acr <- subset(luc,f3<=vitro.cut)
  # ROC enrichment classifier
  tested <- subset(result.OR,id%in%luc$OLFR);tested <- tested[order(tested$FDR,decreasing=F),]
  tp <- fp <- 0; TP <- FP <- DIS <- NULL
  for (i in 1:dim(tested)[1]) {
    if (tested[i,]$id %in% acr$OLFR) {tp <- tp + 1};if (tested[i,]$id %in% no.acr$OLFR) {fp <- fp + 1}
    TP <- c(TP,tp/dim(subset(tested,id%in%acr$OLFR))[1]); FP <- c(FP,fp/dim(subset(tested,id%in%no.acr$OLFR))[1])
    DIS <- c(DIS,tp/length(acr$OLFR)-fp/length(no.acr$OLFR))
  }
  plot(FP,TP,type="l",col="blue",lwd=3,xlab=NA,ylab=NA,frame.plot=F)
  title(xlab="False positive rate",ylab="True positive rate",cex.lab=1.5)
  segments(0,0,1,1,lwd=3)#;abline(v=FP[which.min(abs(tested$FDR-vivo.cut))],lty=2)
  #  legend("bottomright",c("in vivo enrichment classifier","random",paste("in vivo cut-off ( p<",vivo.cut,")")),col=c("blue","black","black"),lwd=c(3,3,1),lty=c(1,1,2),bty="n",cex=0.9)
  #legend("bottomright",c(expression(italic('in vivo')*' enrichment classifier'),"Random"),col=c("blue","black"),lwd=c(3,3),lty=c(1,1),bty="n",cex=1.2)
  legend("bottomright",c("in vivo enrichment classifier","Random"),col=c("blue","black"),lwd=c(3,3),lty=c(1,1),bty="n",cex=1.2)
  DIS <- TP; AUC <- 0; for (i in 2:dim(tested)[1]) {AUC <- AUC + DIS[i]*(FP[i]-FP[i-1])}
  text(0.2,1,paste("AUC:",round(AUC,3)),cex=1.3)
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  text(5,4,paste("AUC:",round(AUC,3)))
  ANS <- subset(tested,id%in%luc$OLFR,select=c("id","FDR"))
  ANS$vitro <- (ANS$id%in%acr$OLFR)
  text(5,5,paste("p=",wilcox.test(FDR~vitro,data=ANS,alternative="greater")$p.value))
  # ROC fold of induction classifier
  tested.luc <- subset(luc,OLFR%in%result.OR$id);tested.luc <- tested.luc[order(tested.luc$f3,decreasing=T),]
  tp <- fp <- 0;TP <- FP <- DIS <- NULL
  for (i in 1:dim(tested.luc)[1]) {
    if (tested.luc[i,]$OLFR %in% subset(result.OR,FDR<vivo.cut&logFC>0)$id) {tp <- tp + 1}
    if (tested.luc[i,]$OLFR %in% subset(result.OR,FDR>=vivo.cut|logFC<=0)$id) {fp <- fp + 1}
    TP <- c(TP,tp/dim(subset(tested.luc,OLFR%in%subset(result.OR,FDR<vivo.cut&logFC>0)$id))[1])
    FP <- c(FP,fp/dim(subset(tested.luc,OLFR%in%subset(result.OR,FDR>=vivo.cut|logFC<=0)$id))[1])
  }
  plot(FP,TP,type="l",col="blue",lwd=3,xlab=NA,ylab=NA,frame.plot=F)
  title(xlab="False positive rate",ylab="True positive rate",cex.lab=1.5)
  segments(0,0,1,1,lwd=3)
  #  abline(v=FP[which.min(abs(tested.luc$f3-vitro.cut))],lty=2)
  #  legend("bottomright",c("in vitro response classifier","random","in vitro cut-off"),col=c("blue","black","black"),lwd=c(3,3,1),lty=c(1,1,2),bty="n",cex=0.9)
  legend("bottomright",c("in vitro response classifier","Random"),col=c("blue","black"),lwd=c(3,3),lty=c(1,1),bty="n",cex=1.2)
  DIS <- TP; AUC <- 0; for (i in 2:dim(tested.luc)[1]) {AUC <- AUC + DIS[i]*(FP[i]-FP[i-1])}
  text(0.2,1,paste("AUC:",round(AUC,3)),cex=1.3)
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  ANS2 <- subset(tested.luc,OLFR%in%result.OR$id,select=c("OLFR","f3"))
  ANS2$vivo <- (ANS2$OLFR%in%subset(result.OR,FDR<vivo.cut&logFC>0)$id)
  text(5,5,paste("p=",wilcox.test(f3~vivo,data=ANS2)$p.value))
  text(5,4,paste("AUC:",round(AUC,3)))
  list<-list(activated=exp(subset(result.OR,id%in%subset(luc,f3>vitro.cut)$OLFR)$logFC),not=exp(subset(result.OR,id%in%subset(luc,f3<=vitro.cut)$OLFR)$logFC))
  q <- data.frame(lapply(list,quantile,probs=c(0.05,0.25,0.5,0.75,0.95),na.rm=T))
  #  par(mfrow=c(1,2))
  boxplot(list,col=c("red","gray"),border=MyTransparent,ylab=NA,ylim=c(0,max(q)*1.5),xlim=c(0.5,5),frame.plot=F,names=NA,cex=0.6,xaxt="n")
  title(ylab="Fold enrichment in vivo",cex.lab=1.5)
  for (j in 1:2) {
    arrows(j,q[2,j],j,q[1,j],angle=90,col=ifelse(j%%2==0,"gray","red"),lwd=3);arrows(j,q[4,j],j,q[5,j],angle=90,col=ifelse(j%%2==0,"gray","red"),lwd=3)
    segments(j-0.4,q[3,j],j+0.4,q[3,j],lwd=3)
  }
  legend("topleft",c(paste("Activated in vitro (n = ",dim(subset(luc,f3>vitro.cut&OLFR%in%result.OR$id))[1],")", sep = ""),paste("Not activated in vitro (n = ",dim(subset(luc,f3<=vitro.cut&OLFR%in%result.OR$id))[1],")", sep = "")),fill=c("red","gray"),border=MyTransparent,bty="n",cex=1.2)
  signif.code <- function (p) {if (p<=0.01) {return ("***")}; if (p>0.01&p<=0.05) {return ("*")}; if (p>0.05) {return ("N.S.")}}
  segments(1,max(q[5,1],q[5,2])+max(q)*0.06,2,max(q[5,1],q[5,2])+max(q)*0.06,lwd=2)
  text(1.5,max(q[5,1],q[5,2])+max(q)*0.12,signif.code(wilcox.test(list$activated,list$not)$p.value),cex=ifelse(signif.code(wilcox.test(list$activated,list$not)$p.value)=="N.S.",1,2))
  
  #  par(mfrow=c(1,1))
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  text(5,5,paste("activated vs not\n, p=",wilcox.test(subset(result.OR,id%in%subset(luc,f3>vitro.cut)$OLFR)$logFC,subset(result.OR,id%in%subset(luc,f3<=vitro.cut)$OLFR)$logFC,alternative="greater")$p.value))
  enriched.confirmed <- intersect(subset(tested,FDR<vivo.cut)$id,subset(luc,f3>vitro.cut)$OLFR)
  return(enriched.confirmed)
  #return(intersect(tested$id,subset(luc,f3>vitro.cut)$OLFR))
}
###
pdf("Ac_ROC.pdf",5,5)
ROC(result.No_Ac.OR,agg.all.luc.dedup.clean,vivo.cut=0.05,vitro.p=0.01)
ROC(result.No_1pAc.OR,agg.all.luc.dedup.clean,vivo.cut=0.05,vitro.p=0.01)
dev.off()

######### plot ACR doses
respond091914
dose.data <- read.table("100614_AcDose.data",header=F)
dose.name <- scan("100614_AcDose.name",what="character")

dose.name.clean <- dose.name[which(dose.name%in%respond091914)]
dose.data.clean <- dose.data[,as.vector(rbind(which(dose.name%in%respond091914)*3-2,which(dose.name%in%respond091914)*3-1,which(dose.name%in%respond091914)*3))]

means <- NULL
for (i in 1:length(dose.name.clean)) {
  means <- c(means,rowMeans(dose.data.clean[8,c(i*3-2,i*3-1,i*3)]))
}
ordered.dose.name <- dose.name.clean[order(means,decreasing=T)]
ordered.dose.data <- dose.data.clean[,as.vector(rbind(order(means,decreasing=T)*3-2,order(means,decreasing=T)*3-1,order(means,decreasing=T)*3))]

ordered.dose.name <- c(ordered.dose.name,"pCI")
ordered.dose.data <- cbind(ordered.dose.data,dose.data[,166:168])


pdf("ACRdose.pdf",5,5)
par(mfrow=c(1,1))
colors <- rep(c("red","orange","green","blue","purple"),11)
colors[ordered.dose.name=="pCI"] <- "gray"
for (i in 1:length(ordered.dose.name)) {
  mod <- i%%5
  sub.dose <- list(ordered.dose.name[i],ordered.dose.data[,(i*3-2):(i*3)])
  dr <- drm(as.vector(as.matrix(sub.dose[[2]]))~rep(c(-30,-9,-8,-7,-6,-5,-4,-3),3),fct = LL2.3(),logDose=10,type="continuous")
  if (mod==1) {
    plot(dr,broken=T,bp=-10,col=colors[i],lwd=2,xlab=NA,ylab=NA,pch=16,ylim=c(min(sub.dose[[2]])*0.9,max(sub.dose[[2]])*1.1),xttrim=F,bty="n",cex=1.5)
    title(ylab="Response",xlab=expression('Log'[10]*'(acetophenone)[M]'),cex.lab=1.5)
    arrows(c(-10,-9,-8,-7,-6,-5,-4,-3),rowMeans(sub.dose[[2]])+apply(sub.dose[[2]],1,std),c(-10,-9,-8,-7,-6,-5,-4,-3),rowMeans(sub.dose[[2]])-apply(sub.dose[[2]],1,std),col=colors[i],angle=90,lwd=2,code=3,length=0.05)
    legend("topleft",ordered.dose.name[i:min(length(ordered.dose.name),(i+4))],pch=16,col=colors[i:min(length(ordered.dose.name),(i+4))],bty="n",cex=1.2)
  }
  if (mod!=1) {
    plot(dr,broken=T,bp=-10,col=colors[i],lwd=2,pch=16,add=T,cex=1.5)
    arrows(c(-10,-9,-8,-7,-6,-5,-4,-3),rowMeans(sub.dose[[2]])+apply(sub.dose[[2]],1,std),c(-10,-9,-8,-7,-6,-5,-4,-3),rowMeans(sub.dose[[2]])-apply(sub.dose[[2]],1,std),col=colors[i],angle=90,lwd=2,code=3,length=0.05) 
  }
}
dev.off()

######## sensititivy in vivo and in vitro
library(drc)
ec50=NULL
maxres=NULL
lower<-upper<-ed<-NULL
for (i in 1:length(ordered.dose.name)) {
  sub.dose <- list(ordered.dose.name[i],ordered.dose.data[,(i*3-2):(i*3)])
  dr <- drm(as.vector(as.matrix(sub.dose[[2]]))~rep(c(-15,-9,-8,-7,-6,-5,-4,-3),3),fct = LL.4(fixed=c(NA,NA,NA,NA)),logDose=10)
  #  plot(dr,broken=T,bp=-10,col="red",lwd=2,xlab="log(acetophenone)[M]",ylab="Response",pch=16,main=sub.dose[[1]],bty="n")
  #  arrows(c(-10,-9,-8,-7,-6,-5,-4,-3),rowMeans(sub.dose[[2]])+apply(sub.dose[[2]],1,std),c(-10,-9,-8,-7,-6,-5,-4,-3),rowMeans(sub.dose[[2]])-apply(sub.dose[[2]],1,std),col="red",angle=90,lwd=2,code=3,length=0.05)
  ec50 <- c(ec50,dr$coefficients[4])
  maxres <- c(maxres,dr$coefficients[3]-dr$coefficients[2])
  ed <- ED(dr,50,interval="delta"); lower <- c(lower, ed[2]); upper <- c(upper, ed[3]) 
}
#ec50 for Olfr736 is off because of high concentration. 
#correct for this individually
i=which(ordered.dose.name=="Olfr736")
sub.dose <- list(ordered.dose.name[i],ordered.dose.data[,(i*3-2):(i*3)])
dr <- drm(as.vector(as.matrix(sub.dose[[2]][-8,]))~rep(c(-15,-9,-8,-7,-6,-5,-4),3),fct = LL.4(fixed=c(NA,NA,NA,NA)),logDose=10)
ec50[i] <- dr$coefficients[4]
ec501p <- subset(result.No_1pAc.OR,id%in%ordered.dose.name)
ec50100p <- subset(result.No_Ac.OR,id%in%ordered.dose.name)
ec50.vivo <- merge(ec501p,ec50100p,by = "id")#x is 1%, y is 100%
ec50.vivo$ec50 <- 0
for (i in 1:dim(ec50.vivo)[1]) {ec50.vivo$ec50[i] <- ec50[which(ordered.dose.name==ec50.vivo$id[i])]}
ec50.vivo$maxres <- 0
for (i in 1:dim(ec50.vivo)[1]) {ec50.vivo$maxres[i] <- maxres[which(ordered.dose.name==ec50.vivo$id[i])]}
ec50.vivo$lower <- 0
for (i in 1:dim(ec50.vivo)[1]) {ec50.vivo$lower[i] <- lower[which(ordered.dose.name==ec50.vivo$id[i])]}
ec50.vivo$upper <- 0
for (i in 1:dim(ec50.vivo)[1]) {ec50.vivo$upper[i] <- upper[which(ordered.dose.name==ec50.vivo$id[i])]}

pdf("EC50.pdf",5.5,5.5)
par(mar=c(4,4.4,2,2)+0.1)
plot(ec50.vivo$logFC.x,log10(ec50.vivo$ec50),pch=16,xlab=NA,ylab=NA,frame.plot = F)
title(xlab=expression('Log'[2]*'(fold change, 1% acetophenone)'), ylab=expression('Log'[10]*'(EC50)') ,cex.lab=1.5)
abline(lm(log10(ec50.vivo$ec50)~ec50.vivo$logFC.x))

plot(ec50.vivo$logFC.y,log10(ec50.vivo$ec50),pch=16,xlab=NA,ylab=NA,frame.plot = F)
title(xlab=expression('Log'[2]*'(fold change, 100% acetophenone)'), ylab=expression('Log'[10]*'(EC50)') ,cex.lab=1.5)
abline(lm(log10(ec50.vivo$ec50)~ec50.vivo$logFC.y))

library(beanplot)
beanplot(log10(subset(ec50.vivo,FDR.x<0.05)$ec50),log10(subset(ec50.vivo,FDR.x>=0.05)$ec50),side="both",col=list("red","orange"),show.names=F,frame.plot=F,ylab=NA)
title(ylab=expression('Log'[10]*'(EC50)') ,cex.lab=1.5)
plot.new();legend("bottomright",c(expression(paste("Enriched by 1% acetophenone (",italic(p)<0.05,")")),expression(paste("Additionally enriched by 100% acetophenone (",italic(p)<0.05,")"))),fill=c("red","orange"),bty="n",cex=0.8)

beanplot(log10(subset(ec50.vivo,FDR.x<0.001)$ec50),log10(subset(ec50.vivo,FDR.x>=0.001&FDR.y<=0.001)$ec50),side="both",col=list("red","orange"),show.names=F,frame.plot=F,ylab=NA)
title(ylab=expression('Log'[10]*'(EC50)') ,cex.lab=1.5)
plot.new();legend("bottomright",c(expression(paste("Enriched by 1% acetophenone (",italic(p)<0.001,")")),expression(paste("Additionally enriched by 100% acetophenone (",italic(p)<0.001,")"))),fill=c("red","orange"),bty="n",cex=0.8)
dev.off()

summary(lm(log10(ec50.vivo$ec50)~ec50.vivo$logFC.x))
summary(lm(log10(ec50.vivo$ec50)~ec50.vivo$logFC.y))
wilcox.test(log(subset(ec50.vivo,FDR.x<0.05)$ec50),log(subset(ec50.vivo,FDR.x>=0.05&FDR.y<0.05)$ec50)) #&id!="Olfr736"


#use all ec50 values to generate contingency table
all.luc.dose.data <- all.luc.dedup[,2:13]
all.luc.dose.name <- as.vector(all.luc.dedup[,1])
all.luc.dose.name[all.luc.dose.name=="Olfr175-ps1"] <- "Olfr174"
#write.table((all.luc.dose.data),"all_luc_dose_data.txt",sep="\t")
library(drc)
ec50=NULL
maxres=NULL
lower<-upper<-ed<-NULL
for (i in 1:length(all.luc.dose.name)) {
  sub.dose <- as.numeric(all.luc.dose.data[i,])
  dr <- tryCatch(drm(sub.dose~c(-15,-15,-15,-5.52,-5.52,-5.52,-4.52,-4.52,-4.52,-3.52,-3.52,-3.52),fct = LL.4(fixed=c(NA,NA,NA,NA)),logDose=10),error=function(err){return("failed")})
  #  plot(dr,broken=T,bp=-10,col="red",lwd=2,xlab="log(acetophenone)[M]",ylab="Response",pch=16,main=sub.dose[[1]],bty="n")
  #  arrows(c(-10,-9,-8,-7,-6,-5,-4,-3),rowMeans(sub.dose[[2]])+apply(sub.dose[[2]],1,std),c(-10,-9,-8,-7,-6,-5,-4,-3),rowMeans(sub.dose[[2]])-apply(sub.dose[[2]],1,std),col="red",angle=90,lwd=2,code=3,length=0.05)
  if (dr=="failed") {
    ec50 <- c(ec50,100)
    maxres <- c(maxres,0)
    lower <- c(lower, 100); upper <- c(upper, 100)     
  }
  if (dr!= "failed"){
    ec50 <- c(ec50,dr$coefficients[4])
    maxres <- c(maxres,dr$coefficients[3]-dr$coefficients[2])
    ed <- ED(dr,50,interval="delta"); lower <- c(lower, ed[2]); upper <- c(upper, ed[3]) 
  }
}
#ec50 substitute with more acurate ones when possible
ec50.new <- ec50
rec = NULL
for (i in 1:length(ec50.new)) {
  if (all.luc.dose.name[i]%in%ec50.vivo$id) {
    ec50.new[i] <- ec50.vivo$ec50[which(ec50.vivo$id==all.luc.dose.name[i])] 
    rec <- c(rec,i)
  }
}
plot(ec50[rec],ec50.new[rec],log="xy")
abline(0,1)



#use a modified ROC function! the return is changed to in vitro positive only

names <- as.vector(agg.all.luc.dedup.clean$OLFR)
names[which(names=="Olfr175-ps1")] <- "Olfr174"
agg.all.luc.dedup.clean$OLFR <- names
responded.invitro <- ROC(result.No_Ac.OR,agg.all.luc.dedup.clean,vivo.cut=0.05,vitro.p=0.01)

a=-4
b=-2
sensitive.names <- intersect(all.luc.dose.name[which(log10(ec50.new)<=a)],responded.invitro)
moderate.names <- intersect(all.luc.dose.name[which(log10(ec50.new)>a&log10(ec50.new)<=b)],responded.invitro)
low.names <- setdiff(intersect(agg.all.luc.dedup.clean$OLFR,result.No_Ac$id),c(sensitive.names,moderate.names))

c1=length(intersect(setdiff(result.No_Ac.OR$id,c(as.character(subset(result.No_Ac.OR,FDR<0.05&logFC>0)$id),as.character(subset(result.No_1pAc.OR,FDR<0.05&logFC>0)$id))),sensitive.names))
c2=length(intersect(setdiff(result.No_Ac.OR$id,c(as.character(subset(result.No_Ac.OR,FDR<0.05&logFC>0)$id),as.character(subset(result.No_1pAc.OR,FDR<0.05&logFC>0)$id))),moderate.names))
c3=length(intersect(setdiff(result.No_Ac.OR$id,c(as.character(subset(result.No_Ac.OR,FDR<0.05&logFC>0)$id),as.character(subset(result.No_1pAc.OR,FDR<0.05&logFC>0)$id))),low.names))

l1=length(setdiff(intersect(subset(result.No_Ac.OR,FDR<0.05&logFC>0)$id,sensitive.names),subset(result.No_1pAc.OR,FDR<0.05&logFC>0)$id))
l2=length(setdiff(intersect(subset(result.No_Ac.OR,FDR<0.05&logFC>0)$id,moderate.names),subset(result.No_1pAc.OR,FDR<0.05&logFC>0)$id))
l3=length(setdiff(intersect(subset(result.No_Ac.OR,FDR<0.05&logFC>0)$id,low.names),subset(result.No_1pAc.OR,FDR<0.05&logFC>0)$id))

h1=length(intersect(subset(result.No_1pAc.OR,FDR<0.05&logFC>0)$id,sensitive.names))
h2=length(intersect(subset(result.No_1pAc.OR,FDR<0.05&logFC>0)$id,moderate.names))
h3=length(intersect(subset(result.No_1pAc.OR,FDR<0.05&logFC>0)$id,low.names))

c(c1,c2,c3,l1,l2,l3,h1,h2,h3)
mtx <- matrix(c(c1,c2,c3,l1,l2,l3,h1,h2,h3),nrow=3)
colnames(mtx) <- c("not enriched","100% only","1%")
rownames(mtx) <- c("sensitive","moderate","low")
#write.table(mtx,file="group.txt",sep="\t",quote=F)
#par(mar=c(5.1,4.1,4.1,2.1)
pdf("barfigure.pdf",5,5)
barplot(mtx/rbind(colSums(mtx),colSums(mtx),colSums(mtx)),border=NA,cex.names=1.3,yaxt="n",col=c("red","orange","gray"),space=1,ylim=c(0,1.1))
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),lab=paste0(c(0,20,40,60,80,100), " %"),las=2,cex.axis=1.2)
text(x=c(1.5,3.5,5.5),y=rep(1.05,3),labels=c(paste("n =",(c1+c2+c3)),paste("n =",(l1+l2+l3)),paste("n =",(h1+h2+h3))),cex=1.2)
#mtext("Percentage of ORs", side=2, line=3, cex=1.5)
#text(x = c(0.8,1.8,3.2), par("usr")[3] - 0.2, labels = c("not enriched","100% only","1%"), srt = 45, pos = 1, xpd = TRUE)
plot.new();legend("topright",c(expression(log[10]*EC50<=-4),expression(log[10]*EC50<=-2),expression(log[10]*EC50>-2~"or no response")),fill=c("red","orange","gray"),bty="n",border=NA)
dev.off()
chisq.test(mtx)
chisq.test(mtx[,-1])

ec50.out <- data.frame(id=all.luc.dose.name,log10EC50=log10(ec50.new))
ec50.out$log10EC50[ec50.out$log10EC50>=2] <- "NA"
ec50.out <- subset(ec50.out,id%in%result.No_Ac.OR$id)
write.table(subset(ec50.out,id%in%result.No_Ac.OR$id),"ec50.xls",quote=F,sep="\t")
write.table(ec50.vivo,"ec50vivo.xls",,quote=F,sep="\t")

########modeling, PCA, SVM and elastic-net
###PCA
dim(All.OR)
AA2prop <- function(aa,AAprop) {
  return(as.numeric(subset(AAprop,AA==aa,select=c(c,p,v))))
}

Prop.OR.raw <- NULL 
for (i in 1:dim(All.OR)[1]) {
  row <- NULL
  for (j in (1:dim(All.OR)[2])) {
    row <- c(row,AA2prop(All.OR[i,j],AAprop))
  }
  Prop.OR.raw <- rbind(Prop.OR.raw,row)
}
#Prop.OR.raw <- Prop.OR
colNam <- NULL
for (i in 1:dim(All.OR)[2]) {colNam <- c(colNam,paste(i,"c",sep="_"),paste(i,"p",sep="_"),paste(i,"v",sep="_"))}
keep.ind <- (colMeans(All.OR=="-") < 0.1)
Prop.OR.raw <- as.data.frame(Prop.OR.raw)
colnames(Prop.OR.raw) <- colNam
row.names(Prop.OR.raw) <- OR.aln$nam

Prop.OR <- Prop.OR.raw[,as.vector(sapply(keep.ind, function (x) rep(x,3)))]
for (i in 1:dim(Prop.OR)[2]) {Prop.OR[is.na(Prop.OR[,i]),i] <- mean(Prop.OR[,i],na.rm=T)}
var.ind <- apply(Prop.OR,2,var)!=0
Prop.OR <- Prop.OR[var.ind]
#OR.pca <- prcomp(~ ., data = prop.OR, na.action = na.omit)
OR.pca <- prcomp(Prop.OR,scale.=T)
plot(OR.pca$x[,1],OR.pca$x[,2],pch=3,cex=0.7)

no.R <- function(result.OR,luc,vivo.cut,vitro.p) {
  set.seed(1);mixmodl.luc <- normalmixEM(luc$f3,mu=c(0,50));plot(mixmodl.luc,which=2,breaks=100)
  vitro.cut <- qnorm(p=(1-vitro.p), mixmodl.luc$mu[1],mixmodl.luc$sigma[1])
  no.acr <- intersect(subset(luc,f3<=vitro.cut)$OLFR,subset(result.OR,FDR>=vivo.cut|logFC<0)$id)
  return(no.acr)
}

no.ACR <- no.R(result.No_Ac.OR,agg.all.luc.dedup.clean,vivo.cut=0.05,vitro.p=0.05)
#no.ACR <- intersect(no.ACR, subset(result.No_Ac.OR,logFC<=0)$id)
no.ACR.1p <- no.R(result.No_1pAc.OR,agg.all.luc.dedup.clean,vivo.cut=0.05,vitro.p=0.05)
#no.ACR.1p <- intersect(no.ACR.1p,subset(result.No_1pAc.OR,logFC<=0)$id)
neither.ACR <- intersect(no.ACR,no.ACR.1p)

ACR.ind <- which(OR.aln$nam%in%respond091914)
notACR.ind <- which(OR.aln$nam%in%neither.ACR)

library(rgl)
col=rep("gray",length(OR.pca$x[,1]))
col[notACR.ind] <- "magenta"
col[ACR.ind] <- "cyan"
plot3d(OR.pca$x[,1],OR.pca$x[,2],OR.pca$x[,3],col=col,pch=16,xlab="PC1",ylab="PC2",zlab="PC3",size=5)
rgl.postscript("pca.pdf","pdf")
pdf("legend.pdf",5,5)
#plot.new();legend("bottomright",pch=16,col=c("cyan","magenta","gray"),c("Activated","Not activated","Not tested"),bty = "n")
plot.new();legend("bottomright",pch=16,col=c("cyan","magenta","gray"),c("Acetophenone ORs","Not acetophenone ORs","Undetermined"),bty = "n")

dev.off()

###SVM
data <- Modata
for (cost in c(0.01,0.1,1,10,100)){
  for (gamma in c(0.001,0.01,0.1)){  
    
    cost = 10
    gamma=0.001
    fold=10
    cycle=100
    test.size=round(dim(data)[1]/fold)
    #BIG.TP<-BIG.FP<-AUC<-NULL
    ANS <- NULL
    set.seed(0)
    for (j in 1:cycle) { 
      repeat{
        test.ind <- sample(1:dim(data)[1],test.size)
        test.data <- data[test.ind,]
        if (var(test.data$res)!=0) {break}
      }
      train.data <- data[-test.ind,]
      col.ind <- apply(train.data,2,var)!=0
      col.ind[1] <- TRUE
      train.data <- train.data[,col.ind]
      test.data <- test.data[,col.ind]
      svm.model <- svm(as.factor(res) ~ ., data=train.data,kernel="radial",cost = cost, gamma = gamma,scale=T,probability=TRUE)
      #  svm.model <- svm(as.factor(res) ~ ., data=train.data,kernel="linear",scale=T,probability=TRUE)
      svm.pred <- predict(svm.model,test.data,na.action = na.exclude, probability=TRUE,decision.values = TRUE)
      ans <- cbind(attr(svm.pred, "decision.values"),test.data$res,attr(svm.pred, "probabilities")[,1])
      ANS <- rbind(ANS,ans)
    }
    ROC.gen(ANS)
    text(5,3,paste("cost:",cost))
    text(5,2,paste("gamma:",gamma))
  }
}
cros.svm.c10g0001 <- ROC.gen(ANS)
svm.model <- svm(as.factor(res) ~ ., data=Modata[apply(Modata,2,var)!=0],kernel="radial",cost=10,gamma=0.001,scale=T,probability=TRUE)
svm.pred <- predict(svm.model,Hudata[apply(Modata,2,var)!=0][-26,],na.action = na.exclude, probability=TRUE,decision.values = TRUE)
#svm.obj <- tune(svm, as.factor(res) ~ ., data=Modata[apply(Modata,2,var)!=0],gamma = 2^(-1:1), cost = 2^(2:4))
pred <- attr(svm.pred, "decision.values")[,1]
ans <- cbind(pred,Hudata$res[-26])
ext.svm.c10g0001 <- ROC.gen(ans)

#elastic-net
## pick alpha
library(glmnet)
set.seed(0)
foldid=sample(1:10,size=length(Modata$res),replace=TRUE)


size=100
alpha <- seq(0,1,1/size)
elnet <- list()
for ( i in 1:length(alpha)) {
  elnet[[i]] <- cv.glmnet(y=as.factor(Modata$res), x=as.matrix(Modata[,-1]),family="binomial",foldid=foldid,alpha=alpha[i])
}

Parameters <- NULL
for (i in 1:length(alpha)) {
  Parameters <- rbind(Parameters,c(alpha[i],elnet[[i]]$lambda.min,min(elnet[[i]]$cvm),elnet[[i]]$cvsd[which.min(elnet[[i]]$cvm)]))
}
colnames(Parameters) <- c("alpha","lambda.min","cvm.min","cvsd")

alpha <- Parameters[which.min(Parameters[,3]),1]

##predict human
set.seed(0)
glmnet.model <- cv.glmnet(y=as.factor(Modata$res), x=as.matrix(Modata[,-1]),family="binomial",alpha=alpha)
glmnet.model$lambda.1se
row.names(coef(glmnet.model,s="lambda.1se"))[which(coef(glmnet.model,s="lambda.1se")!=0)]
coef(glmnet.model,s="lambda.1se")
sum(coef(glmnet.model,s="lambda.1se")!=0)
Hu.elnet.pred <- predict(glmnet.model,as.matrix(Hudata[-26,-1]),type="response",s="lambda.1se")
Hu.ans <- cbind(Hu.elnet.pred,Hudata$res[-26])
ext.elnet <- ROC.gen(Hu.ans)
wilcox.test(Hu.ans[,1]~Hu.ans[,2],alternative="less")
#glm(res~X119_c, data=data, family=binomial(logit))
#Hu.ROC <- roc(Hu.ans[,2],Hu.ans[,1])
#roc.test(Hu.ROC.sim,Hu.ROC)

#fla.small <- paste("factor(res) ~", paste(row.names(coef(glmnet.model,s="lambda.1se"))[which(coef(glmnet.model,s="lambda.1se")!=0)][-1], collapse="+"))

### cross val
fold=10
cycle=100
test.size=round(dim(Modata)[1]/fold)
#BIG.TP<-BIG.FP<-AUC<-NULL
ANS4 <- NULL
set.seed(0)
for (j in 1:cycle) { 
  repeat{
    test.ind <- sample(1:dim(Modata)[1],test.size)
    test.data <- Modata[test.ind,]
    if (var(test.data$res)!=0) {break}
  }
  train.data <- Modata[-test.ind,]
  elnet.model <- cv.glmnet(y=as.factor(train.data$res), x=as.matrix(train.data[,-1]),family="binomial",alpha=alpha)
  #elnet.model <- glm(as.formula(fla.small), data=train.data, family=binomial(logit))
  #elnet.model <- glm(y=as.factor(train.data$res), x=as.matrix(train.data[,-1]),family="binomial",alpha=alpha,lambda=lambda)
  elnet.pred <- predict(elnet.model,as.matrix(test.data[,-1]),type="response")
  ans <- cbind(elnet.pred,test.data$res)
  ANS4 <- rbind(ANS4,ans)
}
cros.elnet.a007 <- ROC.gen(ANS4)
#ROC4 <- roc(ANS4[,2],ANS4[,1])





#similarity
OR.dist <- dist.alignment(OR.aln,matrix="similarity")
OR.dist.mtx <- as.matrix(OR.dist)
similarity.predict <- function(train.OR,train.res,test.OR,OR.dist.mtx) {
  sub.dist.mtx <- OR.dist.mtx[row.names(OR.dist.mtx)%in%train.OR,colnames(OR.dist.mtx)%in%test.OR]
  sub.dist.mtx <- sub.dist.mtx[,match(colnames(sub.dist.mtx),test.OR)]
  most.similar.OR <- row.names(sub.dist.mtx)[apply(sub.dist.mtx,2,which.min)]
  pred.sign <- 2*(train.res[match(most.similar.OR,train.OR)]-0.5)
  pred.val <- 1/apply(sub.dist.mtx,2,min)
  pred <- pred.sign*pred.val
  return(list(pred,sub.dist.mtx,most.similar.OR))
}

fold=10
cycle=100
test.size=round(dim(data)[1]/fold)
#BIG.TP<-BIG.FP<-AUC<-NULL
ANS2 <- NULL
set.seed(1)
for (j in 1:cycle) { 
  repeat{
    test.ind <- sample(1:dim(data)[1],test.size)
    test.data <- data[test.ind,]
    if (var(test.data$res)!=0) {break}
  }
  train.OR <- row.names(data[-test.ind,])
  train.res <- data[-test.ind,1]
  test.OR <- row.names(test.data)
  sim.pred <- similarity.predict(train.OR,train.res,test.OR,OR.dist.mtx)
  #  ans2 <- cbind(sim.pred[[1]],test.data$res)
  ANS2 <- rbind(ANS2,cbind(sim.pred[[1]],test.data$res))#,sim.pred[[3]]))
}
cros.sim <- ROC.gen(ANS2)
MoHudist <- dist.alignment(MoHuOR.aln,matrix="similarity")
MoHudist.mtx <- as.matrix(MoHudist)
ans2 <- cbind(similarity.predict(row.names(Modata),Modata$res,row.names(Hudata[-26,]),MoHudist.mtx)[[1]],Hudata$res[-26])
ext.sim <- ROC.gen(ans2)

pdf("model_ROC.pdf",5,5)
plot(cros.sim[[2]],cros.sim[[3]],type="l",lwd=3,frame.plot=F,xlab=NA,ylab=NA)
title(ylab="True positive rate",xlab="False positive rate",cex.lab=1.5)
points(cros.svm.c10g0001[[2]],cros.svm.c10g0001[[3]],type="l",lwd=3,col="magenta")
points(cros.elnet.a007[[2]],cros.elnet.a007[[3]],type="l",lwd=3,col="green")
legend("bottomright",c("SVM","Elastic net","Similarity"),col=c("magenta","green","black"),lwd=3,bty="n",cex=1.2)

plot(ext.sim[[2]],ext.sim[[3]],type="l",lwd=3,ylab=NA,xlab=NA,frame.plot=F,xlim=c(0,1))
title(ylab="True positive rate",xlab="False positive rate",cex.lab=1.5)
points(ext.svm.c10g0001[[2]],ext.svm.c10g0001[[3]],type="l",lwd=3,col="magenta")
points(ext.elnet[[2]],ext.elnet[[3]],type="l",lwd=3,col="green")
legend("bottomright",c("SVM","Elastic net","Similarity"),col=c("magenta","green","black"),lwd=3,bty="n",cex=1.2)
dev.off()





########## trees
plot.OR.tree <- function(OR.tree,edgecol=MyGray2,edgewid=0.5) { #label2 overrides label1
  #  length <- 1111
  #  plot(OR.tree,cex=1,type="fan",edge.color="azure4",pch=21,show.tip.label=F,tip.color="black",x.lim=c(-0.2,0.4),y.lim=c(-0.41,0.29))
  #  tiplabels(pch=21,cex=0.9,col="azure4",bg="azure3",lwd=0.7)
  plot(OR.tree,cex=1,type="unrooted",edge.color=edgecol,pch=21,show.tip.label=F,tip.color="black",edge.width=edgewid,no.margin=T)#,y.lim=c(-0.41,0.29)) x.lim=c(-0.2,0.4),y.lim=c(0.4,-0.3
  tiplabels(pch=21,cex=0.9,col=MyTransparent,bg=MyTransparent,lwd=0.7)
}

add.OR.tree.label <- function(OR.tree,label,pch=21,pt.col="azure2",bg.col="azure1",lwd=0.7,cex=1,frame="none") {
  length <- length(OR.tree$tip.label)
  pt.color <- rep(MyTransparent,length)
  bg.color <- rep(MyTransparent,length)
  tip.color <- rep(MyTransparent,length)
  for (i in 1:length) {
    if (OR.tree$tip.label[i] %in% label)
    { pt.color[i]<-pt.col;bg.color[i]<-bg.col;tip.color[i]<-bg.col }
  }   
  tiplabels(pch=pch,cex=cex,col=pt.color,bg=bg.color,lwd=lwd)
  #  tiplabels(OR.tree$tip.label,cex=cex,col=tip.color,bg=MyTransparent,frame="none")
}

pdf("trees.pdf",7,7)

edge.100p <- which.edge(OR.tree,intersect(respond091914,subset(result.No_Ac.OR,FDR<0.05&logFC>0)$id))
edge.1p <- which.edge(OR.tree,intersect(respond091914,subset(result.No_1pAc.OR,FDR<0.05&logFC>0)$id))
edge.clr <- rep(MyGray,length(OR.tree$edge))
edge.clr[edge.100p] <- "orange"
edge.clr[edge.1p] <- "red"

edge.wdth <- rep(0.5,length(OR.tree$edge))
edge.wdth[edge.100p] <- 2
edge.wdth[edge.1p] <- 2

plot.OR.tree(OR.tree,edge.clr,edge.wdth)
add.OR.tree.label(OR.tree,c(result.No_1pAc.OR$id,result.No_Ac.OR$id),21,"gray31","azure3",0.5,1)
add.OR.tree.label(OR.tree,intersect(respond091914,subset(result.No_Ac.OR,FDR<0.05&logFC>0)$id),21,"gray31","orange",0.5,1)
add.OR.tree.label(OR.tree,intersect(respond091914,subset(result.No_1pAc.OR,FDR<0.05&logFC>0)$id),21,"gray31","red",0.5,1)
#plot.new();legend("topleft",c("enriched by 1% acetophenone \n and confirmed in vitro","","additionally enriched by 100% acetophenone \n and confirmed in vitro"), pch=c(21,26,21),col=c("gray31",MyTransparent,"gray31"),pt.bg=c("red",MyTransparent,"orange"),bty="n",pt.lwd=0.5,pt.cex=1,cex=1) 
plot.new();legend("topleft",c("Enriched by 1% acetophenone and confirmed in vitro","Additionally enriched by 100% acetophenone and confirmed in vitro"), pch=21,col="gray31",pt.bg=c("red","orange"),bty="n",pt.lwd=0.5,pt.cex=1,cex=1) 

col1 <- "orange"
col2 <- "deepskyblue"
col3 <- "red"
#TMTR092414 <- subset(tmt_092414,RES==1)$OR
edge.Ac <- which.edge(OR.tree,respond091914)
edge.TMT <- which.edge(OR.tree,as.vector(TMTR092414))
edge.AcTMT <- which.edge(OR.tree,intersect(TMTR092414,respond091914))
edge.clr <- rep(MyGray,length(OR.tree$edge))
edge.clr[edge.Ac] <- col1
edge.clr[edge.TMT] <- col2
edge.clr[edge.AcTMT] <- col3

edge.wdth <- rep(0.5,length(OR.tree$edge))
edge.wdth[edge.TMT] <- 2
edge.wdth[edge.Ac] <- 2
plot.OR.tree(OR.tree,edge.clr,edge.wdth)
add.OR.tree.label(OR.tree,c(result.No_1pTMT.OR$id,result.No_TMT.OR$id),21,"gray31","azure3",0.5)
add.OR.tree.label(OR.tree,respond091914,21,"gray31",col1,0.5)
add.OR.tree.label(OR.tree,TMTR092414,21,"gray31",col2,0.5)
add.OR.tree.label(OR.tree,intersect(respond091914,TMTR092414),21,"gray31",col3,0.5)
plot.new();legend("topleft",c("Acetophenone","TMT","Acetophenone and TMT"), pch=21,col="gray31",pt.bg=c(col1,col2,col3),bty="n",pt.lwd=0.5,pt.cex=1,cex=1) 
#plot.new();legend("topleft",c("enriched by 1% TMT\n and confirmed in vitro","","additionally enriched by 100% TMT \n and confirmed in vitro"), pch=c(21,26,21),col=c("gray31",MyTransparent,"gray31"),pt.bg=c("red",MyTransparent,"orange"),bty="n",pt.lwd=0.5,pt.cex=1,cex=1) 
dev.off()



#########
ACR.ind <- which(OR.aln$nam%in%respond091914)
ACR.sub <- NULL
for (i in ACR.ind) {
  ACR.sub <- rbind(ACR.sub, strsplit(OR.aln$seq[[i]],NULL)[[1]])
}
All.OR <- NULL
for (i in 1:length(OR.aln$nam)) {
  All.OR <- rbind(All.OR, strsplit(OR.aln$seq[[i]],NULL)[[1]])
}
quan9.adj <- read.table("quan9adj.txt",header=T)
Hyp.All.TMM <- data.frame(AA.All=Hyp.All.AA,AA.ACR=Hyp.ACR.AA,AA.TMTR=Hyp.TMTR.AA,TMM=OR.TMM.v3,quan.ac=quan9.adj$x,quan.tmt=quan8.adj$x,freq=Hyp.All.AA.freq)
Hyp.All.TMM <- Hyp.All.TMM[Hyp.All.TMM$AA.All!="-",]
#Hyp.All.TMM <- Hyp.All.TMM[keep.ind,]
#Hyp.ACR.TMM <- data.frame(AA=Hyp.ACR.AA,TMM=OR.TMM.v3,quan=quan9.adj$x,All.AA=Hyp.All.AA,All.AA.freq=Hyp.All.AA.freq)
#Hyp.ACR.TMM <- Hyp.ACR.TMM[Hyp.ACR.TMM$AA!="-",]
table(Hyp.All.TMM$TMM)
n <- table(Hyp.All.TMM$TMM)
cord <- rbind(end.curve(x0=20,y0=90,n=n[[1]],direction=1),#EC1
              tmpart(x0=19,y0=89,n=n[[2]],direction=1),#TM1
              smooth.curve(x0=19,x1=29,y0=68,n=n[[3]],direction=-1,n.diff=0),#IC1
              tmpart(x0=29,y0=90,n=n[[4]],direction=-1),#TM2
              smooth.curve(x0=30,x1=40,y0=91,n=n[[5]],direction=1,n.diff=0),#EC2
              tmpart(x0=39,y0=91,n=n[[6]],direction=1),#TM3
              smooth.curve(x0=39,x1=49,y0=64.5,n=n[[7]],direction=-1,n.diff=-2),#IC2
              tmpart(x0=49,y0=87,n=n[[8]],direction=-1),#TM4
              smooth.curve(x0=49,x1=59,y0=86,n=n[[9]],direction=1,n.diff=2),#EC3
              tmpart(x0=59,y0=90,n=n[[10]],direction=1),#TM5
              smooth.curve(x0=59.5,x1=70,y0=66.5,n=n[[11]],direction=-1,n.diff=0),#IC3 #n.diff should be even numbers!?
              tmpart(x0=69,y0=91,n=n[[12]],direction=-1),#TM6
              smooth.curve(x0=71,x1=81,y0=91,n=n[[13]],direction=1,n.diff=0),#EC4
              tmpart(x0=79,y0=91,n=n[[14]],direction=1),#TM7
              end.curve(x0=81,y0=74,n=n[[15]],direction=-1))#IC4  
plot(cord)

pdf("hypo_pocket.pdf",8,8)
plot(cord,pch=21,col=MyGray,bg="gray",cex=2,frame.plot=F,axes=F,xlab=NA,ylab=NA)
points(cord,col=MyGray,bg="gray",type="l",lwd=2)
#points(cord,pch=21,col=MyGray,bg=ifelse(M72.sites$quan<0.1,"blue",MyTransparent),cex=2)
points(cord,pch=21,col=MyGray,bg=ifelse(Hyp.All.TMM$quan.ac<0.05,"orange",MyTransparent),cex=2)
points(cord,pch=21,col=MyGray,bg=ifelse(Hyp.All.TMM$quan.ac<0.01,"red",MyTransparent),cex=2)
points(cord,pch=21,col=ifelse(Hyp.All.TMM$freq>0.9,"magenta",MyTransparent),bg=MyTransparent,cex=2)
#points(cord,pch=21,col=ifelse(Hyp.ACR.TMM$All.AA.freq>0.9&Hyp.ACR.TMM$All.AA!="-","magenta",MyTransparent),bg=MyTransparent,cex=1.8)
#points(cord,pch=21,col=ifelse(Hyp.ACR.TMM$All.AA.freq>0.9&Hyp.ACR.TMM$All.AA!="-","magenta",MyTransparent),bg=MyTransparent,cex=2)
#points(cord,pch=21,col=MyGray,bg=ifelse(M72.sites$quan.4<0.03,"red",MyTransparent),cex=2)
for (i in 1:dim(Hyp.All.TMM)[1]) {text(x=cord[i,1],y=cord[i,2], labels=toupper(Hyp.All.TMM$AA.ACR[i]),col="gray27",cex=0.9)}
#plot.new();legend("bottomleft",c("conserved among acetophenone ORs than random (p<0.01)","conserved among acetophenone ORs than random (p<0.05)", "conserved in all ORs (>90%)"),pch=21,col=c("gray","gray","magenta"),pt.bg=c("red","orange",MyTransparent),bty="n",pt.cex=1,cex=1)
#plot.new();legend("bottomleft",c("conserved among acetophenone ORs than random (p<0.01)","conserved among acetophenone ORs than random (p<0.05)","conserved in all ORs (90%)"),pch=21,col=c("gray","gray","magenta"),pt.bg=c("red","orange","gray"),bty="n",pt.cex=1,cex=1)
plot.new();legend("bottomleft",c(expression(paste("Conserved among acetophenone ORs than random (",italic(p)<0.01,")")),
                                 expression(paste("Conserved among acetophenone ORs than random (",italic(p)<0.05,")")),
                                 "Conserved in all ORs (90%)"),
                  pch=21,col=c("gray","gray","magenta"),pt.bg=c("red","orange","gray"),bty="n",pt.cex=1,cex=1)

#plot.new();legend("bottomleft",c(expression('conserved in acetophenone ORs (p'<='0.00075, FDR adjusted:p<0.01)'),expression('conserved in acetophenone ORs (p'<='0.0085, FDR adjusted:p<0.05)'),"conserved in all ORs (90%)"),pch=21,col=c("gray","gray","magenta"),pt.bg=c("red","orange","gray"),bty="n",pt.cex=1,cex=1)
#expression('log'[10]*'(Mean read counts, no-odor control)')
plot(cord,pch=21,col=MyGray,bg="gray",cex=2,frame.plot=F,axes=F,xlab=NA,ylab=NA)
points(cord,col=MyGray,bg="gray",type="l",lwd=2)
#points(cord,pch=21,col=MyGray,bg=ifelse(M72.sites$quan<0.1,"blue",MyTransparent),cex=2)
points(cord,pch=21,col=MyGray,bg=ifelse(Hyp.All.TMM$quan.tmt<0.05,"orange",MyTransparent),cex=2)
points(cord,pch=21,col=MyGray,bg=ifelse(Hyp.All.TMM$quan.tmt<0.01,"red",MyTransparent),cex=2)
points(cord,pch=21,col=ifelse(Hyp.All.TMM$freq>0.9,"magenta",MyTransparent),bg=MyTransparent,cex=2)
#points(cord,pch=21,col=ifelse(Hyp.ACR.TMM$All.AA.freq>0.9&Hyp.ACR.TMM$All.AA!="-","magenta",MyTransparent),bg=MyTransparent,cex=1.8)
#points(cord,pch=21,col=ifelse(Hyp.ACR.TMM$All.AA.freq>0.9&Hyp.ACR.TMM$All.AA!="-","magenta",MyTransparent),bg=MyTransparent,cex=2)
#points(cord,pch=21,col=MyGray,bg=ifelse(M72.sites$quan.4<0.03,"red",MyTransparent),cex=2)
for (i in 1:dim(Hyp.All.TMM)[1]) {text(x=cord[i,1],y=cord[i,2], labels=toupper(Hyp.All.TMM$AA.TMTR[i]),col="gray27",cex=0.9)}
#plot.new();legend("bottomleft",c("conserved among acetophenone ORs than random (p<0.01)","conserved among acetophenone ORs than random (p<0.05)", "conserved in all ORs (>90%)"),pch=21,col=c("gray","gray","magenta"),pt.bg=c("red","orange",MyTransparent),bty="n",pt.cex=1,cex=1)
plot.new();legend("bottomleft",c(expression(paste("Conserved among TMT ORs than random (",italic(p)<0.01,")")),
                                 expression(paste("Conserved among TMT ORs than random (",italic(p)<0.05,")")),
                                 "Conserved in all ORs (90%)"),
                  pch=21,col=c("gray","gray","magenta"),pt.bg=c("red","orange","gray"),bty="n",pt.cex=1,cex=1)


Hyb.AA <- Hyp.All.TMM$AA.All
Hyb.AA[Hyp.All.TMM$quan.ac<0.05&!is.na(Hyp.All.TMM$quan.ac)] <- Hyp.All.TMM$AA.ACR[Hyp.All.TMM$quan.ac<0.05&!is.na(Hyp.All.TMM$quan.ac)]
Hyb.AA[Hyp.All.TMM$quan.tmt<0.05&!is.na(Hyp.All.TMM$quan.tmt)] <- Hyp.All.TMM$AA.TMTR[Hyp.All.TMM$quan.tmt<0.05&!is.na(Hyp.All.TMM$quan.tmt)]
plot(cord,pch=21,col=MyGray,bg="gray",cex=2,frame.plot=F,axes=F,xlab=NA,ylab=NA)
points(cord,col=MyGray,bg="gray",type="l",lwd=2)
#points(cord,pch=21,col=MyGray,bg=ifelse(M72.sites$quan<0.1,"blue",MyTransparent),cex=2)
points(cord,pch=21,col=MyGray,bg=ifelse(Hyp.All.TMM$quan.ac<0.05,"orange",MyTransparent),cex=2)
points(cord,pch=21,col=MyGray,bg=ifelse(Hyp.All.TMM$quan.tmt<0.05,"deepskyblue",MyTransparent),cex=2)
points(cord,pch=21,col=MyGray,bg=ifelse(Hyp.All.TMM$quan.tmt<0.05&Hyp.All.TMM$quan.ac<0.05,"red",MyTransparent),cex=2)
points(cord,pch=21,col=ifelse(Hyp.All.TMM$freq>0.9,"magenta",MyTransparent),bg=MyTransparent,cex=2)
#points(cord,pch=21,col=ifelse(Hyp.ACR.TMM$All.AA.freq>0.9&Hyp.ACR.TMM$All.AA!="-","magenta",MyTransparent),bg=MyTransparent,cex=1.8)
#points(cord,pch=21,col=ifelse(Hyp.ACR.TMM$All.AA.freq>0.9&Hyp.ACR.TMM$All.AA!="-","magenta",MyTransparent),bg=MyTransparent,cex=2)
#points(cord,pch=21,col=MyGray,bg=ifelse(M72.sites$quan.4<0.03,"red",MyTransparent),cex=2)
for (i in 1:dim(Hyp.All.TMM)[1]) {text(x=cord[i,1],y=cord[i,2], labels=toupper(Hyb.AA[i]),col="gray27",cex=0.9)}
#plot.new();legend("bottomleft",c("conserved among acetophenone ORs than random (p<0.01)","conserved among acetophenone ORs than random (p<0.05)", "conserved in all ORs (>90%)"),pch=21,col=c("gray","gray","magenta"),pt.bg=c("red","orange",MyTransparent),bty="n",pt.cex=1,cex=1)
plot.new();legend("bottomleft",c("Conserved among acetophenone ORs","Conserved among TMT ORs", "Conserved among both acetophenone and TMT ORs","Conserved in all ORs"),pch=21,col=c("gray","gray","gray","magenta"),pt.bg=c("orange","deepskyblue","red","gray"),bty="n",pt.cex=1,cex=1)


dev.off()


### cross-validation with Ivan's result
ACR.ivan <- c("Olfr145","Olfr901","Olfr109","Olfr478","Olfr745","Olfr983","Olfr1484",
              "Olfr736","Olfr922","Olfr1054","Olfr160","Olfr1377","Olfr47","Olfr1501",
              "Olfr476","Olfr19","Olfr143","Olfr229","Olfr30","Olfr556","Olfr1477","Olfr935")
ans.ivan.vitro <- data.frame(vitro=agg.all.luc.dedup.clean$f3,
                       res=ifelse(agg.all.luc.dedup.clean$OLFR%in%ACR.ivan,1,0),
                       id=agg.all.luc.dedup.clean$OLFR)
ans.ivan.vitro <- subset(ans.ivan.vitro,id%in%result.No_Ac.OR$id)
ROC.gen(ans.ivan.vitro)

ans.ivan.vivo100p <- data.frame(vivo=result.No_Ac.OR$FDR*(-1),
                             res=ifelse(result.No_Ac.OR$id%in%ACR.ivan,1,0),
                             id=result.No_Ac.OR$id)
#ans.ivan.vivo <- subset(ans.ivan.vivo,id%in%result.No_Ac.OR$id)
ROC.gen(ans.ivan.vivo100p)

ans.ivan.vivo1p <- data.frame(vivo=result.No_1pAc.OR$FDR*(-1),
                            res=ifelse(result.No_1pAc.OR$id%in%ACR.ivan,1,0),
                            id=result.No_1pAc.OR$id)
#ans.ivan.vivo <- subset(ans.ivan.vivo,id%in%result.No_Ac.OR$id)
ROC.gen(ans.ivan.vivo1p)

pdf("FigS7.pdf",5,5)
ROC.gen(ans.ivan.vitro)
ROC.gen(ans.ivan.vivo100p)
dev.off()
