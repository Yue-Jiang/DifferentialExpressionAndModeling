library(e1071)
library(seqinr)
library(ape)
library(pROC)

AAprop <- read.table("pRIP/OR sequences/AAproperty.txt",header=T,sep="\t")
#AAprop$c <- scale(AAprop$c)
#AAprop$p <- scale(AAprop$p)
#AAprop$v <- scale(AAprop$v)
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

Prop.OR.raw <- as.data.frame(Prop.OR.raw)
colnames(Prop.OR.raw) <- colNam
row.names(Prop.OR.raw) <- OR.aln$nam
keep.ind <- (colMeans(All.OR=="-") < 0.1)
Prop.OR <- Prop.OR.raw[,as.vector(sapply(keep.ind, function (x) rep(x,3)))]
for (i in 1:dim(Prop.OR)[2]) {Prop.OR[is.na(Prop.OR[,i]),i] <- mean(Prop.OR[,i],na.rm=T)}
#var.ind <- apply(Prop.OR,2,var)!=0
#Prop.OR <- Prop.OR[var.ind]

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

data <- data.frame(res=c(rep(1,length(ACR.ind)),rep(0,length(notACR.ind))),Prop.OR[c(ACR.ind,notACR.ind),])
row.names(data) <- row.names(Prop.OR)[c(ACR.ind,notACR.ind)]
#rm.cnst.predcol <- function(data1, data2) ### mark constant colomns except the first in data2, remove them in data1
#  {data1.right <- data1[,-1]; return (cbind(data1[,1],data1.right[apply(data2[,-1],2,var)!=0]))}

ROC.gen <- function(ans,add=F) {  
  ans <- ans[order(ans[,1],decreasing=T),]
  tp <- fp <- 0
  TP <- FP <-  NULL
  tp.n <- sum(ans[,2]==1)
  fp.n <- sum(ans[,2]==0)
  l <- dim(ans)[1]
  for (i in 1:l) {
    if (ans[i,2]==1) {tp <- tp+1}
    if (ans[i,2]==0) {fp <- fp+1}
    TP <- c(TP,tp/tp.n)
    FP <- c(FP,fp/fp.n)
  }
  auc <- TP[1:(l-1)]%*%(FP[2:l]-FP[1:(l-1)])
#  plot(c(0,FP),c(0,TP),type="l",col="blue",lwd=3,xlab="False positive rate",ylab="True positive rate",frame.plot=F,xlim=c(0,1),ylim=c(0,1))
  if (add==F) {
    plot(c(0,FP),c(0,TP),type="l",col="blue",lwd=3,xlab=NA,ylab=NA,frame.plot=F)
    title(xlab="False positive rate",ylab="True positive rate",cex.lab=1.5)
    segments(0,0,1,1,lwd=3)
#    legend("bottomright",c("protein sequence-based classifier","random"),col=c("blue","black"),lwd=c(3,3),lty=c(1,1),bty="n",cex=1)
#    legend("bottomright",c("in vitro response classifier","Random"),col=c("blue","black"),lwd=c(3,3),lty=c(1,1),bty="n",cex=1)
    legend("bottomright",c("in vivo enrichment classifier","Random"),col=c("blue","black"),lwd=c(3,3),lty=c(1,1),bty="n",cex=1)
    text(0.1,1,paste("AUC:",round(auc,3)))
  }
  else {points(c(0,FP),c(0,TP),type="l",col="green")}
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  text(5,5,paste("p=",wilcox.test(ans[,1]~ans[,2],alternative="less")$p.value))
  text(5,4,paste("AUC:",round(auc,3)))
#  plot(sort(ans[,1]),TP,type="l",col="blue",lwd=3,xlab=NA,ylab=NA,frame.plot=F)
#  points(sort(ans[,1]),FP,type="l",col="green")
#  points(sort(ans[,1]),TP/FP,type="l",col="red")
  return (list(auc,FP,TP))
}

library(e1071)
deg=3
fold=10
cycle=10
test.size=round(dim(data)[1]/fold)
#BIG.TP<-BIG.FP<-AUC<-NULL
ANS <- NULL
set.seed(1)
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
  svm.model <- svm(as.factor(res) ~ ., data=train.data,kernel="linear",cost = 100,scale=T,probability=TRUE)
#  svm.model <- svm(as.factor(res) ~ ., data=train.data,kernel="linear",scale=T,probability=TRUE)
  svm.pred <- predict(svm.model,test.data,na.action = na.exclude, probability=TRUE,decision.values = TRUE)
  ans <- cbind(attr(svm.pred, "decision.values"),test.data$res,attr(svm.pred, "probabilities")[,1])
  ANS <- rbind(ANS,ans)
}
  #ans <- ans[order(ans[,2],decreasing=T),]
ROC.gen(ANS)
plot(sort(ANS[,1],decreasing=T),crosval[[3]]/crosval[[2]],type="l")
ROC1 <- roc(ANS[,2],ANS[,1])
pdf("crossvalidation.pdf",5,5)
crosval <- ROC.gen(ANS)
dev.off()

ANS[order(ANS[,1],decreasing=T),][1:80,]


#test.otherOR <- data.frame(res=NA,Prop.OR[row.names(Prop.OR)%in%setdiff(agg.all.luc.dedup.clean$OLFR,c(result.No_1pAc.OR$id,result.No_Ac.OR$id)),])
#test.otherOR <- data.frame(res=c(0,0),Prop.OR[row.names(Prop.OR)%in%setdiff(agg.all.luc.dedup.clean$OLFR,c(subset(result.No_1pAc.OR,baseMean>0)$id,subset(result.No_Ac.OR,baseMean>0)$id)),])

#res <- paste.data()
#test.otherOR$res <- res$res.compute. #this is from testotherORluc.xlsx

#svm.model <- svm(as.factor(res) ~ ., data=data[apply(data,2,var)!=0],
 #                kernel="polynomial",degree=3,
#                 #kernel="linear",
#                 scale=T,probability=TRUE)
#svm.pred <- predict(svm.model,test.otherOR[apply(data,2,var)!=0],na.action = na.exclude, probability=TRUE,decision.values = TRUE)
#ans.ext3 <- cbind(attr(svm.pred, "probabilities")[,1],test.otherOR$res)
#ROC.ext1 <- roc(ans.ext[,2],ans.ext[,1])

#pdf("externalvalidation.pdf",5,5)
#ext <- ROC.gen(ans.ext)
#dev.off()


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
cros.val.sim <- ROC.gen(ANS2)
ROC2 <- roc(ANS2[,2],ANS2[,1])

roc.test(ROC1,ROC2)

#external, using the ORs with less than 1cp as test
#before this first set 1cpm cutoff and regenerate result.No_1pAc.OR and result.No_Ac.OR

result.No_Ac.OR <- filter.OR (result.No_Ac,add=0,cpm=1)
result.No_1pAc.OR <- filter.OR (result.No_1pAc,add=0,cpm=1)

  test.ind <- which(row.names(data)%in%setdiff(row.names(data),c(result.No_1pAc.OR$id,result.No_Ac.OR$id)))
  test.data <- data[test.ind,]
  train.data <- data[-test.ind,]
  col.ind <- apply(train.data,2,var)!=0
  col.ind[1] <- TRUE
  train.data <- train.data[,col.ind]
  test.data <- test.data[,col.ind]
  svm.model <- svm(as.factor(res) ~ ., data=train.data,kernel="polynomial",degree=deg,scale=T,probability=TRUE)
#  svm.model <- svm(as.factor(res) ~ ., data=train.data,kernel="linear",scale=T,probability=TRUE)
  svm.pred <- predict(svm.model,test.data,na.action = na.exclude, probability=TRUE,decision.values = TRUE)
  ans <- cbind(attr(svm.pred, "decision.values"),test.data$res,attr(svm.pred, "probabilities")[,1])
  ROC.ext <- ROC.gen(ans)
  roc.ext <- roc(ans[,2],ans[,1])
  


sim.pred <- similarity.predict(row.names(train.data),train.data$res,row.names(test.data),OR.dist.mtx)
ans.sim <- data.frame(pred=sim.pred[[1]],res=test.data$res,most.similar=sim.pred[[3]])
ROC.gen(ans.sim)
ROC.ext.sim <- ROC.gen(ans.sim)
roc.ext.sim <- roc(ans.sim[,2],ans.sim[,1],direction="<")
roc.test(roc.ext,roc.ext.sim)


#ROC.ext2 <- roc(ans.sim[,2],ans.sim[,1])
#plot(ROC.ext2)
#roc.test(ROC.ext1,ROC.ext2)
#ext.sim <- ROC.gen(ans.sim)


pdf("model.pdf",5,5)
plot(crosval[[2]],crosval[[3]],type="l",col=MyTransparent,xlim=c(-0.1,1.1),ylim=c(-0.1,1.1),frame.plot=F,lwd=3,xlab=NA,ylab=NA)
title(xlab="False positive rate",ylab="True positive rate",cex.lab=1.5)
segments(0,0,1,1,lwd=3,col="gray")
points(cros.val.sim[[2]],cros.val.sim[[3]],type="l",col="green",lwd=3)
points(crosval[[2]],crosval[[3]],type="l",col="blue",lwd=3)
legend("bottomright",c("trained classifier","overall similarity-based","random"),col=c("blue","green","gray"),lwd=3,lty=1,bty="n",cex=1)

plot(ROC.ext[[2]],ROC.ext[[3]],type="l",col=MyTransparent,xlim=c(-0.1,1.1),ylim=c(-0.1,1.1),frame.plot=F,lwd=3,xlab=NA,ylab=NA)
title(xlab="False positive rate",ylab="True positive rate",cex.lab=1.5)
segments(0,0,1,1,lwd=3,col="gray")
points(ROC.ext.sim[[2]],ROC.ext.sim[[3]],type="l",col="green",lwd=3)
points(ROC.ext[[2]],ROC.ext[[3]],type="l",col="blue",lwd=3)
legend("bottomright",c("trained classifier","overall similarity-based","random"),col=c("blue","green","gray"),lwd=3,lty=1,bty="n",cex=1)

dev.off()





#plot(crosval[[2]],crosval[[3]],type="l",col=MyTransparent,xlim=c(-0.1,1.1),ylim=c(-0.1,1.1),frame.plot=F,lwd=3,xlab=NA,ylab=NA)
#title(xlab="False positive rate",ylab="True positive rate",cex.lab=1.5)
#segments(0,0,1,1,lwd=3,col="gray")
#points(cros.val.sim[[2]],cros.val.sim[[3]],type="l",col="green",lwd=3)
#legend("bottomright",c("overall similarity-based","random"),col=c("green","black"),lwd=3,lty=1,bty="n",cex=1)

#plot(ext[[2]],ext[[3]],type="l",col=MyTransparent,xlim=c(-0.1,1.1),ylim=c(-0.1,1.1),frame.plot=F,lwd=3,xlab=NA,ylab=NA)
#title(xlab="False positive rate",ylab="True positive rate",cex.lab=1.5)
#segments(0,0,1,1,lwd=3,col="gray")
#points(ext.sim[[2]],ext.sim[[3]],type="l",col="green",lwd=3)
#points(c(0,ext[[2]]),c(0,ext[[3]]),type="l",col="blue",lwd=3)
#legend("bottomright",c("trained classifier","overall similarity-based","random"),col=c("blue","green","black"),lwd=3,lty=1,bty="n",cex=1)

#plot(ext[[2]],ext[[3]],type="l",col=MyTransparent,xlim=c(-0.1,1.1),ylim=c(-0.1,1.1),frame.plot=F,lwd=3,xlab=NA,ylab=NA)
#title(xlab="False positive rate",ylab="True positive rate",cex.lab=1.5)
#segments(0,0,1,1,lwd=3,col="gray")
#points(ext.sim[[2]],ext.sim[[3]],type="l",col="green",lwd=3)
#legend("bottomright",c("overall similarity-based","random"),col=c("green","black"),lwd=3,lty=1,bty="n",cex=1)




mtx <- as.matrix(1:dim(All.OR)[1]) 
ind.mtx <- mtx[,rep(1,dim(All.OR)[2])]
randomORaa.ind <- apply(ind.mtx,2,sample,500000,replace=T)
randomORprop.ind <- randomORaa.ind[,as.vector(sapply(1:dim(All.OR)[2], function (x) rep(x,3)))]

index2content <- function (index, content.mtx) { return(content.mtx[cbind(index,1:length(index))])}
Prop.pseudo.OR.raw <- NULL
for (i in 1:dim(randomORprop.ind)[1]) {
  Prop.pseudo.OR.raw <- rbind(Prop.pseudo.OR.raw,index2content(randomORprop.ind[i,],Prop.OR.raw))
}

colnames(Prop.pseudo.OR.raw) <- colNam
#row.names(Prop.pseudo.OR.raw) <- OR.aln$nam
#keep.ind <- (colMeans(All.OR=="-") < 0.1)
Prop.pseuso.OR <- Prop.pseudo.OR.raw[,as.vector(sapply(keep.ind, function (x) rep(x,3)))]
for (i in 1:dim(Prop.pseuso.OR)[2]) {Prop.pseuso.OR[is.na(Prop.pseuso.OR[,i]),i] <- mean(Prop.OR[,i],na.rm=T)}


All.OR[cbind(randomORaa.ind[420,],1:558)]
pseudo.OR.data <- data.frame(res=rbinom(dim(Prop.pseuso.OR)[1],1,0.5),Prop.pseuso.OR) #this res just to fill up the space

svm.model <- svm(as.factor(res) ~ ., data=data[apply(data,2,var)!=0],kernel="polynomial",degree=deg,scale=T,probability=TRUE)
svm.pred <- predict(svm.model,pseudo.OR.data[apply(data,2,var)!=0],na.action = na.exclude, probability=TRUE,decision.values = TRUE)
hist(attr(svm.pred, "probabilities")[,1],breaks=100)


##S3VM###

system("~/Desktop/OR_prediction/svmlin-v1.0/svmlin")

system("mkdir temporary")
setwd("temporary/")

mksvmlin.format <- function(Prop) {
  svmlinProp <- NULL
  for (i in 1:dim(Prop)[1]) { #row:i
    row <- NULL
    for (j in 1:dim(Prop)[2]) { #col:j
      if (!is.na(Prop[i,j])) {
          row <- paste(row, paste(j,":",Prop[i,j],sep=""),sep=" ")
        }    
      }
    svmlinProp <- c(svmlinProp,row)
  }
  return (svmlinProp)
}
svmlinProp <- mksvmlin.format(Prop.OR.raw)

labels=rep(0,dim(Prop.OR)[1])
labels[OR.aln$nam%in%respond091914] <- 1
labels[OR.aln$nam%in%neither.ACR] <- -1



fold=10
cycle=100
test.size=round(dim(data)[1]/fold)
#BIG.TP<-BIG.FP<-AUC<-NULL
ANS <- NULL
set.seed(1)
for (j in 1:cycle) { 
  ans <- NULL
  repeat{
    test.ind <- sample((1:dim(Prop.OR)[1])[labels!=0],test.size)
    if (var(labels[test.ind])!=0) {break}
  }
  train.ind <- setdiff(1:dim(Prop.OR)[1],test.ind)
  cat(svmlinProp[train.ind],file="training_examples",sep="\n")
  cat(labels[train.ind],file="training_labels",sep="\n")
  cat(svmlinProp[test.ind],file="test_examples",sep="\n")
  system("~/Desktop/OR_prediction/svmlin-v1.0/svmlin training_examples training_labels")
  system("~/Desktop/OR_prediction/svmlin-v1.0/svmlin -f training_examples.weights test_examples")
  pred <- read.table("test_examples.outputs")$V1
  ans <- cbind(pred,labels[test.ind])
  ANS <- rbind(ANS,ans)
}
ANS[ANS[,2]==-1,2] <- 0
pdf("S3VM.pdf",5,5)
ROC.gen(ANS)
dev.off()

svmlin.test.otherOR <- mksvmlin.format(Prop.OR.raw[row.names(Prop.OR.raw)%in%setdiff(agg.all.luc.dedup.clean$OLFR,c(result.No_1pAc.OR$id,result.No_Ac.OR$id)),])

cat(svmlinProp,file="training_examples",sep="\n")
cat(labels,file="training_labels",sep="\n")
cat(svmlin.test.otherOR,file="test_examples",sep="\n")
system("~/Desktop/OR_prediction/svmlin-v1.0/svmlin -R 0.1 training_examples training_labels")
system("~/Desktop/OR_prediction/svmlin-v1.0/svmlin -f training_examples.weights test_examples")
pred <- read.table("test_examples.outputs")$V1
ans <- cbind(pred,test.otherOR$res)
ans[ans[,2]==-1,2] <- 0
ROC.gen(ans)



###
write.table(intersect(subset(result.No_1pAc.OR,baseMean==0)$id,subset(result.No_Ac.OR,baseMean==0)$id),file="testotherOR.xls",sep="\t",quote=F)


set.seed(1)
muta.start.names <- intersect(neither.ACR,subset(tmt_092414,RES==1)$OR)
svm.model <- svm(as.factor(res) ~ ., data=data[apply(data,2,var)!=0],kernel="polynomial",degree=deg,scale=T,probability=TRUE)
test.tmtOR <- data.frame(res=rep(0,length(muta.start.names)),Prop.OR[row.names(Prop.OR)%in%muta.start.names,])
svm.pred <- predict(svm.model,test.tmtOR[apply(data,2,var)!=0],na.action = na.exclude, probability=TRUE,decision.values = TRUE)
p0.tmtR <- attr(svm.pred, "probabilities")[,1]
p0.tmtR <- p0.tmtR[match(names(p0.tmtR),muta.start.names)]

OR.ind <- 1
keep.aa.ind <- (colMeans(All.OR=="-") < 0.1)
keep.prop.ind <- as.vector(sapply(keep.aa.ind, function (x) rep(x,3)))
p0 <- p0.tmtR[OR.ind]

start.aa <- All.OR[which(OR.aln$nam==muta.start.names[OR.ind]),]
start.prop <- Prop.OR.raw[which(OR.aln$nam==muta.start.names[OR.ind]),] #raw
size <- 5000
P <- NULL
AA <- start.aa

for (i in 1:size) {
#  muta.aa.ind <- sample((1:length(keep.aa.ind))[keep.aa.ind],50)
  muta.aa.ind <- sample(intersect(which(quan7.adj<0.2),(1:length(keep.aa.ind))[keep.aa.ind]),100)
  muta.prop.ind <- c(muta.aa.ind*3-2,muta.aa.ind*3-1,muta.aa.ind*3)
  repeat {
    to.OR.ind <- sample(which(row.names(Prop.OR.raw)%in%respond091914),1) #sample(1:dim(Prop.OR.raw)[1],1)#
    if (!"-"%in%All.OR[to.OR.ind,muta.aa.ind]) break
  }
  to.aa <- start.aa
  to.prop <- start.prop
  to.aa[muta.aa.ind] <- All.OR[to.OR.ind,muta.aa.ind]
  to.prop[muta.prop.ind] <- Prop.OR.raw[to.OR.ind,muta.prop.ind]
  to.data <- data.frame(res=c(0),to.prop[,keep.prop.ind])
  svm.pred <- predict(svm.model,to.data[apply(data,2,var)!=0],na.action = na.exclude, probability=TRUE,decision.values = TRUE)
  p1 <- attr(svm.pred, "probabilities")[,1]
#  if (p1>p0) {
    AA <- rbind(AA,to.aa)
#    start.aa <- to.aa
#    start.prop <- to.prop
#    p0 <- p1
    P <- c(P,p1)
#  }
}
plot(P)
paste(to.aa[to.aa!="-"],collapse="")

###v1, discard###
Olfr1395.mu <- "--------------------------------------------------------------------------mawagnqtli-----------shfvllglf-ddsdlhlflfsifmvmflvalsgnglmillilm----dsrl-htpmyfflswlsfmdlmlistivpkml--anfllg----rgsisfagcgleilffltllgdecfllafmaydryvaiynplryavimsrrvcwlmvigswsfglvdgtihavftlr------------fpycgs-----neidhffcevpavlkls--cadtslyetmiyvcvvlmlllpfsvisasylril---vavlrmrsaegrrkafatcsshmigvslfygaamimymrpqayh------sskqdkvvsafytmitpmlnpliyslrnkdvtgalrkl----------lgkcpcgggtlg----------------------------------------------------------------------------------------------------------"
Olfr31.mu   <- "---------------------------------------------------------------------------meg-dntss-----------tdftlmglf-ndeetqglvfatisvifltalvangilifliht----dahl-htpmyfllshlsfidamysstivpkml--vdylsg----qrtisfvgctaehflfltlvgaeffllglmaydryvaicnplrypvkmsrricwiiiigswfggsldgfiltpitms------------fpfcas-----reinhffceapavlkla--cadtalyetvlyvccvlmllitfsvvivsyaril---aailhmssvegrkkifatcsshmtvvslfygaaiytymvphsyh------spsqdkivsvfytiltpmlnpliysmrnkdvsgglrra-------lgkigssqrvskd-----f------------------------------------------------------------------------------------------------------"
Olfr463.mu  <- "---------------------------------------------------------------------------mepgn-----------ltwvsefvflgf-seiwelqvflfvvffcvylttvvgnlliivtvss----dprl-htpmyfllrnlavldacfssvttpkml--vnfls----ekksisyrgcmveiflfvflgtamcfflsvmaydrlvaisrplhyatimnsqlcmglvvasyvggfahsivqlvlmfs------------lpfcgp-----neldnfycdvpqvlala--cmdtslleflmisnsgmldviwfflllisysvil---vmlrsmsg-earrkiastctthiivvsmifipsiylyarpftpf--------tmdkavsifhtvvtpmlnpiiytlrnqemqaamkrl----------akrlalcnre------------------------------------------------------------------------------------------------------------"
Olfr96.mu   <- "------------------------------------------------------------------------mgilstgnqtvt----------e--fvltgfh-evpglhllffllftvlylsiitgnmllavvvvs----sqrl-htpmyfflvvlsfieivytttvvpkml--vgflq----eaksisvagclleffvfgslatdecfllavmaydrylaichplryphkmgpqwylglvltvwlsgfmvdglvvalmaq------------lrfcgp-----nlvdhfycdfsplmvla--csdtqvaqvltfvlsvvfltvpflivlisyaqiv---vtvlrvpsgtrrtkifstcsshlavvstfygtlmvlyivpssv------hsqllsmviallytvvipifnpviytlrnqevqqalrrl--------lyckptem----------------------------------------------------------------------------------------------------------------"
Olfr1427.mu <- "---------------------------------------------------------------------------melgn-----------htkvtefifcgl-tqsqelqlllffflsivyittvlvnltimvtvtw----esrl-hnpmyfflrnlsvldicfstitvpkvl--vdlls----rdktisfngcvteifffhllggadifslsvmafdrymaifrplhyvvkmsrgrctaliiaswvggfvhsiiqiflllp------------lpfcgp-----nevrsfycdvpqvlkla--ctdtfvlellmisnvglittlwfvlllvsytvil---tmllshtg-egrkkafstctshitvvtlhhvpciyvyarpftal--------pmdravsitlnivvpvlnpmiytlrnqdmksamkrl----------rkrlilseve------------------------------------------------------------------------------------------------------------"
Olfr180.mu  <- "---------------------------------------------------------------------------mektnh-----------slttqfilvgf-sdhpdlktplfllfsviylvtmvgnlglvaviyl----eprl-hnpmyiflgnlalmdsccscavtpkml--enfis----vdrrisyyecvaefyflctaetadcfllaamaydryvaicnplqyhsmmskklyiqmsigtfitgnlhslihvgcllr------------ltfcks-----nridhffcdilpllrls--ctnpninelmiyifsmpiqiftittvlvsyfcil---ltifkmkskdgrgkifstcashlfsvsifh-vcllmyirps--d------egnkdmpvavfytiiipllnpfiyslrnkevvnamkkv----mkthsifknvsas--------mar----------------------------------------------------------------------------------------------------"
Olfr181.mu  <- "---------------------------------------------------------------------------mmkanh-----------sltvefiltgf-sdhtdlktllfllfsaiylvtivgnlglvaliym----eprl-htpmyyflgnlafldsccscaitpkml--enfis----vdrrisypecmvefyflclaetadcfllaamaydryvaicnplqyhtmmskklsiqmsigtfiggnlhslihtvcllr------------lnfcas-----rridhffcdipplykls--ctdpfinelmlfifsmpiqvftittvlvsyscil---ltvlkmkskdgrgkifstcashffsvsifh-icllmyigps--k------nsnkdipvgvfytivipllnpfiyslrnkevvnavkkv----mkthsifknssas--------iah----------------------------------------------------------------------------------------------------"
Olfr449.mu  <- "--------------------------------------------------------------------------mdvdnqt-rv----------t-kfilvgfp-dslsmraavflmflvayiltvaenvgiillvqq----nrpl-hkpmyfflanlsfldacyisvtvpkll--fnfws----msnsisfthcmielyffialmctecvllaamaydryvaicrplhyptimshrlcfrlilgswvigfgistakiyfisr------------lsfcgp-----nvinhvfcdispvlnls--ctdmsiaelvdfvlalviflfpisitvlsygcil---atvlrmp--tgkqkifstcashlvvvtifysatifmyarprai------hafnmnmvisifyaivtpalnpfiyslrnrevkeamkkl---------iycqvirsd--------------------------------------------------------------------------------------------------------------"
Olfr464.mu  <- "---------------------------------------------------------------------------mepqn-----------itwvsefillgf-sqtqelqkllfvlflcvyittvvgnllimltvtf----dsrl-dmpmyfllrnlsvidlsystvtspkml--vdffh----ktktisyqgcmaeifffhtlgtgtvfflsvmaydrfiaiyqpllyvtimntklcvplivacwsggfahstvqtsamlp------------lpfcgp-----nvidnfycdvpqvlrla--ctdtsllqflmiynsgmlvliwfllllvsytvil---vmlrshss-qarrkiastctthiivvsmifipclyiysrpftpf--------pldkavsifytvltpmlnpmiytlrnkemqaamkrl----aehlvltkrnel----------------------------------------------------------------------------------------------------------------"

###v2###
Olfr1395.mu <- "--------------------------------------------------------------------------mawagnqtli-----------shfvllglf-dhsplhlflfsiimviflvalsgnglmillilg----dsrl-htpmyfflswlsfmdlclistitprma--adfllg----rgsisfagcgleillfltllgdecfllafmaydryvaisnplrysvimsrrvcwlmvagsylfglvdgliqavftlr------------fsycas-----qeidhffcevpavlkla--csdtslyetmlyvccvlmlllpfsvisasysril---vavlrmrsaegrrkafatcsshligvslfhgaamitymrpqanh------sseqdkvvsafytmitpmlnpliyslrnkevtgalrki----------lgkcpcgggtlg----------------------------------------------------------------------------------------------------------"
Olfr31.mu   <- "---------------------------------------------------------------------------megnnntss-----------tdftfmglf-nteetsglvfatisviflfalvangimiflihg----dahl-htpmyfllshlsfidamyistitpkml--vdyllg----qrtisfvgctaehflyltlvgaeffllglmaydryvaicnplrypvkmsrriyiiiiagswfggsldgflltpitms------------fsfcrs-----reinhffceapavlkla--cadtnlyetvmyvccvlmllipflvvissyaril---atvyhmssvegrkkifatcsshmtvvtlfygaaiytyvvphsyh------spsqdmifsvfytiltpmlnpliysmrnkdvsgglkra-------lgkigssqrvskd-----f------------------------------------------------------------------------------------------------------"
Olfr463.mu  <- "---------------------------------------------------------------------------mepnn-----------ltwvsefvflgf-seiwelqvflfvvflcvysftvvgnlgiivtvss----dprl-hnpmyyllrnlavldacfstvtapkml--vdfls----ekktisyrgcmveifffhtlggamvfflsvmaydrlvaisrplhyvtimnsqlcmglivaswvggfahstvqlslmlp------------lsfcgp-----nvlrnfycdvppvlala--cmdtsllefllisnvgmldviwfllvlisylvil---vmlrkhsg-earrkifstctthiivvsmifipsiylyvrpftpf--------tmdkavsishtvmtpmlnpmiyslrnqemqaamkrl----------akrlalcnre------------------------------------------------------------------------------------------------------------"
Olfr96.mu   <- "------------------------------------------------------------------------mgilstgnqtvt----------e--fvltgfh-evpglhllffsvftilylsiltgnmgiavvvvs----sqrl-htpmyfflvvlsfieavytstvvpkml--egflq----eat-isvagclleffvfgslgtdecfllavmaydryvaichplryphlmgpqvclgliltvylggfmvaglvvalmaq------------lrfcgp-----nlvdhfycdfsplmvls--csdtqvaqvtlfvcsvvfltvpfglvlisyaqiv---vtvlrvpsgtrrtkifstcgshlavvstfygtlmvlyivpsan------hsqllskviallytvvtpifnpviytlrnqevqqalrrl--------lycnptem----------------------------------------------------------------------------------------------------------------"
Olfr1427.mu <- "---------------------------------------------------------------------------melgn-----------htkvtefifcgl-tdsqelslllffflfivyittvlvnvgimvtvtw----esrl-htpmyfllrnlsfldicfstitvpkvl--vdlis----rrktisfngcfteifffhllggadifllsamafdrymaifrplhyvtimsrgrctaliaaswsggivhsivhtflllp------------lpfcgp-----nvvdsfycdvppvlala--ctdtfvlellmisnnglittlwfvlllvsytvil---tmlrshts-egrkkiistctshitvvtlhfvpciyvyarpftal--------pmdravsiflniivpvlnpmiytlrnqemksamkrl----------rkrlilseve------------------------------------------------------------------------------------------------------------"
Olfr180.mu  <- "---------------------------------------------------------------------------meknnh-----------slttqfilvgf-sdhpdlktplfllfsviylvtmvgnlglvaviyl----eprl-htpmyyflgnlalmdsccscaitpkml--enffs----vdksislyecmaefyflvlaetadcfllaamaydryvaicnplqyhsmmspklyiqmsigtfitsnlhslihvvcllr------------ltfcks-----nridhffcdilplyrls--csdpfinelmiyyfsmpiqvftiltvlvsyfcil---ltifkmkskdgrrkifstcashffsvsify-vcllmyirpf--d------egnkdipvavfytiiipllnpiiyslrnkevvnavkkv----mkthsifknasas--------mar----------------------------------------------------------------------------------------------------"
Olfr181.mu  <- "---------------------------------------------------------------------------mmkanh-----------sltvefiligf-sdhtdlktllfllfsaiylvtlvgnlglvaliym----eprl-htpmyiflgnlafmdsccscaitpkml--enffs----vdrrislyecmvefyllclagttdcfllaamaydryvaicnplqyhtmmskrvyipmiigtfiasnlhatihtgcllr------------lsfcks-----rridhffcdilplyals--ctdtninelmlyifsmpiqvitittvlvsysfil---ltvfkmkskdgrgkifstcashffsvsify-icllmyigps-nk------nsnkdipvsvfytivipllnpfiyslrnkevvnavkkv----mkthsifknssas--------iah----------------------------------------------------------------------------------------------------"
Olfr449.mu  <- "--------------------------------------------------------------------------mdvdnqt-rv----------t-kfilvgfp-gslsmraavflmflvayiftvaenviiillvqq----nrpl-hnpmyyflanlsfletwyssvttpkll--fsfws----msnsisfthcvielflfialmctecvllaamaydryvaicrplhyptimspglcfrlalgswvigfgislakiyfisr------------lsfcgp-----nvinhffcdispvlnls--ctdmsiaelvlfvlalviflfpllitvvsygcil---atvlkmp--tgkqkifstcashlvgvtifysaiifmyarpsai------hafnmnmvisifyaivtpalnpfiyclrnrevkealkkl---------iycqvirsd--------------------------------------------------------------------------------------------------------------"
Olfr464.mu  <- "---------------------------------------------------------------------------mepqn-----------itwvsefillgf-sqtqelqkllfvvflcvyiftvignilimitvtg----dprl-dmpmyfllrnlavidlcystvtspkml--vdffh----ktktisypgcmaeifffvllgggtvfflsvmaydryiaisqplhyvtimntrlcvglvvaswsggfahsivqlslmlp------------lpfcgp-----nvidnfycdvpqvlria--ctdtslleflmisnsgmlvliwfllllisytvil---vmlrshss-earrkaastctshitvvsmyhipciyiysrpftpf--------ehdkavsisytvltpmlnpiiytlrnqemqeamkrl----aehlvltkrnel----------------------------------------------------------------------------------------------------------------"

muta.aa <- rbind(strsplit(Olfr1395.mu,"")[[1]],
                 strsplit(Olfr31.mu,"")[[1]],
                 strsplit(Olfr463.mu,"")[[1]],
                 strsplit(Olfr96.mu,"")[[1]],
                 strsplit(Olfr1427.mu,"")[[1]],
                 strsplit(Olfr180.mu,"")[[1]],
                 strsplit(Olfr181.mu,"")[[1]],
                 strsplit(Olfr449.mu,"")[[1]],
                 strsplit(Olfr464.mu,"")[[1]])
                
Prop.muta.raw <- NULL
for (i in (1:dim(muta.aa)[1])) {
  row <- NULL
  for (j in (1:dim(muta.aa)[2])) {
      row <- c(row,AA2prop(muta.aa[i,j],AAprop))
  }
  Prop.muta.raw <- rbind(Prop.muta.raw,row)
}
Prop.muta.raw <- as.data.frame(Prop.muta.raw)
colnames(Prop.muta.raw) <- colNam

keep.ind <- (colMeans(All.OR=="-") < 0.1)
Prop.muta <- Prop.muta.raw[,as.vector(sapply(keep.ind, function (x) rep(x,3)))]
for (i in 1:dim(Prop.muta)[2]) {Prop.muta[is.na(Prop.muta[,i]),i] <- mean(Prop.OR[,i],na.rm=T)}
muta.data <- data.frame(res=rep(0,dim(muta.aa)[1]),Prop.muta)

svm.model <- svm(as.factor(res) ~ ., data=data[apply(data,2,var)!=0],kernel="polynomial",degree=3,scale=T,probability=TRUE)
predict(svm.model,muta.data[apply(data,2,var)!=0],na.action = na.exclude, probability=TRUE,decision.values = TRUE)
attr(svm.pred, "probabilities")[,1]


MoHuOR.aln <- read.alignment(file="OR_prediction/alignments/Mo+Hu 27.aln",format="clustal")
MoHuOR.aln <- read.alignment(file="pRIP/OR sequences/mouseORprotein_clean_v17_2W1_5P3_1A1_2C1_2J2_2M7_10J5_51E1_51E2_51L1.aln",format="clustal")
MoHuOR.aln$nam[1:1090] <- OR.aln$nam

MoHuOR <- NULL
for (i in 1:length(MoHuOR.aln$nam)) {
  MoHuOR <- rbind(MoHuOR, strsplit(MoHuOR.aln$seq[[i]],NULL)[[1]])
}
Prop.MoHuOR.raw <- NULL 
for (i in 1:dim(MoHuOR)[1]) {
  row <- NULL
  for (j in (1:dim(MoHuOR)[2])) {
    row <- c(row,AA2prop(MoHuOR[i,j],AAprop))
  }
  Prop.MoHuOR.raw <- rbind(Prop.MoHuOR.raw,row)
}
colNam2 <- NULL
for (i in 1:dim(MoHuOR)[2]) {colNam2 <- c(colNam2,paste(i,"c",sep="_"),paste(i,"p",sep="_"),paste(i,"v",sep="_"))}

Prop.MoHuOR.raw <- as.data.frame(Prop.MoHuOR.raw)
colnames(Prop.MoHuOR.raw) <- colNam2
row.names(Prop.MoHuOR.raw) <- c(OR.aln$nam,MoHuOR.aln$nam[1091:1118])
keep.ind2 <- (colMeans(MoHuOR[1:1090,]=="-") < 0.1)
Prop.MoHuOR <- Prop.MoHuOR.raw[,as.vector(sapply(keep.ind2, function (x) rep(x,3)))]
for (i in 1:dim(Prop.MoHuOR)[2]) {Prop.MoHuOR[is.na(Prop.MoHuOR[,i]),i] <- mean(Prop.MoHuOR[1:1090,i],na.rm=T)}

Modata <- data.frame(res=c(rep(1,length(ACR.ind)),rep(0,length(notACR.ind))),Prop.MoHuOR[c(ACR.ind,notACR.ind),])
row.names(Modata) <- row.names(Prop.OR)[c(ACR.ind,notACR.ind)]

#Hudata <- data.frame(res=c(1,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),Prop.MoHuOR[1091:1118,])
Hudata <- data.frame(res=c(1,0,0,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),Prop.MoHuOR[1091:1118,])
#Hudata <- data.frame(res=c(1,1,0,0,0,0,0,0,0,0),Prop.MoHuOR[1091:1100,])

svm.model <- svm(as.factor(res) ~ ., data=Modata[apply(Modata,2,var)!=0],kernel="polynomial",degree=3,scale=T,probability=TRUE)
svm.pred <- predict(svm.model,Hudata[apply(Modata,2,var)!=0],na.action = na.exclude, probability=TRUE,decision.values = TRUE)


pred <- attr(svm.pred, "decision.values")[,1]
ans <- cbind(pred,Hudata$res)

## elastic-net
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
Hu.elnet <- ROC.gen(Hu.ans)
wilcox.test(Hu.ans[,1]~Hu.ans[,2],alternative="less")
#glm(res~X119_c, data=data, family=binomial(logit))
Hu.ROC <- roc(Hu.ans[,2],Hu.ans[,1])
roc.test(Hu.ROC.sim,Hu.ROC)

fla.small <- paste("factor(res) ~", paste(row.names(coef(glmnet.model,s="lambda.1se"))[which(coef(glmnet.model,s="lambda.1se")!=0)][-1], collapse="+"))

### cross-validation on selected model
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
  elnet.model <- glm(as.formula(fla.small), data=train.data, family=binomial(logit))
  #elnet.model <- glm(y=as.factor(train.data$res), x=as.matrix(train.data[,-1]),family="binomial",alpha=alpha,lambda=lambda)
  elnet.pred <- predict(elnet.model,test.data,type="response")
  ans <- cbind(elnet.pred,test.data$res)
  ANS4 <- rbind(ANS4,ans)
}
crosval.4 <- ROC.gen(ANS4)
ROC4 <- roc(ANS4[,2],ANS4[,1])
roc.test(ROC1,ROC4)


## mouse human similarity based
MoHudist <- dist.alignment(MoHuOR.aln,matrix="similarity")
MoHudist.mtx <- as.matrix(MoHudist)
ans2 <- cbind(similarity.predict(row.names(Modata),Modata$res,row.names(Hudata[-26,]),MoHudist.mtx)[[1]],Hudata$res[-26])
Hu.sim <- ROC.gen(ans2)

Hu.ROC.sim <- roc(ans2[,2],ans2[,1])

