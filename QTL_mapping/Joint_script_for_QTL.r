#BiocManager::install("qtl2",type="binary")
#BiocManager::install("qtl",type="binary")
#BiocManager::install("VariantAnnotation",type="binary")
#BiocManager::install('snpStats')

# make snpbinner format data ----------------------------------------------------
library(VariantAnnotation)
library(data.table)
system("python hanzhang_data/vcf_format.py")
vcf=readVcf("hanzhang_data/filter397.final.format.vcf")
snp.matrix <- genotypeToSnpMatrix(vcf, uncertain = FALSE)
snp.matrix.transposed <- t(as(snp.matrix$genotypes, "character"))
snp.matrix.transposed[1:10,1:10]
new_vcf=t(snp.matrix.transposed)
data_for_binner=cbind(row.names(snp.matrix.transposed),str_split_fixed(row.names(snp.matrix.transposed),"_",2)[,2],snp.matrix.transposed)
names(data_for_binner)[1:2]=c("marker","position(bp)")
data_for_binner[1:10,1:10]
data_for_binner=apply(data_for_binner,2,function(x){
  if(x[1]!="1_3999"|x[1]!="3399"){
    x[x=="NA"]="-"
    x[x=="A/A"]="a"
    x[x=="B/B"]="h"
    x[x=="A/B"]="b"
    return(x)
  }else{
    return(x)
  }
})
data_for_binner=as.data.frame(data_for_binner)
names(data_for_binner)[1:2]=c("marker","position(bp)")
data_for_binner[1:10,1:10]
for(i in 1:5){
  sub=data_for_binner[grepl(paste0(i,"_"),data_for_binner[,1]),]
  write.table(sub,file=paste0("hanzhang_data/snp_for_binner_chr",i,".txt"),row.names = F,col.names = T,quote = F,sep="\t")  
}


# make qtl format data ----------------------------------------------------
###load snpbinner result and combined together
for(i in 1:5){
  data=data.frame(fread(paste0("hanzhang_data\\bins_5k_chr",i,".txt")))
  data=data[c(-1,-2),]
  data=rbind(rep(as.character(i),ncol(data)),data)
  data=rbind(paste0(data[1,],"_",data[2,]),data)
  if(i==1){
    data[1,1]=""
    data[2,1]=""
    data[3,1]=""
    snpbin=data  
  }else{
    data=data[,-1]
    snpbin=cbind(snpbin,data)
  }
}
snpbin[1:10,1:10]
table(as.character(snpbin[4,]))
phe=data.frame(fread("hanzhang_data/all_phe.txt"))
row.names(phe)=phe[,1]
phe=phe[as.character(snpbin[4:nrow(snpbin),1]),]
qtl_csv=cbind(rbind(rep("",8),rep("",8),rep("",8),phe[,2:9]),snpbin)
qtl_csv[1,1:8]=names(phe)[2:9]
qtl_csv[1,9]="At_id"
qtl_csv[1:10,1:10]
qtl_csv[1:200,1:10]
write.table(qtl_csv,file="hanzhang_data/qtl_format_snp.csv",row.names = F,col.names = F,quote=F,sep=",")

dim(snpbin)
snpbin[1:10,1:10]
plot(as.numeric(snpbin[3,]))
# QTL mapping -------------------------------------------------------------
library(qtl)
library(qtl2)
library(pheatmap)
#risib means RIL
mydata <- read.cross(format="csv", file="hanzhang_data/qtl_format_snp.csv",
                     na.strings = "NA",genotypes = c("a", "b", "h"),crosstype="risib")
summary(mydata)
plotMap(mydata)
geno.image(mydata,reorder=FALSE,main="Genotype data",alternate.chrid=FALSE)
nm <- est.map(mydata, error.prob=0.001, verbose=FALSE)
plot.map(mydata, nm)
plot(nm)

mydata_qtl2=convert2cross2(mydata)
mydata_qtl2$crosstype


# get kinship -------------------------------------------------------------
pr <- calc_genoprob(mydata_qtl2, error_prob=0.0001, cores=4,quiet=F)
kinship=calc_kinship(pr)
pheatmap(kinship)


# mlm ---------------------------------------------------------------------
out_pg <- scan1(pr, mydata_qtl2$pheno, kinship)
plot(out_pg, mydata_qtl2$gmap, lodcolumn=1, col="tomato")
plot(out_pg, mydata_qtl2$gmap, lodcolumn=2, col="skyblue",add=T)
plot(out_pg, mydata_qtl2$gmap, lodcolumn=3, col="lightgreen",add=T)
plot(out_pg, mydata_qtl2$gmap, lodcolumn=4, col="orange",add=T)
legend("topleft", lwd=2, col=c("tomato","skyblue","lightgreen","orange"), c(names(mydata_qtl2$pheno[1,])[1:4]), bg="gray90", lty=c(1,1,2))

plot(out_pg, mydata_qtl2$gmap, lodcolumn=5, col="tomato")
plot(out_pg, mydata_qtl2$gmap, lodcolumn=6, col="skyblue",add=T)
plot(out_pg, mydata_qtl2$gmap, lodcolumn=7, col="lightgreen",add=T)
plot(out_pg, mydata_qtl2$gmap, lodcolumn=8, col="orange",add=T)
legend("topleft", lwd=2, col=c("tomato","skyblue","lightgreen","orange"), c(names(mydata_qtl2$pheno[1,])[5:8]), bg="gray90", lty=c(1,1,2))

#permutation threshold
operm <- scan1perm(pr, mydata_qtl2$pheno,kinship=kinship, n_perm=1000, cores=10)
summary(operm)

### get the QTL and the effect
find_peaks(out_pg, mydata_qtl2$gmap, threshold=2)

png("chr3_effect.png",width = 500,height=1200)
par(mfrow=c(8,1))
par(mar=c(1,5,1,3))
for(i in 1:8){
  chr3_effect=scan1coef(pr[,"3"],pheno=mydata_qtl2$pheno[,i],kinship=kinship)  
  col <- c("slateblue", "violetred")
  plot(chr3_effect,mydata_qtl2$gmap["3"], columns=1:2, col=col,ylab=paste("QTL effect",names(mydata_qtl2$pheno[1,])[i]))
  last_coef <- unclass(chr3_effect)[nrow(chr3_effect),]
  for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
}
dev.off()

png("chr4_effect.png",width = 500,height=1200)
par(mfrow=c(8,1))
par(mar=c(1,5,1,3))
for(i in 1:8){
  chr3_effect=scan1coef(pr[,"4"],pheno=mydata_qtl2$pheno[,i],kinship=kinship)  
  col <- c("slateblue", "violetred")
  plot(chr3_effect,mydata_qtl2$gmap["4"], columns=1:2, col=col,ylab=paste("QTL effect",names(mydata_qtl2$pheno[1,])[i]))
  last_coef <- unclass(chr3_effect)[nrow(chr3_effect),]
  for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
}
dev.off()


# segregation distortion --------------------------------------------------
gt <- geno.table(mydata,scanone.output=TRUE)
plot(unlist(mydata_qtl2$gmap),-log10(gt$P.value))
par(mfrow=c(2,1))
plot(gt, ylab=expression(paste(-log[10], " P-value")))
plot(gt, lod=2:4, ylab="Genotype frequency")
abline(h=c(0.4, 0.6), lty=2, col="gray")


# epistatic ---------------------------------------------------------------

hyper <- calc.genoprob(mydata, step=3, err=0.001)
out2 <- scantwo(mydata, verbose=FALSE,addcovar = kinship,phe.col=1)
summary(out2,thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6))
summary(out2)



# GbyE --------------------------------------------------------------------
library(stringr)
GbyE=data.frame(fread("hanzhang_data/GbyE/All_allele_info.txt",header = F))
for(i in 10:406){
  GbyE[2:4,i]=sub("\\|","/",str_split_fixed(GbyE[2:4,i],":",2)[,1])
}
gt=data.frame(t(GbyE[,10:406]))
row.names(gt)=gt[,1]
names(gt)=c("","BRC1","GRF8","GA5")
gt=gt[,-1]

phe=data.frame(fread("hanzhang_data/all_phe.txt"))
row.names(phe)=phe[,1]

GbyEdata=cbind(phe[row.names(gt),2:9],gt)
for(i in 1:8){
  GbyEdata[,i]=zscore(GbyEdata[,i])  
}

# PH GbyE -----------------------------------------------------------------
names(GbyEdata)
y=c(GbyEdata[,1],GbyEdata[,3],GbyEdata[,2],GbyEdata[,4])
env=as.factor(c(rep(c("HY","CD"),each=794)))
time=as.factor(c(rep(c("1","2","1","2"),each=397)))

G=as.factor(c(rep(GbyEdata[,9],4)))
lm=lm(y~env+G+time+G:env+G:time+G:env:time)
BRC1_av_ph=anova(lm)

G=as.factor(c(rep(GbyEdata[,10],4)))
lm=lm(y~env+G+G:env)
GRF8_av_ph=anova(lm)

G=as.factor(c(rep(GbyEdata[,11],4)))
lm=lm(y~env+G+G:env)
GA5_av_ph=anova(lm)

BRC1_av_ph
GRF8_av_ph
GA5_av_ph

###exclude time 
y=c(GbyEdata[,1],GbyEdata[,2])
y=c(GbyEdata[,1],GbyEdata[,4])
y=c(GbyEdata[,3],GbyEdata[,4])
y=c(GbyEdata[,3],GbyEdata[,2])
env=as.factor(c(rep(c("HY","CD"),each=397)));G=as.factor(c(rep(GbyEdata[,9],2)));lm=lm(y~env+G+G:env);anova(lm)

y=c(GbyEdata[,1],GbyEdata[,2])
y=c(GbyEdata[,1],GbyEdata[,4])
y=c(GbyEdata[,3],GbyEdata[,4])
y=c(GbyEdata[,3],GbyEdata[,2])
env=as.factor(c(rep(c("HY","CD"),each=397)));G=as.factor(c(rep(GbyEdata[,10],2)));lm=lm(y~env+G+G:env);anova(lm)

###4 env
y=c(GbyEdata[,1],GbyEdata[,3],GbyEdata[,2],GbyEdata[,4]);env=as.factor(c(rep(c("1","2","3","4"),each=397)));G=as.factor(c(rep(GbyEdata[,9],4)));lm=lm(y~env+G+G:env);anova(lm)
y=c(GbyEdata[,1],GbyEdata[,3],GbyEdata[,2],GbyEdata[,4]);env=as.factor(c(rep(c("1","2","3","4"),each=397)));G=as.factor(c(rep(GbyEdata[,10],4)));lm=lm(y~env+G+G:env);anova(lm)
# branch GbyE -----------------------------------------------------------------
names(GbyEdata)
y=c(GbyEdata[,5],GbyEdata[,7],GbyEdata[,6],GbyEdata[,8])
env=as.factor(c(rep(c("HY","CD"),each=794)))

G=as.factor(c(rep(GbyEdata[,9],4)))
lm=lm(y~env+G+G:env)
BRC1_av_br=anova(lm)

G=as.factor(c(rep(GbyEdata[,10],4)))
lm=lm(y~env+G+G:env)
GRF8_av_br=anova(lm)

G=as.factor(c(rep(GbyEdata[,11],4)))
lm=lm(y~env+G+G:env)
GA5_av_br=anova(lm)

BRC1_av_br
GRF8_av_br
GA5_av_br


y=c(GbyEdata[,5],GbyEdata[,6])
y=c(GbyEdata[,5],GbyEdata[,8])
y=c(GbyEdata[,7],GbyEdata[,8])
y=c(GbyEdata[,7],GbyEdata[,6])
env=as.factor(c(rep(c("HY","CD"),each=397)));G=as.factor(c(rep(GbyEdata[,9],2)));lm=lm(y~env+G+G:env);anova(lm)

y=c(GbyEdata[,5],GbyEdata[,6])
y=c(GbyEdata[,5],GbyEdata[,8])
y=c(GbyEdata[,7],GbyEdata[,8])
y=c(GbyEdata[,7],GbyEdata[,6])
env=as.factor(c(rep(c("HY","CD"),each=397)));G=as.factor(c(rep(GbyEdata[,10],2)));lm=lm(y~env+G+G:env);anova(lm)


y=c(GbyEdata[,5],GbyEdata[,7],GbyEdata[,6],GbyEdata[,8]);env=as.factor(c(rep(c("1","2","3","4"),each=397)));G=as.factor(c(rep(GbyEdata[,9],4)));lm=lm(y~env+G+G:env);anova(lm)
y=c(GbyEdata[,5],GbyEdata[,7],GbyEdata[,6],GbyEdata[,8]);env=as.factor(c(rep(c("1","2","3","4"),each=397)));G=as.factor(c(rep(GbyEdata[,10],4)));lm=lm(y~env+G+G:env);anova(lm)
# phe analysis ------------------------------------------------------------
boxplot(GbyEdata[,c(1,3,2,4)],ylim=c(0,80),ylab="PH",frame.plot=F)
mark=1
for(i in c(1,3,2,4)){
  points(x=rep(mark,nrow(GbyEdata)),y=GbyEdata[,i],pch=20)
  mark=mark+1
}
cols=c(rgb(1,0,0,0.2),rgb(0,0,1,0.2),rgb(0,1,0,0.2))
names(cols)=c("0/0","1/1","0/1")
j=10
for(i in 1:nrow(GbyEdata)){
  lines(x=c(1,2),y=c(GbyEdata[i,1],GbyEdata[i,3]),col=cols[GbyEdata[i,j]])
}
for(i in 1:nrow(GbyEdata)){
  lines(x=c(2,3),y=c(GbyEdata[i,3],GbyEdata[i,2]),col=cols[GbyEdata[i,j]])
}
for(i in 1:nrow(GbyEdata)){
  lines(x=c(3,4),y=c(GbyEdata[i,2],GbyEdata[i,4]),col=cols[GbyEdata[i,j]])
}


boxplot(GbyEdata[,c(5,7,6,8)],ylim=c(0,20),ylab="Branch",frame.plot=F)
mark=1
for(i in c(5,7,6,8)){
  points(x=rep(mark,nrow(GbyEdata)),y=GbyEdata[,i],pch=20)
  mark=mark+1
}
cols=c(rgb(1,0,0,0.2),rgb(0,0,1,0.2),rgb(0,1,0,0.2))
names(cols)=c("0/0","1/1","0/1")
j=9
for(i in 1:nrow(GbyEdata)){
  lines(x=c(1,2),y=c(GbyEdata[i,5],GbyEdata[i,7]),col=cols[GbyEdata[i,j]])
}
for(i in 1:nrow(GbyEdata)){
  lines(x=c(2,3),y=c(GbyEdata[i,7],GbyEdata[i,6]),col=cols[GbyEdata[i,j]])
}
for(i in 1:nrow(GbyEdata)){
  lines(x=c(3,4),y=c(GbyEdata[i,6],GbyEdata[i,8]),col=cols[GbyEdata[i,j]])
}


# IIIVmrMLM GbyE ----------------------------------------------------------
install.packages(c("lars","RcppEigen","Rcpp","doParallel","MASS","openxlsx","BEDMatrix","bigmemory","stringr","biglasso","progress","ncvreg","coin","sampling","sbl"),type="binary")
install.packages("IIIVmrMLM_1.0.zip",repos=NULL,type="binary")
#bcftools view filter397.final.format.vcf -Oz -o filter397.final.format.vcf.gz
#bcftools index filter397.final.format.vcf.gz
#bcftools filter filter397.final.format.vcf.gz --regions 3:5379516-7379516 > chr3_sub.vcf
#bcftools filter filter397.final.format.vcf.gz --regions 4:17804118-19804118 > chr4_sub.vcf
#plink -vcf sub_snp.vcf --make-bed --out sub_snp_plink
ids=data.frame(fread("hanzhang_data/qtl_format_snp.csv"))
row.names(kinship)=ids$At_id[3:399]
kinship[1:10,1:10]
for(i in 1:397){
  kinship[i,i]=1
}
write.table(kinship,file="hanzhang_data/GbyE/kin.csv",row.names = T,col.names = F,quote=F,sep=",")

library("IIIVmrMLM")
IIIVmrMLM(fileGen="hanzhang_data/GbyE/sub_snp_plink", 
          filePhe="hanzhang_data/GbyE/all_phe.csv", 
          fileKin="hanzhang_data/GbyE/kin.csv",
          filePS=NULL,fileCov=NULL, 
          method="Multi_env",trait=1:8,n.en=c(2,2,2,2,2,2,2,2),SearchRadius=20, svpal=0.01, 
          DrawPlot=TRUE,Plotformat="pdf",Chr_name_com=NULL,dir="hanzhang_data/GbyE/")
IIIVmrMLM(fileGen="hanzhang_data/GbyE/sub_snp_plink", 
          filePhe="hanzhang_data/GbyE/all_phe_two_trait.csv", 
          fileKin="hanzhang_data/GbyE/kin.csv",
          filePS=NULL,fileCov=NULL, 
          method="Multi_env",trait=1:2,n.en=c(4,4),SearchRadius=20, svpal=0.01, 
          DrawPlot=TRUE,Plotformat="pdf",Chr_name_com=NULL,dir="hanzhang_data/GbyE/")
IIIVmrMLM(fileGen="hanzhang_data/GbyE/three_snp_plink", 
          filePhe="hanzhang_data/GbyE/all_phe.csv", 
          fileKin="hanzhang_data/GbyE/kin.csv",
          filePS=NULL,fileCov=NULL, 
          method="Multi_env",trait=1:8,n.en=c(2,2,2,2,2,2,2,2),SearchRadius=20, svpal=0.01, 
          DrawPlot=TRUE,Plotformat="pdf",Chr_name_com=NULL,dir="hanzhang_data/GbyE/")
IIIVmrMLM(fileGen="hanzhang_data/GbyE/three_snp_plink", 
          filePhe="hanzhang_data/GbyE/all_phe_two_trait.csv", 
          fileKin="hanzhang_data/GbyE/kin.csv",
          filePS=NULL,fileCov=NULL, 
          method="Multi_env",trait=1:2,n.en=c(4,4),SearchRadius=20, svpal=0.01, 
          DrawPlot=TRUE,Plotformat="pdf",Chr_name_com=NULL,dir="hanzhang_data/GbyE/")

# QEI analysis ------------------------------------------------------------
vcf=readVcf("hanzhang_data/filter397.final.format.vcf")
snp.matrix <- genotypeToSnpMatrix(vcf, uncertain = FALSE)
snp.matrix.transposed <- t(as(snp.matrix$genotypes, "character"))
snp.matrix.transposed[1:10,1:10]

### 4 env 2trait
QTN_1=data.frame(fread("hanzhang_data/GbyE/2trait_result-QEI_detection/trait1_midresult.csv"))
QTN_1[QTN_1[,5]<1e-3,]
plot(-log10(QTN_1[,5]))

info=rbind(snp.matrix.transposed["3_6379516",],
           snp.matrix.transposed["4_18594118",])
table(info[1,]==snp.matrix.transposed["3_6181077",])
table(info[2,]==snp.matrix.transposed["4_18951365",])
QTN_1[QTN_1[,1]%in%c("3_6379516","4_18594118"),]


QTN_2=data.frame(fread("hanzhang_data/GbyE/2trait_result-QEI_detection/trait5_midresult.csv"))
QTN_2[QTN_2[,5]<1e-3,]
QTN_2[QTN_2[,1]%in%c("3_6379516","4_18594118"),]
plot(-log10(QTN_2[,5]))
a=snp.matrix.transposed["3_6181077",]
b=snp.matrix.transposed["3_6379516",]
data=cbind(phe[names(snp.matrix.transposed[1,]),],a,b)
par(mfrow=c(2,1))
for(i in 5:8){
  boxplot(data[,6]~data$a);summary(lm(data[,5]~data$a))
  boxplot(data[,6]~data$b);summary(lm(data[,5]~data$b))  
}
convert_info=c(1,2,3)
names(convert_info)=c("A/A","A/B","B/B")

plot(convert_info[a][order(convert_info[a])],pch=20,col="red",yaxt="n")
axis(2,c(1,2,3),c("A/A","A/B","B/B"))
plot(convert_info[b][order(convert_info[a])],pch=20,col="blue",yaxt="n")
axis(2,c(1,2,3),c("A/A","A/B","B/B"))


QTN_2=data.frame(fread("hanzhang_data/GbyE/8trait_result-QEI_detection/trait_midresult.csv"))
QTN_2[QTN_2[,5]<1e-3,]
QTN_2[QTN_2[,1]%in%c("3_6379516","4_18594118"),]


# gwas permutation threshold ----------------------------------------------
phe=data.frame(fread("hanzhang_data/all_phe.txt"))
row.names(phe)=phe[,1]
for(i in 2:9){
  print(i)
  sub=phe[,c(1,1,i)]
  for(j in 1:500){
    sub=cbind(sub,sample(phe[,i],397,replace=T))
  }
  names(sub)=c("family","id",paste0("phe",1:(ncol(sub)-2)))
  write.table(sub,file=paste0("hanzhang_data/permutation/phe_",(i-1),".txt"),col.names = T,row.names = F,quote=F)
}
library(parallel)
permutation_fun=function(x){
  gcta="D:\\OneDrive\\桌面\\搞科研\\软件\\gcta\\exe\\gcta64.exe"
  min_p=c()
  for(j in 1:500){
    system(paste0(gcta," --bfile hanzhang_data/permutation/snp_plink --grm-sparse hanzhang_data/permutation/grm --fastGWA-mlm --pheno ",paste0("hanzhang_data/permutation/phe_",x,".txt")," --mpheno ",j," --thread-num 1 --out F:\\permutation\\",paste0(x,"_",j)))    
    data=data.frame(fread(paste0("F:\\permutation\\",paste0(x,"_",j),".fastGWA")))
    min_p=c(min_p,min(na.omit(data[,10])))
  }
  return(min_p)
}
clnum<-8
cl <- makeCluster(getOption("cl.cores", clnum))
clusterEvalQ(cl,library(data.table))
permutation_p=parLapply(cl,1:8,permutation_fun)
stopCluster(cl)
save.image()

rbind(names(phe)[2:9],unlist(lapply(permutation_p,function(x){-log10(quantile(x,0.05))})))
par(mfrow=c(4,2))
for(i in 1:8){
  num=-log10(unlist(permutation_p[[i]][-1]))
  num[num>100]=100
  hist(num,breaks=seq(0,100,0.25),xlim=c(0,10),main=names(phe)[i+1])
  abline(v=quantile(num,0.95),col="red")
}

i=3
plot(density(-log10(unlist(permutation_p[[i]][-1]))))
abline(v=quantile(-log10(unlist(permutation_p[[i]][-1])),0.95),col="red")


for(x in 1:8){
  data=data.frame(fread(paste0("F:\\permutation\\",paste0(x,"_",j),".fastGWA")))
  names(data)[3]="BP"
  data=data[complete.cases(data),]
  if(x==1){
    sub=data[data$SNP%in%c("3_6379516","4_18594118"),]
  }else{
    sub=rbind(sub,data[data$SNP%in%c("3_6379516","4_18594118"),])
  }
  png(paste0(names(phe)[x+1],"_manhattan.png"),width = 1000,height=500)
  qqman::manhattan(data,col = c("blue4", "orange3"))
  dev.off()
}

save.image()



# use new phe -------------------------------------------------------------
save(permutation_p,file="hanzhang_data/permut_p.RData")
phe=data.frame(fread("hanzhang_data/hy2021.zg.fam"))
phe=phe[,c(1,2,6)]
sub=phe[,c(1,1,3)]
for(j in 1:500){
  sub=cbind(sub,sample(phe[,i],397,replace=T))
}
sub=as.data.frame(sub)
names(sub)=c("family","id",paste0("phe",1:(ncol(sub)-2)))
write.table(sub,file=paste0("hanzhang_data/permutation/phe_ph21.txt"),col.names = T,row.names = F,quote=F)
min_p=c()
x="ph"
for(j in 1:500){
  shell(paste0(gcta," --bfile hanzhang_data/permutation/snp_plink --grm-sparse hanzhang_data/permutation/grm --fastGWA-mlm --pheno hanzhang_data/permutation/phe_ph21.txt --mpheno ",j," --thread-num 2 --out F:\\permutation\\",paste0(x,"_",j)))    
  data=data.frame(fread(paste0("F:\\permutation\\",paste0(x,"_",j),".fastGWA")))
  min_p=c(min_p,min(na.omit(data[,10])))
}
save(min_p,file="hanzhang_data/permut_ph_p.RData")
-log10(quantile(min_p,0.05))
num=-log10(min_p)
num[num>100]=100
hist(num,breaks=seq(0,100,0.25),xlim=c(0,10),main=names(phe)[i+1])
abline(v=quantile(num,0.95),col="red")
