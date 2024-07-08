# Segmental_analysis
 Some segemntal analysis code, which contain the QTL mapping, MTLMM, Chip-seq, ATAC and so on.

 ## Table of contents
|Analysis |Description |
| --- | --- |
|[QTL mapping ](#1)|Inplement SNPbinner and R/qtl(2) to perform QTL mapping|
|[Chip-seq](#2)||
|[ATAC-seq](#3)||
|[ATAC-seq](#4)||
|[Two trait GWAS](#5)||
|[Coloc](#6)||




### <h1 id="1">QTL mapping  </h1>
This section contains the QTL mapping Steps.

To identify the crucial loci in genome which associated with a trait/phenotype in the family data, QTL mapping was more suitable than GWAS (except Magic population). The excistence of the low recombination rate, highly linkaged SNP/QTN and the strong population structure makes GWAS ananlysis ignor some real causual loci.

This analysis need some softwares/packages as below:
- SNPbinner
- R qtl package
- R qtl2 package

The analysis pipeline were come from these webs:
https://github.com/solgenomics/snpbinner
https://www.nature.com/articles/s41588-022-01172-2

## Data preperation
- VCF files for all snps in all individuals.
```txt
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  0-13C   0-24A   0-30A   0-31A   0-33A
1       3999    1_3999  C       T       139073  PASS    *      GT:AD:DP:GQ:PL  0/0     0/0     1/1     0/0     1/1
1       4044    1_4044  C       T       130665  PASS    *       GT:AD:DP:GQ:PL  0/0     0/0     1/1     0/0     1/1
1       4446    1_4446  T       C       158841  PASS    *       GT:AD:DP:GQ:PL  0/0     0/0     1/1     0/0     1/1
1       4996    1_4996  G       A       20340   PASS    *  GT:AD:DP:GQ:PL   0/1     0/1     0/1     0/1     0/1
1       5018    1_5018  A       T       26137.9 PASS    *       GT:AD:DP:GQ:PL  0/1     0/1     0/1     0/1     0/1
1       5020    1_5020  A       T       26130.3 PASS    *       GT:AD:DP:GQ:PL  0/1     0/1     0/1     0/1     0/1
1       5885    1_5885  A       C       164025  PASS    *     GT:AD:DP:GQ:PL  0/0     0/0     1/1     0/0     1/1
1       5963    1_5963  T       G       167931  PASS    *      GT:AD:DP:GQ:PL  0/0     0/0     1/1     0/0     1/1
```
- phenotype files for all individuals.
```txt
id      20HY_ph 21CD_ph 21HY_ph 22CD_ph 20HY_branch     21CD_branch     21HY_branch     22CD_branch
1-18A   23.182  20.455  19.67   18.9775 1       1       0.2     0
1-22A   27.378  21.008  16.29   27.33   1       1       0       0
1-32A   34.91   20.616  426.252 24.255  0       1       1.5     0.5
1-52A   26.836  25.292  15.71   24.38   1.4     1       0.8     0
1-77C   19.374  26.13   23.14   24.525  1.4     1       2       0.25
2-4A    14.876  19.964  15.2675 20.67   1       1       1.5     0
2-5A    24.226  18.15   21.695  20.4525 1.8     1       1.4     0
2-6A    22.25   19.106  19.256  16.065  1.6     1       2.4     0
2-17A   24.734  19.15   11.592  NA      1       1       0       NA
```
- family info
## Data convertion
### Convert SNP genotypes to ABH format
```r
library(VariantAnnotation)
library(data.table)
#load vcf data
vcf=readVcf("vcf.vcf")
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
  write.table(sub,file=paste0("snp_for_binner_chr",i,".txt"),row.names = F,col.names = T,quote = F,sep="\t")  
}
```

After comverting the vcf format in R, the genotypes of SNP were converted to ABH format and seperated into each chromosome.

```txt
marker  position(bp)    0-13C   0-24A   0-30A   0-31A   0-33A   0-43A   0-49C   0-4A
1_3999  3999    a       a       h       a       h       h       h       a
1_4044  4044    a       a       h       a       h       h       h       a
1_4446  4446    a       a       h       a       h       h       h       a
1_4996  4996    b       b       b       b       b       b       a       b
1_5018  5018    b       b       b       b       b       b       a       b
1_5020  5020    b       b       b       b       b       b       a       b
1_5885  5885    a       a       h       a       h       h       h       a
1_5963  5963    a       a       h       a       h       h       h       a
1_6018  6018    a       a       h       a       h       h       h       a
```
### Splited snps into setted bins and genotyped each bin

Parameters see snpbinner github

Get the bin region of each chromosome.
```shell
snpbinner crosspoints --input ../hanzhang_data/snp_for_binner_chr1.txt  --output ../hanzhang_data/crosspoint_chr1.txt --min-ratio 0.01 &
snpbinner crosspoints --input ../hanzhang_data/snp_for_binner_chr2.txt  --output ../hanzhang_data/crosspoint_chr2.txt --min-ratio 0.01 &
snpbinner crosspoints --input ../hanzhang_data/snp_for_binner_chr3.txt  --output ../hanzhang_data/crosspoint_chr3.txt --min-ratio 0.01 &
snpbinner crosspoints --input ../hanzhang_data/snp_for_binner_chr4.txt  --output ../hanzhang_data/crosspoint_chr4.txt --min-ratio 0.01 &
snpbinner crosspoints --input ../hanzhang_data/snp_for_binner_chr5.txt  --output ../hanzhang_data/crosspoint_chr5.txt --min-ratio 0.01 &
```
Genotyped the bin region.

```shell
snpbinner bins --input ../hanzhang_data/crosspoint_chr1.txt --output ../hanzhang_data/bins_5k_chr1.txt --min-bin-size 5000
snpbinner bins --input ../hanzhang_data/crosspoint_chr2.txt --output ../hanzhang_data/bins_5k_chr2.txt --min-bin-size 5000
snpbinner bins --input ../hanzhang_data/crosspoint_chr3.txt --output ../hanzhang_data/bins_5k_chr3.txt --min-bin-size 5000
snpbinner bins --input ../hanzhang_data/crosspoint_chr4.txt --output ../hanzhang_data/bins_5k_chr4.txt --min-bin-size 5000
snpbinner bins --input ../hanzhang_data/crosspoint_chr4.txt --output ../hanzhang_data/bins_5k_chr5.txt --min-bin-size 5000
```

After processing, the format of output file is :
```txt
##bin start,0,366961,427101,438002,452405,462309,489922,754301,762254
##bin end,366961,427101,438002,452405,462309,489922,754301,762254,822248
bin center,183480,397031,432551,445203,457357,476115,622111,758277,792251
0-13C,a,a,a,a,a,a,a,a,a
0-24A,a,a,a,a,a,a,a,a,a
0-30A,h,h,h,h,h,h,h,h,h
0-31A,a,a,a,a,a,a,a,a,a
0-33A,h,h,h,h,h,h,h,h,h
0-43A,h,h,h,h,h,h,h,h,h
0-49C,h,h,h,h,h,h,h,h,h
```

### Combining the binfiles with phenotypes into R package qtls format

```r
#load snpbinner result and combined together
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
```

The format of csv file is :
```txt
X20HY_ph,X21CD_ph,X21HY_ph,X22CD_ph,X20HY_branch,X21CD_branch,X21HY_branch,X22CD_branch,At_id,1_183480,1_397031,1_432551,1_445203,1_457357,1_476115,1_622111,1_758277,1_792251,1_835919,1_855673
,,,,,,,,,1,1,1,1,1,1,1,1,1,1,1
,,,,,,,,,183480,397031,432551,445203,457357,476115,622111,758277,792251,835919,855673
14.0675,15.966,18.898,NA,7,4,10.8,NA,0-13C,a,a,a,a,a,a,a,a,a,a,a
7.502,6.85,12.198,8.56,5,2.2,6,2.75,0-24A,a,a,a,a,a,a,a,a,a,a,a
11.50666667,7.96,9.864,9.515,5,2.666666667,6.4,4.5,0-30A,h,h,h,h,h,h,h,h,h,h,h
21.186,17.734,23.57,18.7825,2.8,1.2,4.6,0.666666667,0-31A,a,a,a,a,a,a,a,a,a,a,a
16.62,21.26,19.992,20.6725,5.333333333,1.6,5.4,0.75,0-33A,h,h,h,h,h,h,h,h,h,h,h
24.4,25.7525,24.274,26.45,3.2,1.2,2.6,1,0-43A,h,h,h,h,h,h,h,h,h,h,h
27.16,25.834,15.346,NA,3,1.4,0,NA,0-49C,h,h,h,h,h,h,h,h,h,h,h
```

## Run QTL mapping

Due to package qtl could not invoke the kinship information,so the loaded data need to be converted into qtl2 format, and perform qtl mapping.

Convert format:
```r
library(qtl)
library(qtl2)
library(pheatmap)
#risib means RIL
mydata <- read.cross(format="csv", file="hanzhang_data/qtl_format_snp.csv",
                     na.strings = "NA",genotypes = c("a", "b", "h"),crosstype="risib")
summary(mydata)
mydata_qtl2=convert2cross2(mydata)
```

Calculate kinship:
```r
pr <- calc_genoprob(mydata_qtl2, error_prob=0.0001, cores=4,quiet=F)
kinship=calc_kinship(pr)
```

Perform QTL mapping:
```r
out_pg <- scan1(pr, mydata_qtl2$pheno, kinship)
plot(out_pg, mydata_qtl2$gmap, lodcolumn=1, col="tomato")
plot(out_pg, mydata_qtl2$gmap, lodcolumn=2, col="skyblue",add=T)
plot(out_pg, mydata_qtl2$gmap, lodcolumn=3, col="lightgreen",add=T)
plot(out_pg, mydata_qtl2$gmap, lodcolumn=4, col="orange",add=T)
legend("topleft", lwd=2, col=c("tomato","skyblue","lightgreen","orange"), c(names(mydata_qtl2$pheno[1,])[1:4]), bg="gray90", lty=c(1,1,2))
```
Get the permutated threshold of each phe:

```r
operm <- scan1perm(pr, mydata_qtl2$pheno,kinship=kinship, n_perm=1000, cores=10)
summary(operm)
```


