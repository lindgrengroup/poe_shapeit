#!/bin/bash
#$ -cwd -V
#$ -j y
#$ -P drong.prjc -q short.qc
#$ -t 740-745
#$ -N miRNAPerm

echo "************************************************************"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "************************************************************"



##maternal 
##paternal 11 266 733 30  85 246 247 248 249 250 251 345 587 618 658 679 925

## for running the script on a head node, provide index as first argument
if [ -z "$SGE_TASK_ID" ]; then SGE_TASK_ID=$1; fi



module load R

cd /well/lindgren/alexd/miRNA/HAPMAP_SHAPEIT

source 0_config.sh

R --vanilla <<EOF

load("/well/lindgren/alexd/miRNA/annotation.RData")
load("/well/lindgren/alexd/miRNA/expression.RData")

b37 <- read.table("/well/lindgren/alexd/miRNA/annotation_b37.bed")
colnames(b37) <- c("chr","posb37","pos2b37","ID")

anno <- cbind(anno,b37)

expression <- mirna

require(data.table)
require(foreach)

## haps <- fread("/well/lindgren/alexd/HAPMAP/SHAPEIT/MERGED_combined.haps", data.table=F) 
## haps.chr <- as.numeric(haps[,1])
## haps.pos <- as.numeric(haps[,3])
## haps.rs <- as.character(haps[,2])
## save(haps.chr,haps.pos,haps.rs,file="haps.RData")

load("haps.RData")

anno <- anno[!(anno[,"chr"] %in% "X"), ]#& !(duplicated(anno[,"probeID"])),]


anno[,"chr"] <- substr(anno[,"chr"],4,6)

require(foreach)

perm=$SGE_TASK_ID


require(data.table)


                                        #load("/well/lindgren/alexd/miRNA/HAPMAP_ALSPAC_METHOD/INPUT/hapmap3_r3_b36_fwd.consensus.qc.poly.CEU.phased.combined.qc.RData")
                                        #load("/well/lindgren/alexd/miRNA/expression.RData")
                                        #load("/well/lindgren/alexd/miRNA/annotation.RData")

                                        #load("${INPUT}.phased.combined.RData")

haps <- fread("${INPUT}.haps",data.table=F)

load("$EXP")
                                        #load("$ANNO")


                                        #haps <- phased.combined
expression <- mirna

                                        #sample <-  fread("/well/lindgren/alexd/miRNA/HAPMAP_ALSPAC_METHOD/INPUT/hapmap3_r3_b36_fwd.consensus.qc.poly.CEU.fam", data.table=F)

sample <-  fread("${INPUT}.sample", data.table=F)
sample <- sample[-1,]

colnames(haps) <- sample[,2]

rownames(sample) <- sample[,2]


colnames(haps)[6:ncol(haps)] <- paste(rep(sample[,2],each=2),c("MT","FT"),sep="_")


children <- sample[sample[,4]!=0,2]

                                        #children <- unique(substr(colnames(haps),1,7))

                                        #expression[,Reduce(intersect,list(children,colnames(expression),sample[children,4], sample[children,5]))]
colnames(expression)
children

sample[children,4]

sum(children%in% colnames(expression))

sum(sample[children,4]%in% colnames(expression))
sum(sample[children,5]%in% colnames(expression))

children <- children[sample[children,4]%in% colnames(expression) & sample[children,5]%in% colnames(expression)]
children



                                        #colnames(expression) <-  link[colnames(expression),1]
colnames(haps)
children
sample

## bim <- foreach(k=1:22,.combine=rbind) %do%{
##   print(k)
##   a <- fread(paste("/well/lindgren/alexd/miRNA/HAPMAP_ALSPAC_METHOD/INPUT/hapmap3_r3_b36_fwd.consensus.qc.poly.CEU_R_",k,".bim",sep=""),data.table=F)
##   a
## }

## head(bim)

##rownames(haps) <- bim[,2]

head(rownames(haps))
                                        #rep.rs <- intersect(c("rs1268538","rs893184","rs16944141","rs1268538","rs893184","rs11064410","rs7955666","rs6033","rs2100700","rs4930081","rs11631662","rs11612312","rs9900527","rs5993624","rs11170164"),rownames(haps))

colnames(expression)

chd.exp <- expression[,sample(children)]
##pat <- haps[,paste(children,"F_A",sep="")]
##mat <- haps[,paste(children,"M_A",sep="")]

colnames(haps)

pat <- haps[,paste(children,"_FT",sep="")]
mat <- haps[,paste(children,"_MT",sep="")]

add <- pat + mat

require(foreach)

probes <- rownames(chd.exp)
snps<- rownames(mat)

dim(chd.exp)
dim(mat)



summary(lm(unlist(chd.exp[1,]) ~ unlist(mat[1,])))


cis.pat <- foreach(i=1:nrow(anno),  .errorhandling ='pass', .combine=rbind) %do% {
  print(i)

  number <- as.character(anno[i,2])

  file <- which(rownames(expression) %in% number)

  chr <- as.numeric(as.character(anno[i,"chr"]))
  pos <- as.numeric(as.character(anno[i,"posb37"]))

  

  
  index <- which(haps.chr==chr & haps.pos > (pos - 500000) & (haps.pos < pos + 500000))

  print(length(index))

  pat.lm <- foreach(j=index,.combine=rbind) %do% {
    if(j%%100==0) {print(j)}
    out <-try(coefficients(summary(lm(unlist(chd.exp[i,]) ~ unlist(pat[j,]))))[2,])
    if(!inherits(out,"try-error")){
      out
    } else {
      rep(NA,4)
    }

  }

  print(head(pat.lm))

  if(!inherits(sub,"try-error")) {
    

    sub <- pat.lm

    rownames(pat.lm) <- haps.rs[index]
    pat.lm <- cbind(haps.ref[index],pat.lm)
  

    sub <- try(cbind(rep(anno[i,1],nrow(sub)),sub))
    sub

  } else {
    rep(NA,7)
  }

}




cis.mat <- foreach(i=1:nrow(anno),  .errorhandling ='pass', .combine=rbind) %do% {
  print(i)

  number <- as.character(anno[i,2])

  file <- which(rownames(expression) %in% number)

  chr <- as.numeric(as.character(anno[i,"chr"]))
  pos <- as.numeric(as.character(anno[i,"posb37"]))

  

  
  index <- which(haps.chr==chr & haps.pos > (pos - 500000) & (haps.pos < pos + 500000))

  print(length(index))

  mat.lm <- foreach(j=index,.combine=rbind) %do% {
    if(j%%100==0) {print(j)}
    out <-try(coefficients(summary(lm(unlist(chd.exp[i,]) ~ unlist(mat[j,]))))[2,])
    if(!inherits(out,"try-error")){
      out
    } else {
      rep(NA,4)
    }

  }

  print(head(mat.lm))

  if(!inherits(sub,"try-error")) {
    

    sub <- mat.lm

    rownames(mat.lm) <- haps.rs[index]
    mat.lm <- cbind(haps.ref[index],mat.lm)
  

    sub <- try(cbind(rep(anno[i,1],nrow(sub)),sub))
    sub

  } else {
    rep(NA,7)
  }

}

head(cis.mat)

cis.add <- foreach(i=1:nrow(anno),  .errorhandling ='pass', .combine=rbind) %do% {
  print(i)

  number <- as.character(anno[i,2])

  file <- which(rownames(expression) %in% number)

  chr <- as.numeric(as.character(anno[i,"chr"]))
  pos <- as.numeric(as.character(anno[i,"posb37"]))

  

  
  index <- which(haps.chr==chr & haps.pos > (pos - 500000) & (haps.pos < pos + 500000))

  print(length(index))

  add.lm <- foreach(j=index,.combine=rbind) %do% {
    if(j%%100==0) {print(j)}
    out <-try(coefficients(summary(lm(unlist(chd.exp[i,]) ~ unlist(add[j,]))))[2,])
    if(!inherits(out,"try-error")){
      out
    } else {
      rep(NA,4)
    }

  }

  print(head(add.lm))

  if(!inherits(sub,"try-error")) {
    

    sub <- add.lm

    rownames(add.lm) <- haps.rs[index]
    add.lm <- cbind(haps.ref[index],add.lm)
  

    sub <- try(cbind(rep(anno[i,1],nrow(sub)),sub))
    sub

  } else {
    rep(NA,7)
  }

}

head(cis.add)

save(cis.mat,cis.pat,cis.add,file=paste("PERM/HAPMAP_SHAPEIT_cis_poe_",perm,".RData",sep=""))


EOF
