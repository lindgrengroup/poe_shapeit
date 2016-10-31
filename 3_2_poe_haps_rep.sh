#!/bin/bash
#$ -cwd -V
#$ -j y
#$ -P drong.prjc -q short.qc
#$ -t 454
#$ -N HM_POE
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

require(foreach)

i=$SGE_TASK_ID


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

chd.exp <- expression[,children]
##pat <- haps[,paste(children,"F_A",sep="")]
##mat <- haps[,paste(children,"M_A",sep="")]

colnames(haps)

pat <- haps[,paste(children,"_FT",sep="")]
mat <- haps[,paste(children,"_MT",sep="")]

require(foreach)

probes <- rownames(chd.exp)
snps<- rownames(mat)

dim(chd.exp)
dim(mat)



summary(lm(unlist(chd.exp[1,]) ~ unlist(mat[1,])))

mat.lm <-   foreach(j=1:nrow(mat),.combine=rbind) %do% {
 if(j%%1000==0) {print(j)}
    out <- try(coefficients(summary(lm(unlist(chd.exp[i,]) ~ unlist(mat[j,]))))[2,])
   if(!inherits(out,"try-error")){
     out
  } else {
    rep(NA,4)
   }
}


save(mat.lm,snps,probes,file=paste("${OUTPUT}/mat_",i,".R.RData",sep=""))

pat.lm <- foreach(j=1:nrow(pat),.combine=rbind) %do% {
 if(j%%100==0) {print(j)}
    out <-try(coefficients(summary(lm(unlist(chd.exp[i,]) ~ unlist(pat[j,]))))[2,])
  if(!inherits(out,"try-error")){
     out
  } else {
    rep(NA,4)
   }

}

save(pat.lm,snps,probes,file=paste("${OUTPUT}/pat_",i,".R.RData",sep=""))


add.lm <-   foreach(j=1:nrow(pat),.combine=rbind) %do% {
 if(j%%100==0) {print(j)}
    out<-try(coefficients(summary(lm(unlist(chd.exp[i,]) ~ (unlist(pat[j,]) + unlist(mat[j,])))))[2,])

  if(!inherits(out,"try-error")){
     out
  } else {
    rep(NA,4)
   }

}

save(add.lm,snps,probes,file=paste("${OUTPUT}/add_",i,".R.RData",sep=""))

EOF
