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

## for running the script on a head node, provide index as first argument
if [ -z "$SGE_TASK_ID" ]; then SGE_TASK_ID=$1; fi

module load R

source 0_config.sh

R --vanilla <<EOF

require(foreach)

i=$SGE_TASK_ID


require(data.table)

haps <- fread("${INPUT}.haps",data.table=F)

load("$EXP")

expression <- mirna


sample <-  fread("${INPUT}.sample", data.table=F)
sample <- sample[-1,]

colnames(haps) <- sample[,2]

rownames(sample) <- sample[,2]


colnames(haps)[6:ncol(haps)] <- paste(rep(sample[,2],each=2),c("MT","FT"),sep="_")


children <- sample[sample[,4]!=0,2]

colnames(expression)
children

sample[children,4]

sum(children%in% colnames(expression))

sum(sample[children,4]%in% colnames(expression))
sum(sample[children,5]%in% colnames(expression))

#children <- children[sample[children,4]%in% colnames(expression) & sample[children,5]%in% colnames(expression)]
children



colnames(haps)
children
sample

head(rownames(haps))

colnames(expression)

chd.exp <- expression[,children]


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
