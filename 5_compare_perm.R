
require(foreach)

dir("PERM")

alpha <- 0.05/1417081

load("HAPMAP_SHAPEIT_cis_poe_1.RData")

cis.mat <- cis.mat[which((unlist(lapply(cis.mat[,4],is.numeric)))),]

cis.mat.org <- cis.mat[which((unlist(lapply(cis.mat[,4],is.numeric)))),]
cis.pat.org <- cis.pat[which((unlist(lapply(cis.pat[,4],is.numeric)))),]
cis.add.org <- cis.add[which((unlist(lapply(cis.add[,4],is.numeric)))),]


cis.obs <- c(sum(as.numeric(cis.mat[,5])<alpha,na.rm=T),sum(as.numeric(cis.pat[,5])<alpha,na.rm=T),sum(as.numeric(cis.add[,5])<alpha,na.rm=T))

cis.mat <- cis.mat[order(unlist(cis.mat[,4])),]
cis.pat <- cis.pat[order(unlist(cis.pat[,4])),]

cis.add <- cis.add[order(unlist(cis.add[,4])),]

cis.perm <- foreach(i=1:1000,.combine=rbind) %do% {
  print(i)
  load(paste("PERM/HAPMAP_SHAPEIT_cis_poe_",i,".RData",sep=""))
  cis.mat <- cis.mat[which((unlist(lapply(cis.mat[,4],is.numeric)))),]
  cis.pat <- cis.pat[which((unlist(lapply(cis.pat[,4],is.numeric)))),]
  cis.add <- cis.add[which((unlist(lapply(cis.add[,4],is.numeric)))),]
  print(head(cis.mat,2))
  c(sum(cis.mat[,5]<alpha,na.rm=T),sum(cis.pat[,5]<alpha,na.rm=T),sum(cis.add[,5]<alpha,na.rm=T))
}
save(cis.perm,file="cis.perm.RData")

cis.maxT <- foreach(i=1:1000,.combine=rbind) %do% {
  print(i)
  load(paste("PERM/HAPMAP_SHAPEIT_cis_poe_",i,".RData",sep=""))
  cis.mat <- cis.mat[which((unlist(lapply(cis.mat[,4],is.numeric)))),]
  cis.pat <- cis.pat[which((unlist(lapply(cis.pat[,4],is.numeric)))),]
  cis.add <- cis.add[which((unlist(lapply(cis.add[,4],is.numeric)))),]
  print(head(cis.mat,2))
  c(max(abs(unlist(cis.mat[,4])),na.rm=T),max(abs(unlist(cis.pat[,4])),na.rm=T),max(unlist(cis.add[,4]),na.rm=T))
}

mat.c.total <- 0
pat.c.total <- 0
add.c.total <- 0

cis.all <- foreach(i=1:1000,.combine=rbind) %do% {
  print(i)
  load(paste("PERM/HAPMAP_SHAPEIT_cis_poe_",i,".RData",sep=""))
  cis.mat <- cis.mat[which((unlist(lapply(cis.mat[,4],is.numeric)))),]
  cis.pat <- cis.pat[which((unlist(lapply(cis.pat[,4],is.numeric)))),]
  cis.add <- cis.add[which((unlist(lapply(cis.add[,4],is.numeric)))),]
  mat.c <- abs(unlist(cis.mat[,4])) >= abs(unlist(cis.mat.org[,4]))
  pat.c <- abs(unlist(cis.pat[,4])) >= abs(unlist(cis.pat.org[,4]))
  add.c <- abs(unlist(cis.add[,4])) >= abs(unlist(cis.add.org[,4]))
  mat.c.total <- mat.c.total + mat.c
  pat.c.total <- pat.c.total + pat.c
  add.c.total <- add.c.total + add.c
  c(max(abs(unlist(cis.mat[,4])),na.rm=T),max(abs(unlist(cis.pat[,4])),na.rm=T),max(unlist(cis.add[,4]),na.rm=T))
}

save(cis.all,file="cis.all.RData")

sum(cis.all[,2]>=cis.obs[2])/1000

sum(max(abs(unlist(cis.pat[,4])),na.rm=T) < cis.maxT[,1])/1000

save(cis.maxT,file="cis.maxT.RData")

