args=(commandArgs(TRUE)) ## retrieve variable from the job call if any
print(args)
if( length(args)>0){
  for(k in 1:length(args)){
    eval(parse(text=args[[k]]))
  }
}
print(CHR)
if( !exists("CHR") ) print("arg missing") ### check if variable CHR exists and if not set it to 22 by default

library(snpStats) # for reading plink files
library(nnet) # for fitting the models

#setwd("/scratch/o/oespinga/kmin940/OAI7/rm")
setwd("/scratch/o/oespinga/kmin940/OAI11/TOPMED/rm")

### read genotype and phenotype data
gen <- read.plink(paste0('rm_chr',CHR)) ### note that variable CHR may come from the script call
phen <- read.table('~/OAI11/gwas/pheno_v11.txt',header=T)
cov = read.table('~/OAI11/gwas/cov_v11.txt',header=T)

### you may need to install  some of these packages
library(parallel)
library(doParallel)
library(iterators)
library(foreach)

cl <- makeCluster(detectCores()-1, type="SOCK", outfile="")

registerDoParallel(cl)

a <- foreach(G=iter(gen$geno, by='column'), i=icount(), .combine='rbind', .multicombine=T, .inorder=F, .verbose=T, .packages = c("nnet")) %dopar% {
  #print(i)

  #fit the alternative model:
  md1 <- multinom(as.factor(phen$Pheno) ~ as.numeric(G)+as.numeric(cov$PC1)+as.factor(cov$SEX)+as.numeric(cov$AGE)+as.numeric(cov$BMI))

  #fit the null model
  md0 <- multinom(as.factor(phen$Pheno) ~ 1+as.numeric(cov$PC1)+as.factor(cov$SEX)+as.numeric(cov$AGE)+as.numeric(cov$BMI))

  #do the LRT
  test <- anova(md1,md0)

  ## extract some additional statistics
  coefs <- coef(md1)[,2]
  se <- summary(md1)$standard.errors[,2]
  z <- coefs/se
  pvals <- pchisq(z^2,1,lower.tail = F)
  maf=mean(3-as.numeric(G))/2

  # write the output (you may want to check that this unlist works here)
  unlist(c(i,exp(coef(md1)[,2]),-0.5*test$"Resid. Dev",test$P[2],se,z,pvals,maf),use.names=F)
}

stopCluster(cl)

a1 = data.frame(gen$map[,-3],a[order(a[,1]), ])
##first 5 columns: map info, 6:index, 7,8,9:odds ratios, 10,11:log likelihood ratio for md1 ,md0, 12:p-value for multi model , 13,14,15:standard errors, 16,17,18:z-values, 19,20,21:p-values, 22:maf

### save results to text file (this could also be saved to R binary file to save space)
write.table(a1,file=paste0("/scratch/o/oespinga/kmin940/OAI11/gwas/nnet_chr",CHR),row.names=F,quote=F,col.names=F)
