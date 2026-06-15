## NOTE: data/Data_extraction_2_v7_LeeAdded.xlsx already includes Lee 2021 in patient_rr.
## This script documents the original sensitivity addition (pre-Lee); the main analysis
## is reproduced via R/nma_analysis.R.

## Sensitivity analysis: add Lee AY 2021 (arbitrated-reference, AI-vs-human) to patient-level DTA-NMA
## mode: base | lee735 | leeimg
args <- commandArgs(trailingOnly=TRUE); mode <- args[1]
N_CHAINS<-4; N_ADAPT<-1000; N_BURN<-4000; N_KEEP<-4000; N_THIN<-1
XL <- "data/Data_extraction_2_v7_LeeAdded.xlsx"
OUT <- "results"
# jags.moddir auto-detected by rjags
suppressMessages({library(readxl); library(rjags); library(coda); library(parallel)})
set.seed(20260615)
model_string <- "
model {
 for(i in 1:nObs){
  tp[i]~dbin(pi[i,1], pos[i]); tn[i]~dbin(pi[i,2], neg[i])
  logit(pi[i,1]) <- mu[i,1]; logit(pi[i,2]) <- mu[i,2]
  MU[i,1] <- threshold.sens[test[i],threshold[i]] + study.re.sens[s[i]] + test.re.sens[s[i],test[i]]
  MU[i,2] <- threshold.spec[test[i],threshold[i]] + study.re.spec[s[i]] + test.re.spec[s[i],test[i]]
  mu[i,1:2] ~ dmnorm(MU[i,], Omega[,]) }
 for(k in 1:ntest){ sens[k,1]<- exp(threshold.sens[k,1])/(1+exp(threshold.sens[k,1]))
  spec[k,1]<- exp(threshold.spec[k,1])/(1+exp(threshold.spec[k,1]))
  threshold.sens[k,1]~dnorm(0,0.05); threshold.spec[k,1]~dnorm(0,0.05) }
 for(k in 1:ns){ study.re.sens[k]~dnorm(0,taustudysens); study.re.spec[k]~dnorm(0,taustudyspec)
  for(l in 1:ntest){ test.re.sens[k,l]~dnorm(0,tautestsens); test.re.spec[k,l]~dnorm(0,tautestspec) } }
 taustudysens<-pow(SDstudysens,-2); SDstudysens~dunif(0,2)
 taustudyspec<-pow(SDstudyspec,-2); SDstudyspec~dunif(0,2)
 tautestsens<-pow(SDtestsens,-2); SDtestsens~dunif(0,2)
 tautestspec<-pow(SDtestspec,-2); SDtestspec~dunif(0,2)
 Omega[1:2,1:2]<-inverse(Sigma.sq[,]); for(m in 1:2){ Sigma.sq[m,m]<-pow(sd[m],2) }
 for(i in 1:2){ for(j in (i+1):2){ Sigma.sq[i,j]<-rho[i,j]*sd[i]*sd[j]; Sigma.sq[j,i]<-Sigma.sq[i,j] } }
 for(m in 1:2){ sd[m]~dunif(0,1) }
 for(i in 1:2){ for(j in (i+1):2){ g[j,i]<-0; a[i,j]~dunif(0,3.1415); rho[i,j]<-inprod(g[,i],g[,j]) } }
 g[1,1]<-1; g[1,2]<-cos(a[1,2]); g[2,2]<-sin(a[1,2])
 for(k in 1:ntest){ dsens[k]<-threshold.sens[k,1]; dspec[k]<-threshold.spec[k,1] }
 for(k in 1:totaltest){ youden[k]<-exp(dsens[k])/(1+exp(dsens[k]))+exp(dspec[k])/(1+exp(dspec[k]))-1 }
 for(k in 1:totaltest){ for(j in 1:totaltest){ comparison[k,j]<-step(youden[j]-youden[k]) } }
 for(k in 1:totaltest){ rkyouden[k]<-sum(comparison[k,]) }
 for(k in 1:totaltest){ for(m in 1:totaltest){ youdensucra[k,m]<-equals(rkyouden[k],m) } } }
"
osa <- as.data.frame(read_excel(XL, sheet="patient_rr"))
osa <- osa[!is.na(osa$tp1)&!is.na(osa$fp1)&!is.na(osa$fn1)&!is.na(osa$tn1),
           c("authoryr","test","dataset_id","tp1","fn1","tn1","fp1")]
add <- function(df, ar, tp,fn,tn,fp,test) rbind(df, data.frame(authoryr="Lee et al. 2021",
        test=test, dataset_id="lee_2021", tp1=tp, fn1=fn, tn1=tn, fp1=fp))
if(mode=="lee735"){            # encounter-scaled N~735, balanced ~50% prevalence (pos367/neg368)
  osa<-add(osa,, 302,65,310,58,"human_only"); osa<-add(osa,, 295,72,299,69,"ai_only")
} else if(mode=="leeimg"){     # image-level back-calc (pos2810/neg2818)
  osa<-add(osa,, 2310,500,2377,441,"human_only"); osa<-add(osa,, 2261,549,2291,527,"ai_only")
}
osa$test_num <- as.numeric(factor(osa$test, levels=c("ai_only","human_only","human_with_ai")))
osa$study_num<- as.numeric(as.factor(osa$dataset_id)); n<-nrow(osa)
data_list<-list(ns=length(unique(osa$study_num)),ntest=3,totaltest=3,nObs=n,id=1:n,
  s=osa$study_num,test=osa$test_num,threshold=rep(1,n),
  tp=osa$tp1,tn=osa$tn1,pos=osa$tp1+osa$fn1,neg=osa$fp1+osa$tn1)
cat(sprintf("[%s] tables=%d groups=%d  arms: %s\n",mode,n,data_list$ns,
   paste(names(table(osa$test)),table(osa$test),sep="=",collapse=" ")))
rngs<-c("base::Wichmann-Hill","base::Marsaglia-Multicarry","base::Super-Duper","base::Mersenne-Twister")
mk<-function(sd0,c) list(SDstudysens=sd0,SDstudyspec=sd0,SDtestsens=sd0,SDtestspec=sd0,
   a=structure(c(NA,1),dim=c(1,2)),sd=c(0.5,0.5),.RNG.name=rngs[c],.RNG.seed=1000+c)
inits<-lapply(1:N_CHAINS,function(c) mk(c(1.0,0.8,1.2,0.5)[c],c))
params<-c("sens","spec","youden","youdensucra")
run<-function(c){ # jags.moddir auto-detected by rjags
  suppressMessages(library(rjags))
  jm<-jags.model(textConnection(model_string),data=data_list,inits=inits[[c]],n.chains=1,n.adapt=N_ADAPT,quiet=TRUE)
  update(jm,N_BURN); coda.samples(jm,params,n.iter=N_KEEP,thin=N_THIN)[[1]] }
cl<-mclapply(1:N_CHAINS,run,mc.cores=min(N_CHAINS,detectCores()))
ok<-sapply(cl,function(x) inherits(x,"mcmc")); if(!all(ok)){print(cl[!ok]);stop("chain fail")}
M<-as.matrix(as.mcmc.list(cl)); arms<-c("DL only","Human only","Human + DL")
per<-do.call(rbind,lapply(1:3,function(k) data.frame(arm=arms[k],
  sens=round(mean(M[,sprintf("sens[%d,1]",k)]),3),
  spec=round(mean(M[,sprintf("spec[%d,1]",k)]),3),
  youden=round(mean(M[,sprintf("youden[%d]",k)]),3),
  youden_lo=round(quantile(M[,sprintf("youden[%d]",k)],.025),3),
  youden_hi=round(quantile(M[,sprintf("youden[%d]",k)],.975),3))))
K<-3; Prank<-sapply(1:3,function(k) sapply(1:3,function(m) mean(M[,sprintf("youdensucra[%d,%d]",k,m)])))
rk<-do.call(rbind,lapply(1:3,function(k){pr<-Prank[,k]; data.frame(arm=arms[k],
  P_best=round(pr[1],3), SUCRA=round(sum(sapply(1:K,function(m)(K-m)/(K-1)*pr[m])),3),
  mean_rank=round(sum((1:K)*pr),3))}))
cat("\n== PER-ARM ==\n"); print(per,row.names=FALSE)
cat("\n== RANKING ==\n"); print(rk,row.names=FALSE)
write.csv(per,file.path(OUT,paste0(mode,"_per_arm.csv")),row.names=FALSE)
write.csv(rk, file.path(OUT,paste0(mode,"_ranking.csv")),row.names=FALSE)
