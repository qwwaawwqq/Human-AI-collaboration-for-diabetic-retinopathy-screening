#!/usr/bin/env Rscript
## ============================================================================
##  Deep learning vs human readers vs human-AI collaboration for diabetic
##  retinopathy screening: Bayesian bivariate DTA network meta-analysis.
##
##  SINGLE INTEGRATED ANALYSIS SCRIPT.  Run from the repository root:
##      Rscript R/dr_nma_analysis.R
##  Quick smoke-test (tiny MCMC, ~1 min) instead of the full 50k run:
##      DR_QUICK=1 Rscript R/dr_nma_analysis.R
##
##  Reproduces, writing all output to results/:
##    * pooled Se/Sp/Youden, SUCRA/ranking, pairwise contrasts, heterogeneity,
##      convergence              (patient & eye)                  -> *_per_arm/_pairwise/_ranking/_het/_conv.csv
##    * unit-of-analysis + leave-one-out (drop Li 2024b) sensitivity (eye)  -> eye_sensitivity.csv
##    * expected outcomes per 10,000 screened                     -> *_prevalence.csv
##    * univariate categorical meta-regression (26 covariates)    -> metareg_*.csv
##    * Figure 2 (forest), Figure 3 (HSROC), Figure 4 (multivariate MR)  -> Figure*.tiff
##  Requires: R >=4.3 with rjags, coda, readxl, parallel; JAGS >=4.3.
## ============================================================================
suppressMessages({library(readxl); library(rjags); library(coda); library(parallel)})
set.seed(20260612)
XL  <- "data/Data_extraction_2_v7_LeeAdded.xlsx"
OUT <- "results"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
ARMS <- c("DL only","Human only","Human + DL")
QUICK <- nzchar(Sys.getenv("DR_QUICK"))
N_CHAINS <- if (QUICK) 2 else 4
N_ADAPT  <- if (QUICK) 200 else 1000
N_BURN   <- if (QUICK) 300 else 10000
N_KEEP   <- if (QUICK) 500 else 10000
N_THIN   <- if (QUICK) 1   else 10

## ---- (1) JAGS bivariate NMA model (no covariates) --------------------------
model_string <- "
model {
 for(i in 1:nObs){
  tp[i]~dbin(pi[i,1], pos[i]); tn[i]~dbin(pi[i,2], neg[i])
  logit(pi[i,1]) <- mu[i,1]; logit(pi[i,2]) <- mu[i,2]
  MU[i,1] <- threshold.sens[test[i],threshold[i]] + study.re.sens[s[i]] + test.re.sens[s[i],test[i]]
  MU[i,2] <- threshold.spec[test[i],threshold[i]] + study.re.spec[s[i]] + test.re.spec[s[i],test[i]]
  mu[i,1:2] ~ dmnorm(MU[i,], Omega[,])
 }
 for(k in 1:ntest){
  sens[k,1]<- exp(threshold.sens[k,1])/(1+exp(threshold.sens[k,1]))
  spec[k,1]<- exp(threshold.spec[k,1])/(1+exp(threshold.spec[k,1]))
  threshold.sens[k,1]~dnorm(0,0.05); threshold.spec[k,1]~dnorm(0,0.05)
 }
 for(k in 1:ns){
  study.re.sens[k]~dnorm(0,taustudysens); study.re.spec[k]~dnorm(0,taustudyspec)
  for(l in 1:ntest){ test.re.sens[k,l]~dnorm(0,tautestsens); test.re.spec[k,l]~dnorm(0,tautestspec) }
 }
 taustudysens <- pow(SDstudysens,-2); SDstudysens ~ dunif(0,2)
 taustudyspec <- pow(SDstudyspec,-2); SDstudyspec ~ dunif(0,2)
 tautestsens <- pow(SDtestsens,-2);  SDtestsens ~ dunif(0,2)
 tautestspec <- pow(SDtestspec,-2);  SDtestspec ~ dunif(0,2)
 Omega[1:2,1:2] <- inverse(Sigma.sq[,])
 for(m in 1:2){ Sigma.sq[m,m] <- pow(sd[m],2) }
 for(i in 1:2){ for(j in (i+1):2){ Sigma.sq[i,j] <- rho[i,j]*sd[i]*sd[j]; Sigma.sq[j,i] <- Sigma.sq[i,j] } }
 for(m in 1:2){ sd[m] ~ dunif(0,1) }
 for(i in 1:2){ for(j in (i+1):2){ g[j,i] <- 0; a[i,j] ~ dunif(0,3.1415); rho[i,j] <- inprod(g[,i], g[,j]) } }
 g[1,1] <- 1; g[1,2] <- cos(a[1,2]); g[2,2] <- sin(a[1,2])
 for(k in 1:ntest){ dsens[k] <- threshold.sens[k,1]; dspec[k] <- threshold.spec[k,1] }
 for(k in 1:totaltest){ youden[k] <- exp(dsens[k])/(1+exp(dsens[k])) + exp(dspec[k])/(1+exp(dspec[k])) - 1 }
 for(k in 1:totaltest){ for(j in 1:totaltest){ comparison[k,j] <- step(youden[j]-youden[k]) } }
 for(k in 1:totaltest){ rkyouden[k] <- sum(comparison[k,]) }
 for(k in 1:totaltest){ for(m in 1:totaltest){ youdensucra[k,m] <- equals(rkyouden[k], m) } }
}
"

## ---- (2) data prep ---------------------------------------------------------
prep <- function(sheet, group_col="dataset_id", one_per=FALSE, drop_author=NULL){
  osa <- as.data.frame(read_excel(XL, sheet=sheet))
  osa <- osa[!is.na(osa$tp1) & !is.na(osa$fp1) & !is.na(osa$fn1) & !is.na(osa$tn1), ]
  if (!is.null(drop_author)) osa <- osa[osa$authoryr != drop_author, ]
  if (one_per){                          # keep the largest 2x2 table per group x arm
    osa$.n <- osa$tp1+osa$fn1+osa$fp1+osa$tn1
    key <- paste(osa[[group_col]], osa$test, sep="||")
    ord <- order(key, -osa$.n); osa <- osa[ord, ]; key <- key[ord]   # reorder key WITH osa
    osa <- osa[!duplicated(key), ]
  }
  osa$test_num  <- as.numeric(factor(osa$test, levels=c("ai_only","human_only","human_with_ai")))
  osa$study_num <- as.numeric(as.factor(osa[[group_col]])); n <- nrow(osa)
  list(osa=osa, data=list(ns=length(unique(osa$study_num)), ntest=3, totaltest=3, nObs=n,
       s=osa$study_num, test=osa$test_num, threshold=rep(1,n),
       tp=osa$tp1, tn=osa$tn1, pos=osa$tp1+osa$fn1, neg=osa$fp1+osa$tn1))
}

## ---- (3) fit (parallel chains) ---------------------------------------------
fit_nma <- function(data_list){
  rngs <- c("base::Wichmann-Hill","base::Marsaglia-Multicarry","base::Super-Duper","base::Mersenne-Twister")
  mk <- function(sd0,c) list(SDstudysens=sd0,SDstudyspec=sd0,SDtestsens=sd0,SDtestspec=sd0,
                             a=structure(c(NA,1),dim=c(1,2)), sd=c(0.5,0.5), .RNG.name=rngs[c], .RNG.seed=1000+c)
  inits <- lapply(1:N_CHAINS, function(c) mk(c(1.0,0.8,1.2,0.5)[c], c))
  params <- c("sens","spec","youden","youdensucra","SDstudysens","SDstudyspec","SDtestsens","SDtestspec","sd","rho")
  one <- function(c){ suppressMessages(library(rjags))
    jm <- jags.model(textConnection(model_string), data=data_list, inits=inits[[c]], n.chains=1, n.adapt=N_ADAPT, quiet=TRUE)
    update(jm, N_BURN); coda.samples(jm, params, n.iter=N_KEEP, thin=N_THIN)[[1]] }
  cl <- mclapply(1:N_CHAINS, one, mc.cores=min(N_CHAINS, detectCores()))
  if(!all(sapply(cl, inherits, "mcmc"))) stop("chain failure")
  samp <- as.mcmc.list(cl); list(samp=samp, M=as.matrix(samp))
}

## ---- (4) summaries (per-arm, pairwise, ranking, heterogeneity, draws) ------
ranks_of <- function(M){  # SUCRA + P(best) from youdensucra
  Pr <- sapply(1:3, function(k) sapply(1:3, function(mm) mean(M[,sprintf("youdensucra[%d,%d]",k,mm)])))
  data.frame(arm=ARMS, P_best=Pr[1,], SUCRA=sapply(1:3, function(k) sum(((3-(1:3))/2)*Pr[,k])),
             mean_rank=sapply(1:3, function(k) sum((1:3)*Pr[,k])), row.names=NULL)
}
summarise <- function(fit, tag, data_list){
  M <- fit$M
  per <- do.call(rbind, lapply(1:3, function(k) data.frame(arm=ARMS[k],
    sens=mean(M[,sprintf("sens[%d,1]",k)]), sens_lo=quantile(M[,sprintf("sens[%d,1]",k)],.025), sens_hi=quantile(M[,sprintf("sens[%d,1]",k)],.975),
    spec=mean(M[,sprintf("spec[%d,1]",k)]), spec_lo=quantile(M[,sprintf("spec[%d,1]",k)],.025), spec_hi=quantile(M[,sprintf("spec[%d,1]",k)],.975),
    youden=mean(M[,sprintf("youden[%d]",k)]), youden_lo=quantile(M[,sprintf("youden[%d]",k)],.025), youden_hi=quantile(M[,sprintf("youden[%d]",k)],.975), row.names=NULL)))
  write.csv(per, file.path(OUT,paste0(tag,"_per_arm.csv")), row.names=FALSE)
  prs <- list(c(3,1),c(3,2),c(2,1))
  pw <- do.call(rbind, lapply(prs, function(p){ i<-p[1]; j<-p[2]
    dY<-M[,sprintf("youden[%d]",i)]-M[,sprintf("youden[%d]",j)]
    dSp<-M[,sprintf("spec[%d,1]",i)]-M[,sprintf("spec[%d,1]",j)]
    data.frame(comparison=paste(ARMS[i],"vs",ARMS[j]),
      dYouden=mean(dY), dY_lo=quantile(dY,.025), dY_hi=quantile(dY,.975), P_sup=mean(dY>0),
      dSpec=mean(dSp), dSp_lo=quantile(dSp,.025), dSp_hi=quantile(dSp,.975), row.names=NULL)}))
  write.csv(pw, file.path(OUT,paste0(tag,"_pairwise.csv")), row.names=FALSE)
  rk <- ranks_of(M); write.csv(rk, file.path(OUT,paste0(tag,"_ranking.csv")), row.names=FALSE)
  het <- data.frame(param=c("SDstudy_sens","SDstudy_spec","SDtest_sens","SDtest_spec"),
    mean=sapply(c("SDstudysens","SDstudyspec","SDtestsens","SDtestspec"), function(p) mean(M[,p])))
  write.csv(het, file.path(OUT,paste0(tag,"_het.csv")), row.names=FALSE)
  kc <- c(sprintf("sens[%d,1]",1:3),sprintf("spec[%d,1]",1:3),sprintf("youden[%d]",1:3))
  gd <- tryCatch(gelman.diag(fit$samp[,kc],autoburnin=FALSE,multivariate=FALSE)$psrf[,1], error=function(e) NA)
  write.csv(data.frame(param=kc, Rhat=if(length(gd)>1) round(gd,3) else NA), file.path(OUT,paste0(tag,"_conv.csv")), row.names=FALSE)
  saveRDS(list(sens=M[,sprintf("sens[%d,1]",1:3)], spec=M[,sprintf("spec[%d,1]",1:3)],
               youden=M[,sprintf("youden[%d]",1:3)], arms=ARMS, tag=tag), file.path(OUT,paste0(tag,"_draws.rds")))
  cat(sprintf("  [%s] collab Youden=%.3f  SUCRA=%.2f  P(best)=%.2f  (max Rhat=%.3f)\n",
      tag, per$youden[3], rk$SUCRA[3], rk$P_best[3], if(length(gd)>1) max(gd,na.rm=TRUE) else NA))
  invisible(list(per=per, pw=pw, rk=rk))
}

## ---- (5) per-10,000-screened outcomes -------------------------------------
prevalence <- function(tag, prevs=c(0.02,0.05,0.10,0.20)){
  d <- readRDS(file.path(OUT,paste0(tag,"_draws.rds"))); rows <- list()
  for (p in prevs) for (k in 1:3){ Se<-d$sens[,k]; Sp<-d$spec[,k]
    miss<-10000*p*(1-Se); fr<-10000*(1-p)*(1-Sp)
    rows[[length(rows)+1]] <- data.frame(tag=tag, prevalence=sprintf("%d%%",round(p*100)), arm=d$arms[k],
      missed=round(mean(miss)), missed_lo=round(quantile(miss,.025)), missed_hi=round(quantile(miss,.975)),
      false_ref=round(mean(fr)), fr_lo=round(quantile(fr,.025)), fr_hi=round(quantile(fr,.975))) }
  write.csv(do.call(rbind,rows), file.path(OUT,paste0(tag,"_prevalence.csv")), row.names=FALSE)
}

## ---- (6) univariate categorical meta-regression ---------------------------
mr_model <- "model{
 for(i in 1:nObs){ tp[i]~dbin(pi[i,1],pos[i]); tn[i]~dbin(pi[i,2],neg[i])
  logit(pi[i,1])<-mu[i,1]; logit(pi[i,2])<-mu[i,2]
  MU[i,1]<-threshold.sens[test[i]]+study.re.sens[s[i]]+test.re.sens[s[i],test[i]]+beta_Se[cov_x[i]]
  MU[i,2]<-threshold.spec[test[i]]+study.re.spec[s[i]]+test.re.spec[s[i],test[i]]+beta_Sp[cov_x[i]]
  mu[i,1:2]~dmnorm(MU[i,],Omega[,]) }
 beta_Se[1]<-0; beta_Sp[1]<-0
 for(k in 2:n_lev){ beta_Se[k]~dnorm(0,0.05); beta_Sp[k]~dnorm(0,0.05) }
 for(k in 1:ntest){ threshold.sens[k]~dnorm(0,0.05); threshold.spec[k]~dnorm(0,0.05) }
 for(k in 1:ns){ study.re.sens[k]~dnorm(0,taus1); study.re.spec[k]~dnorm(0,taus2)
  for(l in 1:ntest){ test.re.sens[k,l]~dnorm(0,taut1); test.re.spec[k,l]~dnorm(0,taut2) } }
 taus1<-pow(SDs1,-2);SDs1~dunif(0,2); taus2<-pow(SDs2,-2);SDs2~dunif(0,2)
 taut1<-pow(SDt1,-2);SDt1~dunif(0,2); taut2<-pow(SDt2,-2);SDt2~dunif(0,2)
 Omega[1:2,1:2]<-inverse(Sig[,]); for(j in 1:2){Sig[j,j]<-pow(sd[j],2)}
 for(i in 1:2){for(j in (i+1):2){Sig[i,j]<-rho[i,j]*sd[i]*sd[j]; Sig[j,i]<-Sig[i,j]}}
 for(j in 1:2){sd[j]~dunif(0,1)}
 for(i in 1:2){for(j in (i+1):2){g[j,i]<-0; ang[i,j]~dunif(0,3.1415); rho[i,j]<-inprod(g[,i],g[,j])}}
 g[1,1]<-1; g[1,2]<-cos(ang[1,2]); g[2,2]<-sin(ang[1,2]) }"
COVS <- c("readerexpertise","study_design","multicenter","economic","healthcare","external","dr",
 "othereyetarget","pupil","ctype","camera","criteria","analysis_type","rstandard","disagree",
 "rimage","field","vendor","commercial","certified","quality_check","architecture","region","aigradable","hgradable")
metareg <- function(sheet, suffix){
  osa <- as.data.frame(read_excel(XL,sheet=sheet))
  osa <- osa[!is.na(osa$tp1)&!is.na(osa$fp1)&!is.na(osa$fn1)&!is.na(osa$tn1),]
  osa$test_num<-as.numeric(factor(osa$test,levels=c("ai_only","human_only","human_with_ai")))
  osa$study_num<-as.numeric(as.factor(osa$dataset_id))
  base<-list(ns=length(unique(osa$study_num)),ntest=3,nObs=nrow(osa),s=osa$study_num,test=osa$test_num,
    tp=osa$tp1,tn=osa$tn1,pos=osa$tp1+osa$fn1,neg=osa$fp1+osa$tn1)
  mkinit<-function(c) list(SDs1=1,SDs2=0.8,SDt1=1,SDt2=0.5,sd=c(0.5,0.5),
    ang=structure(c(NA,1),dim=c(1,2)),.RNG.name=c("base::Wichmann-Hill","base::Super-Duper")[c],.RNG.seed=99+c)
  ni <- if (QUICK) 500 else 1500
  res<-data.frame(); cv_list <- if (QUICK) head(COVS,3) else COVS
  for(cv in cv_list){ r<-osa[[cv]]; if(is.null(r)||all(is.na(r))) next
    r[is.na(r)]<-"NR"; f<-as.factor(as.character(r)); ft<-sort(table(f),decreasing=TRUE)
    f<-relevel(f,ref=names(ft)[1]); lev<-levels(f); nl<-length(lev); if(nl<2) next
    dl<-c(base,list(cov_x=as.numeric(f),n_lev=nl))
    ok<-tryCatch({ jm<-jags.model(textConnection(mr_model),data=dl,inits=list(mkinit(1),mkinit(2)),n.chains=2,n.adapt=500,quiet=TRUE)
      update(jm,ni); sm<-coda.samples(jm,c("beta_Se","beta_Sp"),n.iter=ni); TRUE},error=function(e) FALSE)
    if(!ok) next; M<-as.matrix(sm)
    for(k in 2:nl){ bs<-M[,sprintf("beta_Se[%d]",k)]; bp<-M[,sprintf("beta_Sp[%d]",k)]
      res<-rbind(res,data.frame(level=suffix,cov=cv,contrast=paste0(lev[k]," vs ",lev[1]),
        Se_OR=round(exp(mean(bs)),3),Se_lo=round(exp(quantile(bs,.025)),3),Se_hi=round(exp(quantile(bs,.975)),3),
        Se_sig=(quantile(bs,.025)>0|quantile(bs,.975)<0),
        Sp_OR=round(exp(mean(bp)),3),Sp_lo=round(exp(quantile(bp,.025)),3),Sp_hi=round(exp(quantile(bp,.975)),3),
        Sp_sig=(quantile(bp,.025)>0|quantile(bp,.975)<0))) } }
  write.csv(res, file.path(OUT,paste0("metareg_",suffix,".csv")), row.names=FALSE); res
}

## ---- (6b) multivariate meta-regression (design matrix of significant covariates) ----
mv_model <- "model{
 for(i in 1:nObs){ tp[i]~dbin(pi[i,1],pos[i]); tn[i]~dbin(pi[i,2],neg[i])
  logit(pi[i,1])<-mu[i,1]; logit(pi[i,2])<-mu[i,2]
  MU[i,1]<-threshold.sens[test[i]]+study.re.sens[s[i]]+test.re.sens[s[i],test[i]]+inprod(X[i,],bSe[])
  MU[i,2]<-threshold.spec[test[i]]+study.re.spec[s[i]]+test.re.spec[s[i],test[i]]+inprod(X[i,],bSp[])
  mu[i,1:2]~dmnorm(MU[i,],Omega[,]) }
 for(p in 1:P){ bSe[p]~dnorm(0,0.05); bSp[p]~dnorm(0,0.05) }
 for(k in 1:ntest){ threshold.sens[k]~dnorm(0,0.05); threshold.spec[k]~dnorm(0,0.05) }
 for(k in 1:ns){ study.re.sens[k]~dnorm(0,taus1); study.re.spec[k]~dnorm(0,taus2)
  for(l in 1:ntest){ test.re.sens[k,l]~dnorm(0,taut1); test.re.spec[k,l]~dnorm(0,taut2) } }
 taus1<-pow(SDs1,-2);SDs1~dunif(0,2); taus2<-pow(SDs2,-2);SDs2~dunif(0,2)
 taut1<-pow(SDt1,-2);SDt1~dunif(0,2); taut2<-pow(SDt2,-2);SDt2~dunif(0,2)
 Omega[1:2,1:2]<-inverse(Sig[,]); for(j in 1:2){Sig[j,j]<-pow(sd[j],2)}
 for(i in 1:2){for(j in (i+1):2){Sig[i,j]<-rho[i,j]*sd[i]*sd[j]; Sig[j,i]<-Sig[i,j]}}
 for(j in 1:2){sd[j]~dunif(0,1)}
 for(i in 1:2){for(j in (i+1):2){g[j,i]<-0; ang[i,j]~dunif(0,3.1415); rho[i,j]<-inprod(g[,i],g[,j])}}
 g[1,1]<-1; g[1,2]<-cos(ang[1,2]); g[2,2]<-sin(ang[1,2]) }"
multivar <- function(sheet, suffix, uni){
  if(is.null(uni)||nrow(uni)==0) return(invisible())
  sig <- unique(uni$cov[uni$Se_sig | uni$Sp_sig]); if(length(sig)==0){ cat("  (no sig covariates; multivariate skipped)\n"); return(invisible()) }
  osa<-as.data.frame(read_excel(XL,sheet=sheet)); osa<-osa[!is.na(osa$tp1)&!is.na(osa$fp1)&!is.na(osa$fn1)&!is.na(osa$tn1),]
  osa$test_num<-as.numeric(factor(osa$test,levels=c("ai_only","human_only","human_with_ai"))); osa$study_num<-as.numeric(as.factor(osa$dataset_id))
  for(cv in sig){ r<-as.character(osa[[cv]]); r[is.na(r)]<-"NR"; f<-as.factor(r); ft<-sort(table(f),decreasing=TRUE); osa[[cv]]<-relevel(f,ref=names(ft)[1]) }
  X<-model.matrix(as.formula(paste("~",paste(sig,collapse="+"))), data=osa)[,-1,drop=FALSE]; P<-ncol(X)
  dl<-list(ns=length(unique(osa$study_num)),ntest=3,nObs=nrow(osa),s=osa$study_num,test=osa$test_num,
    tp=osa$tp1,tn=osa$tn1,pos=osa$tp1+osa$fn1,neg=osa$fp1+osa$tn1,X=X,P=P)
  mkinit<-function(c) list(SDs1=1,SDs2=0.8,SDt1=1,SDt2=0.5,sd=c(0.5,0.5),ang=structure(c(NA,1),dim=c(1,2)),
    .RNG.name=c("base::Wichmann-Hill","base::Super-Duper")[c],.RNG.seed=77+c)
  ni<-if(QUICK)500 else 4000
  jm<-jags.model(textConnection(mv_model),data=dl,inits=list(mkinit(1),mkinit(2)),n.chains=2,n.adapt=500,quiet=TRUE)
  update(jm,ni); M<-as.matrix(coda.samples(jm,c("bSe","bSp"),n.iter=ni))
  out<-do.call(rbind,lapply(1:P,function(p){ bs<-M[,sprintf("bSe[%d]",p)]; bp<-M[,sprintf("bSp[%d]",p)]
    data.frame(contrast=colnames(X)[p], Se_b=mean(bs),Se_lo=quantile(bs,.025),Se_hi=quantile(bs,.975),Se_sig=(quantile(bs,.025)>0|quantile(bs,.975)<0),
      Sp_b=mean(bp),Sp_lo=quantile(bp,.025),Sp_hi=quantile(bp,.975),Sp_sig=(quantile(bp,.025)>0|quantile(bp,.975)<0),row.names=NULL)}))
  write.csv(out, file.path(OUT,paste0("metareg_",suffix,"_multivar.csv")), row.names=FALSE); cat(sprintf("  multivariate %s: %d terms\n",suffix,P)); out
}

## ---- (7) figures ----------------------------------------------------------
TIF <- function(f,w,h) tiff(file.path(OUT,f), width=w, height=h, units="in", res=300, compression="lzw", type="cairo")
fig_forest <- function(){              # Figure 2
  dP<-readRDS(file.path(OUT,"patient_draws.rds")); dE<-readRDS(file.path(OUT,"eye_draws.rds"))
  TIF("Figure2_forest.tiff",10,6); par(mfrow=c(2,2), mar=c(4,8,3,1), family="Times", las=1)
  panel<-function(d,mat,xlab,main){ m<-colMeans(d[[mat]]); lo<-apply(d[[mat]],2,quantile,.025); hi<-apply(d[[mat]],2,quantile,.975)
    plot(m,3:1,xlim=range(c(lo,hi)),ylim=c(.5,3.5),pch=19,yaxt="n",xlab=xlab,ylab="",main=main,cex=1.3)
    segments(lo,3:1,hi,3:1,lwd=2); axis(2,at=3:1,labels=ARMS,tick=FALSE)
    text(m,(3:1)+.25,sprintf("%.3f [%.3f, %.3f]",m,lo,hi),cex=.8) }
  panel(dP,"sens","Sensitivity","(A) Patient level — Sensitivity"); panel(dP,"spec","Specificity","(A) Patient level — Specificity")
  panel(dE,"sens","Sensitivity","(B) Eye level — Sensitivity");     panel(dE,"spec","Specificity","(B) Eye level — Specificity")
  dev.off()
}
ellipse <- function(x,y,col){ mu<-c(mean(x),mean(y)); S<-cov(cbind(x,y)); e<-eigen(S); th<-seq(0,2*pi,len=100)
  r<-sqrt(qchisq(.95,2)); pts<-t(mu+ r*(e$vectors %*% diag(sqrt(e$values)) %*% rbind(cos(th),sin(th)))); lines(pts,col=col,lwd=1.5) }
fig_hsroc <- function(){               # Figure 3
  cols<-c("#1b9e77","#d95f02","#7570b3")
  dr<-function(d,main){ plot(NA,xlim=c(0,.4),ylim=c(.6,1),xlab="1 - Specificity",ylab="Sensitivity",main=main)
    for(k in 1:3){ x<-1-d$spec[,k]; y<-d$sens[,k]; points(mean(x),mean(y),pch=19,col=cols[k],cex=1.6); ellipse(x,y,cols[k]) }
    legend("bottomright",ARMS,col=cols,pch=19,bty="n",cex=.9) }
  dP<-readRDS(file.path(OUT,"patient_draws.rds")); dE<-readRDS(file.path(OUT,"eye_draws.rds"))
  TIF("Figure3_HSROC.tiff",11,5.5); par(mfrow=c(1,2),family="Times"); dr(dP,"(A) Patient level"); dr(dE,"(B) Eye level (primary)"); dev.off()
}
fig_multivar_mr <- function(){         # Figure 4: multivariate meta-regression forest (sensitivity scale)
  fp<-file.path(OUT,"metareg_patient_multivar.csv"); fe<-file.path(OUT,"metareg_eye_multivar.csv")
  if(!file.exists(fe)){ cat("  (multivariate CSV absent — Figure 4 skipped)\n"); return(invisible()) }
  HEAD<-1.6
  draw<-function(d,main,cx){ col<-ifelse(d$Se_sig,"red","#666666"); r<-range(c(d$Se_lo,d$Se_hi,0)); xl<-r+c(-1,1)*diff(r)*.08; y<-nrow(d):1
    plot(d$Se_b,y,pch=19,col=col,xlim=xl,ylim=c(.5,nrow(d)+HEAD),yaxt="n",xlab="Beta coefficient (log-odds, sensitivity)",ylab="",main=main,cex.main=1.05)
    segments(d$Se_lo,y,d$Se_hi,y,col=col,lwd=1.8); abline(v=0,lty=2,col="grey50")
    axis(2,at=y,labels=d$contrast,las=2,cex.axis=cx,tick=FALSE)
    text(xl[1],nrow(d)+HEAD*.78,"Red = significant",adj=c(0,.5),cex=.75,col="red") }
  P<-if(file.exists(fp)) read.csv(fp) else NULL; E<-read.csv(fe)
  np<-if(!is.null(P)) nrow(P) else 0
  TIF("Figure4_multivariate_MR.tiff",10,max(6,(np+nrow(E))*0.5+2)); par(mfrow=c(if(np>0)2 else 1,1), mar=c(5,12,3,2), family="Times", las=1)
  if(np>0) draw(P,"(A) Patient level — multivariate meta-regression",0.75)
  draw(E,"(B) Eye level — multivariate meta-regression",0.70); dev.off()
}

## ============================ RUN ==========================================
cat("== Primary NMA (patient, eye) ==\n")
fitP <- fit_nma(prep("patient_rr")$data); summarise(fitP, "patient", prep("patient_rr")$data)
fitE <- fit_nma(prep("eye_rr")$data);     summarise(fitE, "eye", prep("eye_rr")$data)
prevalence("patient"); prevalence("eye")

cat("== Sensitivity (eye): unit-of-analysis + leave-one-out ==\n")
specs <- list(
  list(name="Base (dataset_id, all)",     args=list("eye_rr","dataset_id",FALSE,NULL)),
  list(name="By publication (all)",        args=list("eye_rr","authoryr", FALSE,NULL)),
  list(name="One per dataset x arm",       args=list("eye_rr","dataset_id",TRUE, NULL)),
  list(name="One per publication x arm",   args=list("eye_rr","authoryr", TRUE, NULL)),
  list(name="Leave-one-out (drop Li 2024b)", args=list("eye_rr","dataset_id",FALSE,"Li et al. 2024b")))
sens_rows <- do.call(rbind, lapply(specs, function(s){
  p <- do.call(prep, s$args); f <- fit_nma(p$data); M <- f$M
  rk <- ranks_of(M); dY <- M[,"youden[3]"]-M[,"youden[1]"]
  data.frame(specification=s$name, collab_tables=sum(p$osa$test=="human_with_ai"),
    collab_SUCRA=round(rk$SUCRA[3],2), collab_Pbest=round(rk$P_best[3],2),
    collab_vs_AI_dYouden=round(mean(dY),3), lo=round(quantile(dY,.025),3), hi=round(quantile(dY,.975),3),
    rank=paste(ARMS[order(-rk$mean_rank)],collapse=" > "), row.names=NULL) }))
write.csv(sens_rows, file.path(OUT,"eye_sensitivity.csv"), row.names=FALSE); print(sens_rows, row.names=FALSE)

cat("== Univariate meta-regression ==\n")
uniP <- metareg("patient_rr","patient"); uniE <- metareg("eye_rr","eye")

cat("== Multivariate meta-regression ==\n")
multivar("patient_rr","patient",uniP); multivar("eye_rr","eye",uniE)

cat("== Figures 2-4 ==\n")
fig_forest(); fig_hsroc(); fig_multivar_mr()
cat("\nDONE. Outputs in results/.\n")
