## Faithful rjags port of the bivariate DTA-NMA in DR_Pathway_NMA.R
## + pairwise posterior contrasts, P(superiority), SUCRA/ranking, unit-of-analysis sensitivity.
## Usage: R -f nma_analysis.R --args <sheet> <tag> <group_col> <one_per:0|1>
args <- commandArgs(trailingOnly=TRUE)
sheet     <- args[1]
tag       <- args[2]
group_col <- ifelse(length(args)>=3 && nchar(args[3])>0, args[3], "dataset_id")
one_per   <- (length(args)>=4 && args[4]=="1")

N_CHAINS<-4; N_ADAPT<-1000; N_BURN<-4000; N_KEEP<-4000; N_THIN<-1
XL <- "data/Data_extraction_2_v7_LeeAdded.xlsx"
OUT <- "results"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# jags.moddir auto-detected by rjags
suppressMessages({library(readxl); library(rjags); library(coda)})
set.seed(20260612)

## ---- exact model string from DR_Pathway_NMA.R ----
model_string <- "
model {
 for(i in 1:nObs){
  tp[i]~dbin(pi[i,1], pos[i])
  tn[i]~dbin(pi[i,2], neg[i])
  logit(pi[i,1]) <- mu[i,1]
  logit(pi[i,2]) <- mu[i,2]
  MU[i,1] <- threshold.sens[test[i],threshold[i]] + study.re.sens[s[i]] + test.re.sens[s[i],test[i]]
  MU[i,2] <- threshold.spec[test[i],threshold[i]] + study.re.spec[s[i]] + test.re.spec[s[i],test[i]]
  mu[i,1:2] ~ dmnorm(MU[i,], Omega[,])
 }
 for(k in 1:ntest){
  sens[k,1]<- exp(threshold.sens[k,1])/(1+exp(threshold.sens[k,1]))
  spec[k,1]<- exp(threshold.spec[k,1])/(1+exp(threshold.spec[k,1]))
  threshold.sens[k,1]~dnorm(0,0.05)
  threshold.spec[k,1]~dnorm(0,0.05)
 }
 for(k in 1:ns){
  study.re.sens[k]~dnorm(0,taustudysens)
  study.re.spec[k]~dnorm(0,taustudyspec)
  for(l in 1:ntest){
   test.re.sens[k,l]~dnorm(0,tautestsens)
   test.re.spec[k,l]~dnorm(0,tautestspec)
  }
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

## ---- data prep (mirrors the original script) ----
osa <- as.data.frame(read_excel(XL, sheet=sheet))
osa <- osa[!is.na(osa$tp1) & !is.na(osa$fp1) & !is.na(osa$fn1) & !is.na(osa$tn1), ]
if (one_per) {
  osa$.n <- osa$tp1+osa$fn1+osa$fp1+osa$tn1
  key <- paste(osa[[group_col]], osa$test, sep="||")
  osa <- osa[order(key, -osa$.n), ]
  osa <- osa[!duplicated(key), ]
}
osa$test_num  <- as.numeric(factor(osa$test, levels=c("ai_only","human_only","human_with_ai")))
osa$study_num <- as.numeric(as.factor(osa[[group_col]]))
n <- nrow(osa)
data_list <- list(ns=length(unique(osa$study_num)), ntest=3, totaltest=3, nObs=n,
                  id=1:n, s=osa$study_num, test=osa$test_num, threshold=rep(1,n),
                  tp=osa$tp1, tn=osa$tn1, pos=osa$tp1+osa$fn1, neg=osa$fp1+osa$tn1)
cat(sprintf("[%s] sheet=%s tables=%d groups=%d (group=%s, one_per=%s)\n",
            tag, sheet, n, data_list$ns, group_col, one_per))
cat("  tables per arm: ", paste(names(table(osa$test)), table(osa$test), sep="="), "\n")

rngs <- c("base::Wichmann-Hill","base::Marsaglia-Multicarry","base::Super-Duper","base::Mersenne-Twister")
mk_init <- function(sd0, c) list(SDstudysens=sd0,SDstudyspec=sd0,SDtestsens=sd0,SDtestspec=sd0,
                                 a=structure(c(NA,1),dim=c(1,2)), sd=c(0.5,0.5),
                                 .RNG.name=rngs[c], .RNG.seed=1000+c)
inits <- lapply(1:N_CHAINS, function(c) mk_init(c(1.0,0.8,1.2,0.5)[c], c))

params <- c("sens","spec","youden","youdensucra","SDstudysens","SDstudyspec",
            "SDtestsens","SDtestspec","sd","rho")
library(parallel)
run_chain <- function(c){
  # jags.moddir auto-detected by rjags
  suppressMessages(library(rjags))
  jm <- jags.model(textConnection(model_string), data=data_list, inits=inits[[c]],
                   n.chains=1, n.adapt=N_ADAPT, quiet=TRUE)
  update(jm, N_BURN)
  coda.samples(jm, params, n.iter=N_KEEP, thin=N_THIN)[[1]]
}
chain_list <- mclapply(1:N_CHAINS, run_chain, mc.cores=min(N_CHAINS, detectCores()))
ok <- sapply(chain_list, function(x) inherits(x,"mcmc"))
if(!all(ok)){ cat("CHAIN ERROR:\n"); print(chain_list[!ok]); stop("parallel chain failure") }
samp <- as.mcmc.list(chain_list)
M <- as.matrix(samp)
arms <- c("DL only","Human only","Human + DL")
qs <- function(x) c(mean=mean(x), lo=quantile(x,.025,names=FALSE), hi=quantile(x,.975,names=FALSE))

## per-arm
per <- do.call(rbind, lapply(1:3, function(k){
  data.frame(arm=arms[k],
    sens_mean=mean(M[,sprintf("sens[%d,1]",k)]), sens_lo=quantile(M[,sprintf("sens[%d,1]",k)],.025),
    sens_hi=quantile(M[,sprintf("sens[%d,1]",k)],.975),
    spec_mean=mean(M[,sprintf("spec[%d,1]",k)]), spec_lo=quantile(M[,sprintf("spec[%d,1]",k)],.025),
    spec_hi=quantile(M[,sprintf("spec[%d,1]",k)],.975),
    youden_mean=mean(M[,sprintf("youden[%d]",k)]), youden_lo=quantile(M[,sprintf("youden[%d]",k)],.025),
    youden_hi=quantile(M[,sprintf("youden[%d]",k)],.975), row.names=NULL)}))
write.csv(per, file.path(OUT,paste0(tag,"_per_arm.csv")), row.names=FALSE)

## pairwise contrasts: (Human+DL vs DL), (Human+DL vs Human), (Human vs DL)
prs <- list(c(3,1),c(3,2),c(2,1))
pw <- do.call(rbind, lapply(prs, function(p){
  i<-p[1]; j<-p[2]
  dSe<-M[,sprintf("sens[%d,1]",i)]-M[,sprintf("sens[%d,1]",j)]
  dSp<-M[,sprintf("spec[%d,1]",i)]-M[,sprintf("spec[%d,1]",j)]
  dY <-M[,sprintf("youden[%d]",i)]-M[,sprintf("youden[%d]",j)]
  data.frame(comparison=paste(arms[i],"-",arms[j]),
    dSens_mean=mean(dSe),dSens_lo=quantile(dSe,.025),dSens_hi=quantile(dSe,.975),P_sens_sup=mean(dSe>0),
    dSpec_mean=mean(dSp),dSpec_lo=quantile(dSp,.025),dSpec_hi=quantile(dSp,.975),P_spec_sup=mean(dSp>0),
    dYouden_mean=mean(dY),dYouden_lo=quantile(dY,.025),dYouden_hi=quantile(dY,.975),P_youden_sup=mean(dY>0),
    row.names=NULL)}))
write.csv(pw, file.path(OUT,paste0(tag,"_pairwise.csv")), row.names=FALSE)

## ranking
K<-3
Prank <- sapply(1:3, function(k) sapply(1:3, function(mm) mean(M[,sprintf("youdensucra[%d,%d]",k,mm)])))
# Prank[m,k] = P(arm k has rank m); rank 1 = best
rank_tab <- do.call(rbind, lapply(1:3, function(k){
  pr <- Prank[,k]
  data.frame(arm=arms[k], P_best=pr[1],
             SUCRA=sum(sapply(1:K, function(m) (K-m)/(K-1)*pr[m])),
             mean_rank=sum((1:K)*pr), row.names=NULL)}))
write.csv(rank_tab, file.path(OUT,paste0(tag,"_ranking.csv")), row.names=FALSE)

## heterogeneity
het <- data.frame(
  param=c("SDstudy_sens","SDstudy_spec","SDtest_sens","SDtest_spec"),
  mean=c(mean(M[,"SDstudysens"]),mean(M[,"SDstudyspec"]),mean(M[,"SDtestsens"]),mean(M[,"SDtestspec"])),
  lo=c(quantile(M[,"SDstudysens"],.025),quantile(M[,"SDstudyspec"],.025),quantile(M[,"SDtestsens"],.025),quantile(M[,"SDtestspec"],.025)),
  hi=c(quantile(M[,"SDstudysens"],.975),quantile(M[,"SDstudyspec"],.975),quantile(M[,"SDtestsens"],.975),quantile(M[,"SDtestspec"],.975)))
write.csv(het, file.path(OUT,paste0(tag,"_het.csv")), row.names=FALSE)

## convergence: Gelman-Rubin on sens/spec/youden
key_cols <- c(sprintf("sens[%d,1]",1:3),sprintf("spec[%d,1]",1:3),sprintf("youden[%d]",1:3))
gd <- tryCatch(gelman.diag(samp[,key_cols], autoburnin=FALSE, multivariate=FALSE)$psrf, error=function(e) NA)
ess <- tryCatch(effectiveSize(samp[,key_cols]), error=function(e) NA)
conv <- data.frame(param=key_cols, Rhat=if(is.matrix(gd)) round(gd[,1],3) else NA,
                   ESS=if(length(ess)>1) round(ess) else NA)
write.csv(conv, file.path(OUT,paste0(tag,"_conv.csv")), row.names=FALSE)

## save draws for prevalence calc
saveRDS(list(sens=M[,sprintf("sens[%d,1]",1:3)], spec=M[,sprintf("spec[%d,1]",1:3)],
             youden=M[,sprintf("youden[%d]",1:3)], arms=arms, meta=list(tag=tag,n=n,ns=data_list$ns)),
        file.path(OUT,paste0(tag,"_draws.rds")))

cat("\n==== ",tag," ====\n"); print(per, digits=3)
cat("\n-- pairwise --\n"); print(pw, digits=3)
cat("\n-- ranking --\n"); print(rank_tab, digits=3)
cat("\n-- convergence (max Rhat=",round(max(conv$Rhat,na.rm=TRUE),3),") --\n")
cat("DONE_",tag,"\n",sep="")
writeLines("ok", file.path(OUT,paste0(tag,"_DONE.flag")))
