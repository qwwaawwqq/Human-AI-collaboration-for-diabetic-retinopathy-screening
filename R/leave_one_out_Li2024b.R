## Leave-one-out: drop Li et al. 2024b (DeepDR-LLM) entirely; eye-level. Run from repo root.
suppressMessages({library(readxl); library(rjags); library(coda)})
set.seed(20260612)
# run from repository root
src <- readLines("R/nma_analysis.R")
i0 <- grep('^model_string <- "', src); rest <- src[(i0+1):length(src)]; i1 <- i0 + which(rest == '"')[1]
eval(parse(text=paste(src[i0:i1], collapse="\n")))
osa0 <- as.data.frame(read_excel("data/Data_extraction_2_v7_LeeAdded.xlsx", sheet="eye_rr"))
osa0 <- osa0[!is.na(osa0$tp1)&!is.na(osa0$fp1)&!is.na(osa0$fn1)&!is.na(osa0$tn1),]
fit <- function(osa,grp,label){
  osa$test_num<-as.numeric(factor(osa$test,levels=c("ai_only","human_only","human_with_ai")))
  osa$study_num<-as.numeric(as.factor(osa[[grp]])); n<-nrow(osa)
  dl<-list(ns=length(unique(osa$study_num)),ntest=3,totaltest=3,nObs=n,id=1:n,s=osa$study_num,
           test=osa$test_num,threshold=rep(1,n),tp=osa$tp1,tn=osa$tn1,pos=osa$tp1+osa$fn1,neg=osa$fp1+osa$tn1)
  jm<-jags.model(textConnection(model_string),data=dl,n.chains=4,n.adapt=1000,quiet=TRUE); update(jm,10000)
  M<-as.matrix(coda.samples(jm,c("sens","spec","youden"),n.iter=40000,thin=40))
  Y<-cbind(M[,"youden[1]"],M[,"youden[2]"],M[,"youden[3]"]); arms<-c("AI","Human","Collab")
  ranks<-t(apply(Y,1,function(r) rank(-r,ties.method="min"))); sucra<-(3-colMeans(ranks))/2
  pb<-tabulate(max.col(Y),3)/nrow(Y); dCA<-Y[,3]-Y[,1]
  cat(sprintf("\n### %s | tables=%d studies=%d | collab tables=%d\n",label,n,dl$ns,sum(osa$test=="human_with_ai")))
  for(k in 1:3) cat(sprintf("  %-7s Se=%.3f Sp=%.3f Youden=%.3f [%.3f,%.3f] SUCRA=%.2f P(best)=%.2f\n",
     arms[k],mean(M[,sprintf("sens[%d,1]",k)]),mean(M[,sprintf("spec[%d,1]",k)]),
     mean(Y[,k]),quantile(Y[,k],.025),quantile(Y[,k],.975),sucra[k],pb[k]))
  cat(sprintf("  Collab vs AI dYouden = %+.3f [%.3f, %.3f] | RANK %s\n",
     mean(dCA),quantile(dCA,.025),quantile(dCA,.975),paste(arms[order(-colMeans(Y))],collapse=">")))
}
# leave-one-out: drop Li et al. 2024b entirely
loo <- osa0[osa0$authoryr != "Li et al. 2024b", ]
cat(sprintf("Dropped Li 2024b: %d -> %d eye tables\n", nrow(osa0), nrow(loo)))
fit(loo,"dataset_id","LEAVE-ONE-OUT: drop Li 2024b (DeepDR-LLM), all remaining tables")
cat("\nDONE\n")
