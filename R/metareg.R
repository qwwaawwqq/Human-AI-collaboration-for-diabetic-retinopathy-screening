a<-commandArgs(trailingOnly=TRUE); sheet<-a[1]; covarg<-a[2]; suf<-a[3]
XL<-"data/Data_extraction_2_v7_LeeAdded.xlsx"
OUT<-"results"
# jags.moddir auto-detected by rjags
suppressMessages({library(readxl);library(rjags);library(coda)}); set.seed(2026)
osa<-as.data.frame(read_excel(XL,sheet=sheet))
osa<-osa[!is.na(osa$tp1)&!is.na(osa$fp1)&!is.na(osa$fn1)&!is.na(osa$tn1),]
osa$test_num<-as.numeric(factor(osa$test,levels=c("ai_only","human_only","human_with_ai")))
osa$study_num<-as.numeric(as.factor(osa$dataset_id)); n<-nrow(osa)
base<-list(ns=length(unique(osa$study_num)),ntest=3,nObs=n,s=osa$study_num,test=osa$test_num,
  tp=osa$tp1,tn=osa$tn1,pos=osa$tp1+osa$fn1,neg=osa$fp1+osa$tn1)
m<-"model{
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
mkinit<-function(c) list(SDs1=1,SDs2=0.8,SDt1=1,SDt2=0.5,sd=c(0.5,0.5),
  ang=structure(c(NA,1),dim=c(1,2)),.RNG.name=c("base::Wichmann-Hill","base::Super-Duper")[c],.RNG.seed=99+c)
covs<-c("readerexpertise","study_design","multicenter","economic","healthcare","external","dr",
 "othereyetarget","pupil","ctype","camera","criteria","analysis_type","rstandard","disagree",
 "rimage","field","vendor","commercial","certified","quality_check","architecture","region",
 "aigradable","hgradable")
if(covarg!="ALL") covs<-strsplit(covarg,",")[[1]]
res<-data.frame()
for(cv in covs){ r<-osa[[cv]]; if(is.null(r)||all(is.na(r))){cat("skip(missing)",cv,"\n");next}
 r[is.na(r)]<-"NR"; f<-as.factor(as.character(r)); ft<-sort(table(f),decreasing=TRUE)
 f<-relevel(f,ref=names(ft)[1]); lev<-levels(f); nl<-length(lev); if(nl<2){cat("skip(const)",cv,"\n");next}
 dl<-c(base,list(cov_x=as.numeric(f),n_lev=nl))
 ok<-tryCatch({ jm<-jags.model(textConnection(m),data=dl,inits=list(mkinit(1),mkinit(2)),n.chains=2,n.adapt=500,quiet=TRUE)
  update(jm,1500); sm<-coda.samples(jm,c("beta_Se","beta_Sp"),n.iter=1500); TRUE},error=function(e){cat("ERR",cv,"::",e$message,"\n");FALSE})
 if(!ok)next; M<-as.matrix(sm)
 for(k in 2:nl){ bs<-M[,sprintf("beta_Se[%d]",k)]; bp<-M[,sprintf("beta_Sp[%d]",k)]
  res<-rbind(res,data.frame(level=ifelse(sheet=="patient_rr","patient","eye"),cov=cv,
   contrast=paste0(lev[k]," vs ",lev[1]),
   Se_OR=round(exp(mean(bs)),3),Se_lo=round(exp(quantile(bs,.025)),3),Se_hi=round(exp(quantile(bs,.975)),3),
   Se_sig=(quantile(bs,.025)>0|quantile(bs,.975)<0),
   Sp_OR=round(exp(mean(bp)),3),Sp_lo=round(exp(quantile(bp,.025)),3),Sp_hi=round(exp(quantile(bp,.975)),3),
   Sp_sig=(quantile(bp,.025)>0|quantile(bp,.975)<0))) }
 cat("done",cv,"(",nl,"lev)\n") }
fn<-file.path(OUT,paste0("metareg_",suf,".csv")); write.csv(res,fn,row.names=FALSE)
cat("SAVED",fn,"rows=",nrow(res),"\n")
