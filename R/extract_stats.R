## Extract canonical Se/Sp/Youden/SUCRA/per-10,000 from the fitted models. Run from repo root.

arms <- c("DL only","Human only","Human + DL")   # test order 1,2,3
getcol <- function(s, name, k){
  x <- s[[name]]; d <- dim(x)
  if(length(d)==2) return(x[,k])
  if(length(d)==3) return(x[,k,1])
  stop("unexpected dims for ", name)
}
q <- function(v) sprintf("%.3f [%.3f, %.3f]", mean(v), quantile(v,.025), quantile(v,.975))
analyze <- function(f, label){
  if(!file.exists(f)){ cat("\n#####", label, ": MODEL NOT READY\n"); return(invisible()) }
  m <- readRDS(f); s <- m$BUGSoutput$sims.list; K <- 3
  cat("\n##########", label, "##########\n-- per-arm Se / Sp / Youden --\n")
  for(k in 1:3) cat(sprintf("%-11s Se %s  Sp %s  Y %s\n", arms[k],
      q(getcol(s,"sens",k)), q(getcol(s,"spec",k)), q(getcol(s,"youden",k))))
  cat("-- SUCRA & P(best) (Youden) --\n")
  for(k in 1:3){ pr <- sapply(1:3, function(mm) mean(s$youdensucra[,k,mm]))
    cat(sprintf("%-11s SUCRA %.3f  P(best) %.3f\n", arms[k], sum(pr*((K-(1:K))/(K-1))), pr[1])) }
  cat("-- pairwise (dYouden / dSe / dSp) --\n")
  for(p in list(c(3,1),c(3,2),c(1,2)))
    cat(sprintf("%s vs %s : dY %s | dSe %s | dSp %s\n", arms[p[1]], arms[p[2]],
      q(getcol(s,"youden",p[1])-getcol(s,"youden",p[2])),
      q(getcol(s,"sens",p[1])-getcol(s,"sens",p[2])),
      q(getcol(s,"spec",p[1])-getcol(s,"spec",p[2]))))
  cat("-- per 10,000 screened (missed / false referrals) --\n")
  for(pv in c(.02,.05,.10,.20)) for(k in 1:3){
    se<-mean(getcol(s,"sens",k)); sp<-mean(getcol(s,"spec",k))
    cat(sprintf("  prev %2.0f%% %-11s missed %3.0f  false %4.0f\n", pv*100, arms[k],
      pv*10000*(1-se), (1-pv)*10000*(1-sp))) }
}
analyze("results/AI_NMA_model.rds","PATIENT-LEVEL (40 studies, +Lee)")
analyze("results/eye_rr_NMA_model.rds","EYE-LEVEL (PRIMARY)")
