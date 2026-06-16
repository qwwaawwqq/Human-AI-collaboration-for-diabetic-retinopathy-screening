## Figure 4: combined patient + eye multivariate meta-regression forest plot (run from repo root).
# run from repository root
pat <- read.csv("results/AI_NMA_Multivariate_MR.csv", check.names=FALSE)
eye <- read.csv("results/eye_rr_NMA_Multivariate_MR.csv", check.names=FALSE)
HEAD <- 1.6   # fixed, small headroom (was n+2.2 -> big gap)
draw <- function(means,los,his,cols,labels,n,y,main,lab_cex){
  r <- range(c(los,his,0),na.rm=TRUE); xlim <- r + c(-1,1)*diff(r)*0.08
  plot(means,y,pch=19,col=cols,xlim=xlim,ylim=c(0.5,n+HEAD),yaxt="n",
       xlab="Beta coefficient",ylab="",main=main,cex=1.1,cex.axis=1,cex.lab=1.1,cex.main=1.2)
  segments(los,y,his,y,col=cols,lwd=1.8); abline(v=0,lty=2,col="grey50",lwd=1.2)
  axis(2,at=y,labels=labels,las=2,cex.axis=0.85,tick=FALSE)
  par(xpd=NA)
  for(i in seq_len(n)) text(means[i],y[i]-0.28,sprintf("%.2f [%.2f, %.2f]",means[i],los[i],his[i]),
                            adj=c(0.5,0.5),cex=lab_cex,col=cols[i])
  par(xpd=FALSE)
  # legend INSIDE the headroom, left-aligned, clearly below the title (no overlap)
  text(xlim[1], n+HEAD*0.80, "Red = Significant",            adj=c(0,0.5), cex=0.8, col="red")
  text(xlim[1], n+HEAD*0.45, "Line = 95% Credible Interval", adj=c(0,0.5), cex=0.8, col="black")
}
mk <- function(d){ list(lab=paste0(d$Covariate,": ",d$Level), n=nrow(d), y=nrow(d):1,
  cse=ifelse(d$Se_Sig,"red","#666666"), csp=ifelse(d$Sp_Sig,"red","#666666")) }
P<-mk(pat); E<-mk(eye)
pad <- 1.8                       # was 4 -> inflated the 3-row panel
h_top <- P$n + pad; h_bot <- E$n + pad
render <- function(file,dev){
  dev(file, width=14, height=max(7,(h_top+h_bot)*0.55+2.5), units="in", res=300,
      compression=if(identical(dev,tiff))"lzw" else NULL, type="cairo")
  layout(matrix(1:4,nrow=2,byrow=TRUE), heights=c(h_top,h_bot))
  par(mar=c(5,11,3.0,4.5), family="Times", las=1)
  draw(pat$Se_Mean,pat$Se_Lo,pat$Se_Hi,P$cse,P$lab,P$n,P$y,"Multivariate MR – Sensitivity",0.75)
  mtext("(A)",side=3,line=1.2,adj=0,font=2,cex=1.2,family="Times")
  draw(pat$Sp_Mean,pat$Sp_Lo,pat$Sp_Hi,P$csp,P$lab,P$n,P$y,"Multivariate MR – Specificity",0.75)
  draw(eye$Se_Mean,eye$Se_Lo,eye$Se_Hi,E$cse,E$lab,E$n,E$y,"Multivariate MR – Sensitivity",0.70)
  mtext("(B)",side=3,line=1.2,adj=0,font=2,cex=1.2,family="Times")
  draw(eye$Sp_Mean,eye$Sp_Lo,eye$Sp_Hi,E$csp,E$lab,E$n,E$y,"Multivariate MR – Specificity",0.70)
  layout(1); dev.off()
}
render("results/Figure4_MultivariateMR.tiff", tiff)
png("results/Figure4_MultivariateMR.png", width=1400, height=max(700,(h_top+h_bot)*55+250), res=130)
layout(matrix(1:4,nrow=2,byrow=TRUE), heights=c(h_top,h_bot)); par(mar=c(5,11,3.0,4.5),family="Times",las=1)
draw(pat$Se_Mean,pat$Se_Lo,pat$Se_Hi,P$cse,P$lab,P$n,P$y,"Multivariate MR – Sensitivity",0.75); mtext("(A)",side=3,line=1.2,adj=0,font=2,cex=1.2)
draw(pat$Sp_Mean,pat$Sp_Lo,pat$Sp_Hi,P$csp,P$lab,P$n,P$y,"Multivariate MR – Specificity",0.75)
draw(eye$Se_Mean,eye$Se_Lo,eye$Se_Hi,E$cse,E$lab,E$n,E$y,"Multivariate MR – Sensitivity",0.70); mtext("(B)",side=3,line=1.2,adj=0,font=2,cex=1.2)
draw(eye$Sp_Mean,eye$Sp_Lo,eye$Sp_Hi,E$csp,E$lab,E$n,E$y,"Multivariate MR – Specificity",0.70); dev.off()
cat("rendered. h_top=",h_top," h_bot=",h_bot,"\n")
