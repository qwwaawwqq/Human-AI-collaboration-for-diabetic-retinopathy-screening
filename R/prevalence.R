## Per-10,000-screened outcomes across a prevalence grid, from saved posterior draws.
## Illustrative scenario estimates (addresses R1 minor 3 / R4 h).
OUT <- "results"
args <- commandArgs(trailingOnly=TRUE)
tag <- args[1]               # e.g. eye_base
prevs <- c(0.02,0.05,0.10,0.20)
d <- readRDS(file.path(OUT,paste0(tag,"_draws.rds")))
arms <- d$arms
rows <- list()
for (p in prevs){
  dis <- 10000*p; hea <- 10000*(1-p)
  for (k in 1:3){
    Se <- d$sens[,k]; Sp <- d$spec[,k]
    missed <- dis*(1-Se)          # false negatives
    falsref <- hea*(1-Sp)         # false positives = unnecessary referrals
    rows[[length(rows)+1]] <- data.frame(
      tag=tag, prevalence=sprintf("%d%%",round(p*100)), arm=arms[k],
      missed_mean=round(mean(missed)), missed_lo=round(quantile(missed,.025)), missed_hi=round(quantile(missed,.975)),
      false_referrals_mean=round(mean(falsref)), false_lo=round(quantile(falsref,.025)), false_hi=round(quantile(falsref,.975)))
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT,paste0(tag,"_prevalence.csv")), row.names=FALSE)
print(res, row.names=FALSE)
