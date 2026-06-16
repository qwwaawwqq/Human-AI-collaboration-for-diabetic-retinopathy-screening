## Bivariate DTA-NMA full pipeline (50,000 iter / 10,000 burn-in, 4 chains; JAGS via R2jags/rjags).
## Generates Figures 2-4 (forest, HSROC, multivariate MR), supplementary convergence/funnel/MR plots,
## the fitted models (results/*.rds) and the result CSVs. Run from the repository root.
## Pre-computed models + CSVs are provided in results/ so figures reproduce without re-fitting.
dir.create("results", showWarnings=FALSE, recursive=TRUE)

## Data
# pos & neg are included as data 
# because jags does not allow creating new data in Model
## ranking treatments
# rank() in jags is different from that in winbugs
# wrote a loop to calculate the ranking
## Constraints
# add constraints on estimations of sensitivity & specificity
# Include bsens & bspec in Data
# csens & cspec as stochastic nodes 
# jags only allow stochastic nodes on the left of <- 

# Call the library "rjags"
# Load the necessary libraries
library(R2jags)
set.seed(20260615)  # reproducibility for 50k re-run

# Define the JAGS model as a string
model_string <- "
model {

temp.var<- id[1]

# Bivariate model
 	for(i in 1:nObs){																	#Loop through observations

# pos: cases; neg: noncases		
#	pos[i] <- tp[i] + fn[i]
#	neg[i] <- fp[i] + tn[i]
	
# pi[,1]: true sensitivity; pi[,2]: true specificity	
  tp[i]~dbin(pi[i,1], pos[i])
	tn[i]~dbin(pi[i,2], neg[i])

# mu[,1]: logit true sensitivity; mu[,2]: logit true specificity				
	logit(pi[i,1]) <- mu[i,1] 
	logit(pi[i,2]) <- mu[i,2]
			
# Model for linear predictor
# threshold.sens[test[i],threshold[i]]: fixed test & threshold effect for sensitivity
# threshold.spec[test[i],threshold[i]]: fixed test & threshold effect for specifcity
# study.re.sens[s[i]]: study random effect for sensitivity
# study.re.spec[s[i]]: study random effect for specificity 
# test.re.sens[s[i],test[i]]: within study sensitivity heterogeneity for tests
# test.re.spec[s[i],test[i]]: within study specificity heterogeneity for tests	
MU[i,1] <-  threshold.sens[test[i],threshold[i]] + study.re.sens[s[i]] + test.re.sens[s[i],test[i]]
MU[i,2] <-  threshold.spec[test[i],threshold[i]] + study.re.spec[s[i]] + test.re.spec[s[i],test[i]]
	
# For models assuming a common correlation and heterogeneity parameter across tests
# Correlation between sensitivity and specificity across tests
mu[i,1:2] ~ dmnorm(MU[i,], Omega[,])

}
			
# Back transform on to sensitivity and specificity scale					
# We have 5 tests and each has different thresholds									
# Loop through each test
# test 1 to ntest
	for(k in 1:ntest){									
# Loop through the number of thresholds							
	sens[k,1]<- exp(threshold.sens[k,1])/(1+exp(threshold.sens[k,1]))
	spec[k,1]<- exp(threshold.spec[k,1])/(1+exp(threshold.spec[k,1]))	
			
# Priors on the fixed test and threshold effects		
	threshold.sens[k,1]~dnorm(0,0.05)
	threshold.spec[k,1]~dnorm(0,0.05)				
}
			
# Priors on the random study and test effects 			
	for(k in 1:ns){										
# Loop through the number of studies
	study.re.sens[k]~dnorm(0,taustudysens)
	study.re.spec[k]~dnorm(0,taustudyspec)

for(l in 1:ntest){										
# Loop through the number of tests
	test.re.sens[k,l]~dnorm(0,tautestsens)
	test.re.spec[k,l]~dnorm(0,tautestspec)
		}
}

# Priors on the study-level random effects for sensitivity				
	taustudysens <- pow(SDstudysens,-2)
	SDstudysens ~ dunif(0,2)
# Priors on the study-level random effects for specificity
# Different tau for sens & spec			
	taustudyspec <- pow(SDstudyspec,-2)
	SDstudyspec ~ dunif(0,2)
# Priors on the within-study random effects for sensitivity across tests		
	tautestsens <- pow(SDtestsens,-2)
	SDtestsens ~ dunif(0,2)
# Priors on the within-study random effects for specificity across tests				
	tautestspec <- pow(SDtestspec,-2)
	SDtestspec ~ dunif(0,2)

# For models assuming common correlation and heterogeneity parameter
# Specifying covariance matrix

	 Omega[1:2,1:2] <- inverse(Sigma.sq[,])         
	 for(m in 1:2) { 
      	 Sigma.sq[m,m] <- pow(sd[m],2)  # Diagonal elements
	 } 		   
	
  	 for(i in 1:2) {
       	 for(j in (i+1):2) {
       	 Sigma.sq[i,j] <- rho[i,j]*sd[i]*sd[j]  # off-diagonal elements
          	 Sigma.sq[j,i] <- Sigma.sq[i,j]    # off-diagonal elements
      	 }
 	 }
	
# Prior on the common between-study standard deviation
	 for(m in 1:2) {    
      	 sd[m] ~ dunif(0, 1)                 
  	 } 

# Spherical parameterisation to estimate correlation parameter
  	for(i in 1:2) {
       	for(j in (i+1):2) {
         	 g[j,i] <- 0
	       a[i,j] ~ dunif(0, 3.1415)
# inprod: inner product of two vectors	
	       rho[i,j] <- inprod(g[,i], g[,j])  
      		}
  } 
  
  	g[1,1] <- 1
  	g[1,2] <- cos(a[1,2])
 	  g[2,2] <- sin(a[1,2])

# Calculate logit sensitivity and specificity for each test/threshold combination			
	for (k in 1:ntest) {
			dsens[k] <- threshold.sens[k,1]
			dspec[k] <- threshold.spec[k,1]
	}
# Alternative code		
#for (i in 1:5) {
#    for (j in 1:3) {
#    dsens[(i-1)*3 + j] <- threshold.sens[i, j]
#    dspec[(i-1)*3 + j] <- threshold.spec[i, j]
#    }
#  }		

# Youden's J index
  # Step 1: Calculate the youden values
  for (k in 1:totaltest) {
    youden[k] <- exp(dsens[k]) / (1 + exp(dsens[k])) + exp(dspec[k]) / (1 + exp(dspec[k])) - 1
  }

  # Step 2: Initialize rank comparison matrix
  for (k in 1:totaltest) {
    for (j in 1:totaltest) {
      comparison[k, j] <- step(youden[j] - youden[k])
    }
  }

  # Step 3: Calculate ranks
  for (k in 1:totaltest) {
    rkyouden[k] <- sum(comparison[k, ])
  }

  # Step 4: Compute youdensucra
  for (k in 1:totaltest) {
    for (m in 1:totaltest) {
      youdensucra[k, m] <- equals(rkyouden[k], m)
    }
  }

}
"
## Import data from Excel file
## AI vs Human vs AI+Human -- Network Meta-Analysis
library("readxl")
library(dplyr)
library(stringr)

osa <- as.data.frame(read_excel("Data_extraction_2_v7_LeeAdded.xlsx", sheet = "patient_rr"))

# Keep only complete rows
osa <- osa[!is.na(osa$tp1) & !is.na(osa$fp1) & !is.na(osa$fn1) & !is.na(osa$tn1), ]

# Map 3 test arms: 1 = ai_only, 2 = human_only, 3 = human_with_ai
osa$test_num  <- as.numeric(factor(osa$test, levels = c("ai_only", "human_only", "human_with_ai")))
osa$study_num <- as.numeric(as.factor(osa$dataset_id))   # group by dataset_id

cat("Rows per test arm:\n"); print(table(osa$test))

ntest_arms <- 3   # AI, Human, Human+AI

data_list <- list(
  ns        = length(unique(osa$study_num)),
  ntest     = ntest_arms,
  totaltest = ntest_arms,
  nObs      = nrow(osa),
  id        = 1:nrow(osa),
  s         = osa$study_num,
  test      = osa$test_num,
  threshold = rep(1, nrow(osa)),
  tp        = osa$tp1,
  tn        = osa$tn1,
  pos       = osa$tp1 + osa$fn1,
  neg       = osa$fp1 + osa$tn1
)
cat("Data list summary:\n"); str(data_list)

# Initial values
ns_art <- length(unique(osa$study_num))

inits1 <- list(
  SDstudysens   = 1, SDstudyspec = 1, SDtestsens = 1, SDtestspec = 1,
  a             = matrix(c(NA, 1), nrow = 1, ncol = 2),
  mu            = matrix(rep(0.2, nrow(osa)*2), nrow = nrow(osa), ncol = 2),
  sd            = c(1, 1),
  study.re.sens = rep(1, ns_art),
  study.re.spec = rep(1, ns_art)
)

inits2 <- list(
  SDstudysens   = 0.8, SDstudyspec = 0.8, SDtestsens = 0.8, SDtestspec = 0.8,
  a             = matrix(c(NA, 1), nrow = 1, ncol = 2),
  mu            = matrix(rep(0.1, nrow(osa)*2), nrow = nrow(osa), ncol = 2),
  sd            = c(1, 1),
  study.re.sens = rep(1, ns_art),
  study.re.spec = rep(1, ns_art)
)

inits3 <- list(
  SDstudysens   = 1.2, SDstudyspec = 1.2, SDtestsens = 1.2, SDtestspec = 1.2,
  a             = matrix(c(NA, 1), nrow = 1, ncol = 2),
  mu            = matrix(rep(0.3, nrow(osa)*2), nrow = nrow(osa), ncol = 2),
  sd            = c(1, 1),
  study.re.sens = rep(1, ns_art),
  study.re.spec = rep(1, ns_art)
)

inits4 <- list(
  SDstudysens   = 0.5, SDstudyspec = 0.5, SDtestsens = 0.5, SDtestspec = 0.5,
  a             = matrix(c(NA, 1), nrow = 1, ncol = 2),
  mu            = matrix(rep(0.0, nrow(osa)*2), nrow = nrow(osa), ncol = 2),
  sd            = c(1, 1),
  study.re.sens = rep(1, ns_art),
  study.re.spec = rep(1, ns_art)
)

# Parameters to be monitored
parameters <- c("sens", "spec", "threshold.sens", "threshold.spec", "SDstudysens",
                "SDstudyspec", "SDtestsens", "SDtestspec", "sd", "rho",
                "youden", "youdensucra")

## running jags
# DEBUG: n.iter=1000, n.burnin=500  |  PRODUCTION: n.iter=10000, n.burnin=5000
DEBUG_MODE <- FALSE   # <-- set FALSE for production run
n_iter   <- ifelse(DEBUG_MODE, 1000,  50000)
n_burnin <- ifelse(DEBUG_MODE,  500,   10000)

if (file.exists("AI_NMA_model.rds")) {
  cat("Loading existing AI_NMA_model.rds...\n")
  mod.fit.R2jags <- readRDS("results/AI_NMA_model.rds")
} else {
  mod.fit.R2jags <- jags(data = data_list, inits = list(inits1, inits2, inits3, inits4),
                    parameters.to.save = parameters, n.chains = 4,
                    n.iter = n_iter, n.burnin = n_burnin,
                    model.file = textConnection(model_string))
  saveRDS(mod.fit.R2jags, "results/AI_NMA_model.rds")
}

features <- mod.fit.R2jags$BUGSoutput

library(stringr); library(dplyr); library(ggplot2); library(gridExtra); library(patchwork)

# ── Summarise NMA posterior per arm ──────────────────────────────────────────
selected_params <- features$summary[grep("^sens|^spec", rownames(features$summary)), ]
plot_data <- as.data.frame(selected_params) %>%
  mutate(Parameter = rownames(selected_params),
         indices = str_extract_all(Parameter, "\\d+"),
         i = as.numeric(sapply(indices, `[`, 1)),
         j = as.numeric(sapply(indices, `[`, 2))) %>%
  arrange(i, j)

arm_labels <- c("DL only", "Human only", "Human + DL")
arm_cols   <- c("#2166AC","#D6604D","#1A9641")

sens_plot_data <- plot_data[grep("^sens", rownames(plot_data)), ]
spec_plot_data <- plot_data[grep("^spec", rownames(plot_data)), ]
sens_plot_data$test <- factor(arm_labels, levels = rev(arm_labels))
spec_plot_data$test <- factor(arm_labels, levels = rev(arm_labels))

# ── Compute observed Se / Sp per row from raw cell counts ────────────────────
osa$obs_Se <- osa$tp1 / (osa$tp1 + osa$fn1)
osa$obs_Sp <- osa$tn1 / (osa$fp1 + osa$tn1)
se_lo_w <- function(p, n) mapply(function(pp,nn){
  z <- qnorm(0.975); ((pp+z^2/(2*nn))-z*sqrt((pp*(1-pp)/nn)+z^2/(4*nn^2)))/(1+z^2/nn)}, p, n)
se_hi_w <- function(p, n) mapply(function(pp,nn){
  z <- qnorm(0.975); ((pp+z^2/(2*nn))+z*sqrt((pp*(1-pp)/nn)+z^2/(4*nn^2)))/(1+z^2/nn)}, p, n)
osa$obs_Se_lo <- se_lo_w(osa$obs_Se, osa$tp1+osa$fn1)
osa$obs_Se_hi <- se_hi_w(osa$obs_Se, osa$tp1+osa$fn1)
osa$obs_Sp_lo <- se_lo_w(osa$obs_Sp, osa$fp1+osa$tn1)
osa$obs_Sp_hi <- se_hi_w(osa$obs_Sp, osa$fp1+osa$tn1)
osa$raw_counts <- sprintf(" (%d, %d, %d, %d)", osa$tp1, osa$fp1, osa$fn1, osa$tn1)
osa$arm_col  <- arm_cols[osa$test_num]
osa$arm_label <- factor(arm_labels[osa$test_num], levels = arm_labels)
# Build study label: pure authoryr
osa$study_label <- osa$authoryr

# Sort by test arm then year (newest at top) BEFORE de-duplication
osa$year_num <- as.numeric(stringr::str_extract(osa$authoryr, "\\d{4}"))
osa_sorted <- osa[order(osa$test_num, -osa$year_num, osa$authoryr), ]

  # (User requested: do NOT deduplicate author years, ensure every row has a label)

# Append raw counts to whatever label string remains (even if blank)
# (REMOVED: User wants raw numbers evenly spaced on right, not pasted into the label)
# We still need unique internal IDs so lines don't collapse in sorting/plotting
osa_sorted$plot_id <- make.unique(paste0(osa_sorted$authoryr, "_", osa_sorted$index))
osa_sorted$plot_id <- factor(osa_sorted$plot_id, levels = rev(unique(osa_sorted$plot_id)))

# ── Build full COMBINED forest plot (individual + summary) ──────────────────
make_full_forest <- function(osa_df, sens_sum, spec_sum,
                             arm_labs, arm_colors, title_prefix="") {
  n_studies <- nrow(osa_df)
  n_arms    <- length(arm_labs)
  # Individual study points
  df_ind <- data.frame(
    plot_id= osa_df$plot_id,
    label  = osa_df$study_label,
    TP     = osa_df$tp1, FP = osa_df$fp1, FN = osa_df$fn1, TN = osa_df$tn1,
    Se_m   = osa_df$obs_Se, Se_lo = osa_df$obs_Se_lo, Se_hi = osa_df$obs_Se_hi,
    Sp_m   = osa_df$obs_Sp, Sp_lo = osa_df$obs_Sp_lo, Sp_hi = osa_df$obs_Sp_hi,
    arm    = factor(arm_labs[osa_df$test_num], levels=arm_labs),
    col    = arm_colors[osa_df$test_num],
    type   = "Study",
    stringsAsFactors = FALSE
  )

  # Summary rows (one per arm, with posterior mean + 95% CrI)
  df_sum <- data.frame(
    plot_id= paste0("Summary_", arm_labs),
    label  = paste0("\u25BC Summary: ", arm_labs),
    TP     = "", FP = "", FN = "", TN = "",
    Se_m   = sens_sum[,"mean"], Se_lo = sens_sum[,"2.5%"], Se_hi = sens_sum[,"97.5%"],
    Sp_m   = spec_sum[,"mean"], Sp_lo = spec_sum[,"2.5%"], Sp_hi = spec_sum[,"97.5%"],
    arm    = factor(arm_labs, levels=arm_labs),
    col    = arm_colors,
    type   = "Summary",
    stringsAsFactors = FALSE
  )
  df_all <- rbind(df_ind, df_sum)

  # y-positions: studies first (top), then summary rows (bottom, separated)
  y <- nrow(df_all):1
  y_ind <- y[(n_studies+1):length(y)]          # individual rows from top
  y_sum <- y[1:n_arms]                          # summary rows at bottom

  # Put summary ABOVE, individuals below (standard forest plot layout)
  # Actually standard: summary at bottom. Use y order = summary first (larger y = top)
  y_all <- c(y_sum + n_studies + 2, y_ind)    # summary above gap, individuals below

  # Simpler: just reverse — summary at top as pooled rows
  y_pos <- seq(nrow(df_all)+n_arms+1, n_arms+2)  # gap of n_arms between summary&studies
  y_pos_sum <- (nrow(df_all)+n_arms+2):(nrow(df_all)+2)
  # Use straightforward y positions
  y_plot_ind <- n_arms:1 + n_studies + n_arms
  y_plot_sum <- (n_arms-1):0 + n_arms

  # Actually, simplest readable approach:
  y_final <- c((n_studies + n_arms + 1):(n_arms + 2),   # individual (top)
               (n_arms):1)                                # summary (bottom)

  list(df=df_all, y=y_final,
       n_sum=n_arms, n_ind=n_studies,
       arm_labs=arm_labs, arm_colors=arm_colors)
}

fp_data <- make_full_forest(osa_sorted, sens_plot_data[,c("mean","2.5%","97.5%")],
                             spec_plot_data[,c("mean","2.5%","97.5%")],
                             arm_labels, arm_cols)

plot_forest_panel <- function(fp, metric="Se", xlab="Sensitivity", title="", show_left_text=TRUE) {
  df   <- fp$df
  
  # Map y coordinates explicitly using the pre-calculated gap vectors from make_full_forest
  # This preserves the spacing above the summary rows, while entirely preventing identical studies
  # (like "Li 2021") from overlapping or generating un-plotted gap lines.
  df$y <- fp$y
  
  m    <- df[[paste0(metric,"_m")]]
  lo   <- df[[paste0(metric,"_lo")]]
  hi   <- df[[paste0(metric,"_hi")]]
  is_s <- df$type == "Summary"

  xlim_use <- range(c(lo,hi,0,1), na.rm=TRUE)
  
  # Expand left margin only for the first panel
  if (show_left_text) {
    par(mar=c(5, 35, 2, 11), las=1)
  } else {
    par(mar=c(5, 1, 2, 11), las=1)
  }
  
  plot(m, df$y, pch=ifelse(is_s, 18, 19),
       cex=ifelse(is_s, 1.8, 0.85),
       col=as.character(df$col),
       xlim=xlim_use, ylim=c(0.0, max(df$y, na.rm=TRUE) + 2.0),
       yaxs="i",
       yaxt="n", xlab=xlab, ylab="",
       main=title, cex.axis=0.9, cex.lab=1.0, cex.main=1.1)
  segments(lo, df$y, hi, df$y, col=as.character(df$col),
           lwd=ifelse(is_s, 2.2, 0.9))
  abline(v=0.5, lty=3, col="grey70")
  # Horizontal rule above summaries
  abline(h=fp$n_sum + 1.5, lty=1, col="grey40", lwd=0.8)
  
  # Draw custom explicit table data on the left!
  par(xpd=NA)
  if (show_left_text) {
    # Push the left-most columns out much further to prevent overlap
    col_names <- c("Study"=-1.50, "TP"=-0.60, "FP"=-0.45, "FN"=-0.30, "TN"=-0.15)
    
    # Print headers above the very top (highest Y) row
    top_y <- max(df$y, na.rm=TRUE)
    text(x = col_names["Study"], y = top_y + 1.2, labels = "Study", adj=0, font=2, cex=1.1)
    text(x = col_names["TP"],    y = top_y + 1.2, labels = "TP",    adj=1, font=2, cex=1.1)
    text(x = col_names["FP"],    y = top_y + 1.2, labels = "FP",    adj=1, font=2, cex=1.1)
    text(x = col_names["FN"],    y = top_y + 1.2, labels = "FN",    adj=1, font=2, cex=1.1)
    text(x = col_names["TN"],    y = top_y + 1.2, labels = "TN",    adj=1, font=2, cex=1.1)
  }
  
  # Print the Mean header explicitly for both panels on the right side
  top_y <- max(df$y, na.rm=TRUE)
  text(x = 1.05, y = top_y + 1.2, labels = "Mean [95% CrI]", adj=0, font=2, cex=1.1)

  # For each row, print the author/counts explicitly aligned
  # The axis() yaxt is replaced entirely by this text generation.
  for(i in seq_along(df$y)) {
    yi <- df$y[i]
    if(!is.na(yi)) {
      if (show_left_text) {
        text(x = col_names["Study"], y = yi, labels = as.character(df$label[i]), adj=0, cex=1.1)
        text(x = col_names["TP"],    y = yi, labels = df$TP[i], adj=1, cex=1.1)
        text(x = col_names["FP"],    y = yi, labels = df$FP[i], adj=1, cex=1.1)
        text(x = col_names["FN"],    y = yi, labels = df$FN[i], adj=1, cex=1.1)
        text(x = col_names["TN"],    y = yi, labels = df$TN[i], adj=1, cex=1.1)
      }
      
      # Add Mean [95% CI] numbers to the right of the plot line
      right_text <- sprintf("%.2f [%.2f, %.2f]", m[i], lo[i], hi[i])
      text(x = 1.05, y = yi, labels = right_text, pos=4, cex=1.1, col=as.character(df$col[i]))
    }
  }
  par(xpd=FALSE)
  
  # Colour legend
  legend("bottomleft", legend=fp$arm_labs, col=fp$arm_colors,
         pch=19, bty="n", cex=0.85, pt.cex=1.1)
}

tiff(paste0("AI_NMA_ForestPlot_Full.tiff"),
     width=20.0, height=max(10, nrow(fp_data$df)*0.24 + 2.5),
     units="in", res=300, compression="lzw", type="cairo")
layout(matrix(c(1,2), nrow=1), widths=c(13.4, 6.6))
plot_forest_panel(fp_data, "Se", "Sensitivity (observed / posterior)", "Sensitivity", show_left_text=TRUE)
plot_forest_panel(fp_data, "Sp", "Specificity (observed / posterior)", "Specificity", show_left_text=FALSE)
dev.off()
cat("Full forest plot saved: AI_NMA_ForestPlot_Full.tiff\n")

# Keep simple arm-level forest plot too
sens_forest_plot <- ggplot(sens_plot_data,
    aes(x=test, y=mean, ymin=`2.5%`, ymax=`97.5%`, colour=test)) +
  geom_pointrange(size=0.9, linewidth=1.2, fatten=2) +
  scale_colour_manual(values=c("DL only"=arm_cols[1],
                               "Human only"=arm_cols[2],"Human + DL"=arm_cols[3]),
                      guide="none") +
  coord_flip() + theme_bw(base_size=13, base_family="Times") + ylim(0.5,1) +
  labs(title="Sensitivity (Pooled NMA)", x=NULL, y="Posterior Mean \u00b1 95% CrI") +
  theme(plot.title=element_text(hjust=0.5,face="bold"))

spec_forest_plot <- ggplot(spec_plot_data,
    aes(x=test, y=mean, ymin=`2.5%`, ymax=`97.5%`, colour=test)) +
  geom_pointrange(size=0.9, linewidth=1.2, fatten=2) +
  scale_colour_manual(values=c("DL only"=arm_cols[1],
                               "Human only"=arm_cols[2],"Human + DL"=arm_cols[3]),
                      guide="none") +
  coord_flip() + theme_bw(base_size=13, base_family="Times") + ylim(0.5,1) +
  labs(title="Specificity (Pooled NMA)", x=NULL, y=NULL) +
  theme(plot.title=element_text(hjust=0.5,face="bold"))

tiff("AI_NMA_ForestPlot_Pooled.tiff", width=10, height=4, units="in", res=300, compression="lzw", type="cairo")
print(sens_forest_plot + spec_forest_plot + plot_layout(ncol=2))
dev.off()
cat("Pooled forest plot saved: AI_NMA_ForestPlot_Pooled.tiff\n")

# Save model object
# saveRDS already handled during model fitting

cat("Base NMA model saved: AI_NMA_model.rds\n")

################################################################################
# PART A: SUMMARY TABLE -- Sensitivity, Specificity, Youden by Test Arm
################################################################################
cat("\n=== PART A: Accuracy Summary ===\n")

test_labels <- c("DL only", "Human only", "Human + DL")
feat <- mod.fit.R2jags$BUGSoutput$summary

sens_rows <- feat[grep("^sens\\[",  rownames(feat)), ]
spec_rows <- feat[grep("^spec\\[",  rownames(feat)), ]
youd_rows <- feat[grep("^youden\\[",rownames(feat)), ]

accuracy_table <- data.frame(
  Test        = test_labels,
  Sensitivity = sprintf("%.3f [%.3f, %.3f]", sens_rows[,"mean"], sens_rows[,"2.5%"], sens_rows[,"97.5%"]),
  Specificity = sprintf("%.3f [%.3f, %.3f]", spec_rows[,"mean"], spec_rows[,"2.5%"], spec_rows[,"97.5%"]),
  Youden      = sprintf("%.3f [%.3f, %.3f]", youd_rows[,"mean"], youd_rows[,"2.5%"], youd_rows[,"97.5%"])
)
print(accuracy_table)
write.csv(accuracy_table, "AI_NMA_Accuracy_Summary.csv", row.names = FALSE)

# SUCRA
sucra_sims <- mod.fit.R2jags$BUGSoutput$sims.list$youdensucra
sucra_vals <- sapply(1:ntest_arms, function(k)
  mean(sapply(1:(ntest_arms-1), function(m) sucra_sims[,k,m])))
cat("SUCRA (Youden):", round(sucra_vals,3), "\n")

# Heterogeneity — SD, tau, and I² (logit-normal approximation: I² = τ²/(τ² + π²/3))
sims <- mod.fit.R2jags$BUGSoutput$sims.list
tau_se  <- sims$SDstudysens;  tau_sp  <- sims$SDstudyspec
tau2_se <- tau_se^2;           tau2_sp <- tau_sp^2
I2_se   <- tau2_se / (tau2_se + (pi^2/3))
I2_sp   <- tau2_sp / (tau2_sp + (pi^2/3))
het <- data.frame(
  Metric = c("Tau (between-study SD) Se", "Tau (between-study SD) Sp",
             "SD within-study Se",        "SD within-study Sp",
             "I² (between-study) Se",      "I² (between-study) Sp"),
  Mean   = c(mean(tau_se),  mean(tau_sp),
             mean(sims$SDtestsens), mean(sims$SDtestspec),
             mean(I2_se),   mean(I2_sp)),
  Lo95   = c(quantile(tau_se,0.025),  quantile(tau_sp,0.025),
             quantile(sims$SDtestsens,0.025), quantile(sims$SDtestspec,0.025),
             quantile(I2_se,0.025),  quantile(I2_sp,0.025)),
  Hi95   = c(quantile(tau_se,0.975),  quantile(tau_sp,0.975),
             quantile(sims$SDtestsens,0.975), quantile(sims$SDtestspec,0.975),
             quantile(I2_se,0.975),  quantile(I2_sp,0.975))
)
print(het)
write.csv(het, "AI_NMA_Heterogeneity.csv", row.names = FALSE)

################################################################################
# PART B: HSROC CURVE (publication quality)
################################################################################
cat("\n=== PART B: HSROC Curve ===\n")

inv_logit <- function(x) 1/(1+exp(-x))

make_ellipse <- function(muSe, muSp, Sigma, level = 0.95, n = 300) {
  cc  <- sqrt(qchisq(level, 2))
  th  <- seq(0, 2*pi, length.out = n)
  eig <- eigen(Sigma, symmetric = TRUE)
  A   <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)))
  XY  <- sweep(t(A %*% (cc * rbind(cos(th), sin(th)))), 2, c(muSe, muSp), "+")
  data.frame(FPR = 1 - inv_logit(XY[,2]), TPR = inv_logit(XY[,1]))
}

thr_se <- sims$threshold.sens   # iterations x ntest x 1
thr_sp <- sims$threshold.spec

arm_cols  <- c("#2166AC","#D6604D","#1A9641")   # blue / red / green
arm_lty   <- c(1L, 2L, 4L)

osa_total_n <- osa$tp1 + osa$fn1 + osa$fp1 + osa$tn1   # study size
osa_obs_fpr <- 1 - osa$obs_Sp   # observed FPR per row
osa_obs_tpr <- osa$obs_Se        # observed TPR per row
# Scale point size: sqrt scaling so large studies aren't overwhelming
# Point-size cap (Fig3 fix #10): smaller largest point so big studies don't swamp the cluster
osa_cex <- 0.4 + 1.0 * sqrt(osa_total_n / max(osa_total_n, na.rm=TRUE))
# Capture patient test_num as a distinct vector (osa is reassigned to eye data later at L~1038)
osa_test_num <- osa$test_num

## -----------------------------------------------------------------------------
## Reusable HSROC panel drawer (single source of truth for both patient & eye
## panels and the combined Fig3 device). All inputs passed explicitly so it is
## immune to later reassignment of `osa`. Encodes Fig3 SCI fixes:
##   #9 fill alpha 0.20 + outline alpha 0.55; #11 outline lwd 0.5;
##   #12 padded xlim/ylim (no xaxs/yaxs="i" clipping); #13 grey gridlines;
##   #16 50% credible ellipse for parity.
## -----------------------------------------------------------------------------
draw_hsroc_panel <- function(test_num, obs_fpr, obs_tpr, pt_cex,
                             thr_se, thr_sp, ntest, arm_cols, arm_lty,
                             test_labels, main_title) {
  par(mar=c(5,5,3,2), family="Times")
  plot(NA, xlim=c(-0.02,1.02), ylim=c(-0.02,1.02), las=1,
       xaxt="n", yaxt="n",
       xlab="1 \u2212 Specificity  (FPR)", ylab="Sensitivity  (TPR)",
       main=main_title,
       cex.axis=1.1, cex.lab=1.2, cex.main=1.3)
  axis(1, at=seq(0,1,0.2), cex.axis=1.1)
  axis(2, at=seq(0,1,0.2), las=1, cex.axis=1.1)
  abline(0,1, lty=3, col="grey70", lwd=1.2)
  abline(h=c(0.2,0.4,0.6,0.8), v=c(0.2,0.4,0.6,0.8), col="grey90", lwd=0.6)
  # --- Individual study data points (behind ellipses) ---
  for(k in 1:ntest) {
    idx_k <- which(test_num == k)
    points(obs_fpr[idx_k], obs_tpr[idx_k],
           pch=21,
           bg  = adjustcolor(arm_cols[k], alpha.f=0.20),
           col = adjustcolor(arm_cols[k], alpha.f=0.55),
           cex = pt_cex[idx_k],
           lwd = 0.5)
  }
  # --- Posterior confidence ellipses + summary points (on top) ---
  for(k in 1:ntest) {
    se_draw <- thr_se[,k,1]; sp_draw <- thr_sp[,k,1]
    muSe <- mean(se_draw); muSp <- mean(sp_draw)
    Sig  <- cov(cbind(se_draw, sp_draw))
    ell  <- make_ellipse(muSe, muSp, Sig, level=0.95)
    ell2 <- make_ellipse(muSe, muSp, Sig, level=0.50)
    lines(ell$FPR,  ell$TPR,  col=arm_cols[k], lty=arm_lty[k], lwd=2)
    lines(ell2$FPR, ell2$TPR, col=arm_cols[k], lty=arm_lty[k], lwd=1)
    points(1-inv_logit(muSp), inv_logit(muSe),
           pch=21, bg=arm_cols[k], col="white", cex=2.5, lwd=2)
  }
  legend("bottomright",
         legend=c(test_labels, "Individual study (size \u221a N)"),
         col=c(arm_cols, "grey40"),
         lty=c(arm_lty, NA), pch=c(21,21,21,21),
         pt.bg=c(arm_cols, "grey80"),
         pt.cex=c(1.5,1.5,1.5,1.2),
         lwd=2, bty="n", cex=0.95, title="Test arm")
}

## -----------------------------------------------------------------------------
## Reusable meta-regression (Fig4) single-subpanel drawer. One source of truth
## for the patient & eye standalone panels and the combined Fig4 device. Encodes
## SCI fixes: #2 8% x-axis padding; #3 label y-offset -0.28 + cex param;
## #4 top headroom ylim +2.2; #5 corner mtext line -0.6/-1.5 cex 0.85;
## #6 left covariate-label cex 0.85. Caller sets par(mar=...) (#7) and layout.
## -----------------------------------------------------------------------------
draw_mr_panel <- function(means, los, his, cols, labels, n, y,
                          main_title, lab_cex=0.75) {
  r <- range(c(los, his, 0), na.rm=TRUE)
  xlim <- r + c(-1, 1) * diff(r) * 0.08
  plot(means, y, pch=19, col=cols,
       xlim=xlim, ylim=c(0.5, n + 2.2),
       yaxt="n", xlab="Beta coefficient", ylab="",
       main=main_title, cex=1.1, cex.axis=1, cex.lab=1.1, cex.main=1.2)
  segments(los, y, his, y, col=cols, lwd=1.8)
  abline(v=0, lty=2, col="grey50", lwd=1.2)
  axis(2, at=y, labels=labels, las=2, cex.axis=0.85, tick=FALSE)
  par(xpd=NA)
  for(i in seq_len(n)) {
    est_str <- sprintf("%.2f [%.2f, %.2f]", means[i], los[i], his[i])
    text(means[i], y[i] - 0.28, est_str, adj=c(0.5, 0.5), cex=lab_cex, col=cols[i])
  }
  par(xpd=FALSE)
  mtext("Red = Significant", side=3, line=-0.6, adj=0.02, cex=0.85, col="red")
  mtext("Line = 95% Credible Interval", side=3, line=-1.5, adj=0.02, cex=0.85, col="black")
}

tiff("AI_NMA_HSROC.tiff", width=8, height=8, units="in", res=300, compression="lzw", type="cairo")
draw_hsroc_panel(osa_test_num, osa_obs_fpr, osa_obs_tpr, osa_cex,
                 thr_se, thr_sp, ntest_arms, arm_cols, arm_lty, test_labels,
                 "HSROC \u2013 DL only  vs  Human only  vs  Human\u2009+\u2009DL")
dev.off()
cat("HSROC saved: AI_NMA_HSROC.tiff\n")

################################################################################
# PART B2: MCMC DIAGNOSTICS (trace + density)
################################################################################
cat("\n=== PART B2: MCMC Diagnostics ===\n")

library(lattice)

mod.fit.R2jags.mcmc <- as.mcmc(mod.fit.R2jags)

key_params <- c("sens[1,1]","sens[2,1]","sens[3,1]",
                "spec[1,1]","spec[2,1]","spec[3,1]",
                "SDstudysens","SDstudyspec")

# Trace plots
cairo_pdf("AI_NMA_Traceplot.pdf", width=14, height=10)
print(xyplot(mod.fit.R2jags.mcmc[, key_params],
             xlab="Iteration", ylab="Value",
             main="Trace plots \u2013 Key Parameters",
             par.strip.text=list(cex=0.85),
             scales=list(y=list(relation="free"))))
dev.off()
cat("Trace plot saved: AI_NMA_Traceplot.pdf\n")

# Density plots
cairo_pdf("AI_NMA_Densityplot.pdf", width=14, height=10)
print(densityplot(mod.fit.R2jags.mcmc[, key_params],
                  main="Posterior Density \u2013 Key Parameters",
                  par.strip.text=list(cex=0.85),
                  scales=list(y=list(relation="free"))))
dev.off()
cat("Density plot saved: AI_NMA_Densityplot.pdf\n")

# Gelman-Rubin convergence plot
cairo_pdf("AI_NMA_GelmanPlot.pdf", width=8, height=10)
mcmc_sub <- mod.fit.R2jags.mcmc[, key_params, drop=FALSE]
gelman.plot(mcmc_sub, main="", auto.layout=TRUE)
dev.off()
cat("Gelman-Rubin plot saved: AI_NMA_GelmanPlot.pdf\n")

################################################################################
# PART C: UNIVARIATE META-REGRESSION  (categorical, vs baseline level 1)
################################################################################
# PART C: UNIVARIATE META-REGRESSION  (categorical, vs baseline level 1)
################################################################################
cat("\n=== PART C: Univariate Meta-Regression (categorical) ===\n")

candidate_covs <- c("reader", "study", "study_design", "multicenter", 
                    "economic", "healthcare", "external", "dr", "othereyetarget",
                    "pupil", "ctype", "camera", "criteria", "analysis_type", 
                    "rstandard", "disagree", "rimage", "field", "vendor", 
                    "commercial", "certified", "quality_check", "architecture", 
                    "region", "aigradable", "hgradable")

# Categorical NMA meta-regression model
# beta_Se[1] = 0 (reference),  beta_Se[k] ~ dnorm(0,0.05) for k>=2
univ_cat_model <- "
model {
  for(i in 1:nObs) {
    tp[i] ~ dbin(pi[i,1], pos[i])
    tn[i] ~ dbin(pi[i,2], neg[i])
    logit(pi[i,1]) <- mu[i,1]
    logit(pi[i,2]) <- mu[i,2]
    MU[i,1] <- threshold.sens[test[i],1] + study.re.sens[s[i]] +
                test.re.sens[s[i],test[i]] + beta_Se[cov_x[i]]
    MU[i,2] <- threshold.spec[test[i],1] + study.re.spec[s[i]] +
                test.re.spec[s[i],test[i]] + beta_Sp[cov_x[i]]
    mu[i,1:2] ~ dmnorm(MU[i,], Omega[,])
  }
  beta_Se[1] <- 0;  beta_Sp[1] <- 0
  for(k in 2:n_lev){ beta_Se[k]~dnorm(0,0.05); beta_Sp[k]~dnorm(0,0.05) }
  for(k in 1:ntest){ threshold.sens[k,1]~dnorm(0,0.05); threshold.spec[k,1]~dnorm(0,0.05) }
  for(k in 1:ns){
    study.re.sens[k]~dnorm(0,taustudysens); study.re.spec[k]~dnorm(0,taustudyspec)
    for(l in 1:ntest){ test.re.sens[k,l]~dnorm(0,tautestsens); test.re.spec[k,l]~dnorm(0,tautestspec) }
  }
  taustudysens<-pow(SDstudysens,-2); SDstudysens~dunif(0,2)
  taustudyspec<-pow(SDstudyspec,-2); SDstudyspec~dunif(0,2)
  tautestsens <-pow(SDtestsens,-2);  SDtestsens ~dunif(0,2)
  tautestspec <-pow(SDtestspec,-2);  SDtestspec ~dunif(0,2)
  Omega[1:2,1:2]<-inverse(Sigma.sq[,])
  for(m in 1:2){Sigma.sq[m,m]<-pow(sd[m],2)}
  for(i in 1:2){for(j in (i+1):2){Sigma.sq[i,j]<-rho[i,j]*sd[i]*sd[j]; Sigma.sq[j,i]<-Sigma.sq[i,j]}}
  for(m in 1:2){sd[m]~dunif(0,1)}
  for(i in 1:2){for(j in (i+1):2){g[j,i]<-0; a[i,j]~dunif(0,3.1415); rho[i,j]<-inprod(g[,i],g[,j])}}
  g[1,1]<-1; g[1,2]<-cos(a[1,2]); g[2,2]<-sin(a[1,2])
}
"

univ_results <- data.frame()
for(cov_name in candidate_covs) {
  cat("  Testing:", cov_name, "\n")
  cov_raw <- osa[[cov_name]]
  if(is.null(cov_raw) || all(is.na(cov_raw))) { cat("  -> missing, skip\n"); next }
  # Convert to factor, fill NA with "NR" (own category)
  cov_raw[is.na(cov_raw)] <- "NR"
  cov_factor <- as.factor(cov_raw)
  freq_table <- sort(table(cov_factor), decreasing = TRUE)
  ref_level <- names(freq_table)[1]
  cov_factor <- relevel(cov_factor, ref = ref_level)
  cov_levels <- levels(cov_factor)
  n_lev <- length(cov_levels)
  if(n_lev < 2) { cat("  -> constant, skip\n"); next }

  dl <- c(data_list, list(cov_x = as.numeric(cov_factor), n_lev = n_lev))
  in1 <- c(inits1); in2 <- c(inits2)

  tryCatch({
    fit <- jags(data=dl, inits=list(in1, in2),
                parameters.to.save=c("beta_Se","beta_Sp"),
                n.chains=2, n.iter=n_iter, n.burnin=n_burnin,
                model.file=textConnection(univ_cat_model))
    s <- fit$BUGSoutput$summary
    # Extract k>=2 (vs reference level 1)
    for(k in 2:n_lev) {
      rn_se <- paste0("beta_Se[",k,"]"); rn_sp <- paste0("beta_Sp[",k,"]")
      bs  <- s[rn_se,]; bsp <- s[rn_sp,]
      univ_results <- rbind(univ_results, data.frame(
        Covariate = cov_name,
        Level     = paste0(cov_levels[k], " vs ", cov_levels[1]),
        Se_Mean   = bs["mean"],  Se_Lo = bs["2.5%"],  Se_Hi = bs["97.5%"],
        Se_Sig    = (bs["2.5%"] > 0 | bs["97.5%"] < 0),
        Se_Pval   = 2 * min(mean(fit$BUGSoutput$sims.list$beta_Se[,k] > 0),
                            mean(fit$BUGSoutput$sims.list$beta_Se[,k] < 0)),
        Sp_Mean   = bsp["mean"], Sp_Lo = bsp["2.5%"], Sp_Hi = bsp["97.5%"],
        Sp_Sig    = (bsp["2.5%"] > 0 | bsp["97.5%"] < 0),
        Sp_Pval   = 2 * min(mean(fit$BUGSoutput$sims.list$beta_Sp[,k] > 0),
                            mean(fit$BUGSoutput$sims.list$beta_Sp[,k] < 0))
      ))
    }
  }, error = function(e) cat("  -> Error:", e$message, "\n"))
}

write.csv(univ_results, "AI_NMA_Univariate_MR.csv", row.names=FALSE)
cat("Univariate MR saved: AI_NMA_Univariate_MR.csv\n")
print(univ_results)

# --- Publication-quality univariate forest plot ---
if(nrow(univ_results) > 0) {
  label_u <- paste0(univ_results$Covariate, ": ", univ_results$Level)
  n_u <- nrow(univ_results); y_u <- n_u:1

  tiff("AI_NMA_Univariate_MR.tiff", width=14, height=max(6, n_u*0.38+2),
       units="in", res=300, compression="lzw", type="cairo")
  par(mfrow=c(1,2), mar=c(5,max(nchar(label_u))/2.5+1, 3.5, 4.5), family="Times", las=1)

  # Red for significant, Grey for non-significant
  col_u_se <- ifelse(univ_results$Se_Sig, "red", "#666666")
  col_u_sp <- ifelse(univ_results$Sp_Sig, "red", "#666666")

  # Sensitivity panel
  xlim_se <- range(c(univ_results$Se_Lo, univ_results$Se_Hi, 0), na.rm=TRUE)
  plot(univ_results$Se_Mean, y_u, pch=19, col=col_u_se,
       xlim=xlim_se, ylim=c(0.5, n_u+1.5),
       yaxt="n", xlab="Beta coefficient", ylab="",
       main="Univariate MR \u2013 Sensitivity", cex=1.1, cex.axis=1, cex.lab=1.1, cex.main=1.2)
  segments(univ_results$Se_Lo, y_u, univ_results$Se_Hi, y_u, col=col_u_se, lwd=1.8)
  abline(v=0, lty=2, col="grey50", lwd=1.2)
  axis(2, at=y_u, labels=label_u, las=2, cex.axis=0.78, tick=FALSE)
  # p-value annotation with expression for italic p
  par(xpd=NA)
  for(i in seq_len(n_u)) {
    est_str <- sprintf("%.2f [%.2f, %.2f]", univ_results$Se_Mean[i], univ_results$Se_Lo[i], univ_results$Se_Hi[i])
    text(univ_results$Se_Mean[i], y_u[i] - 0.35, est_str,
         adj=c(0.5, 0.5), cex=0.75, col=col_u_se[i])
  }
  par(xpd=FALSE)
  mtext("Red = Significant", side=3, line=-1.2, adj=0.02, cex=0.8, col="red")
  mtext("Line = 95% Credible Interval", side=3, line=-2.2, adj=0.02, cex=0.8, col="black")

  # Specificity panel
  xlim_sp <- range(c(univ_results$Sp_Lo, univ_results$Sp_Hi, 0), na.rm=TRUE)
  plot(univ_results$Sp_Mean, y_u, pch=19, col=col_u_sp,
       xlim=xlim_sp, ylim=c(0.5, n_u+1.5),
       yaxt="n", xlab="Beta coefficient", ylab="",
       main="Univariate MR \u2013 Specificity", cex=1.1, cex.axis=1, cex.lab=1.1, cex.main=1.2)
  segments(univ_results$Sp_Lo, y_u, univ_results$Sp_Hi, y_u, col=col_u_sp, lwd=1.8)
  abline(v=0, lty=2, col="grey50", lwd=1.2)
  axis(2, at=y_u, labels=label_u, las=2, cex.axis=0.78, tick=FALSE)
  # p-value annotation with expression for italic p
  par(xpd=NA)
  for(i in seq_len(n_u)) {
    est_str <- sprintf("%.2f [%.2f, %.2f]", univ_results$Sp_Mean[i], univ_results$Sp_Lo[i], univ_results$Sp_Hi[i])
    text(univ_results$Sp_Mean[i], y_u[i] - 0.35, est_str,
         adj=c(0.5, 0.5), cex=0.75, col=col_u_sp[i])
  }
  par(xpd=FALSE)
  mtext("Red = Significant", side=3, line=-1.2, adj=0.02, cex=0.8, col="red")
  mtext("Line = 95% Credible Interval", side=3, line=-2.2, adj=0.02, cex=0.8, col="black")
  dev.off()
  cat("Univariate MR plot saved: AI_NMA_Univariate_MR.tiff\n")
}

################################################################################
# PART D: MULTIVARIATE META-REGRESSION  (categorical, significant only)
################################################################################
cat("\n=== PART D: Multivariate Meta-Regression (categorical) ===\n")

sig_covs <- unique(univ_results$Covariate[univ_results$Se_Sig | univ_results$Sp_Sig])

if(length(sig_covs) == 0) {
  cat("No significant covariates – skipping multivariate step.\n")
} else {
  cat("Fitting categorical multivariate model for:", paste(sig_covs, collapse=", "), "\n")

  # Build factor indices and level counts for each sig covariate
  dl_mv <- data_list
  cov_factor_list <- list()
  n_lev_vec      <- integer()
  for(v in sig_covs) {
    x <- osa[[v]]; x[is.na(x)] <- "NR"
    fv <- as.factor(x)
    cov_factor_list[[v]] <- fv
    n_lev_vec[v] <- length(levels(fv))
    dl_mv[[paste0("cov_", v)]]   <- as.numeric(fv)
    dl_mv[[paste0("nlev_", v)]]  <- n_lev_vec[v]
  }

  # Dynamically build JAGS model string for categorical MV MR
  body_mu1 <- "    MU[i,1] <- threshold.sens[test[i],1]+study.re.sens[s[i]]+test.re.sens[s[i],test[i]]"
  body_mu2 <- "    MU[i,2] <- threshold.spec[test[i],1]+study.re.spec[s[i]]+test.re.spec[s[i],test[i]]"
  for(v in sig_covs) {
    body_mu1 <- paste0(body_mu1, "+bSe_", v, "[cov_", v, "[i]]")
    body_mu2 <- paste0(body_mu2, "+bSp_", v, "[cov_", v, "[i]]")
  }
  prior_blocks <- ""
  for(v in sig_covs) {
    prior_blocks <- paste0(prior_blocks,
      "  bSe_",v,"[1]<-0; bSp_",v,"[1]<-0\n",
      "  for(k_ in 2:nlev_",v,"){bSe_",v,"[k_]~dnorm(0,0.05);bSp_",v,"[k_]~dnorm(0,0.05)}\n"
    )
  }
  params_mv <- c(paste0("bSe_", sig_covs), paste0("bSp_", sig_covs), "SDstudysens", "SDstudyspec")

  mv_cat_model <- paste0(
    "model {\n",
    "  for(i in 1:nObs) {\n",
    "    tp[i]~dbin(pi[i,1],pos[i]); tn[i]~dbin(pi[i,2],neg[i])\n",
    "    logit(pi[i,1])<-mu[i,1]; logit(pi[i,2])<-mu[i,2]\n",
    body_mu1, "\n",
    body_mu2, "\n",
    "    mu[i,1:2]~dmnorm(MU[i,],Omega[,])\n  }\n",
    prior_blocks,
    "  for(k in 1:ntest){threshold.sens[k,1]~dnorm(0,0.05);threshold.spec[k,1]~dnorm(0,0.05)}\n",
    "  for(k in 1:ns){study.re.sens[k]~dnorm(0,taustudysens);study.re.spec[k]~dnorm(0,taustudyspec)\n",
    "    for(l in 1:ntest){test.re.sens[k,l]~dnorm(0,tautestsens);test.re.spec[k,l]~dnorm(0,tautestspec)}}\n",
    "  taustudysens<-pow(SDstudysens,-2);SDstudysens~dunif(0,2)\n",
    "  taustudyspec<-pow(SDstudyspec,-2);SDstudyspec~dunif(0,2)\n",
    "  tautestsens<-pow(SDtestsens,-2);SDtestsens~dunif(0,2)\n",
    "  tautestspec<-pow(SDtestspec,-2);SDtestspec~dunif(0,2)\n",
    "  Omega[1:2,1:2]<-inverse(Sigma.sq[,])\n",
    "  for(m in 1:2){Sigma.sq[m,m]<-pow(sd[m],2)}\n",
    "  for(i in 1:2){for(j in (i+1):2){Sigma.sq[i,j]<-rho[i,j]*sd[i]*sd[j];Sigma.sq[j,i]<-Sigma.sq[i,j]}}\n",
    "  for(m in 1:2){sd[m]~dunif(0,1)}\n",
    "  for(i in 1:2){for(j in (i+1):2){g[j,i]<-0;a[i,j]~dunif(0,3.1415);rho[i,j]<-inprod(g[,i],g[,j])}}\n",
    "  g[1,1]<-1;g[1,2]<-cos(a[1,2]);g[2,2]<-sin(a[1,2])\n}\n"
  )

  fit_mv <- jags(data=dl_mv, inits=list(inits1, inits2, inits3, inits4),
                 parameters.to.save=params_mv,
                 n.chains=4, n.iter=n_iter, n.burnin=n_burnin,
                 model.file=textConnection(mv_cat_model))

  s_mv    <- fit_mv$BUGSoutput$summary
  sims_mv <- fit_mv$BUGSoutput$sims.list
  mv_res  <- data.frame()

  for(v in sig_covs) {
    lev <- levels(cov_factor_list[[v]])
    n_k <- n_lev_vec[v]
    for(k in 2:n_k) {
      se_d <- sims_mv[[paste0("bSe_",v)]][,k]
      sp_d <- sims_mv[[paste0("bSp_",v)]][,k]
      mv_res <- rbind(mv_res, data.frame(
        Covariate = v,
        Level     = paste0(lev[k], " vs ", lev[1]),
        Se_Mean   = mean(se_d), Se_Lo=quantile(se_d,0.025), Se_Hi=quantile(se_d,0.975),
        Se_Sig    = (quantile(se_d,0.025)>0 | quantile(se_d,0.975)<0),
        Se_Pval   = 2*min(mean(se_d>0), mean(se_d<0)),
        Sp_Mean   = mean(sp_d), Sp_Lo=quantile(sp_d,0.025), Sp_Hi=quantile(sp_d,0.975),
        Sp_Sig    = (quantile(sp_d,0.025)>0 | quantile(sp_d,0.975)<0),
        Sp_Pval   = 2*min(mean(sp_d>0), mean(sp_d<0))
      ))
    }
  }

  write.csv(mv_res, "AI_NMA_Multivariate_MR.csv", row.names=FALSE)
  cat("Multivariate MR saved: AI_NMA_Multivariate_MR.csv\n")
  print(mv_res)

  # --- Extract DIC, I2, and Tau for Multivariate model ---
  dic_mv <- fit_mv$BUGSoutput$DIC
  tau_se_mv  <- sims_mv$SDstudysens;  tau_sp_mv  <- sims_mv$SDstudyspec
  tau2_se_mv <- tau_se_mv^2;          tau2_sp_mv <- tau_sp_mv^2
  I2_se_mv   <- tau2_se_mv / (tau2_se_mv + (pi^2/3))
  I2_sp_mv   <- tau2_sp_mv / (tau2_sp_mv + (pi^2/3))

  cat("\nPatient-Level Multivariate Model Metrics:\n")
  cat(sprintf("DIC: %.2f\n", dic_mv))
  cat(sprintf("Tau (Se): %.3f [%.3f, %.3f]\n", mean(tau_se_mv), quantile(tau_se_mv, 0.025), quantile(tau_se_mv, 0.975)))
  cat(sprintf("Tau (Sp): %.3f [%.3f, %.3f]\n", mean(tau_sp_mv), quantile(tau_sp_mv, 0.025), quantile(tau_sp_mv, 0.975)))
  cat(sprintf("I2 (Se): %.1f%% [%.1f%%, %.1f%%]\n", mean(I2_se_mv)*100, quantile(I2_se_mv, 0.025)*100, quantile(I2_se_mv, 0.975)*100))
  cat(sprintf("I2 (Sp): %.1f%% [%.1f%%, %.1f%%]\n", mean(I2_sp_mv)*100, quantile(I2_sp_mv, 0.025)*100, quantile(I2_sp_mv, 0.975)*100))
  
  # Base Model metrics for comparison
  dic_base <- mod.fit.R2jags$BUGSoutput$DIC
  
  comp_df <- data.frame(
    Model = c("Base NMA", "Multivariate MR"),
    DIC = c(dic_base, dic_mv),
    Tau_Se = c(sprintf("%.3f [%.3f, %.3f]", mean(tau_se), quantile(tau_se,0.025), quantile(tau_se,0.975)),
               sprintf("%.3f [%.3f, %.3f]", mean(tau_se_mv), quantile(tau_se_mv,0.025), quantile(tau_se_mv,0.975))),
    Tau_Sp = c(sprintf("%.3f [%.3f, %.3f]", mean(tau_sp), quantile(tau_sp,0.025), quantile(tau_sp,0.975)),
               sprintf("%.3f [%.3f, %.3f]", mean(tau_sp_mv), quantile(tau_sp_mv,0.025), quantile(tau_sp_mv,0.975))),
    I2_Se = c(sprintf("%.1f%% [%.1f%%, %.1f%%]", mean(I2_se)*100, quantile(I2_se,0.025)*100, quantile(I2_se,0.975)*100),
              sprintf("%.1f%% [%.1f%%, %.1f%%]", mean(I2_se_mv)*100, quantile(I2_se_mv,0.025)*100, quantile(I2_se_mv,0.975)*100)),
    I2_Sp = c(sprintf("%.1f%% [%.1f%%, %.1f%%]", mean(I2_sp)*100, quantile(I2_sp,0.025)*100, quantile(I2_sp,0.975)*100),
              sprintf("%.1f%% [%.1f%%, %.1f%%]", mean(I2_sp_mv)*100, quantile(I2_sp_mv,0.025)*100, quantile(I2_sp_mv,0.975)*100))
  )
  write.csv(comp_df, "AI_NMA_Model_Comparison.csv", row.names=FALSE)
  cat("\nModel comparison saved: AI_NMA_Model_Comparison.csv\n")


  # --- Publication-quality multivariate forest plot (draw_mr_panel defined at top level) ---
  if(nrow(mv_res) > 0) {
    label_mv <- paste0(mv_res$Covariate, ": ", mv_res$Level)
    n_mv <- nrow(mv_res); y_mv <- n_mv:1
    # Red for significant, Grey for non-significant
    col_mv_se <- ifelse(mv_res$Se_Sig, "red", "#666666")
    col_mv_sp <- ifelse(mv_res$Sp_Sig, "red", "#666666")

    tiff("AI_NMA_Multivariate_MR.tiff", width=14, height=max(7, n_mv*0.55+2.5),
         units="in", res=300, compression="lzw", type="cairo")
    par(mfrow=c(1,2), mar=c(5, 11, 3.5, 4.5), family="Times", las=1)
    draw_mr_panel(mv_res$Se_Mean, mv_res$Se_Lo, mv_res$Se_Hi, col_mv_se,
                  label_mv, n_mv, y_mv, "Multivariate MR \u2013 Sensitivity", lab_cex=0.75)
    draw_mr_panel(mv_res$Sp_Mean, mv_res$Sp_Lo, mv_res$Sp_Hi, col_mv_sp,
                  label_mv, n_mv, y_mv, "Multivariate MR \u2013 Specificity", lab_cex=0.75)
    dev.off()
    cat("Multivariate MR plot saved: AI_NMA_Multivariate_MR.tiff\n")
  }
}

################################################################################
# PART E: DEEKS' FUNNEL PLOT  (publication quality)
################################################################################
cat("\n=== PART E: Deeks' Funnel Plot ===\n")

osa$tp_c <- osa$tp1+0.5; osa$fp_c <- osa$fp1+0.5
osa$fn_c <- osa$fn1+0.5; osa$tn_c <- osa$tn1+0.5

log_DOR <- log((osa$tp_c*osa$tn_c)/(osa$fp_c*osa$fn_c))
ess     <- 4*(osa$tp_c+osa$fn_c)*(osa$fp_c+osa$tn_c) /
             (osa$tp_c+osa$fn_c+osa$fp_c+osa$tn_c)

arm_cols3 <- c(ai_only="#2166AC", human_only="#D6604D", human_with_ai="#1A9641")
pt_cols   <- arm_cols3[osa$test]

fit_dk  <- lm(log_DOR ~ I(1/sqrt(ess)), weights=ess)
pval_dk <- summary(fit_dk)$coefficients[2,4]

tiff("AI_NMA_Deeks_Funnel.tiff", width=8, height=8, units="in", res=300, compression="lzw", type="cairo")
par(mar=c(5,5,3,2), family="Times")
plot(1/sqrt(ess), log_DOR, pch=21, bg=pt_cols, col="white", cex=1.4, lwd=1.2,
     xlab="1 / \u221a(Effective Sample Size)", ylab="Log Diagnostic Odds Ratio",
     main="Deeks' Funnel Plot \u2013 DL vs Human vs Human+DL",
     cex.axis=1.1, cex.lab=1.2, cex.main=1.3, las=1)
abline(fit_dk, col="#444444", lwd=2)
abline(h=0, lty=3, col="grey60")
legend("topright", legend=c("DL only","Human only","Human + DL"),
       pch=21, pt.bg=c("#2166AC","#D6604D","#1A9641"), col="white",
       pt.cex=1.4, pt.lwd=1.2, bty="n", cex=1.05)
legend("topleft", bty="n", cex=1.05,
       legend=bquote("Deeks' test:  "*italic(p)*"\u2009=\u2009"*.(format(round(pval_dk,3), nsmall=3))))
dev.off()
cat("Deeks' Funnel Plot saved: AI_NMA_Deeks_Funnel.tiff\n")
while(!is.null(dev.list())) dev.off()

cat("\n=== patient_rr PIPELINE COMPLETE ===\n")

################################################################################
# ██████████████  eye_rr PIPELINE  ██████████████████████████████████████████
# Identical analysis on the eye-level dataset
################################################################################
cat("\n\n=== Running pipeline on eye_rr dataset ===\n")

# Output goes to the same working directory
OUTDIR <- "results/"

osa <- as.data.frame(read_excel(paste0(OUTDIR, "Data_extraction_2_v7_LeeAdded.xlsx"), sheet = "eye_rr"))
osa <- osa[!is.na(osa$tp1) & !is.na(osa$fp1) & !is.na(osa$fn1) & !is.na(osa$tn1), ]
osa$test_num  <- as.numeric(factor(osa$test, levels = c("ai_only", "human_only", "human_with_ai")))
osa$study_num <- as.numeric(as.factor(osa$dataset_id))   # group by dataset_id
cat("eye_rr rows:\n"); print(table(osa$test))

ntest_arms <- 3
data_list <- list(
  ns = length(unique(osa$study_num)), ntest = ntest_arms, totaltest = ntest_arms,
  nObs = nrow(osa), id = 1:nrow(osa),
  s = osa$study_num, test = osa$test_num, threshold = rep(1, nrow(osa)),
  tp = osa$tp1, tn = osa$tn1,
  pos = osa$tp1 + osa$fn1, neg = osa$fp1 + osa$tn1
)
ns_art <- length(unique(osa$study_num))
inits1 <- list(SDstudysens=1, SDstudyspec=1, SDtestsens=1, SDtestspec=1,
               a=matrix(c(NA,1),1,2), mu=matrix(rep(0.2,nrow(osa)*2),nrow(osa),2),
               sd=c(1,1), study.re.sens=rep(1,ns_art), study.re.spec=rep(1,ns_art))
inits2 <- list(SDstudysens=0.8, SDstudyspec=0.8, SDtestsens=0.8, SDtestspec=0.8,
               a=matrix(c(NA,1),1,2), mu=matrix(rep(0.1,nrow(osa)*2),nrow(osa),2),
               sd=c(1,1), study.re.sens=rep(1,ns_art), study.re.spec=rep(1,ns_art))
inits3 <- list(SDstudysens=1.2, SDstudyspec=1.2, SDtestsens=1.2, SDtestspec=1.2,
               a=matrix(c(NA,1),1,2), mu=matrix(rep(0.3,nrow(osa)*2),nrow(osa),2),
               sd=c(1,1), study.re.sens=rep(1,ns_art), study.re.spec=rep(1,ns_art))
inits4 <- list(SDstudysens=0.5, SDstudyspec=0.5, SDtestsens=0.5, SDtestspec=0.5,
               a=matrix(c(NA,1),1,2), mu=matrix(rep(0.0,nrow(osa)*2),nrow(osa),2),
               sd=c(1,1), study.re.sens=rep(1,ns_art), study.re.spec=rep(1,ns_art))

cat("eye_rr: Fitting NMA model ...\n")
if (file.exists(paste0(OUTDIR, "eye_rr_NMA_model.rds"))) {
  cat("Loading existing eye_rr_NMA_model.rds...\n")
  mod_eye <- readRDS(paste0(OUTDIR, "eye_rr_NMA_model.rds"))
} else {
  mod_eye <- jags(data=data_list, inits=list(inits1,inits2,inits3,inits4),
                  parameters.to.save=parameters, n.chains=4,
                  n.iter=n_iter, n.burnin=n_burnin,
                  model.file=textConnection(model_string))
  saveRDS(mod_eye, paste0(OUTDIR, "eye_rr_NMA_model.rds"))
}
cat("eye_rr NMA model saved.\n")

# --- accuracy summary ---
feat_e  <- mod_eye$BUGSoutput$summary
sims_e  <- mod_eye$BUGSoutput$sims.list
se_rows <- feat_e[grep("^sens\\[", rownames(feat_e)), ]
sp_rows <- feat_e[grep("^spec\\[", rownames(feat_e)), ]
yd_rows <- feat_e[grep("^youden\\[", rownames(feat_e)), ]
acc_eye <- data.frame(
  Test        = test_labels,
  Sensitivity = sprintf("%.3f [%.3f, %.3f]", se_rows[,"mean"], se_rows[,"2.5%"], se_rows[,"97.5%"]),
  Specificity = sprintf("%.3f [%.3f, %.3f]", sp_rows[,"mean"], sp_rows[,"2.5%"], sp_rows[,"97.5%"]),
  Youden      = sprintf("%.3f [%.3f, %.3f]", yd_rows[,"mean"], yd_rows[,"2.5%"], yd_rows[,"97.5%"])
)
write.csv(acc_eye, paste0(OUTDIR, "eye_rr_NMA_Accuracy_Summary.csv"), row.names=FALSE)
print(acc_eye)
# Heterogeneity
tau_se_e  <- sims_e$SDstudysens;  tau_sp_e  <- sims_e$SDstudyspec
tau2_se_e <- tau_se_e^2;          tau2_sp_e <- tau_sp_e^2
I2_se_e   <- tau2_se_e / (tau2_se_e + (pi^2/3))
I2_sp_e   <- tau2_sp_e / (tau2_sp_e + (pi^2/3))

het_eye <- data.frame(
  Metric = c("Tau (between-study SD) Se", "Tau (between-study SD) Sp",
             "SD within-study Se",        "SD within-study Sp",
             "I² (between-study) Se",      "I² (between-study) Sp"),
  Mean   = c(mean(tau_se_e),  mean(tau_sp_e),
             mean(sims_e$SDtestsens), mean(sims_e$SDtestspec),
             mean(I2_se_e),   mean(I2_sp_e)),
  Lo95   = c(quantile(tau_se_e,0.025),  quantile(tau_sp_e,0.025),
             quantile(sims_e$SDtestsens,0.025), quantile(sims_e$SDtestspec,0.025),
             quantile(I2_se_e,0.025),  quantile(I2_sp_e,0.025)),
  Hi95   = c(quantile(tau_se_e,0.975),  quantile(tau_sp_e,0.975),
             quantile(sims_e$SDtestsens,0.975), quantile(sims_e$SDtestspec,0.975),
             quantile(I2_se_e,0.975),  quantile(I2_sp_e,0.975))
)
print(het_eye)
write.csv(het_eye, paste0(OUTDIR, "eye_rr_NMA_Heterogeneity.csv"), row.names=FALSE)

# --- Forest plot (individual + summary) ---
plot_data_e <- as.data.frame(feat_e[grep("^sens|^spec", rownames(feat_e)),]) %>%
  mutate(Parameter=rownames(.),
         indices=str_extract_all(Parameter,"\\d+"),
         i=as.numeric(sapply(indices,`[`,1)),
         j=as.numeric(sapply(indices,`[`,2))) %>% arrange(i,j)
sens_pd_e <- plot_data_e[grep("^sens",rownames(plot_data_e)),]
spec_pd_e <- plot_data_e[grep("^spec",rownames(plot_data_e)),]
sens_pd_e$test <- factor(arm_labels, levels=rev(arm_labels))
spec_pd_e$test <- factor(arm_labels, levels=rev(arm_labels))

osa$obs_Se <- osa$tp1/(osa$tp1+osa$fn1); osa$obs_Sp <- osa$tn1/(osa$fp1+osa$tn1)
osa$obs_Se_lo <- se_lo_w(osa$obs_Se, osa$tp1+osa$fn1); osa$obs_Se_hi <- se_hi_w(osa$obs_Se, osa$tp1+osa$fn1)
osa$obs_Sp_lo <- se_lo_w(osa$obs_Sp, osa$fp1+osa$tn1); osa$obs_Sp_hi <- se_hi_w(osa$obs_Sp, osa$fp1+osa$tn1)
osa$arm_col <- arm_cols[osa$test_num]
osa$study_label <- osa$authoryr
osa$plot_id <- paste0(osa$authoryr, "_", osa$index)
osa$year_num <- as.numeric(stringr::str_extract(osa$authoryr, "\\d{4}"))
osa_sorted_e <- osa[order(osa$test_num, -osa$year_num, osa$authoryr), ]
osa_sorted_e$plot_id <- factor(osa_sorted_e$plot_id, levels=rev(unique(osa_sorted_e$plot_id)))

fp_eye <- make_full_forest(osa_sorted_e,
                            sens_pd_e[,c("mean","2.5%","97.5%")],
                            spec_pd_e[,c("mean","2.5%","97.5%")],
                            arm_labels, arm_cols)

tiff(paste0(OUTDIR,"eye_rr_NMA_ForestPlot_Full.tiff"),
     width=20.0, height=max(6, nrow(fp_eye$df)*0.20 + 1.2), units="in", res=300, compression="lzw", type="cairo")
layout(matrix(c(1,2), nrow=1), widths=c(13.4, 6.6))
plot_forest_panel(fp_eye, "Se", "Sensitivity (observed / posterior)", "Sensitivity", show_left_text=TRUE)
plot_forest_panel(fp_eye, "Sp", "Specificity (observed / posterior)", "Specificity", show_left_text=FALSE)
dev.off()

# Also pooled plot
se_fp_e <- ggplot(sens_pd_e, aes(x=test,y=mean,ymin=`2.5%`,ymax=`97.5%`,colour=test)) +
  geom_pointrange(size=0.9,linewidth=1.2, fatten=2) +
  scale_colour_manual(values=setNames(arm_cols,arm_labels),guide="none") +
  coord_flip() + theme_bw(base_size=13, base_family="Times") + ylim(0.5,1) +
  labs(title="Sensitivity (Pooled NMA)", x=NULL, y="Posterior Mean \u00b1 95% CrI") +
  theme(plot.title=element_text(hjust=0.5,face="bold"))
sp_fp_e <- ggplot(spec_pd_e, aes(x=test,y=mean,ymin=`2.5%`,ymax=`97.5%`,colour=test)) +
  geom_pointrange(size=0.9,linewidth=1.2, fatten=2) +
  scale_colour_manual(values=setNames(arm_cols,arm_labels),guide="none") +
  coord_flip() + theme_bw(base_size=13, base_family="Times") + ylim(0.5,1) +
  labs(title="Specificity (Pooled NMA)", x=NULL, y=NULL) +
  theme(plot.title=element_text(hjust=0.5,face="bold"))
tiff(paste0(OUTDIR,"eye_rr_NMA_ForestPlot_Pooled.tiff"), width=10, height=4, units="in", res=300, compression="lzw", type="cairo")
print(se_fp_e + sp_fp_e + plot_layout(ncol=2)); dev.off()

## --- COMBINED Fig2 (single stacked patchwork): (A) patient-level, (B) eye-level ---
## Generated in-R so the submitted composite carries the Times font, ylim(0.5,1),
## fatten, and real (A)/(B) tags with no external downscale/assembly.
if(exists("sens_forest_plot") && exists("spec_forest_plot") &&
   exists("se_fp_e") && exists("sp_fp_e")) {
  fig2 <- (sens_forest_plot + spec_forest_plot) / (se_fp_e + sp_fp_e) +
          plot_layout(heights=c(1,1)) +
          plot_annotation(tag_levels='A',
                          theme=theme(plot.tag=element_text(face='bold', size=14, family='Times')))
  tiff("AI_NMA_Fig2.tiff", width=10, height=8, units="in",
       res=600, compression="lzw", type="cairo")
  print(fig2); dev.off()
  cat("Combined pooled forest (Fig2) saved: AI_NMA_Fig2.tiff\n")
}

# --- HSROC ---
thr_se_e <- sims_e$threshold.sens; thr_sp_e <- sims_e$threshold.spec

osa_total_n_e <- osa$tp1 + osa$fn1 + osa$fp1 + osa$tn1
osa_obs_fpr_e <- 1 - osa$obs_Sp
osa_obs_tpr_e <- osa$obs_Se
# Point-size cap (Fig3 fix #10) to match the patient panel
osa_cex_e     <- 0.4 + 1.0 * sqrt(osa_total_n_e / max(osa_total_n_e, na.rm=TRUE))
osa_test_num_e <- osa$test_num   # eye test_num (osa currently holds eye data)

# Standalone eye-level HSROC (via shared drawer = same SCI fixes as patient panel)
tiff(paste0(OUTDIR,"eye_rr_NMA_HSROC.tiff"), width=8, height=8, units="in", res=300, compression="lzw", type="cairo")
draw_hsroc_panel(osa_test_num_e, osa_obs_fpr_e, osa_obs_tpr_e, osa_cex_e,
                 thr_se_e, thr_sp_e, ntest_arms, arm_cols, arm_lty, test_labels,
                 "HSROC \u2013 DL only  vs  Human only  vs  Human\u2009+\u2009DL")
dev.off(); cat("eye_rr HSROC saved.\n")

## --- COMBINED Fig3 (single stacked device): (A) patient-level, (B) eye-level ---
## Single in-R device so the submitted composite carries every SCI fix and the
## (A)/(B) tags, with no externally-assembled divergence.
tiff("AI_NMA_HSROC_combined.tiff", width=8, height=16, units="in",
     res=300, compression="lzw", type="cairo")
layout(matrix(1:2, ncol=1))
draw_hsroc_panel(osa_test_num, osa_obs_fpr, osa_obs_tpr, osa_cex,
                 thr_se, thr_sp, ntest_arms, arm_cols, arm_lty, test_labels,
                 "HSROC \u2013 DL only  vs  Human only  vs  Human\u2009+\u2009DL")
mtext("(A)", side=3, line=1, adj=0, font=2, cex=1.3, family="Times")
draw_hsroc_panel(osa_test_num_e, osa_obs_fpr_e, osa_obs_tpr_e, osa_cex_e,
                 thr_se_e, thr_sp_e, ntest_arms, arm_cols, arm_lty, test_labels,
                 "HSROC \u2013 DL only  vs  Human only  vs  Human\u2009+\u2009DL")
mtext("(B)", side=3, line=1, adj=0, font=2, cex=1.3, family="Times")
layout(1)
dev.off(); cat("Combined HSROC (Fig3) saved: AI_NMA_HSROC_combined.tiff\n")

# --- MCMC diagnostics ---
mod_eye_mcmc <- as.mcmc(mod_eye)
cairo_pdf(paste0(OUTDIR,"eye_rr_NMA_Traceplot.pdf"), width=14, height=10)
print(xyplot(mod_eye_mcmc[,key_params], main="Trace plots \u2013 Key Parameters",
             scales=list(y=list(relation="free")))); dev.off()
cairo_pdf(paste0(OUTDIR,"eye_rr_NMA_Densityplot.pdf"), width=14, height=10)
print(densityplot(mod_eye_mcmc[,key_params], main="Posterior Density \u2013 Key Parameters",
                  scales=list(y=list(relation="free")))); dev.off()
cairo_pdf(paste0(OUTDIR,"eye_rr_NMA_GelmanPlot.pdf"), width=8, height=10)
mcmc_sub_eye <- mod_eye_mcmc[, key_params, drop=FALSE]
gelman.plot(mcmc_sub_eye, main="", auto.layout=TRUE); dev.off()
cat("eye_rr Gelman-Rubin plot saved.\n")

# --- Univariate MR (categorical) ---
cat("\neye_rr Univariate Meta-Regression...\n")
univ_eye <- data.frame()
for(cov_name in candidate_covs) {
  cov_raw <- osa[[cov_name]]
  if(is.null(cov_raw)||all(is.na(cov_raw))){next}
  cov_raw[is.na(cov_raw)] <- "NR"
  cov_factor <- as.factor(cov_raw)
  freq_table <- sort(table(cov_factor), decreasing = TRUE)
  ref_level <- names(freq_table)[1]
  cov_factor <- relevel(cov_factor, ref = ref_level)
  cov_levels <- levels(cov_factor)
  n_lev <- length(cov_levels)
  if(n_lev < 2) next
  dl <- c(data_list, list(cov_x=as.numeric(cov_factor), n_lev=n_lev))
  tryCatch({
    fit <- jags(data=dl, inits=list(inits1,inits2,inits3,inits4),
                parameters.to.save=c("beta_Se","beta_Sp"),
                n.chains=4, n.iter=n_iter, n.burnin=n_burnin,
                model.file=textConnection(univ_cat_model))
    s <- fit$BUGSoutput$summary
    for(k in 2:n_lev) {
      bs <- s[paste0("beta_Se[",k,"]"),]; bsp <- s[paste0("beta_Sp[",k,"]"),]
      univ_eye <- rbind(univ_eye, data.frame(
        Covariate=cov_name, Level=paste0(cov_levels[k]," vs ",cov_levels[1]),
        Se_Mean=bs["mean"], Se_Lo=bs["2.5%"], Se_Hi=bs["97.5%"],
        Se_Sig=(bs["2.5%"]>0|bs["97.5%"]<0),
        Se_Pval=2*min(mean(fit$BUGSoutput$sims.list$beta_Se[,k]>0),
                      mean(fit$BUGSoutput$sims.list$beta_Se[,k]<0)),
        Sp_Mean=bsp["mean"], Sp_Lo=bsp["2.5%"], Sp_Hi=bsp["97.5%"],
        Sp_Sig=(bsp["2.5%"]>0|bsp["97.5%"]<0),
        Sp_Pval=2*min(mean(fit$BUGSoutput$sims.list$beta_Sp[,k]>0),
                      mean(fit$BUGSoutput$sims.list$beta_Sp[,k]<0))
      ))
    }
  }, error=function(e) cat("  ->",cov_name,"error:",e$message,"\n"))
}
    write.csv(univ_eye, paste0(OUTDIR,"eye_rr_NMA_Univariate_MR.csv"), row.names=FALSE)

if(nrow(univ_eye)>0){
  label_e <- paste0(univ_eye$Covariate,": ",univ_eye$Level)
  n_e <- nrow(univ_eye); y_e <- n_e:1
  col_e_se <- ifelse(univ_eye$Se_Sig,"red","#666666")
  col_e_sp <- ifelse(univ_eye$Sp_Sig,"red","#666666")
  tiff(paste0(OUTDIR,"eye_rr_NMA_Univariate_MR.tiff"),
       width=14, height=max(6,n_e*0.38+2), units="in", res=300, compression="lzw", type="cairo")
  par(mfrow=c(1,2), mar=c(5,max(nchar(label_e))/2.5+1,3.5,1), family="Times", las=1)
  
  xlim_se <- range(c(univ_eye$Se_Lo, univ_eye$Se_Hi, 0), na.rm=TRUE)
  plot(univ_eye$Se_Mean,y_e,pch=19,col=col_e_se,xlim=xlim_se,ylim=c(0.5,n_e+1.5),
       yaxt="n",xlab="Beta coefficient",ylab="",main="Univariate MR \u2013 Sensitivity",
       cex=1.1,cex.axis=1,cex.lab=1.1,cex.main=1.2)
  segments(univ_eye$Se_Lo,y_e,univ_eye$Se_Hi,y_e,col=col_e_se,lwd=1.8)
  abline(v=0,lty=2,col="grey50",lwd=1.2)
  axis(2,at=y_e,labels=label_e,las=2,cex.axis=0.78,tick=FALSE)
  par(xpd=NA)
  for(i in seq_len(n_e)) {
    est_str <- sprintf("%.2f [%.2f, %.2f]", univ_eye$Se_Mean[i], univ_eye$Se_Lo[i], univ_eye$Se_Hi[i])
    text(univ_eye$Se_Mean[i], y_e[i] - 0.35, est_str,
         adj=c(0.5, 0.5), cex=0.75, col=col_e_se[i])
  }
  par(xpd=FALSE)
  mtext("Red = Significant", side=3, line=-1.2, adj=0.02, cex=0.8, col="red")
  mtext("Line = 95% Credible Interval", side=3, line=-2.2, adj=0.02, cex=0.8, col="black")

  xlim_sp <- range(c(univ_eye$Sp_Lo, univ_eye$Sp_Hi, 0), na.rm=TRUE)
  plot(univ_eye$Sp_Mean,y_e,pch=19,col=col_e_sp,xlim=xlim_sp,ylim=c(0.5,n_e+1.5),
       yaxt="n",xlab="Beta coefficient",ylab="",main="Univariate MR \u2013 Specificity",
       cex=1.1,cex.axis=1,cex.lab=1.1,cex.main=1.2)
  segments(univ_eye$Sp_Lo,y_e,univ_eye$Sp_Hi,y_e,col=col_e_sp,lwd=1.8)
  abline(v=0,lty=2,col="grey50",lwd=1.2)
  axis(2,at=y_e,labels=label_e,las=2,cex.axis=0.78,tick=FALSE)
  par(xpd=NA)
  for(i in seq_len(n_e)) {
    est_str <- sprintf("%.2f [%.2f, %.2f]", univ_eye$Sp_Mean[i], univ_eye$Sp_Lo[i], univ_eye$Sp_Hi[i])
    text(univ_eye$Sp_Mean[i], y_e[i] - 0.35, est_str,
         adj=c(0.5, 0.5), cex=0.75, col=col_e_sp[i])
  }
  par(xpd=FALSE)
  mtext("Red = Significant", side=3, line=-1.2, adj=0.02, cex=0.8, col="red")
  mtext("Line = 95% Credible Interval", side=3, line=-2.2, adj=0.02, cex=0.8, col="black")
  
  dev.off()
  cat("eye_rr Univariate MR plot saved.\n")
}

################################################################################
# PART D: MULTIVARIATE META-REGRESSION (categorical, eye-level)
################################################################################
cat("\neye_rr Multivariate Meta-Regression (categorical)...\n")

sig_covs_e <- unique(univ_eye$Covariate[univ_eye$Se_Sig | univ_eye$Sp_Sig])

if(length(sig_covs_e) == 0) {
  cat("No significant covariates in eye_rr – skipping multivariate step.\n")
} else {
  cat("Fitting categorical multivariate model for:", paste(sig_covs_e, collapse=", "), "\n")

  # Build factor indices and level counts for each sig covariate
  dl_mv_e <- data_list
  cov_factor_list_e <- list()
  n_lev_vec_e      <- integer()
  for(v in sig_covs_e) {
    x <- osa[[v]]; x[is.na(x)] <- "NR"
    fv <- as.factor(x)
    cov_factor_list_e[[v]] <- fv
    n_lev_vec_e[v] <- length(levels(fv))
    dl_mv_e[[paste0("cov_", v)]]   <- as.numeric(fv)
    dl_mv_e[[paste0("nlev_", v)]]  <- n_lev_vec_e[v]
  }

  body_mu1_e <- "    MU[i,1] <- threshold.sens[test[i],1]+study.re.sens[s[i]]+test.re.sens[s[i],test[i]]"
  body_mu2_e <- "    MU[i,2] <- threshold.spec[test[i],1]+study.re.spec[s[i]]+test.re.spec[s[i],test[i]]"
  for(v in sig_covs_e) {
    body_mu1_e <- paste0(body_mu1_e, "+bSe_", v, "[cov_", v, "[i]]")
    body_mu2_e <- paste0(body_mu2_e, "+bSp_", v, "[cov_", v, "[i]]")
  }
  prior_blocks_e <- ""
  for(v in sig_covs_e) {
    prior_blocks_e <- paste0(prior_blocks_e,
      "  bSe_",v,"[1]<-0; bSp_",v,"[1]<-0\n",
      "  for(k_ in 2:nlev_",v,"){bSe_",v,"[k_]~dnorm(0,0.05);bSp_",v,"[k_]~dnorm(0,0.05)}\n"
    )
  }
  params_mv_e <- c(paste0("bSe_", sig_covs_e), paste0("bSp_", sig_covs_e), "SDstudysens", "SDstudyspec")

  mv_cat_model_e <- paste0(
    "model {\n",
    "  for(i in 1:nObs) {\n",
    "    tp[i]~dbin(pi[i,1],pos[i]); tn[i]~dbin(pi[i,2],neg[i])\n",
    "    logit(pi[i,1])<-mu[i,1]; logit(pi[i,2])<-mu[i,2]\n",
    body_mu1_e, "\n",
    body_mu2_e, "\n",
    "    mu[i,1:2]~dmnorm(MU[i,],Omega[,])\n  }\n",
    prior_blocks_e,
    "  for(k in 1:ntest){threshold.sens[k,1]~dnorm(0,0.05);threshold.spec[k,1]~dnorm(0,0.05)}\n",
    "  for(k in 1:ns){study.re.sens[k]~dnorm(0,taustudysens);study.re.spec[k]~dnorm(0,taustudyspec)\n",
    "    for(l in 1:ntest){test.re.sens[k,l]~dnorm(0,tautestsens);test.re.spec[k,l]~dnorm(0,tautestspec)}}\n",
    "  taustudysens<-pow(SDstudysens,-2);SDstudysens~dunif(0,2)\n",
    "  taustudyspec<-pow(SDstudyspec,-2);SDstudyspec~dunif(0,2)\n",
    "  tautestsens<-pow(SDtestsens,-2);SDtestsens~dunif(0,2)\n",
    "  tautestspec<-pow(SDtestspec,-2);SDtestspec~dunif(0,2)\n",
    "  Omega[1:2,1:2]<-inverse(Sigma.sq[,])\n",
    "  for(m in 1:2){Sigma.sq[m,m]<-pow(sd[m],2)}\n",
    "  for(i in 1:2){for(j in (i+1):2){Sigma.sq[i,j]<-rho[i,j]*sd[i]*sd[j];Sigma.sq[j,i]<-Sigma.sq[i,j]}}\n",
    "  for(m in 1:2){sd[m]~dunif(0,1)}\n",
    "  for(i in 1:2){for(j in (i+1):2){g[j,i]<-0;a[i,j]~dunif(0,3.1415);rho[i,j]<-inprod(g[,i],g[,j])}}\n",
    "  g[1,1]<-1;g[1,2]<-cos(a[1,2]);g[2,2]<-sin(a[1,2])\n}\n"
  )

   fit_mv_e <- jags(data=dl_mv_e, inits=list(inits1, inits2, inits3, inits4),
                   parameters.to.save=params_mv_e,
                   n.chains=4, n.iter=n_iter, n.burnin=n_burnin,
                   model.file=textConnection(mv_cat_model_e))

  s_mv_e    <- fit_mv_e$BUGSoutput$summary
  sims_mv_e <- fit_mv_e$BUGSoutput$sims.list
  mv_res_e  <- data.frame()

  for(v in sig_covs_e) {
    lev <- levels(cov_factor_list_e[[v]])
    n_k <- n_lev_vec_e[v]
    for(k in 2:n_k) {
      se_d <- sims_mv_e[[paste0("bSe_",v)]][,k]
      sp_d <- sims_mv_e[[paste0("bSp_",v)]][,k]
      mv_res_e <- rbind(mv_res_e, data.frame(
        Covariate = v,
        Level     = paste0(lev[k], " vs ", lev[1]),
        Se_Mean   = mean(se_d), Se_Lo=quantile(se_d,0.025), Se_Hi=quantile(se_d,0.975),
        Se_Sig    = (quantile(se_d,0.025)>0 | quantile(se_d,0.975)<0),
        Se_Pval   = 2*min(mean(se_d>0), mean(se_d<0)),
        Sp_Mean   = mean(sp_d), Sp_Lo=quantile(sp_d,0.025), Sp_Hi=quantile(sp_d,0.975),
        Sp_Sig    = (quantile(sp_d,0.025)>0 | quantile(sp_d,0.975)<0),
        Sp_Pval   = 2*min(mean(sp_d>0), mean(sp_d<0))
      ))
    }
  }

  write.csv(mv_res_e, paste0(OUTDIR,"eye_rr_NMA_Multivariate_MR.csv"), row.names=FALSE)
  
  # --- Extract DIC, I2, and Tau for Eye-level Multivariate model ---
  dic_mv_e <- fit_mv_e$BUGSoutput$DIC
  tau_se_mv_e  <- sims_mv_e$SDstudysens;  tau_sp_mv_e  <- sims_mv_e$SDstudyspec
  tau2_se_mv_e <- tau_se_mv_e^2;          tau2_sp_mv_e <- tau_sp_mv_e^2
  I2_se_mv_e   <- tau2_se_mv_e / (tau2_se_mv_e + (pi^2/3))
  I2_sp_mv_e   <- tau2_sp_mv_e / (tau2_sp_mv_e + (pi^2/3))

  cat("\nEye-Level Multivariate Model Metrics:\n")
  cat(sprintf("DIC: %.2f\n", dic_mv_e))
  cat(sprintf("Tau (Se): %.3f [%.3f, %.3f]\n", mean(tau_se_mv_e), quantile(tau_se_mv_e, 0.025), quantile(tau_se_mv_e, 0.975)))
  cat(sprintf("Tau (Sp): %.3f [%.3f, %.3f]\n", mean(tau_sp_mv_e), quantile(tau_sp_mv_e, 0.025), quantile(tau_sp_mv_e, 0.975)))
  cat(sprintf("I2 (Se): %.1f%% [%.1f%%, %.1f%%]\n", mean(I2_se_mv_e)*100, quantile(I2_se_mv_e, 0.025)*100, quantile(I2_se_mv_e, 0.975)*100))
  cat(sprintf("I2 (Sp): %.1f%% [%.1f%%, %.1f%%]\n", mean(I2_sp_mv_e)*100, quantile(I2_sp_mv_e, 0.025)*100, quantile(I2_sp_mv_e, 0.975)*100))
  
  # Base Model metrics for comparison
  dic_base_e <- mod_eye$BUGSoutput$DIC
  
  comp_df_e <- data.frame(
    Model = c("Base NMA", "Multivariate MR"),
    DIC = c(dic_base_e, dic_mv_e),
    Tau_Se = c(sprintf("%.3f [%.3f, %.3f]", mean(tau_se_e), quantile(tau_se_e,0.025), quantile(tau_se_e,0.975)),
               sprintf("%.3f [%.3f, %.3f]", mean(tau_se_mv_e), quantile(tau_se_mv_e,0.025), quantile(tau_se_mv_e,0.975))),
    Tau_Sp = c(sprintf("%.3f [%.3f, %.3f]", mean(tau_sp_e), quantile(tau_sp_e,0.025), quantile(tau_sp_e,0.975)),
               sprintf("%.3f [%.3f, %.3f]", mean(tau_sp_mv_e), quantile(tau_sp_mv_e,0.025), quantile(tau_sp_mv_e,0.975))),
    I2_Se = c(sprintf("%.1f%% [%.1f%%, %.1f%%]", mean(I2_se_e)*100, quantile(I2_se_e,0.025)*100, quantile(I2_se_e,0.975)*100),
              sprintf("%.1f%% [%.1f%%, %.1f%%]", mean(I2_se_mv_e)*100, quantile(I2_se_mv_e,0.025)*100, quantile(I2_se_mv_e,0.975)*100)),
    I2_Sp = c(sprintf("%.1f%% [%.1f%%, %.1f%%]", mean(I2_sp_e)*100, quantile(I2_sp_e,0.025)*100, quantile(I2_sp_e,0.975)*100),
              sprintf("%.1f%% [%.1f%%, %.1f%%]", mean(I2_sp_mv_e)*100, quantile(I2_sp_mv_e,0.025)*100, quantile(I2_sp_mv_e,0.975)*100))
  )
  write.csv(comp_df_e, paste0(OUTDIR,"eye_rr_NMA_Model_Comparison.csv"), row.names=FALSE)
  cat("\nEye-level Model comparison saved: eye_rr_NMA_Model_Comparison.csv\n")

  
  if(nrow(mv_res_e) > 0) {
    label_mv_e <- paste0(mv_res_e$Covariate, ": ", mv_res_e$Level)
    n_mv_e <- nrow(mv_res_e); y_mv_e <- n_mv_e:1
    col_mv_se_e <- ifelse(mv_res_e$Se_Sig, "red", "#666666")
    col_mv_sp_e <- ifelse(mv_res_e$Sp_Sig, "red", "#666666")

    # Standalone eye-level Multivariate MR (via shared drawer; dense -> lab_cex 0.70)
    tiff(paste0(OUTDIR,"eye_rr_NMA_Multivariate_MR.tiff"), width=14, height=max(7, n_mv_e*0.55+2.5),
         units="in", res=300, compression="lzw", type="cairo")
    par(mfrow=c(1,2), mar=c(5, 11, 3.5, 2.5), family="Times", las=1)
    draw_mr_panel(mv_res_e$Se_Mean, mv_res_e$Se_Lo, mv_res_e$Se_Hi, col_mv_se_e,
                  label_mv_e, n_mv_e, y_mv_e, "Multivariate MR \u2013 Sensitivity", lab_cex=0.70)
    draw_mr_panel(mv_res_e$Sp_Mean, mv_res_e$Sp_Lo, mv_res_e$Sp_Hi, col_mv_sp_e,
                  label_mv_e, n_mv_e, y_mv_e, "Multivariate MR \u2013 Specificity", lab_cex=0.70)
    dev.off()
    cat("eye_rr Multivariate MR plot saved.\n")

    ## --- COMBINED Fig4 (single stacked device): (A) patient-level, (B) eye-level ---
    ## Proportional row heights so the 24-row eye panel is not squashed against
    ## the ~4-row patient panel. Cells 1,2 = patient Se/Sp; 3,4 = eye Se/Sp.
    if(exists("mv_res") && is.data.frame(mv_res) && nrow(mv_res) > 0) {
      pad <- 4
      h_top <- n_mv   + pad
      h_bot <- n_mv_e + pad
      tiff("AI_NMA_Multivariate_MR_combined.tiff", width=14,
           height=max(7, (h_top + h_bot) * 0.55 + 2.5),
           units="in", res=300, compression="lzw", type="cairo")
      layout(matrix(c(1,2,3,4), nrow=2, byrow=TRUE), heights=c(h_top, h_bot))
      par(mar=c(5, 11, 3.5, 4.5), family="Times", las=1)
      # (A) patient-level
      draw_mr_panel(mv_res$Se_Mean, mv_res$Se_Lo, mv_res$Se_Hi, col_mv_se,
                    label_mv, n_mv, y_mv, "Multivariate MR \u2013 Sensitivity", lab_cex=0.75)
      mtext("(A)", side=3, line=1.4, adj=0, font=2, cex=1.2, family="Times")
      draw_mr_panel(mv_res$Sp_Mean, mv_res$Sp_Lo, mv_res$Sp_Hi, col_mv_sp,
                    label_mv, n_mv, y_mv, "Multivariate MR \u2013 Specificity", lab_cex=0.75)
      # (B) eye-level
      draw_mr_panel(mv_res_e$Se_Mean, mv_res_e$Se_Lo, mv_res_e$Se_Hi, col_mv_se_e,
                    label_mv_e, n_mv_e, y_mv_e, "Multivariate MR \u2013 Sensitivity", lab_cex=0.70)
      mtext("(B)", side=3, line=1.4, adj=0, font=2, cex=1.2, family="Times")
      draw_mr_panel(mv_res_e$Sp_Mean, mv_res_e$Sp_Lo, mv_res_e$Sp_Hi, col_mv_sp_e,
                    label_mv_e, n_mv_e, y_mv_e, "Multivariate MR \u2013 Specificity", lab_cex=0.70)
      layout(1)
      dev.off()
      cat("Combined Multivariate MR (Fig4) saved: AI_NMA_Multivariate_MR_combined.tiff\n")
    } else {
      cat("Combined Fig4 skipped: patient-level mv_res not available.\n")
    }
  }
}


# --- Deeks Funnel ---
osa$tp_c <- osa$tp1+0.5; osa$fp_c <- osa$fp1+0.5
osa$fn_c <- osa$fn1+0.5; osa$tn_c <- osa$tn1+0.5
log_DOR_e <- log((osa$tp_c*osa$tn_c)/(osa$fp_c*osa$fn_c))
ess_e     <- 4*(osa$tp_c+osa$fn_c)*(osa$fp_c+osa$tn_c)/(osa$tp_c+osa$fn_c+osa$fp_c+osa$tn_c)
pt_cols_e <- arm_cols3[osa$test]
fit_dk_e  <- lm(log_DOR_e ~ I(1/sqrt(ess_e)), weights=ess_e)
pval_dk_e <- summary(fit_dk_e)$coefficients[2,4]
tiff(paste0(OUTDIR,"eye_rr_NMA_Deeks_Funnel.tiff"), width=8, height=8, units="in", res=300, compression="lzw", type="cairo")
par(mar=c(5,5,3,2), family="Times")
plot(1/sqrt(ess_e), log_DOR_e, pch=21, bg=pt_cols_e, col="white", cex=1.4, lwd=1.2,
     xlab="1 / \u221a(ESS)", ylab="Log Diagnostic Odds Ratio",
     main="Deeks' Funnel Plot \u2013 DL vs Human vs Human+DL",
     cex.axis=1.1, cex.lab=1.2, cex.main=1.2, las=1)
abline(fit_dk_e, col="#444444", lwd=2); abline(h=0,lty=3,col="grey60")
legend("topright", legend=c("DL only","Human only","Human + DL"),
       pch=21, pt.bg=c("#2166AC","#D6604D","#1A9641"), col="white", pt.cex=1.4, bty="n")
legend("topleft", bty="n", cex=1.05, legend=bquote("Deeks' test: "*italic(p)*"\u2009=\u2009"*.(format(round(pval_dk_e,3),nsmall=3))))
dev.off()
cat("eye_rr Deeks Funnel saved.\n")

while(!is.null(dev.list())) dev.off()

cat("\n=== ALL DONE (patient_rr + eye_rr) ===\n")
cat("patient_rr outputs: AI_NMA_*\n")
cat("eye_rr outputs:     eye_rr_NMA_*\n")

