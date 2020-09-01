p_sim <- function(xa_ast, xb_ast, n, rab, disstest = c("RSDT", "BSDT", "UDT"), alternative = c("two.sided", "less", "greater")) {
  
  alternative <- match.arg(alternative)
  
  disstest <- match.arg(disstest)
  
  Sigma <- matrix(c(1, rab, rab, 1), ncol = 2)
  
  con <- MASS::mvrnorm(n + 1, mu = c(0, 0), Sigma = Sigma)
  
  xa_ast <- con[1, 1] + xa_ast
  xb_ast <- con[1, 2] + xb_ast
  
  
  pval <- switch(disstest,
                 RSDT = singcar::RSDT(xa_ast, xb_ast, con[ -1, 1], con[-1, 2], alternative = alternative)[["p.value"]],
                 BSDT = singcar::BSDT(xa_ast, xb_ast, con[ -1, 1], con[-1, 2], iter = 1000, alternative = alternative)[["p.value"]],
                 UDT  = singcar::UDT(xa_ast, xb_ast, con[ -1, 1], con[-1, 2], alternative = alternative)[["p.value"]]
  )
  
  pval
}


pwr_sim <- function(nsim, xa_ast, xb_ast, n, rab, disstest = c("RSDT", "BSDT", "UDT"), alternative = c("two.sided", "less", "greater"), alpha = 0.05) {
  
  alternative <- match.arg(alternative)
  
  disstest <- match.arg(disstest)
  
  pval <- vector(length = nsim)
  
  for(i in 1:nsim) {
    
    pval[i] <- p_sim(xa_ast = xa_ast, xb_ast = xb_ast, n = n, rab = rab, disstest = disstest, alternative = alternative)
    
  }
  
  
  sum(pval < alpha)/length(pval)
  
}


mc_data <- function(nsim, xa_ast, xb_ast, n, rab, disstest, pb, par, alternative, alpha = 0.05){
  
  mc_power_data <- data.frame(xa_ast = double(), xb_ast = double(),
                              ncon = double(), rab = double(),
                              testnum = integer(), alternative = character())
  
  
  
  
  
  # SETUP OF DATA STRUCTURE
  for (tst in 1:length(disstest)) {
    for (r in rab) {
      for (n_ in n) {
        for (xb in xb_ast) {
          for (xa in xa_ast) {
            for (alt in 1:length(alternative)) {
              for (a in alpha) {
                
                pd <- data.frame(xa_ast = xa, xb_ast = xb, ncon = n_, rab = r, testnum = tst, altnum = alt, alpha = a)
                
                mc_power_data <- rbind(mc_power_data, pd)
                
              }
            }
          }
        }
      }
    }
  }
  
  if (par == TRUE) {
    
    library(parallel)
    # FOR PARALLELL PROCESSING
    cl = makeCluster(12, type = "PSOCK")
    clusterExport(cl, varlist = c("pwr_sim", "p_sim"), envir=environment())
    on.exit(stopCluster(cl))
    
    
    if(pb == FALSE) {
      
      
      
      mc_power_data$power <- parApply(cl, mc_power_data, 1, function(x) pwr_sim(nsim,
                                                                                xa_ast = x[1],
                                                                                xb_ast = x[2],
                                                                                n = x[3],
                                                                                rab = x[4],
                                                                                disstest = disstest[x[5]],
                                                                                alternative = alternative[x[6]],
                                                                                alpha = x[7])) 
    }
    
    if(pb == TRUE) {
      
      library(pbapply)
      # FOR PARALLELL PROCESSING
      
      mc_power_data$power <- pbapply(mc_power_data, 1, function(x) pwr_sim(nsim,
                                                                           xa_ast = x[1],
                                                                           xb_ast = x[2],
                                                                           n = x[3],
                                                                           rab = x[4],
                                                                           disstest = disstest[x[5]],
                                                                           alternative = alternative[x[6]],
                                                                           alpha = x[7]),
                                     cl=cl) 
    }
    
  }
  

  if (par == FALSE) {
    
    mc_power_data$power <- apply(mc_power_data, 1, function(x) pwr_sim(nsim,
                                                                         xa_ast = x[1],
                                                                         xb_ast = x[2],
                                                                         n = x[3],
                                                                         rab = x[4],
                                                                         disstest = disstest[x[5]],
                                                                         alternative = alternative[x[6]],
                                                                         alpha = x[7])) 
    
  }
  
  
  
  mc_power_data$test <- ifelse(mc_power_data$testnum == 1, disstest[1], 
                               ifelse(mc_power_data$testnum == 2, disstest[2], 
                                      disstest[3]))
  mc_power_data$tail <- ifelse(mc_power_data$altnum == 1, alternative[1], 
                               ifelse(mc_power_data$altnum == 2, alternative[2], 
                                      alternative[3]))
  
  mc_power_data <- mc_power_data[ , -c(5, 6)]
  
  mc_power_data
}
pwr_RSDT_plot <- mc_data(nsim = 1000000, xa_ast= c(-1, -1.5, -2, -2.5, -3, -4), xb_ast = 0,
                        n = seq(4, 50, by = 2), rab = 0.5,
                        disstest = c("RSDT"), par = FALSE, alternative = c("two.sided"), alpha = c(0.05))


# save(pwr_RSDT_plot, file = "RSDT_plot.Rdata")

pwr_RSDT_tab <- mc_data(nsim = 1000000, xa_ast= seq(-4, -0.5, by = 0.5), xb_ast = 0,
                         n = seq(4, 24, by = 2), rab = 0.5,
                         disstest = c("RSDT"), par = FALSE, alternative = c("two.sided", "less"), alpha = c(0.05, 0.1))


# save(pwr_RSDT_tab, file = "RSDT_tab1M.Rdata")


pwr_RSDT_cor <- mc_data(nsim = 1000000, xa_ast= seq(-4, -1, by = 1), xb_ast = 0,
                        n = 100, rab = seq(0, 0.9, by = 0.05),
                        disstest = c("RSDT"), par = FALSE, alternative = c("two.sided"), alpha = c(0.05))


# save(pwr_RSDT_cor, file = "RSDT_plot_cor_1M.Rdata")



