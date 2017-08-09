# title         Ant Colony Optimization (ACO)
# GitHub        https://github.com/ulrich-schroeders/syntax-publications
# date          2016-03-01
# version       1.0.0
# reference     Schroeders, U., Wilhelm, O., & Olaru, G. (2016). Meta-heuristics in short scale construction: Ant Colony Optimization and Genetic Algorithm. PLOS ONE, 11, e0167110. http://doi.org/10.1371/journal.pone.0167110
# Additional information: This script is a revised and adopted version of the script provided by Leite (2007). [Leite, W. L. (2007). Ant Colony Optmization (ACO) Algorithm [Computer software]. Retrieved January 1, 2016, from http://education.ufl.edu/leite/code/]

# set workingDir, read data, etc.
library(lavaan)
library(psych)

dat <- read.table("03_Daten/ACO_ppvt_2016-03-21.dat", sep=" ", header=TRUE, na.strings=c(-99))
items.ppvt <- colnames(dat)[grep("vog900", colnames(dat))]

covariates <- c("read.speed", "math", "read", "mental.speed", 
            "reas", "intGer", "motGer", "age.2", "sex", "native.lang", 
            "grade.germ", "grade.math")
            
nitems <- 15 # number of items in the short scale
iter <- 60   # number of iterations
ants <- 150  # number of iterations (start with a iter/ants ratio of 2/3)
evaporation <- 0.90

# summary files
summaryfile <- paste0("ACO_ppvt_run_60150_", nitems, ".csv")
summaryfile2 <- paste0("ACO_ppvt_60150_", nitems, ".csv")

#FUNCTION START:
antcolony <- function(evaporation, items.ppvt, nitems, iter, ants, summaryfile, summaryfile2) {
                               
  best.pheromone <- 0
  best.so.far.pheromone <- 0
  item.vector <- items.ppvt
  
  #creates the table of initial pheromone levels.
  include <- rep(2, length(items.ppvt))
  
  #puts initial best solution (all items selected).
  best.so.far.solution <- include
    
  #creates a list to store factors.
  selected.items <- items.ppvt
  
  #starts counting the iterations
  count <- 1
  
  #starts counting continuous runs regardless of result.
  run <- 1
  
  #defines initial solutions.
  previous.solution <- include
  set.seed(789)
  
  #starts loop through iterations.
  while (count <= iter) { 
    
    #sends a number of ants per time.
    ant  <- 1
    while (ant <= ants) {
      
      mod.1dim <- NULL

      #selects the items for a short form for the factor
      positions <- is.element(item.vector, items.ppvt)
      prob <- include[positions]/ sum(include[positions])
      items <- sample(items.ppvt, size = nitems, replace = F, prob)
        
      #stores selected items 
      selected.items <- items
            
      # specifies CFA model
      mod.1dim <-  paste("lv =~", paste(selected.items, collapse = " + "))
     
      #creates a 0/1 vector of the same length of the long form indicating
      #whether an item was selected or not for the short form.
      select.indicator <- is.element(item.vector, selected.items)
      notselect.indicator <- (select.indicator == FALSE)
      
      # estimates CFA model
      fit.1dim <- cfa(model = mod.1dim, 
                      data = dat, 
                      ordered = items.ppvt,
                      estimator = 'WLSMV',
                      std.lv = TRUE)

      # save the complete lavaan ouput
      out <- capture.output(summary(fit.1dim, fit.measures=TRUE, standardized=TRUE))

      # reads the fit statistics (CFI, RMSEA, Chi^2)
      CFI <- fitMeasures(fit.1dim, "cfi.scaled")[[1]]
      RMSEA <- fitMeasures(fit.1dim, "rmsea.scaled")[[1]]
      Chi <- fitMeasures(fit.1dim, "chisq.scaled")[[1]]

      # optimization #1: model fit
      phi.CFI <- 1/(1+exp(95-100*CFI))
      phi.RMSEA <- 1-(1/(1+exp(5-100*RMSEA)))
      phi.fit <- (phi.CFI + phi.RMSEA)/2
        
      # optimization #2: factor saturation
      std <- standardizedSolution(fit.1dim)
      min.lam <- min(std[1:nitems, 4])
      mean.lam <- mean(std[1:nitems, 4])
      all.lam <- paste(std[1:nitems, 4], collapse="/")
      sum.lam <- sum(std[1:nitems, 4])^2
      sum.the <- sum(1-std[1:nitems, 4]^2)
      omega <- sum.lam/(sum.lam + sum.the)
      phi.lam <- 1/(1+exp(9-10*omega))
        
      # optimization #3: sensitivity
      mean.diff <- mean(describe(dat[, selected.items])[[3]])
      phi.sens <- -5*(mean.diff - .625)^2+1
        
      # optimization #4: differences in correlation to covariates
      dat$ppvt.short <- rowSums(dat[, selected.items])
      cor.covariates <- cor(dat[, c("ppvt", "ppvt.short", covariates)], use="p")
      cor.diff.cova.1 <- cor.covariates[1,3] - cor.covariates[2,3]
      cor.diff.cova.2 <- cor.covariates[1,4] - cor.covariates[2,4]
      cor.diff.cova.3 <- cor.covariates[1,5] - cor.covariates[2,5]
      cor.diff.cova.4 <- cor.covariates[1,6] - cor.covariates[2,6]
      cor.diff.cova.5 <- cor.covariates[1,7] - cor.covariates[2,7]
      cor.diff.cova.6 <- cor.covariates[1,8] - cor.covariates[2,8]
      cor.diff.cova.7 <- cor.covariates[1,9] - cor.covariates[2,9]
      cor.diff.cova.8 <- cor.covariates[1,10] - cor.covariates[2,10]
      cor.diff.cova.9 <- cor.covariates[1,11] - cor.covariates[2,11]
      cor.diff.cova.10 <- cor.covariates[1,12] - cor.covariates[2,12]
      cor.diff.cova.11 <- cor.covariates[1,13] - cor.covariates[2,13]
      cor.diff.cova.12 <- cor.covariates[1,14] - cor.covariates[2,14]

      max.cor.diff <- max(abs(c(cor.diff.cova.1, cor.diff.cova.2, cor.diff.cova.3,
                          cor.diff.cova.4, cor.diff.cova.5, cor.diff.cova.6,
                          cor.diff.cova.7, cor.diff.cova.8, cor.diff.cova.9,
                          cor.diff.cova.10, cor.diff.cova.11, cor.diff.cova.12)))
      phi.cor <- 1/-(1+exp(3-100*max.cor.diff))+1

      # max f(x), optimize = maximize
      pheromone <- phi.fit + phi.lam + phi.cor + phi.sens
        
      #saves information about the selected items and the model fit they generated.
      fit.info <- matrix(c(select.indicator, run, count, ant, 
        Chi, CFI, RMSEA, phi.CFI, phi.RMSEA, phi.fit, 
        min.lam, mean.lam, all.lam, phi.lam, 
        cor.diff.cova.1, cor.diff.cova.2, cor.diff.cova.3, 
        cor.diff.cova.4, cor.diff.cova.5, cor.diff.cova.6,
        cor.diff.cova.7, cor.diff.cova.8, cor.diff.cova.9,
        cor.diff.cova.10, cor.diff.cova.11, cor.diff.cova.12,
        max.cor.diff, phi.cor,
        mean.diff, phi.sens,
        pheromone, round(include,2)), 1, )
            
      write.table(fit.info, file = summaryfile, append = T,
                  quote = F, sep = ";", row.names = F, col.names = F)
        
      #adjusts count based on outcomes and selects best solution.
      if (pheromone >= best.pheromone) {
          
        # updates solution.
        best.solution <- select.indicator
        best.pheromone <- pheromone

        # updates best model fit
        best.Chi <- Chi
        best.RMSEA <- RMSEA
        best.CFI <- CFI
        best.phi.CFI <- phi.CFI
        best.phi.RMSEA <- phi.RMSEA
        best.phi.fit <- phi.fit

        # updates best factor loadings
        best.min.lam <- min.lam
        best.mean.lam <- mean.lam
        best.all.lam <- all.lam
        best.phi.lam <- phi.lam
              
        # updates best correlation differences
        best.cor.diff.cova.1 <- cor.diff.cova.1
        best.cor.diff.cova.2 <- cor.diff.cova.2
        best.cor.diff.cova.3 <- cor.diff.cova.3
        best.cor.diff.cova.4 <- cor.diff.cova.4
        best.cor.diff.cova.5 <- cor.diff.cova.5
        best.cor.diff.cova.6 <- cor.diff.cova.6
        best.cor.diff.cova.7 <- cor.diff.cova.7
        best.cor.diff.cova.8 <- cor.diff.cova.8
        best.cor.diff.cova.9 <- cor.diff.cova.9
        best.cor.diff.cova.10 <- cor.diff.cova.10
        best.cor.diff.cova.11 <- cor.diff.cova.11
        best.cor.diff.cova.12 <- cor.diff.cova.12
        best.max.cor.diff <- max.cor.diff
        best.phi.cor <- phi.cor
              
        # updates best distribution overlap
        best.mean.diff <- mean.diff
        best.phi.sens <- phi.sens
      } 
        
      #Move to next ant.
      ant <- ant + 1
        
      #ends loop through ants.
    }
    
    #adjusts pheromone only if the current pheromone is better than the previous.
    if (best.pheromone > best.so.far.pheromone) {
      
      #implements pheromone evaporation.
      include <- include * evaporation
      
      #adjusts the pheromone levels.
      include.pheromone <- best.solution * best.pheromone * run * 0.2
      
      #updates pheromone.
      include <- include + include.pheromone
      
      # updates best so far solution and pheromone.
      best.so.far.solution <- best.solution
      best.so.far.pheromone <- best.pheromone
      best.so.far.Chi <- best.Chi
      best.so.far.RMSEA <- best.RMSEA
      best.so.far.CFI <- best.CFI
      best.so.far.phi.CFI <- best.phi.CFI
      best.so.far.phi.RMSEA <- best.phi.RMSEA
      best.so.far.phi.fit <- best.phi.fit

      best.so.far.min.lam <- best.min.lam
      best.so.far.mean.lam <- best.mean.lam
      best.so.far.all.lam <- best.all.lam
      best.so.far.phi.lam <- best.phi.lam

      best.so.far.cor.diff.cova.1 <- best.cor.diff.cova.1
      best.so.far.cor.diff.cova.2 <- best.cor.diff.cova.2
      best.so.far.cor.diff.cova.3 <- best.cor.diff.cova.3
      best.so.far.cor.diff.cova.4 <- best.cor.diff.cova.4
      best.so.far.cor.diff.cova.5 <- best.cor.diff.cova.5
      best.so.far.cor.diff.cova.6 <- best.cor.diff.cova.6
      best.so.far.cor.diff.cova.7 <- best.cor.diff.cova.7
      best.so.far.cor.diff.cova.8 <- best.cor.diff.cova.8
      best.so.far.cor.diff.cova.9 <- best.cor.diff.cova.9
      best.so.far.cor.diff.cova.10 <- best.cor.diff.cova.10
      best.so.far.cor.diff.cova.11 <- best.cor.diff.cova.11
      best.so.far.cor.diff.cova.12 <- best.cor.diff.cova.12
      best.so.far.max.cor.diff <- best.max.cor.diff
      best.so.far.phi.cor <- best.phi.cor

      best.so.far.mean.diff <- best.mean.diff
      best.so.far.phi.sens <- best.phi.sens
        
      #re-starts count.
      count <- 1
      
      #end if clause for pheromone adjustment.
    } else {
      
      #advances count.
      count <- count + 1
    }
    
    #ends loop.
    run <- run + 1
  }
  
  title.final.solution <- matrix(c("Items", "Chi", "CFI", "RMSEA", "phi.CFI", "phi.RMSEA", "phi.fit", 
  "min.lam", "mean.lam", "all.lam", "phi.lam", 
  "cor.diff.cova.1", "cor.diff.cova.2", "cor.diff.cova.3", 
  "cor.diff.cova.4", "cor.diff.cova.5", "cor.diff.cova.6", 
  "cor.diff.cova.7", "cor.diff.cova.8", "cor.diff.cova.9", 
  "cor.diff.cova.10", "cor.diff.cova.11", "cor.diff.cova.12", 
  "max.cor.diff", "phi.cor",
  "mean.diff", "phi.sens",
  "pheromone", item.vector), 1)
  
  write.table(titel.final.solution, file = summaryfile2, append = T, 
              quote = F, sep = ";", row.names = F, col.names = F)
  
  # Compile a matrix with the final solution.
  final.solution <- matrix(c(sum(nitems), best.so.far.Chi, best.so.far.CFI, best.so.far.RMSEA, best.so.far.phi.CFI, best.so.far.phi.RMSEA, best.so.far.phi.fit, 
  best.so.far.min.lam, best.so.far.mean.lam, best.so.far.all.lam, best.so.far.phi.lam, 
  best.so.far.cor.diff.cova.1, best.so.far.cor.diff.cova.2, best.so.far.cor.diff.cova.3,
  best.so.far.cor.diff.cova.4, best.so.far.cor.diff.cova.5, best.so.far.cor.diff.cova.6,
  best.so.far.cor.diff.cova.7, best.so.far.cor.diff.cova.8, best.so.far.cor.diff.cova.9,
  best.so.far.cor.diff.cova.10, best.so.far.cor.diff.cova.11, best.so.far.cor.diff.cova.12,
  best.so.far.max.cor.diff, best.so.far.phi.cor,
  best.so.far.mean.diff, best.so.far.phi.sens,
  best.so.far.pheromone, best.so.far.solution), 1)
  
  write.table(final.solution, file = summaryfile2, append = T,
              quote = F, sep = ";", row.names = F, col.names = F)
  return(best.so.far.solution)
}

# run the function
short <- antcolony(evaporation, items.ppvt, nitems, iter, ants, summaryfile, summaryfile2) 


# title         Genetic Algorithm (GA)
# GitHub        https://github.com/ulrich-schroeders/syntax-publications
# date          2016-03-01
# version       1.0.0
# reference     Schroeders, U., Wilhelm, O., & Olaru, G. (2016). Meta-heuristics in short scale construction: Ant Colony Optimization and Genetic Algorithm. PLOS ONE, 11, e0167110. http://doi.org/10.1371/journal.pone.0167110
 
library(GAabbreviate)

# set workingDir, read data, etc.
library(lavaan)
library(psych)

dat <- read.table("03_Daten/ACO_ppvt_2016-03-21.dat", sep=" ", header=TRUE, na.strings=c(-99))
items.ppvt <- colnames(dat)[grep("vog900", colnames(dat))]

GAA.15 = GAabbreviate(dat[, items.ppvt], rowSums(dat[, items.ppvt]), itemCost = 0.01, 
  maxItems = 30, popSize = 500, maxiter = 300, run = 100, verbose = TRUE)
