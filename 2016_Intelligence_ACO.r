# title         Ant Colony Optimization (ACO)
# GitHub        https://github.com/ulrich-schroeders/syntax-publications
# date          2016-03-01
# version       1.0.0
# reference     Schroeders, U., Wilhelm, O., & Olaru, G. (2016). The influence of item sampling on sex differences in knowledge tests. Intelligence, 58(3), 22–32. https://doi.org/10.1016/j.intell.2016.06.003
# Additional information: This script is a revised and adopted version of the script provided by Leite (2015). [Leite, W. L. (2015). Ant Colony Optmization (ACO) Algorithm [Computer software]. Retrieved January 1, 2016, from http://education.ufl.edu/leite/code/]

# set workingDir, read data, etc.
library(lavaan)

# list with items; three fators with 20 items each...
list.items <- list(c(paste0("item", rep(1:20))), c(paste0("item", rep(21:40))), c(paste0("item", rep(41:60))))

# ... should be reduced to 10 items per factor
i.per.f <- c(10, 10, 10)
factors <- c("f1", "f2", "f3")

iter <- 40 # number of iterations
ants <- 60 # number of iterations (start with a iter/ants ratio of 2/3)
evaporation <- 0.90

# summary files
summaryfile <- paste0("ACO_sex_maxF_run_", sum(i.per.f), "_", iter, "_", ants, ".csv")
summaryfile2 <- paste0("ACO_sex_maxF_", sum(i.per.f), "_", iter, "_", ants, ".csv")

#FUNCTION START:
antcolony <- function(evaporation, list.items, i.per.f, iter, ants, factors, summaryfile, summaryfile2) {

  best.pheromone <- -2
  best.so.far.pheromone <- -2
 
  item.vector <- unlist(list.items, use.names = F)
  
  #creates the table of initial pheromone levels.
  include <- rep(2, length(unlist(list.items)))
  
  #puts initial best solution (all items selected).
  best.so.far.solution <- include
    
  #creates a list to store factors.
  selected.items <- list.items
  
  #starts counting the iterations
  count <- 1
  
  #starts counting continuous runs regardless of result.
  run <- 1
  
  #defines initial solutions.
  previous.solution <- include
  
  #starts loop through iterations.
  while (count <= iter) { 
    
    #sends a number of ants per time.
    ant  <- 1
    while (ant <= ants) {
      
      mod.3dim <- NULL
      #selects items for all factors.
      for (fac in 1:length(list.items)) {
        
        #selects the items for a short form for the factor
        positions <- is.element(item.vector, list.items[[fac]])
        prob <- include[positions]/ sum(include[positions])
        
        items <- sample(list.items[[fac]], size = i.per.f[fac], replace = F, prob)
        
        #stores selected items. 
        selected.items[[fac]] <- items
        
        # specifies CFA model
        mod.3dim <- paste(mod.3dim, "\n", paste0(factors[fac], " =~ ", paste(selected.items[[fac]], collapse=" + ")))

      } #finishes factor loop
      print(mod.3dim)
      
      #creates a vector of selected items.
      selected.items <- lapply(selected.items,sort)
      selected.vector <- unlist(selected.items, use.names = F)
      
      #creates a 0/1 vector of the same length of the full form indicating
      #whether an item was selected or not for the short form.
      select.indicator <- is.element(item.vector,selected.vector)
      notselect.indicator <- (select.indicator == FALSE)
      
      # estimates configural MI model
      fit.3dim.config <- cfa(model = mod.3dim, 
                      data = dat, 
                      ordered = unlist(list.items),
                      estimator = 'WLSMV',
                      group.label = c(1,2),
                      group = "sex",
                      std.lv = TRUE)

      # estimates strong MI model
      fit.3dim <- cfa(model = mod.3dim, 
                      data = dat, 
                      ordered = unlist(list.items),
                      estimator = 'WLSMV',
                      group.label = c(1,2),
                      group = "sex",
                      group.equal = c("loadings", "thresholds"),
                      std.lv = TRUE)
      
      # checks if convergence problems exist
      out <- capture.output(summary(fit.3dim))
      out2 <- capture.output(summary(fit.3dim.config))
      if(length(grep("did NOT converge", c(out, out2), ignore.case=T)) == 0) {
      
        # optimization #1: model fit & measurement invariance
        CFI <- fitMeasures(fit.3dim, "cfi.scaled")[[1]]
        RMSEA <- fitMeasures(fit.3dim, "rmsea.scaled")[[1]]
        Chi <- fitMeasures(fit.3dim, "chisq.scaled")[[1]]
      
        phi.CFI <- 1/(1+exp(95-100*CFI))
        phi.RMSEA <- 1-(1/(1+exp(5-100*RMSEA)))
        
        CFI.config <- fitMeasures(fit.3dim.config, "cfi.scaled")[[1]]
        phi.mi <- 1-(1/(1+exp(5-500*(abs(CFI.config - CFI)))))
        phi.fit <- (phi.CFI + phi.RMSEA + phi.mi)/3
      
        # optimization #2: reliability
        std <- standardizedSolution(fit.3dim)
      
        # group 1
        sum.lam <- sum(std[std$group==1 & std$op=="=~" & std$lhs==factors[1], 5])^2
        sum.the <- sum(1-std[std$group==1 & std$op=="=~" & std$lhs==factors[1], 5]^2)
        omega.1sci <- sum.lam/(sum.lam + sum.the)

        sum.lam <- sum(std[std$group==1 & std$op=="=~" & std$lhs==factors[2], 5])^2
        sum.the <- sum(1-std[std$group==1 & std$op=="=~" & std$lhs==factors[2], 5]^2)
        omega.1hum <- sum.lam/(sum.lam + sum.the)

        sum.lam <- sum(std[std$group==1 & std$op=="=~" & std$lhs==factors[3], 5])^2
        sum.the <- sum(1-std[std$group==1 & std$op=="=~" & std$lhs==factors[3], 5]^2)
        omega.1soc <- sum.lam/(sum.lam + sum.the)      
      
        # group 2
        sum.lam <- sum(std[std$group==2 & std$op=="=~" & std$lhs==factors[1], 5])^2
        sum.the <- sum(1-std[std$group==2 & std$op=="=~" & std$lhs==factors[1], 5]^2)
        omega.2sci <- sum.lam/(sum.lam + sum.the)

        sum.lam <- sum(std[std$group==2 & std$op=="=~" & std$lhs==factors[2], 5])^2
        sum.the <- sum(1-std[std$group==2 & std$op=="=~" & std$lhs==factors[2], 5]^2)
        omega.2hum <- sum.lam/(sum.lam + sum.the)

        sum.lam <- sum(std[std$group==2 & std$op=="=~" & std$lhs==factors[3], 5])^2
        sum.the <- sum(1-std[std$group==2 & std$op=="=~" & std$lhs==factors[3], 5]^2)
        omega.2soc <- sum.lam/(sum.lam + sum.the)      

        omega.sci <- (omega.1sci + omega.2sci)/2
        omega.hum <- (omega.1hum + omega.2hum)/2
        omega.soc <- (omega.1soc + omega.2soc)/2
        phi.rel <- 1/(1+exp(7-10*((omega.sci + omega.hum + omega.soc)/3)))
      
        # optimization #3: sex diff (negative values = girl favor; starting -.080)
        sci.diff <- std[std$group==2 & std$op=="~1" & std$lhs==factors[1], 5]
        hum.diff <- std[std$group==2 & std$op=="~1" & std$lhs==factors[2], 5]
        soc.diff <- std[std$group==2 & std$op=="~1" & std$lhs==factors[3], 5]
        phi.diff <- 1 + 3*((sci.diff + hum.diff + soc.diff)/3 + .080)

        # max f(x), optimize = maximize 
        pheromone <- phi.fit + phi.rel + phi.diff
      
        # saves information about the selected items and the model fit they generated.
        fit.info <- matrix(c(select.indicator, run, count, ant, Chi, CFI, RMSEA, phi.CFI, phi.RMSEA, phi.fit, CFI.config, phi.mi, omega.sci, omega.hum, omega.soc, phi.rel, sci.diff, hum.diff, soc.diff, phi.diff, pheromone, round(include,2)), 1, )    
    
        write.table(fit.info, file = summaryfile, append = T,
                  quote = F, sep = ";", row.names = F, col.names = F)          
      }  # end check convergence problems
      
      else {
        pheromone = -2.2
        write.table("convergence problems", file = summaryfile, append = T,
                    quote = F, sep = ";", row.names = F, col.names = F)
      }
      
      # adjusts count based on outcomes and selects best solution.
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
          
          best.CFI.config <- CFI.config
          best.phi.mi <- phi.mi          

          # updates best reliability
          best.omega.sci <- omega.sci
          best.omega.hum <- omega.hum
          best.omega.soc <- omega.soc
          best.phi.rel <- phi.rel 
          
          # sex differences
          best.sci.diff <- sci.diff
          best.hum.diff <- hum.diff
          best.soc.diff <- soc.diff
          best.phi.diff <- phi.diff
        } 
        
        #Move to next ant.
        ant <- ant + 1
        
      #ends loop through ants.
    }
    
    #adjusts pheromone only if the current pheromone is better than the previous.
    if (best.pheromone > best.so.far.pheromone) {
      
      #implements pheromone evaporation.
      include <- include * evaporation
      
      #Adjusts the pheromone levels.
      include.pheromone <- best.solution * best.pheromone * run * 0.2
      
      #updates pheromone.
      include <- include + include.pheromone
      
      #updates best so far solution and pheromone
      best.so.far.solution <- best.solution
      best.so.far.pheromone <- best.pheromone
      
      best.so.far.Chi <- best.Chi
      best.so.far.RMSEA <- best.RMSEA
      best.so.far.CFI <- best.CFI
      best.so.far.phi.CFI <- best.phi.CFI
      best.so.far.phi.RMSEA <- best.phi.RMSEA
      best.so.far.phi.fit <- best.phi.fit

      best.so.far.CFI.config <- best.CFI.config
      best.so.far.phi.mi <- best.phi.mi
      
      best.so.far.omega.sci <- best.omega.sci
      best.so.far.omega.hum <- best.omega.hum
      best.so.far.omega.soc <- best.omega.soc
      best.so.far.phi.rel <- best.phi.rel 
      
      best.so.far.sci.diff <- best.sci.diff
      best.so.far.hum.diff <- best.hum.diff
      best.so.far.soc.diff <- best.soc.diff
      best.so.far.phi.diff <- best.phi.diff
      
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
  "CFI.config", "phi.mi",
  "omega.sci", "omega.hum", "omega.soc", "phi.rel",
  "sci_diff","hum_diff","soc_diff", "phi.diff", "pheromone", item.vector), 1)
  write.table(title.final.solution, file = summaryfile2, append = T,
              quote = F, sep = ";", row.names = F, col.names = F)
  
  # Compile a matrix with the final solution.
  final.solution <- matrix(c(sum(i.per.f), best.so.far.Chi, best.so.far.CFI, best.so.far.RMSEA, best.so.far.phi.CFI, best.so.far.phi.RMSEA, best.so.far.phi.fit, 
  best.so.far.CFI.config, best.so.far.phi.mi,
  best.so.far.omega.sci, best.so.far.omega.hum, best.so.far.omega.soc,best.so.far.phi.rel, best.so.far.sci.diff, best.so.far.hum.diff, best.so.far.soc.diff, best.so.far.phi.diff, 
  best.so.far.pheromone, best.so.far.solution),1)
  
  write.table(final.solution, file = summaryfile2, append = T,
              quote = F, sep = ";", row.names = F, col.names = F)
}

# run the function
short <- antcolony(evaporation, list.items, i.per.f, iter, ants, factors, summaryfile, summaryfile2)
