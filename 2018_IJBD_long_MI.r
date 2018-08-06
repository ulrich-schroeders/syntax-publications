# Please cite as:
#   Edossa, A. K., Schroeders, U., Weinert, S., & Artelt, C. (2018). 
#   The development of emotional and behavioral self-regulation and their effects on academic 
#   in childhood. International Journal of Behavioral Development, 42(2), 192–202. 
#   https://doi.org/10.1177/0165025416687412
#
# contact    Ashenafi Kassahun Edossa <ashenafiedossa@gmail.com>
#            Ulrich Schroeders <ulrich.schroeders@gmail.com>
#
# title      Longitudinal Measurement Invariance Testing with Categorical Data
# GitHub     https://github.com/ulrich-schroeders/syntax-publications
# date       2015-08-04
# version    1.0.0


# Data availability statement: 
# The manuscript used data from the Millennium Cohort Study, which is a longitudinal birth cohort study that follows the lives of children born in 2000 and 2001 in the United Kingdom (more information: Hansen, K. (2014). Millennium cohort study first, second, third, fourth and fifth surveys: A guide to the datasets, (8th ed.). London: Centre for Longitudinal Studies, Institute of Education, University of London.).
# Information and data set: https://discover.ukdataservice.ac.uk/series/?sn=2000031

# 
#  Testing for Longitudinal Measurement Invariance with Continous and Categorical Data.
#  (more information in Edossa et al., 2018)
#  ---------------------------------------------------------------------------
#  Continuous                Factor                 Residual      Factor 
#  variables                 Loadings   Intercepts  Variances     Means
#  ---------------------------------------------------------------------------
#    Configural invariance   *          *           *             Fixed at 0
#    Weak invariance         Fixed      *           *             Fixed at 0
#    Strong invariance       Fixed      *           *             Fixed at 0
#    Strict invariance       Fixed      Fixed       Fixed         Fixed at 0/*
# 
#  ---------------------------------------------------------------------------
#  Categorical variables     Factor     Tresholds   Residual      Factor
#                            Loadings               Variances     Means
#  ---------------------------------------------------------------------------
#    Configural invariance   (*         *)          Fixed at 1    Fixed at 0
#    Strong invariance       (Fixed     Fixed)      Fixed at 1/*  Fixed at 0/*
#    Strict invariance       (Fixed     Fixed)      Fixed at 1    Fixed at 0/*
#  ---------------------------------------------------------------------------
#  Note. The asterisk (*) indicates that the parameter is freely estimated.
#  Fixed = the parameter is fixed to equity over time points; Fixed at 1 =
#  the residual variances are fixed to 1 at all time points; Fixed at 0 = 
#  factor means are fixed at 0 at all time points. Fixed at 0/* = factor
#  means are fixed at 0 at the first time point and freely estimated at the
#  other time points. Fixed at 1/* = the residual variances are fixed to 1 
#  at the first time point and freely estimated at the other time points.
#  Parameters in parentheses need to be varied in tandem.

# Explanation of variable labels
# ------------------------------
# er = emotional regulation
# br = behavioral regulation
# .x = age at which assessment took place

library(lavaan)
 
# configrual measurement invariance
-----------------------------------
mod.configural <- '
# measurement part
er.3 =~ er1.3 + er2.3 + er3.3 + er5.3
er.5 =~ er1.5 + er2.5 + er3.5 + er5.5
er.7 =~ er1.7 + er2.7 + er3.7 + er5.7

br.3 =~ isr4.3 + isr5.3 + task.3 + think.3
br.5 =~ isr4.5 + isr5.5 + task.5 + think.5
br.7 =~ isr4.7 + isr5.7 + task.7 + think.7

# mean structure
er.3 ~ 0*1
er.5 ~ 0*1
er.7 ~ 0*1
br.3 ~ 0*1
br.5 ~ 0*1
br.7 ~ 0*1

# residual covariances
er1.3 ~~ er1.5 + er1.7
er1.5 ~~ er1.7

er2.3 ~~ er2.5 + er2.7
er2.5 ~~ er2.7

er3.3 ~~ er3.5 + er3.7
er3.5 ~~ er3.7

er5.3 ~~ er5.5 + er5.7
er5.5 ~~ er5.7

isr4.3 ~~ isr4.5 + isr4.7
isr4.5 ~~ isr4.7

isr5.3 ~~ isr5.5 + isr5.7
isr5.5 ~~ isr5.7

think.3 ~~ think.5 + think.7
think.5 ~~ think.7

task.3 ~~ task.5 + task.7
task.5 ~~ task.7 '

fit.configural <- cfa(mod.configural, data=mcs, missing="pairwise", 
                      estimator="WLSMV", parameterization="theta",
                      ordered =c("isr4.3", "isr4.5", "isr4.7",
                                 "isr5.3", "isr5.5", "isr5.7",  
                                 "er1.3",  "er1.5",  "er1.7",  
                                 "er2.3",  "er2.5",  "er2.7", 
                                 "er3.3",  "er3.5",  "er3.7", 
                                 "er5.3",  "er5.5",  "er5.7",
                                 "task.3", "task.5", "task.7",  
                                 "think.3","think.5","think.7"))
summary(fit.configural, fit.measures=TRUE, standardized=TRUE)


# strong measurement invariance
-------------------------------
mod.strong <- '
# measurement part
er.3 =~ er1.3 + a2*er2.3 + a3*er3.3 + a4*er5.3 
er.5 =~ er1.5 + a2*er2.5 + a3*er3.5 + a4*er5.5 
er.7 =~ er1.7 + a2*er2.7 + a3*er3.7 + a4*er5.7 

br.3 =~ isr4.3 + b2*isr5.3 + b3*task.3 + b4*think.3 
br.5 =~ isr4.5 + b2*isr5.5 + b3*task.5 + b4*think.5 
br.7 =~ isr4.7 + b2*isr5.7 + b3*task.7 + b4*think.7 

# mean structure
er.3 ~ 0*1
er.5 ~ NA*1 
er.7 ~ NA*1 

br.3 ~ 0*1
br.5 ~ NA*1 
br.7 ~ NA*1 

# thresholds
er1.3 | c1*t1 + c2*t2; er1.5 | c1*t1 + c2*t2; er1.7 | c1*t1 + c2*t2
er2.3 | d1*t1 + d2*t2; er2.5 | d1*t1 + d2*t2; er2.7 | d1*t1 + d2*t2 
er3.3 | e1*t1 + e2*t2; er3.5 | e1*t1 + e2*t2; er3.7 | e1*t1 + e2*t2 
er5.3 | f1*t1 + f2*t2; er5.5 | f1*t1 + f2*t2; er5.7 | f1*t1 + f2*t2 

isr4.3  | g1*t1 + g2*t2; isr4.5  | g1*t1 + g2*t2; isr4.7  | g1*t1 + g2*t2 
isr5.3  | h1*t1 + h2*t2; isr5.5  | h1*t1 + h2*t2; isr5.7  | h1*t1 + h2*t2 
task.3  | i1*t1 + i2*t2; task.5  | i1*t1 + i2*t2; task.7  | i1*t1 + i2*t2 
think.3 | j1*t1 + j2*t2; think.5 | j1*t1 + j2*t2; think.7 | j1*t1 + j2*t2 

# residual covariances
er1.3 ~~ er1.5 + er1.7
er1.5 ~~ er1.7

er2.3 ~~ er2.5 + er2.7
er2.5 ~~ er2.7 

er3.3 ~~ er3.5 + er3.7
er3.5 ~~ er3.7 

er5.3 ~~ er5.5 + er5.7
er5.5 ~~ er5.7

isr4.3 ~~ isr4.5 + isr4.7
isr4.5 ~~ isr4.7

isr5.3 ~~ isr5.5 + isr5.7
isr5.5 ~~ isr5.7 

think.3 ~~ think.5 + think.7
think.5 ~~ think.7

task.3 ~~ task.5 + task.7
task.5 ~~ task.7 

# residual variances
er1.3 ~~ 1*er1.3 
er2.3 ~~ 1*er2.3 
er3.3 ~~ 1*er3.3 
er5.3 ~~ 1*er5.3 
er1.5 ~~ NA*er1.5 
er2.5 ~~ NA*er2.5 
er3.5 ~~ NA*er3.5 
er5.5 ~~ NA*er5.5 
er1.7 ~~ NA*er1.7 
er2.7 ~~ NA*er2.7 
er3.7 ~~ NA*er3.7 
er5.7 ~~ NA*er5.7 
isr4.3 ~~ 1*isr4.3 
isr5.3 ~~ 1*isr5.3 
task.3 ~~ 1*task.3 
think.3 ~~ 1*think.3
isr4.5 ~~ NA*isr4.5 
isr5.5 ~~ NA*isr5.5 
task.5 ~~ NA*task.5 
think.5 ~~ NA*think.5
isr4.7 ~~ NA*isr4.7 
isr5.7 ~~ NA*isr5.7 
task.7 ~~ NA*task.7 
think.7 ~~ NA*think.7  '

fit.strong<- cfa(mod.strong, data=mcs, missing="pairwise", 
                 estimator="WLSMV", parameterization="theta",
                 ordered =c("isr4.3", "isr4.5", "isr4.7",
                            "isr5.3", "isr5.5", "isr5.7",  
                            "er1.3",  "er1.5",  "er1.7",  
                            "er2.3",  "er2.5",  "er2.7", 
                            "er3.3",  "er3.5",  "er3.7", 
                            "er5.3",  "er5.5",  "er5.7",
                            "task.3", "task.5", "task.7",  
                            "think.3","think.5","think.7"))
summary(fit.strong, fit.measures=TRUE, standardized=TRUE)


# strict measurement invariance
-------------------------------
mod.strict <- '
# measurement part
er.3 =~ er1.3 + a2*er2.3 + a3*er3.3 + a4*er5.3 
er.5 =~ er1.5 + a2*er2.5 + a3*er3.5 + a4*er5.5 
er.7 =~ er1.7 + a2*er2.7 + a3*er3.7 + a4*er5.7 

br.3 =~ isr4.3 + b2*isr5.3 + b3*task.3 + b4*think.3 
br.5 =~ isr4.5 + b2*isr5.5 + b3*task.5 + b4*think.5 
br.7 =~ isr4.7 + b2*isr5.7 + b3*task.7 + b4*think.7 

# mean structure
er.3 ~ 0*1
er.5 ~ NA*1 
er.7 ~ NA*1 

br.3 ~ 0*1
br.5 ~ NA*1 
br.7 ~ NA*1 

# thresholds
er1.3 | c1*t1 + c2*t2; er1.5| c1*t1 + c2*t2; er1.7 | c1*t1 + c2*t2
er2.3 | d1*t1 + d2*t2; er2.5| d1*t1 + d2*t2; er2.7 | d1*t1 + d2*t2 
er3.3 | e1*t1 + e2*t2; er3.5| e1*t1 + e2*t2; er3.7 | e1*t1 + e2*t2 
er5.3 | f1*t1 + f2*t2; er5.5| f1*t1 + f2*t2; er5.7 | f1*t1 + f2*t2 

isr4.3  | g1*t1 + g2*t2; isr4.5  | g1*t1 + g2*t2; isr4.7  | g1*t1 + g2*t2 
isr5.3  | h1*t1 + h2*t2; isr5.5  | h1*t1 + h2*t2; isr5.7  | h1*t1 + h2*t2 
task.3  | i1*t1 + i2*t2; task.5  | i1*t1 + i2*t2; task.7  | i1*t1 + i2*t2 
think.3 | j1*t1 + j2*t2; think.5 | j1*t1 + j2*t2; think.7 | j1*t1 + j2*t2

# residual covariances
er1.3 ~~ er1.5 + er1.7
er1.5 ~~ er1.7

er2.3 ~~ er2.5 + er2.7
er2.5 ~~ er2.7 

er3.3 ~~ er3.5 + er3.7
er3.5 ~~ er3.7 

er5.3 ~~ er5.5 + er5.7
er5.5 ~~ er5.7

isr4.3 ~~ isr4.5 + isr4.7
isr4.5 ~~ isr4.7

isr5.3 ~~ isr5.5 + isr5.7
isr5.5 ~~ isr5.7 

think.3 ~~ think.5 + think.7
think.5 ~~ think.7

task.3 ~~ task.5 + task.7
task.5 ~~ task.7 '

fit.strict<- cfa(mod.strict, data=mcs, missing="pairwise", 
                 estimator="WLSMV", parameterization="theta",
                 ordered =c("isr4.3", "isr4.5", "isr4.7",
                            "isr5.3", "isr5.5", "isr5.7",  
                            "er1.3",  "er1.5",  "er1.7",  
                            "er2.3",  "er2.5",  "er2.7", 
                            "er3.3",  "er3.5",  "er3.7", 
                            "er5.3",  "er5.5",  "er5.7",
                            "task.3", "task.5", "task.7",  
                            "think.3","think.5","think.7"))
summary(fit.strict, fit.measures=TRUE, standardized=TRUE)

# testing the difference between the configural vs. strong MI model
-------------------------------------------------------------------
anova(fit.configural, fit.strong, method="satorra.bentler.2010") 
