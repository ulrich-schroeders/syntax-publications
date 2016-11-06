# title         Bivariate Latent Change Score Model
# GitHub        https://github.com/ulrich-schroeders/syntax-publications
# date          2016-11-06
# version       1.0.0
# reference     Schroeders, U., Schipolowski, S., Zettler, I., Golle, J., & Wilhelm, O. (2016). Do the smart get smarter? Development of fluid and crystallized intelligence in 3rd grade. Intelligence. https://doi.org/10.1016/j.intell.2016.08.003

library(lavaan)

mod.lc <- '
# -----------
# gf
# -----------
# measurement part, equal factor loadings across time
  gf_pre =~ gfv1 + lam2*gfn1 + lam3*gff1
  gf_post =~ gfv2 + lam2*gfn2 + lam3*gff2

  gf_change =~ 0    # introduction of the latent difference variable

# (co)variances
  gf_post ~~ 0*gf_post  # no residual, perfect regression
  gf_pre ~~ NA*gf_pre   # freely estimate gf_pre variance
  gf_change ~~ NA*gf_change   # freely estimate change gf variance

# prediction part, difference factor
  gf_post ~ 1*gf_pre    # perfect regression of gf_pre on gf_post 
  gf_post ~ 1*gf_change    # perfect regression of gf_change on gf_post
  gf_change ~~ gf_pre       # gf_pre covaries with gf_change

# mean structure
  gf_pre ~ 0*1            # fixing the mean to zero
  gf_post ~ 0*1           # fixing the mean to zero
  gf_change ~ NA*1           # estimating the mean for the gf_change factor freely

# equal intercepts across time
  gfv1 ~ int1*1
  gfv2 ~ int1*1
  gfn1 ~ int2*1
  gfn2 ~ int2*1
  gff1 ~ int3*1
  gff2 ~ int3*1
      	
# residual variances across time
  gfv1 ~~ res1*gfv1
  gfn1 ~~ res2*gfn1
  gff1 ~~ res3*gff1
  gfv2 ~~ res1*gfv2
  gfn2 ~~ res2*gfn2
  gff2 ~~ res3*gff2
	
# residual covariances
  gfv1 ~~ gfv2
  gfn1 ~~ gfn2
  gff1 ~~ gff2 

# -----------
# gc
# -----------
# measurement part, equal factor loadings across time
  gc_pre =~ gc1Sci + lam5*gc1Hum + lam6*gc1Soc
  gc_post =~ gc2Sci + lam5*gc2Hum + lam6*gc2Soc

  gc_change =~ 0    # introduction of the latent difference variable

# (co)variances
  gc_post ~~ 0*gc_post  # no residual, perfect regression
  gc_pre ~~ NA*gc_pre   # freely estimate gc_pre variance
  gc_change ~~ NA*gc_change   # freely estimate change gc variance

# prediction part, difference factor
  gc_post ~ 1*gc_pre    # perfect regression of gc_pre on gc_post 
  gc_post ~ 1*gc_change    # perfect regression of gc_change on gc_post
  gc_change ~~ gc_pre       # gc_pre covaries with gc_change

# mean structure
  gc_pre ~ 0*1            # fixing the mean to zero
  gc_post ~ 0*1           # fixing the mean to zero
  gc_change ~ NA*1           # estimating the mean for the gc_change factor freely

# equal intercepts across time
  gc1Sci ~ int5*1
  gc2Sci ~ int5*1
  gc1Hum ~ int6*1
  gc2Hum ~ int6*1
  gc1Soc ~ int7*1
  gc2Soc ~ int7*1
      	
# residual variances across time
  gc1Sci ~~ res5*gc1Sci
  gc1Hum ~~ res6*gc1Hum
  gc1Soc ~~ res7*gc1Soc
  gc2Sci ~~ res5*gc2Sci
  gc2Hum ~~ res6*gc2Hum
  gc2Soc ~~ res7*gc2Soc
	
# residual covariances
  gc1Sci ~~ gc2Sci
  gc1Hum ~~ gc2Hum
  gc1Soc ~~ gc2Soc 
  
# ------
# gf gc paths
# ------
  gf_change ~~ gc_pre
  gc_change ~~ gf_pre
  gc_pre ~~ gf_pre
  gc_change ~~ gf_change
  gf_post ~~ 0*gc_post
'

fit.lc <- cfa(model=mod.lc,             # model specification
              data=dat,                 # data
			  estimator= "ML",          # estimator
			  missing="ML",             # FIML
			  orthogonal=TRUE)          # all latent variables are uncorrelated
summary(fit.lc, fit.measures=TRUE, standardized=TRUE)
