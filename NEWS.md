# surv2sampleComp 1.0-5
### Minor improvements
* Added `conf.int` argument in `surv2sample()`, `rmstaug()`, and `rmstreg()` to specify the confidence coefficient of confidence intervals.
* Added `tau_start` argument in `rmst2sample()` to calculate the restricted mean survival time beyond a specific time point (to calculate AUC between `tau_start` and `tau`).
* The integrated restricted mean survival time (Zhao et al., 2012) is now calculated in `surv2sample()`.
* No default value for the truncation time point, `tau`, in `rmst2sample()`. `tau` needs to be specified by users.

# surv2sampleComp 1.0-4
* Fixed bug

# surv2sampleComp 1.0-3
* Added average of t-years and percentiles
* Fixed median survival time (quantile) 

# surv2sampleComp 1.0-2
* Added pbc sample data

# surv2sampleComp 1.0-1
* Initial release
