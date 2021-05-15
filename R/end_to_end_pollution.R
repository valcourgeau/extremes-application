library(gammaextremes)

max_length <- 20000
pollution_data <- read.csv('data/clean_pollution_data.csv')
pollution_data <- pollution_data[,-1]

cols_to_test <- 1:6
mode <- 'two-stage'
type <- 'exp'
bounds <- 'multiplier'
sample.length <- 5000
trials <- 200

# 2-stage

# GMM
method <- 'GMM'
depth <- 4

# # full length
gmm_2_stage_full <- lapply(cols_to_test,
  function(test_column){
    data <- pollution_data[1:max_length, test_column]
    return(SubSampleFit(
      data=data, depth=depth, trials=1,
      sample.length=max_length-1, method=method,
      mode=mode, type=type,
      bounds=bounds,
      parallel=F,
      seed=42))
  }
)
rlist::list.save(gmm_2_stage_full, 'results/gmm_2_stage_full.yaml')

# subs
gmm_2_stage_sub <- lapply(1:6,
       function(test_column){
         data <- pollution_data[1:max_length, test_column]
         return(SubSampleFit(
           data=data, depth=depth, trials=trials,
           sample.length=sample.length, method=method,
           mode=mode, type=type,
           bounds=bounds,
           parallel=T,
           seed=42))
       }
)
rlist::list.save(gmm_2_stage_sub, 'results/gmm_2_stage_sub.yaml')

# variance
gmm_2_stage_variance <- lapply(cols_to_test,
                           function(test_column){
                             data <- pollution_data[1:max_length, test_column]
                             i_guess <- gmm_2_stage_full[[test_column]]$estimators[1,4] #PairwiseLikelihood$InitGuess(data=data, depth=depth, n_trials=20)
                             i_guess_model <- CompositeMarginalMLE(data)
                             # init <- c(-0.009792636, 0.3141497, 19.96388, 0.220771)
                             ts_var <- gammaextremes::TrawlGMM$TwoStageVariance(
                               data=data, params=c(i_guess_model, i_guess),
                               depth=depth, type='exp', max_length=200)
                             return(ts_var)
                           }
)

rlist::list.save(gmm_2_stage_variance, 'results/gmm_2_stage_variance.yaml')

# PL
method <- 'PL'
depth <- 4

# full length
pl_2_stage_full <- lapply(cols_to_test,
                     function(test_column){
                       data <- pollution_data[1:max_length, test_column]
                       return(gammaextremes::SubSampleFit(
                         data=data, depth=depth, trials=1,
                         sample.length=max_length-1, method=method,
                         mode=mode, type=type,
                         bounds=bounds,
                         parallel=F,
                         seed=42))
                     }
)
rlist::list.save(pl_2_stage_full, 'results/pl_2_stage_full.yaml')

# subs
pl_2_stage_sub <- lapply(cols_to_test,
             function(test_column){
               data <- pollution_data[1:max_length, test_column]
               return(gammaextremes::SubSampleFit(
                 data=data, depth=depth, trials=trials,
                 sample.length=sample.length, method=method,
                 mode=mode, type=type,
                 bounds=bounds,
                 parallel=T,
                 seed=42))
             }
)
rlist::list.save(pl_2_stage_sub, 'results/pl_2_stage_sub.yaml')

# variance
pl_2_stage_variance <- lapply(cols_to_test,
                               function(test_column){
                                 data <- pollution_data[1:max_length, test_column]
                                 i_guess <- gmm_2_stage_full[[test_column]]$estimators[1,4] #PairwiseLikelihood$InitGuess(data=data, depth=depth, n_trials=20)
                                 i_guess_model <- CompositeMarginalMLE(data)
                                 # init <- c(-0.009792636, 0.3141497, 19.96388, 0.220771)
                                 ts_var <- gammaextremes::PairwiseLikelihood$TwoStageVariance(
                                   data=data, params=c(i_guess_model, i_guess),
                                   depth=depth, type='exp', max_length=200)
                                 return(ts_var)
                               }
)

rlist::list.save(pl_2_stage_variance, 'results/pl_2_stage_variance.yaml')

# full parameter
mode <- 'full'
# GMM
method <- 'GMM'
bounds <- 'multiplier'
depth <- 4

# full length
gmm_full_parameter_full <- lapply(cols_to_test,
                           function(test_column){
                             data <- pollution_data[1:max_length, test_column]
                             return(SubSampleFit(
                               data=data, depth=depth, trials=1,
                               sample.length=max_length-1, method=method,
                               mode=mode, type=type,
                               bounds=bounds,
                               parallel=F,
                               seed=42))
                           }
)
rlist::list.save(gmm_full_parameter_full, 'results/gmm_full_parameter_full.yaml')

# subs
gmm_full_parameter_sub <- lapply(cols_to_test,
                          function(test_column){
                            data <- pollution_data[1:max_length, test_column]
                            return(SubSampleFit(
                              data=data, depth=depth, trials=trials,
                              sample.length=sample.length, method=method,
                              mode=mode, type=type,
                              bounds=bounds,
                              parallel=T,
                              seed=42))
                          }
)
rlist::list.save(gmm_full_parameter_sub, 'results/gmm_full_parameter_sub.yaml')

# PL
method <- 'PL'
depth <- 4
bounds <- 'multiplier'
mode <- 'full'

# full length
pl_full_parameter_full <- lapply(cols_to_test,
                          function(test_column){
                            data <- pollution_data[1:max_length, test_column]
                            return(gammaextremes::SubSampleFit(
                              data=data, depth=depth, trials=1,
                              sample.length=max_length-1, method=method,
                              mode=mode, type=type,
                              bounds=bounds,
                              parallel=T,
                              seed=42))
                          }
)
rlist::list.save(pl_full_parameter_full, 'results/pl_full_parameter_full.yaml')

# subs
pl_full_parameter_sub <- lapply(cols_to_test,
                         function(test_column){
                           data <- pollution_data[1:max_length, test_column]
                           return(gammaextremes::SubSampleFit(
                             data=data, depth=depth, trials=trials,
                             sample.length=sample.length, method=method,
                             mode=mode, type=type,
                             bounds=bounds,
                             parallel=T,
                             seed=42))
                         }
)
rlist::list.save(pl_full_parameter_sub, 'results/pl_full_parameter_sub.yaml')