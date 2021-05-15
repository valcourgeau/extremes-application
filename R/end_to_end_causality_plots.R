library(xvine)
library(dplyr)

max_length <- 20000
max_length_val <- 25000

traffic_files <- c(
  'data/MaryleboneRoadTrafficCount2000.xls',
  'data/MaryleboneRoadTrafficCount2001.xls',
  'data/MaryleboneRoadTrafficCount2002.xls'
)
traffic_data <- lapply(traffic_files, function(tf){
  mrt <- readxl::read_excel(tf, sheet = 2)
  mrt$Date <- stringr::str_remove(mrt$Date, "UTC")
  mrt$Date <- lubridate::parse_date_time(mrt$Date, "%Y-%m-%d %H:%M:%S", tz = 'UTC')
  mrt_tibble <- dplyr::tibble(mrt)
  traffic_tibble <- dplyr::group_by(mrt_tibble, Date) %>% 
    summarize(traffic = sum(CLS1 + CLS2 + CLS3 + CLS4 + CLS5 + CLS6, na.rm = TRUE))
})
traffic_data <- dplyr::bind_rows(traffic_data)$traffic
# traffic_data_tmp <- stats::stl(ts(traffic_data, frequency = 24), s.window = 24, t.window = 24*7)$time.series[, 3] # [1:max_length, 3]
traffic_data_val <- traffic_data[(max_length+1):max_length_val]
traffic_data <- traffic_data[1:max_length]

# pollution_data <- read.csv('data/non_gpd_clean_data.csv')

pollution_files <- c(
  'data/pollution_2000.csv',
  'data/pollution_2001.csv',
  'data/pollution_2002.csv'
)
pollution_data <- lapply(pollution_files, function(pf){
  mpf <- read.csv(pf)
  mpf <- mpf[,-c(1,2,4,6,8,10,12,14)]
  mpf[mpf == 'No data'] <- 0.0
  colnames(mpf) <- c('O3', 'NO', 'NO2', 'SO2', 'CO', 'PM10')
  mpf <- apply(mpf, MARGIN = 2, zoo::na.locf)
  mpf
})
pollution_data <- do.call(rbind, pollution_data)
pollution_data <-  apply(pollution_data, c(1,2), as.numeric)
pollution_data <- apply(pollution_data, MARGIN = 2, zoo::na.locf)[1:max_length,]
pollution_data <- cbind(pollution_data, traffic_data)
pollution_data <- apply(pollution_data, 2, function(x){x + rnorm(length(x), 0, 0.05*sd(x))})
save_core_pollution_data <- pollution_data
pollution_data <- t(t(pollution_data) / apply(pollution_data, MARGIN = 2, sd))

# d_values <- apply(pollution_data, MARGIN = 2, FUN = function(x){qq <- fracdiff::fdGPH(x); qq$data <- x; qq})
# all_data <- lapply(d_values, function(x){x$fd <- fracdiff::diffseries(x$data, d = x$d); x})
# p_data <- do.call(cbind, lapply(all_data, function(x){x$fd}))
# p_data[,7] <- pollution_data[,7]

old_pollution_data <- pollution_data
# pollution_data <- t(t(p_data) / apply(p_data, MARGIN = 2, sd))
colnames(pollution_data) <- c(colnames(pollution_data)[1:6], 'v/hr')



set.seed(42)
# thres_selection <- list()
# for(i in seq_len(ncol(pollution_data))){
#   if(i < 7){
#     qts <- quantile(pollution_data[,i], seq(from=0.8, to=.96, by=.02))
#   }else{
#     qts <- quantile(pollution_data[,i], seq(from=0.92, to=0.98, by=.01))
#   }
#   tryCatch({
#     thres_selection[[i]] <- eva::gpdSeqTests(
#       data = pollution_data[,i],
#       thresholds = qts,
#       method = 'ad',
#       nsim = 500,
#       allowParallel = T,
#       numCores = parallel::detectCores()
#     )
#   })
# }
# rlist::list.save(thres_selection, 'results/plots/causality/files/thres_selection.RData')
thres_selection <- rlist::list.load('results/plots/causality/files/thres_selection.RData')
u_all <- lapply(
  thres_selection,
  function(selec){
    sol_idx <- which(selec$ForwardStop > .05)
    if(length(sol_idx) == 0){
      print('no idx found.')
      sol <- tail(selec$threshold, 1)
    }else{
      sol <- selec$threshold[min(sol_idx)]
    }
    return(sol)
})
u_all <- unlist(u_all)

# Hyperparamters
u0_single <- .78
k_markov <- 7
conditional_on_col <- 7

# u0s_fixed <- rep(u0_single, ncol(pollution_data))
u0s_fixed <- vapply(seq_len(ncol(pollution_data)), function(i){mean(pollution_data[,i] < u_all[i])}, 1.0)
u0s_fixed
# u0s_fixed[7] <- .75
pollution_data_it <- xvine::apply_integral_transform(pollution_data, u0s_fixed)
#unfold
params_gpd_ests <- pollution_data_it$par.ests
params_gpd_ses <- pollution_data_it$par.ses
pollution_data_it <- pollution_data_it$data

# EMPIRICAL RESULTS
it_stack <- lapply(
  data.frame(t(build_stack(pollution_data_it, k = k_markov))), 
  function(x) matrix(x, byrow = T, ncol=ncol(pollution_data_it))
) # returns a list
it_stack <- abind::abind(it_stack, along = 3)

rit_exp_stack <- lapply(
  data.frame(t(build_stack(xvine::apply_reverse_exponential(pollution_data_it)$data, k = k_markov))), 
  function(x) matrix(x, byrow = T, ncol=ncol(pollution_data_it))
) # returns a list
rit_exp_stack <- abind::abind(rit_exp_stack, along = 3)

rit_stack <- lapply(
  data.frame(t(build_stack(pollution_data, k = k_markov))), 
  function(x) matrix(x, byrow = T, ncol=ncol(pollution_data_it))
) # returns a list
rit_stack <- abind::abind(rit_stack, along = 3)

# F/CF PARAMS
stack_pollution_data <- xvine::build_stack(pollution_data, k_markov)
f_cf_params <- lapply(
  seq_len(ncol(pollution_data)),
  function(target){
    it_data <- xvine::integral_transform(pollution_data[,target], u0=u0s_fixed[target])
    return(list('it'=it_data, 'u'=u_all[[target]], 'data'=pollution_data[,target], 'conditional_on_col'=conditional_on_col, 'target'=target))
  }
)
f_cf_params_matrix <- vapply(
  f_cf_params, 
  function(f_cf_target){
    prms <- as.vector(
      c(c(rbind(f_cf_target$it$par.ests, f_cf_target$it$par.ses), f_cf_target$u), f_cf_target$it$u0)
    )
  },
  FUN.VALUE = rep(0, 6)
)
f_cf_params_matrix <- t(f_cf_params_matrix)
f_cf_params_matrix
# f_cf_params_matrix <- cbind(f_cf_params_matrix, vapply(d_values, function(x){x$d}, 1.0))
# f_cf_params_matrix <- cbind(f_cf_params_matrix, vapply(d_values, function(x){x$sd.reg}, 1.0))
f_cf_params_matrix <- cbind(f_cf_params_matrix, apply(pollution_data, 2, function(x){unlist(tseries::adf.test(x, k = 12)$statistic)}))
# colnames(f_cf_params_matrix) <- c('gamma', 'gamma_sd', 'sigma', 'sigma_sd', 'u', 'u0', 'd', 'd_sd', 'ADF statistic')
colnames(f_cf_params_matrix) <- c('gamma', 'gamma_sd', 'sigma', 'sigma_sd', 'u', 'u0', 'ADF statistic')
rownames(f_cf_params_matrix) <- colnames(pollution_data) 
f_cf_params_matrix
write.csv(round(f_cf_params_matrix, 3), file = 'results/plots/causality/f_cf_params_matrix.csv')



# LOAD ABOVE FOR DATA PREP





pdf('results/plots/causality/files/marginals.pdf', height = 10, width = 11)
# png('results/plots/causality/files/marginals.png', height = 700, width = 600)
par(mfrow=c(7,2), mar=c(3.9,2.2,0,0))
plotting_colours <- RColorBrewer::brewer.pal(n = ncol(pollution_data), name = 'Dark2')
max_quant <- 99
# Exploration plots

layout(
  matrix(seq_len((1+ncol(pollution_data)) * 5), nrow=1+ncol(pollution_data), ncol=5, byrow = T),
  width=c(2,4,4,4,4,4), heights = c(1, rep(2, ncol(pollution_data)), mar=c(4,5,3,2))
) 
{
  plotting_colours <- RColorBrewer::brewer.pal(n = ncol(pollution_data), name = 'Dark2')
  plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
  plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
  title(main='A', cex.main=3, line = -2)
  plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
  title(main='B', cex.main=3, line = -2)
  plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
  title(main='C', cex.main=3, line = -2)
  plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
  title(main='D', cex.main=3, line = -2)
  
  for(i in seq_len(ncol(pollution_data))){
    plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
    title(colnames(pollution_data)[i], cex.main=1.6, line=-4)
    
    plot(old_pollution_data[1:3000,i], xlab='Index', bty= 'n',
         main='', ylab='', cex.main=1.8, cex.lab=1.5, cex.axis=1.4, col=plotting_colours[i], type='l')
    # plot(pollution_data[1:10000,i], xlab='Values', bty= 'n',
    #      main='', ylab='', cex.main=1.8, cex.lab=1.5, cex.axis=1.4, col=plotting_colours[i], type='l')
    hist(pmin(pollution_data[, i], 8), probability = T, xlab='Values', breaks = 25,
         main='', ylab='', cex.main=1.8, cex.lab=1.5, cex.axis=1.4, col=plotting_colours[i])
    hist_lims_exc <- c(f_cf_params[[i]]$u*.9, 5)
    hist(pmin(f_cf_params[[i]]$data[f_cf_params[[i]]$data > f_cf_params[[i]]$u], 8), probability = T, xlab='Values',
         breaks = 10,  main='', ylab='', cex.main=1.8, cex.lab=1.5, cex.axis=1.4, col=plotting_colours[i])
    th_quantiles <- evir::qgpd(seq_len(max_quant)/100, mu = f_cf_params[[i]]$u,
                               xi = f_cf_params[[i]]$it$par.ests[1], beta = f_cf_params[[i]]$it$par.ests[2])
    empi_quantiles <- quantile(
      f_cf_params[[i]]$data[f_cf_params[[i]]$data > f_cf_params[[i]]$u],
      seq_len(max_quant)/100)
    plot(th_quantiles, empi_quantiles, ylab='',
         xlab='GPD quantiles', bty='n', cex.main=1.5, cex.lab=1.3, cex=1.4, col=plotting_colours[i], pch=16, cex.axis=1.5)
    lines(c(min(th_quantiles), 5*max(empi_quantiles)), c(min(th_quantiles), 5*max(empi_quantiles)), cex=1.5, lwd=2)
  }
}
dev.off()
graphics.off()

pdf('results/plots/causality/files/acf.pdf', height = 10, width = 11)
png('results/plots/causality/files/acf.png', height = 700, width = 600)
par(mfrow=c(7,2), mar=c(3.9,2.2,0,0))
layout(
  matrix(seq_len((1+ncol(pollution_data)) * 5), nrow=1+ncol(pollution_data), ncol=5, byrow = T),
  width=c(2,4,4,4,4,4), heights = c(1, rep(2, ncol(pollution_data)), mar=c(4,5,3,2))
)
{
  plotting_colours <- RColorBrewer::brewer.pal(n = ncol(pollution_data), name = 'Dark2')
  plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
  plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
  title(main='A', cex.main=3, line = -2)
  plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
  title(main='B', cex.main=3, line = -2)
  plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
  title(main='C', cex.main=3, line = -2)
  plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
  title(main='D', cex.main=3, line = -2)
  
  acf_vals <- acf(pollution_data, plot=F)
  pacf_vals <- pacf(pollution_data, plot=F)
  
  for(i in seq_len(ncol(pollution_data))){
    plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
    title(colnames(pollution_data)[i], cex.main=1.6, line=-4)
    barplot(acf(pollution_data[,i], lag.max = 15, plot=F)$acf[,,1], names.arg=0:15,  bty='n', lwd=3, cex.lab=1.5, cex.axis=1.3, cex.names = 1.4,
            width = .3, space = .8, col=plotting_colours[i], xlab='', border=1)
    abline(h=0, col='black')
    abline(h=qnorm((1 + .975)/2)/sqrt(nrow(pollution_data)), lty=2, col='blue')
    abline(h=-qnorm((1 + .975)/2)/sqrt(nrow(pollution_data)), lty=2, col='blue')
    
    barplot(pacf(pollution_data[,i], lag.max = 15, plot=F)$acf[,,1], names.arg=1:15,  bty='n', lwd=3, cex.lab=1.5, cex.axis=1.3, cex.names = 1.4,
            width = .3, space = .8, col=plotting_colours[i], xlab='', border=1)
    abline(h=0, col='black')
    abline(h=qnorm((1 + .975)/2)/sqrt(nrow(pollution_data)), lty=2, col='blue')
    abline(h=-qnorm((1 + .975)/2)/sqrt(nrow(pollution_data)), lty=2, col='blue')
    
    barplot(acf_vals$acf[1:15, i, 7],names.arg=1:15,  bty='n', lwd=3, cex.lab=1.5, cex.axis=1.3, cex.names = 1.4,
            width = .3, space = .8, col=plotting_colours[i], xlab='', border=1)
    abline(h=0, col='black')
    abline(h=qnorm((1 + .975)/2)/sqrt(nrow(pollution_data)), lty=2, col='blue')
    abline(h=-qnorm((1 + .975)/2)/sqrt(nrow(pollution_data)), lty=2, col='blue')
    
    barplot(pacf_vals$acf[1:15, i, 7],names.arg=1:15,  bty='n', lwd=3, cex.lab=1.5, cex.axis=1.3, cex.names = 1.4,
            width = .3, space = .8, col=plotting_colours[i], xlab='', border=1)
    abline(h=0, col='black')
    abline(h=qnorm((1 + .975)/2)/sqrt(nrow(pollution_data)), lty=2, col='blue')
    abline(h=-qnorm((1 + .975)/2)/sqrt(nrow(pollution_data)), lty=2, col='blue')
  }
}
dev.off()

aic_markov_level <- 10
# svine_fit_aic <- xvine::fit_markov_svines(
#   data = pollution_data_it,
#   k.markov = aic_markov_level,
#   family_set="archimedean",
#   selcrit="aic",
#   tree_crit="tau",
#   threshold=0.05
# )
svine_fit_aic <- rlist::list.load('results/plots/causality/files/svine_fit_aic.RData')

plot(svine_fit_aic$copula, 1, var_names = 'use')

# above_idx_data <- which(it_stack[1,conditional_on_col,] > u0s_fixed[conditional_on_col])
# rlist::list.save(list(above_idx_data), 'results/plots/causality/files/above_idx_data.RData')
# below_idx_data <- sample(which(it_stack[1,conditional_on_col,] <= u0s_fixed[conditional_on_col]), size=length(above_idx_data), replace=FALSE)
# rlist::list.save(list(below_idx_data), 'results/plots/causality/files/below_idx_data.RData')
above_idx_data <- rlist::list.load('results/plots/causality/files/above_idx_data.RData')[[1]]
below_idx_data <- rlist::list.load('results/plots/causality/files/below_idx_data.RData')[[1]]
# vine-generated data
{
# split_idx <- seq(from=1, to=length(above_idx_data)+1, by=50)
# 
# above_cond_sims_svines_aic <- lapply(
#   seq_len(length(split_idx)-1),
#    function(j){lapply(
#       split_idx[j]:(split_idx[j+1]-1),
#       function(i){
#         cond_sims_svines <- svines::svinecop_sim(
#           n = k_markov-1, rep = 1, model = svine_fit_aic$copula,
#           past = rbind(matrix(runif((aic_markov_level-1)*ncol(pollution_data)), nrow=aic_markov_level-1, ncol=ncol(pollution_data)), it_stack[1,,above_idx_data[i]]),
#           cores = 1)
#         return(rbind(it_stack[1,,above_idx_data[i]], cond_sims_svines))
#       }
#   )
#   }
# )
# combined_above_cond_sims_svines_aic <- abind::abind(
#   lapply(above_cond_sims_svines_aic, function(splt)abind::abind(splt, along = 3)),
#   along = 3
# )
# 
# below_cond_sims_svines_aic <- lapply(
#   seq_len(length(split_idx)-1),
#  function(j){lapply(
#    split_idx[j]:(split_idx[j+1]-1),
#    function(i){
#      cond_sims_svines <- svines::svinecop_sim(
#        n = k_markov-1, rep = 1, model = svine_fit_aic$copula,
#        past = rbind(matrix(runif((aic_markov_level-1)*ncol(pollution_data)), nrow=aic_markov_level-1, ncol=ncol(pollution_data)), it_stack[1,,below_idx_data[i]]),
#        cores = 1)
#      return(rbind(it_stack[1,,below_idx_data[i]], cond_sims_svines))
#    }
#  )
#  }
# )
# combined_below_cond_sims_svines_aic <- abind::abind(
#   lapply(below_cond_sims_svines_aic, function(splt)abind::abind(splt, along = 3)),
#   along = 3
# )
# 
# combined_cond_sims_svines_aic <- abind::abind(
#   list(combined_above_cond_sims_svines_aic, combined_below_cond_sims_svines_aic), along=3
# )
# rlist::list.save(list(combined_cond_sims_svines_aic), 'results/plots/causality/files/cond_sims_svines_aic.RData')
}
combined_cond_sims_svines_aic <- rlist::list.load('results/plots/causality/files/cond_sims_svines_aic.RData')[[1]]
cond_sims_svines_aic_exp <- xvine::apply_reverse_exponential(combined_cond_sims_svines_aic)$data

# subsetting data
sub_rit_exp_stack <- rit_exp_stack[,, c(above_idx_data[1:1500], below_idx_data[1:1500])]
sub_it_stack <- it_stack[,, c(above_idx_data[1:1500], below_idx_data[1:1500])]

sub_it_stack[1,conditional_on_col,1:10]
combined_cond_sims_svines_aic[1,conditional_on_col,1:10]


corr_empirical <- apply(sub_it_stack, MARGIN = c(1,2), function(x){cor(x, sub_it_stack[1,conditional_on_col,])})
corr_sims <- apply(combined_cond_sims_svines_aic, MARGIN = c(1,2), function(x){cor(x, combined_cond_sims_svines_aic[1,conditional_on_col,])})
colnames(corr_empirical) <- colnames(pollution_data)
rownames(corr_empirical) <- c('t', vapply(1:6, function(i) paste('t +', i), '2'))
colnames(corr_sims) <- colnames(pollution_data)
rownames(corr_sims) <- c('t', vapply(1:6, function(i) paste('t +', i), '2'))
corr_error <- ((corr_sims - corr_empirical) )
colnames(corr_error) <- colnames(pollution_data)
rownames(corr_error) <- c('t', vapply(1:6, function(i) paste('t +', i), '2'))

pdf('results/plots/causality/files/correlation.pdf', height = 3.5, width = 11)
par(mfrow=c(1,3), mar=c(2,2.2,3,10))
corrplot::corrplot(
  corr_empirical, method = "color",
  cl.lim = c(min(corr_empirical), max(corr_empirical)), 
  number.cex = 1,
  addCoef.col = "black", # Add coefficient of correlation
  tl.cex=1.7, tl.col='black', tl.srt = 45, cl.cex=1.15
)
corrplot::corrplot(
  corr_sims,  method = "color",
  cl.lim = c(min(corr_sims), max(corr_sims)), 
  number.cex = 1,
  addCoef.col = "black", # Add coefficient of correlation
  tl.cex=1.7, tl.col='black', tl.srt = 45, cl.cex=1.15
)
corrplot::corrplot(
  corr_error,  method = "color",
  cl.lim = c(min(corr_error), max(corr_error)), 
  number.cex = 1,
  addCoef.col = "black", # Add coefficient of correlation
   tl.cex=1.7, tl.col='black', tl.srt = 45, cl.cex=1.15,
)
dev.off()
graphics.off()

w_t <- rep(.5, (ncol(pollution_data)) * (k_markov - 1))
w_t <- w_t / sum(w_t)
v <- qexp(.7)
q_source <- qexp(u0s_fixed[conditional_on_col])

{
svines_rit_exp_varying_v <- vapply(
  seq_len(100) / 100,
  function(v){
    causal_f <- xvine::wrapper_all(
      data = cond_sims_svines_aic_exp,
      col_source = conditional_on_col,
      u0_target = qexp(v), 
      u0_source = q_source
    ) 
    w2 <- rep(0, ncol(pollution_data)*(k_markov-1))
    idx_to_add <- (1+(conditional_on_col-1)*(k_markov-1)):(conditional_on_col*(k_markov-1))
    w2[idx_to_add] <- rep(0, k_markov-1)
    w2[setdiff(1:length(w2), idx_to_add)] <- w_t[setdiff(1:length(w2), idx_to_add)]
    c_p <- causal_f(w2)
    return(c(c_p$factual, c_p$counterfactual))
  }, c(1, 1)
)
  w_t <- rep(.5, (ncol(pollution_data)) * (k_markov - 1))
  w_t <- w_t / sum(w_t)
svines_rit_exp_pns_varying_v <- vapply(
  seq_len(100) / 100,
  function(v){
    poc <- xvine::proba_necessary_sufficient_causation
    wrap_pn_no_reg <-  xvine::wrapper_pn_all(
      data = cond_sims_svines_aic_exp,
      col_source = conditional_on_col, 
      u0_target = v,
      u0_source = q_source,
      poc = poc,
      lambda = NULL, p=1
    )
    
    return(wrap_pn_no_reg(w_t / sum(w_t)))
  }, c(1)
)
svines_rit_exp_pn_varying_v <- vapply(
  seq_len(100) / 100,
  function(v){
    poc <- xvine::proba_necessary_causation
    wrap_pn_no_reg <-  xvine::wrapper_pn_all(
      data = cond_sims_svines_aic_exp,
      col_source = conditional_on_col, 
      u0_target = v,
      u0_source = q_source,
      poc = poc,
      lambda = NULL, p=1
    )
    
    return(wrap_pn_no_reg(w_t / sum(w_t)))
  }, c(1)
)
svines_rit_exp_ps_varying_v <- vapply(
  seq_len(100) / 100,
  function(v){
    poc <- xvine::proba_sufficient_causation
    wrap_pn_no_reg <-  xvine::wrapper_pn_all(
      data = cond_sims_svines_aic_exp,
      col_source = conditional_on_col, 
      u0_target = v,
      u0_source = q_source,
      poc = poc,
      lambda = NULL, p=1
    )
    
    return(wrap_pn_no_reg(w_t / sum(w_t)))
  }, c(1)
)
}

# subset empirical
{
rit_exp_varying_v <- vapply(
  seq_len(100) / 100,
  function(v){
    causal_f <- xvine::wrapper_all(
      data = sub_rit_exp_stack,
      col_source = conditional_on_col,
      u0_target = qexp(v), 
      u0_source = q_source
    ) 
    w2 <- rep(0, ncol(pollution_data)*(k_markov-1))
    idx_to_add <- (1+(conditional_on_col-1)*(k_markov-1)):(conditional_on_col*(k_markov-1))
    w2[idx_to_add] <- rep(0, k_markov-1)
    w2[setdiff(1:length(w2), idx_to_add)] <- w_t[setdiff(1:length(w2), idx_to_add)]
    c_p <- causal_f(w2)
    return(c(c_p$factual, c_p$counterfactual))
  }, c(1, 1)
)
rit_exp_pns_varying_v <- vapply(
  seq_len(100) / 100,
  function(v){
    poc <- xvine::proba_necessary_sufficient_causation
    wrap_pn_no_reg <-  xvine::wrapper_pn_all(
      data = sub_rit_exp_stack,
      col_source = conditional_on_col, 
      u0_target = v,
      u0_source = q_source,
      poc = poc,
      lambda = NULL, p=1
    )
    
    return(wrap_pn_no_reg(w_t / sum(w_t)))
  }, c(1)
)
rit_exp_pn_varying_v <- vapply(
  seq_len(100) / 100,
  function(v){
    poc <- xvine::proba_necessary_causation
    wrap_pn_no_reg <-  xvine::wrapper_pn_all(
      data = sub_rit_exp_stack,
      col_source = conditional_on_col, 
      u0_target = v,
      u0_source = q_source,
      poc = poc,
      lambda = NULL, p=1
    )
    
    return(wrap_pn_no_reg(w_t / sum(w_t)))
  }, c(1)
)
rit_exp_ps_varying_v <- vapply(
  seq_len(100) / 100,
  function(v){
    poc <- xvine::proba_sufficient_causation
    wrap_pn_no_reg <-  xvine::wrapper_pn_all(
      data = sub_rit_exp_stack,
      col_source = conditional_on_col, 
      u0_target = v,
      u0_source = q_source,
      poc = poc,
      lambda = NULL, p=1
    )
    
    return(wrap_pn_no_reg(w_t / sum(w_t)))
  }, c(1)
)
}

# without source
{
  w_t <- rep(.5, (ncol(pollution_data) -1) * (k_markov - 1))
  w_t <- w_t / sum(w_t)
  rit_exp_varying_v_wo_source <- vapply(
    seq_len(100) / 100,
    function(v){
      causal_f <- xvine::wrapper_all(
        data = sub_rit_exp_stack,
        col_source = conditional_on_col,
        u0_target = qexp(v), 
        u0_source = q_source
      ) 
      w2 <- rep(0, ncol(pollution_data)*(k_markov-1))
      idx_to_add <- (1+(conditional_on_col-1)*(k_markov-1)):(conditional_on_col*(k_markov-1))
      w2[idx_to_add] <- rep(0, k_markov-1)
      w2[setdiff(1:length(w2), idx_to_add)] <- w_t[setdiff(1:length(w2), idx_to_add)]
      c_p <- causal_f(w2)
      return(c(c_p$factual, c_p$counterfactual))
    }, c(1, 1)
  )
  rit_exp_pns_varying_v_wo_source <- vapply(
    seq_len(100) / 100,
    function(v){
      poc <- xvine::proba_necessary_sufficient_causation
      wrap_pn_no_reg <-  xvine::wrapper_pn_all_but_source(
        data = sub_rit_exp_stack,
        col_source = conditional_on_col, 
        u0_target = v,
        u0_source = q_source,
        poc = poc,
        lambda = NULL, p=1
      )
      
      return(wrap_pn_no_reg(w_t / sum(w_t)))
    }, c(1)
  )
  rit_exp_pn_varying_v_wo_source <- vapply(
    seq_len(100) / 100,
    function(v){
      poc <- xvine::proba_necessary_causation
      wrap_pn_no_reg <-  xvine::wrapper_pn_all_but_source(
        data = sub_rit_exp_stack,
        col_source = conditional_on_col, 
        u0_target = v,
        u0_source = q_source,
        poc = poc,
        lambda = NULL, p=1
      )
      
      return(wrap_pn_no_reg(w_t / sum(w_t)))
    }, c(1)
  )
  rit_exp_ps_varying_v_wo_source <- vapply(
    seq_len(100) / 100,
    function(v){
      poc <- xvine::proba_sufficient_causation
      wrap_pn_no_reg <-  xvine::wrapper_pn_all_but_source(
        data = sub_rit_exp_stack,
        col_source = conditional_on_col, 
        u0_target = v,
        u0_source = q_source,
        poc = poc,
        lambda = NULL, p=1
      )
      
      return(wrap_pn_no_reg(w_t / sum(w_t)))
    }, c(1)
  )
}

par(mfrow=c(1, 3))
plotting_colours <- RColorBrewer::brewer.pal(n = 6, name = 'Dark2')
lwd_val <- 5
cex_lab <- 2.8
cex_axis <- 2
pdf('results/plots/causality/files/p_f_cf.pdf', height = 5, width = 11)
{
par(mfrow=c(2,3), mar=c(4.5,5,.5,.5))
matplot(qexp(seq_len(100)/100), t(rit_exp_varying_v), lty=1:2, lwd=lwd_val, type="l", pch=15:16, cex=1.5, xlab='v',
        ylab='Probability', bty='n', cex.axis=cex_axis, cex.lab=cex_lab, col=plotting_colours)
legend("topright", legend=c(expression(paste( p[f])), expression(paste( p[cf]))),
       lty=1:3, lwd=lwd_val,  col=plotting_colours, cex=2, bty = "n", pt.cex = 2)

matplot(qexp(seq_len(25)/26), rit_exp_varying_v[1,seq(from=1, to=100, by=4)]-rit_exp_varying_v[2,seq(from=1, to=100, by=4)], lty=1, lwd=lwd_val, type="l", col=plotting_colours[3],
        cex=1.5,pch=15:16, xlab='v', ylab=expression(paste( p[f] - p[cf])), bty='n', cex.axis=cex_axis, cex.lab=cex_lab, ylim=c(0, .25))

rit_exp_poc <- rbind(exp(rit_exp_pns_varying_v), exp(rit_exp_pn_varying_v), exp(rit_exp_ps_varying_v))
matplot(qexp(seq_len(25)/26), t(rit_exp_poc[,seq(from=1, to=100, by=4)]), type="l", lwd=lwd_val, lty=1:3, pch=15:17, xlab='v', ylab='PC',
        bty='n', cex.axis=cex_axis, cex.lab=cex_lab, col=plotting_colours[4:6])
legend("topright", legend=c('PNS', 'PN', 'PS'), lty=1:3, lwd=lwd_val, col=plotting_colours[4:6], cex=1.5, bty = "n")

# without source
# matplot(qexp(seq_len(100)/100), t(rit_exp_varying_v_wo_source), lty=1, lwd=lwd_val, type="l", pch=15:16, cex=1.5, xlab='v',
#         ylab='Probability', bty='n', cex.axis=cex_axis, cex.lab=cex_lab, col=plotting_colours)
# legend("topright", legend=c(expression(paste( p[f])), expression(paste( p[cf]))), lty=1, lwd=lwd_val,
#        col=plotting_colours, cex=2, bty = "n", pt.cex = 2)
# 
# matplot(qexp(seq_len(100)/100), rit_exp_varying_v_wo_source[1,]-rit_exp_varying_v_wo_source[2,], lty=1, lwd=lwd_val, type="l", col=plotting_colours[3],
#         cex=1.5,pch=15:16, xlab='v', ylab=expression(paste( p[f] - p[cf])), bty='n', cex.axis=cex_axis, cex.lab=cex_lab)
# 
# svines_rit_exp_poc <- rbind(exp(rit_exp_pns_varying_v_wo_source), exp(rit_exp_pn_varying_v_wo_source), exp(rit_exp_ps_varying_v_wo_source))
# matplot(qexp(seq_len(25)/26), t(svines_rit_exp_poc[,seq(from=1, to=100, by=4)]), type="l", lwd=lwd_val, lty=1, pch=15:17, xlab='v', ylab='PC',
#         bty='n', cex.axis=cex_axis, cex.lab=cex_lab, col=plotting_colours[4:6])
# legend("topright", legend=c('PNS', 'PN', 'PS'), lty=1, lwd=lwd_val, col=plotting_colours[4:6], cex=1.5, bty = "n")

# svine data
matplot(qexp(seq_len(100)/100), t(svines_rit_exp_varying_v), lty=1:2, lwd=lwd_val, type="l", pch=15:16, cex=1.5, xlab='v',
        ylab='Probability', bty='n', cex.axis=cex_axis, cex.lab=cex_lab, col=plotting_colours)
legend("topright", legend=c(expression(paste( p[f])), expression(paste( p[cf]))), lty=1:2, lwd=lwd_val, 
       col=plotting_colours, cex=2, bty = "n", pt.cex = 2)

matplot(qexp(seq_len(25)/26), svines_rit_exp_varying_v[1,seq(from=1, to=100, by=4)]-svines_rit_exp_varying_v[2,seq(from=1, to=100, by=4)],
        lty=1, lwd=lwd_val, type="l", col=plotting_colours[3],
        cex=1.5,pch=15:16, xlab='v', ylab=expression(paste( p[f] - p[cf])), bty='n', cex.axis=cex_axis, cex.lab=cex_lab)

svines_rit_exp_poc <- rbind(exp(svines_rit_exp_pns_varying_v), exp(svines_rit_exp_pn_varying_v), exp(svines_rit_exp_ps_varying_v))
matplot(qexp(seq_len(25)/26), t(svines_rit_exp_poc[,seq(from=1, to=100, by=4)]), type="l", lwd=lwd_val, pch=15:17, xlab='v', ylab='PC',
        bty='n', cex.axis=cex_axis, cex.lab=cex_lab, col=plotting_colours[4:6], lty=1:3)
legend("topright", legend=c('PNS', 'PN', 'PS'), lty=1:3, lwd=lwd_val, col=plotting_colours[4:6], cex=1.5, bty = "n")
}
dev.off()
graphics.off()



######### COMPUTATION ###########
q_target <- qexp(.8)

lmbds <- list(NULL, .001, .01, .1, .5, 1, 5, 10, 100, 200)
# PNS
{
  # start values
  penalty_weights <- 0
  poc <- xvine::proba_necessary_sufficient_causation
  wrap_pn_no_reg <-  xvine::wrapper_pn_all_but_source(
    data = cond_sims_svines_aic_exp,
    col_source = conditional_on_col, 
    u0_target = q_target,
    u0_source = q_source,
    poc = poc,
    lambda = NULL, p=1
  )
  w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
  w_t <- w_t / sum(w_t)
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c('wrap_pn_no_reg', 'wrap_pn', 'cond_sims_svines_aic_exp', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
  system.time({
    start_w_deoptim <- DEoptim::DEoptim(
      lower = rep(0, length(w_t)),
      upper = rep(1, length(w_t)),
      fn = function(w){-wrap_pn_no_reg(w) + (sum(w)-1)^2*penalty_weights},
      control = DEoptim::DEoptim.control(
        trace = T, parallelType=1, strategy=6, cl=cl, reltol = 1e-1, NP=20*length(w_t), p=.5, steptol=50, itermax=20)
    )
    print(start_w_deoptim$optim$bestmem)
  })
  tryCatch(parallel::stopCluster(cl)) 
  wrap_pn_no_reg(start_w_deoptim$optim$bestmem)
  print(sum(start_w_deoptim$optim$bestmem))
  
  #p1 p2
  penalty_weights <- 0
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  aws_pns_p1 <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all_but_source(
        data = cond_sims_svines_aic_exp,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=1
      )
      fn_opt <- function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights}
      parallel::clusterExport(cl, c('wrap_pn', 'cond_sims_svines_aic_exp', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
      system.time({
        all_p1_hot_start <- optimParallel::optimParallel(
          par=start_w_deoptim$optim$bestmem, fn = fn_opt,
          lower = rep(0, length(w_t)), upper = rep(1, length(w_t)),
          control = list(trace=1), parallel = list(cl=cl, forward=F, loginfo=F)  
        )
      })
      all_p1_hot_start$value_no_reg <- wrap_pn_no_reg(all_p1_hot_start$par)
      return(all_p1_hot_start)
    }
  )
  rlist::list.save(x = aws_pns_p1, file = 'results/plots/causality/poc/aws_pns_p1.RData')
  
  aws_pns_p2 <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all_but_source(
        data = cond_sims_svines_aic_exp,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=2
      )
      w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
      parallel::clusterExport(cl, c('wrap_pn', 'cond_sims_svines_aic_exp', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
      system.time({
        all_p2_hot_start <- optimParallel::optimParallel(
          par=start_w_deoptim$optim$bestmem,
          fn = function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights},
          lower = rep(0, length(w_t)),
          upper = rep(1, length(w_t)),
          control = list(trace=0),
          parallel = list(cl=cl, forward=F, loginfo=F)  
        )
      })
      all_p2_hot_start$value_no_reg <- wrap_pn_no_reg(all_p2_hot_start$par)
      return(all_p2_hot_start)
    }
  )
  rlist::list.save(x = aws_pns_p2, file = 'results/plots/causality/poc/aws_pns_p2.RData')
  tryCatch(parallel::stopCluster(cl)) 
}

# PNS all
{
  # start values
  penalty_weights <- 0
  poc <- xvine::proba_necessary_sufficient_causation
  wrap_pn_no_reg <-  xvine::wrapper_pn_all(
    data = cond_sims_svines_aic_exp,
    col_source = conditional_on_col, 
    u0_target = q_target,
    u0_source = q_source,
    poc = poc,
    lambda = NULL, p=1
  )
  w_t <- rep(.5, (ncol(pollution_data)) * (k_markov - 1))
  w_t <- w_t / sum(w_t)
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c('wrap_pn_no_reg', 'wrap_pn', 'cond_sims_svines_aic_exp', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
  system.time({
    start_w_deoptim_a_pns <- DEoptim::DEoptim(
      lower = rep(0, length(w_t)),
      upper = rep(1, length(w_t)),
      fn = function(w){-wrap_pn_no_reg(w) + (sum(w)-1)^2*penalty_weights},
      control = DEoptim::DEoptim.control(
        trace = T, parallelType=1, strategy=6, cl=cl, reltol = 1e-1, NP=20*length(w_t), p=.5, steptol=50, itermax=20)
    )
    print(start_w_deoptim_a_pns$optim$bestmem)
  })
  tryCatch(parallel::stopCluster(cl)) 
  wrap_pn_no_reg(start_w_deoptim_a_pns$optim$bestmem)
  print(sum(start_w_deoptim_a_pns$optim$bestmem))
  
  #p1 p2
  penalty_weights <- 0
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  a_pns_p1 <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all(
        data = cond_sims_svines_aic_exp,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=1
      )
      fn_opt <- function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights}
      parallel::clusterExport(cl, c('wrap_pn', 'cond_sims_svines_aic_exp', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
      system.time({
        all_p1_hot_start <- optimParallel::optimParallel(
          par=start_w_deoptim_a_pns$optim$bestmem, fn = fn_opt,
          lower = rep(0, length(w_t)), upper = rep(1, length(w_t)),
          control = list(trace=1), parallel = list(cl=cl, forward=F, loginfo=F)  
        )
      })
      all_p1_hot_start$value_no_reg <- wrap_pn_no_reg(all_p1_hot_start$par)
      return(all_p1_hot_start)
    }
  )
  rlist::list.save(x = a_pns_p1, file = 'results/plots/causality/poc/a_pns_p1.RData')
  
  a_pns_p2 <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all(
        data = cond_sims_svines_aic_exp,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=2
      )
      w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
      parallel::clusterExport(cl, c('wrap_pn', 'cond_sims_svines_aic_exp', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
      system.time({
        all_p2_hot_start <- optimParallel::optimParallel(
          par=start_w_deoptim_a_pns$optim$bestmem,
          fn = function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights},
          lower = rep(0, length(w_t)),
          upper = rep(1, length(w_t)),
          control = list(trace=0),
          parallel = list(cl=cl, forward=F, loginfo=F)  
        )
      })
      all_p2_hot_start$value_no_reg <- wrap_pn_no_reg(all_p2_hot_start$par)
      return(all_p2_hot_start)
    }
  )
  rlist::list.save(x = a_pns_p2, file = 'results/plots/causality/poc/a_pns_p2.RData')
  tryCatch(parallel::stopCluster(cl)) 
}

# PN all
{
  # start values
  penalty_weights <- 0
  poc <- xvine::proba_necessary_causation
  wrap_pn_no_reg <-  xvine::wrapper_pn_all(
    data = cond_sims_svines_aic_exp,
    col_source = conditional_on_col, 
    u0_target = q_target,
    u0_source = q_source,
    poc = poc,
    lambda = NULL, p=1
  )
  w_t <- rep(.5, (ncol(pollution_data)) * (k_markov - 1))
  w_t <- w_t / sum(w_t)
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c('wrap_pn_no_reg', 'wrap_pn', 'cond_sims_svines_aic_exp', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
  system.time({
    start_w_deoptim_a_pn <- DEoptim::DEoptim(
      lower = rep(0, length(w_t)),
      upper = rep(1, length(w_t)),
      fn = function(w){-wrap_pn_no_reg(w) + (sum(w)-1)^2*penalty_weights},
      control = DEoptim::DEoptim.control(
        trace = T, parallelType=1, strategy=6, cl=cl, reltol = 1e-1, NP=20*length(w_t), p=.5, steptol=50, itermax=20)
    )
    print(start_w_deoptim_a_pn$optim$bestmem)
  })
  tryCatch(parallel::stopCluster(cl)) 
  wrap_pn_no_reg(start_w_deoptim_a_pn$optim$bestmem)
  print(sum(start_w_deoptim_a_pn$optim$bestmem))
  
  #p1 p2
  penalty_weights <- 0
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  a_pn_p1 <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all(
        data = cond_sims_svines_aic_exp,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=1
      )
      fn_opt <- function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights}
      parallel::clusterExport(cl, c('wrap_pn', 'cond_sims_svines_aic_exp', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
      system.time({
        all_p1_hot_start <- optimParallel::optimParallel(
          par=start_w_deoptim_a_pns$optim$bestmem, fn = fn_opt,
          lower = rep(0, length(w_t)), upper = rep(1, length(w_t)),
          control = list(trace=1), parallel = list(cl=cl, forward=F, loginfo=F)  
        )
      })
      all_p1_hot_start$value_no_reg <- wrap_pn_no_reg(all_p1_hot_start$par)
      return(all_p1_hot_start)
    }
  )
  rlist::list.save(x = a_pn_p1, file = 'results/plots/causality/poc/a_pn_p1.RData')
  
  a_pn_p2 <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all(
        data = cond_sims_svines_aic_exp,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=2
      )
      w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
      parallel::clusterExport(cl, c('wrap_pn', 'cond_sims_svines_aic_exp', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
      system.time({
        all_p2_hot_start <- optimParallel::optimParallel(
          par=start_w_deoptim_a_pns$optim$bestmem,
          fn = function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights},
          lower = rep(0, length(w_t)),
          upper = rep(1, length(w_t)),
          control = list(trace=0),
          parallel = list(cl=cl, forward=F, loginfo=F)  
        )
      })
      all_p2_hot_start$value_no_reg <- wrap_pn_no_reg(all_p2_hot_start$par)
      return(all_p2_hot_start)
    }
  )
  rlist::list.save(x = a_pn_p2, file = 'results/plots/causality/poc/a_pn_p2.RData')
  tryCatch(parallel::stopCluster(cl)) 
}

# PS all
{
  # start values
  penalty_weights <- 0
  poc <- xvine::proba_sufficient_causation
  wrap_pn_no_reg <-  xvine::wrapper_pn_all(
    data = cond_sims_svines_aic_exp,
    col_source = conditional_on_col, 
    u0_target = q_target,
    u0_source = q_source,
    poc = poc,
    lambda = NULL, p=1
  )
  w_t <- rep(.5, (ncol(pollution_data)) * (k_markov - 1))
  w_t <- w_t / sum(w_t)
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c('wrap_pn_no_reg', 'wrap_pn', 'cond_sims_svines_aic_exp', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
  system.time({
    start_w_deoptim_a_ps <- DEoptim::DEoptim(
      lower = rep(0, length(w_t)),
      upper = rep(1, length(w_t)),
      fn = function(w){-wrap_pn_no_reg(w) + (sum(w)-1)^2*penalty_weights},
      control = DEoptim::DEoptim.control(
        trace = T, parallelType=1, strategy=6, cl=cl, reltol = 1e-1, NP=20*length(w_t), p=.5, steptol=50, itermax=20)
    )
    print(start_w_deoptim_a_ps$optim$bestmem)
  })
  tryCatch(parallel::stopCluster(cl)) 
  wrap_pn_no_reg(start_w_deoptim_a_ps$optim$bestmem)
  print(sum(start_w_deoptim_a_ps$optim$bestmem))
  
  #p1 p2
  penalty_weights <- 0
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  a_ps_p1 <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all(
        data = cond_sims_svines_aic_exp,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=1
      )
      fn_opt <- function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights}
      parallel::clusterExport(cl, c('wrap_pn', 'cond_sims_svines_aic_exp', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
      system.time({
        all_p1_hot_start <- optimParallel::optimParallel(
          par=start_w_deoptim_a_pns$optim$bestmem, fn = fn_opt,
          lower = rep(0, length(w_t)), upper = rep(1, length(w_t)),
          control = list(trace=1), parallel = list(cl=cl, forward=F, loginfo=F)  
        )
      })
      all_p1_hot_start$value_no_reg <- wrap_pn_no_reg(all_p1_hot_start$par)
      return(all_p1_hot_start)
    }
  )
  rlist::list.save(x = a_ps_p1, file = 'results/plots/causality/poc/a_ps_p1.RData')
  
  a_ps_p2 <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all(
        data = cond_sims_svines_aic_exp,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=2
      )
      w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
      parallel::clusterExport(cl, c('wrap_pn', 'cond_sims_svines_aic_exp', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
      system.time({
        all_p2_hot_start <- optimParallel::optimParallel(
          par=start_w_deoptim_a_pns$optim$bestmem,
          fn = function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights},
          lower = rep(0, length(w_t)),
          upper = rep(1, length(w_t)),
          control = list(trace=0),
          parallel = list(cl=cl, forward=F, loginfo=F)  
        )
      })
      all_p2_hot_start$value_no_reg <- wrap_pn_no_reg(all_p2_hot_start$par)
      return(all_p2_hot_start)
    }
  )
  rlist::list.save(x = a_ps_p2, file = 'results/plots/causality/poc/a_ps_p2.RData')
  tryCatch(parallel::stopCluster(cl)) 
}



# LOAD sparsity results
{
  aws_pns_p1 <- rlist::list.load('results/plots/causality/poc/aws_pns_p1.RData')
  aws_pns_p2 <- rlist::list.load('results/plots/causality/poc/aws_pns_p2.RData')
  a_pns_p1 <- rlist::list.load('results/plots/causality/poc/a_pns_p1.RData')
  a_pns_p2 <- rlist::list.load('results/plots/causality/poc/a_pns_p2.RData')
  
  aws_pn_p1 <- rlist::list.load('results/plots/causality/poc/aws_pn_p1.RData')
  aws_pn_p2 <-  rlist::list.load('results/plots/causality/poc/aws_pn_p2.RData')
  a_pn_p1 <- rlist::list.load('results/plots/causality/poc/a_pn_p1.RData')
  a_pn_p2 <-  rlist::list.load('results/plots/causality/poc/a_pn_p2.RData')
  
  aws_ps_p1 <- rlist::list.load('results/plots/causality/poc/aws_ps_p1.RData')
  aws_ps_p2 <- rlist::list.load('results/plots/causality/poc/aws_ps_p2.RData')
  a_ps_p1 <- rlist::list.load('results/plots/causality/poc/a_ps_p1.RData')
  a_ps_p2 <-  rlist::list.load('results/plots/causality/poc/a_ps_p2.RData')

  all_pns <- rlist::list.load('results/plots/causality/poc/all_pns.RData')
}

lmbds <- list(.0001, .001, .01, .1, .5, 1, 5, 10, 100, 200)

values_no_reg_pns_p1 <- vapply(aws_pns_p1, function(res){
  exp(res$value_no_reg)
}, 1.0)
values_no_reg_pns_p2 <- vapply(aws_pns_p2, function(res){
  exp(res$value_no_reg)
}, 1.0)
values_no_reg_pns_p1_all <- vapply(a_pns_p1, function(res){
  exp(res$value_no_reg)
}, 1.0)
values_no_reg_pns_p2_all <- vapply(a_pns_p2, function(res){
  exp(res$value_no_reg)
}, 1.0)
values_no_reg_pn_p1 <- vapply(aws_pn_p1, function(res){
  exp(res$value_no_reg)
}, 1.0)
values_no_reg_pn_p2 <- vapply(aws_pn_p2, function(res){
  exp(res$value_no_reg)
}, 1.0)
values_no_reg_pn_p1_all <- vapply(a_pn_p1, function(res){
  exp(res$value_no_reg)
}, 1.0)
values_no_reg_pn_p2_all <- vapply(a_pn_p2, function(res){
  exp(res$value_no_reg)
}, 1.0)
values_no_reg_ps_p1 <- vapply(aws_ps_p1, function(res){
  exp(res$value_no_reg)
}, 1.0)
values_no_reg_ps_p2 <- vapply(aws_ps_p2, function(res){
  exp(res$value_no_reg)
}, 1.0)
values_no_reg_ps_p1_all <- vapply(a_ps_p1, function(res){
  exp(res$value_no_reg)
}, 1.0)
values_no_reg_ps_p2_all <- vapply(a_ps_p2, function(res){
  exp(res$value_no_reg)
}, 1.0)


pns_p1_mat <- lapply(aws_pns_p1, function(res){
  ww <- matrix(res$par/sum(res$par), ncol = ncol(pollution_data) - 1, k_markov-1)
  colnames(ww) <- head(colnames(pollution_data), -1)
  rownames(ww) <- paste('t +', seq_len(nrow(ww)))
  ww
})
pns_p2_mat <- lapply(aws_pns_p2, function(res){
  ww <- matrix(res$par/sum(res$par), ncol(pollution_data) - 1, k_markov-1)
  colnames(ww) <- head(colnames(pollution_data), -1)
  rownames(ww) <- paste('t +', seq_len(nrow(ww)))
  ww
})
pns_p1_mat_all <- lapply(a_pns_p1, function(res){
  ww <- matrix(res$par/sum(res$par), ncol = ncol(pollution_data) , k_markov-1)
  colnames(ww) <- colnames(pollution_data)
  rownames(ww) <- paste('t +', seq_len(nrow(ww)))
  ww
})
pns_p2_mat_all <- lapply(a_pns_p2, function(res){
  ww <- matrix(res$par/sum(res$par), ncol = ncol(pollution_data) , k_markov-1)
  colnames(ww) <- colnames(pollution_data)
  rownames(ww) <- paste('t +', seq_len(nrow(ww)))
  ww
})
pn_p1_mat <- lapply(aws_pn_p1, function(res){
  ww <- matrix(res$par/sum(res$par), ncol = ncol(pollution_data) - 1, k_markov-1)
  colnames(ww) <- head(colnames(pollution_data), -1)
  rownames(ww) <- paste('t +', seq_len(nrow(ww)))
  ww
})
pn_p2_mat <- lapply(aws_pn_p2, function(res){
  ww <- matrix(res$par/sum(res$par), ncol = ncol(pollution_data) - 1, k_markov-1)
  colnames(ww) <- head(colnames(pollution_data), -1)
  rownames(ww) <- paste('t +', seq_len(nrow(ww)))
  ww
})
pn_p1_mat_all <- lapply(a_pn_p1, function(res){
  ww <- matrix(res$par/sum(res$par), ncol = ncol(pollution_data) , k_markov-1)
  colnames(ww) <- colnames(pollution_data)
  rownames(ww) <- paste('t +', seq_len(nrow(ww)))
  ww
})
pn_p2_mat_all <- lapply(a_pn_p2, function(res){
  ww <- matrix(res$par/sum(res$par), ncol = ncol(pollution_data) , k_markov-1)
  colnames(ww) <- colnames(pollution_data)
  rownames(ww) <- paste('t +', seq_len(nrow(ww)))
  ww
})
ps_p1_mat <- lapply(aws_ps_p1, function(res){
  ww <- matrix(res$par/sum(res$par), ncol = ncol(pollution_data) - 1, k_markov-1)
  colnames(ww) <- head(colnames(pollution_data), -1)
  rownames(ww) <- paste('t +', seq_len(nrow(ww)))
  ww
})
ps_p2_mat <- lapply(aws_ps_p2, function(res){
  ww <- matrix(res$par/sum(res$par), ncol = ncol(pollution_data) - 1, k_markov-1)
  colnames(ww) <- head(colnames(pollution_data), -1)
  rownames(ww) <- paste('t +', seq_len(nrow(ww)))
  ww
})
ps_p1_mat_all <- lapply(a_ps_p1, function(res){
  ww <- matrix(res$par/sum(res$par), ncol = ncol(pollution_data) , k_markov-1)
  colnames(ww) <- colnames(pollution_data)
  rownames(ww) <- paste('t +', seq_len(nrow(ww)))
  ww
})
ps_p2_mat_all <- lapply(a_ps_p2, function(res){
  ww <- matrix(res$par/sum(res$par), ncol = ncol(pollution_data) , k_markov-1)
  colnames(ww) <- colnames(pollution_data)
  rownames(ww) <- paste('t +', seq_len(nrow(ww)))
  ww
})
all_pns_p1_mat <- lapply(all_pns, function(res){
  ww <- matrix(res$par/sum(res$par), ncol = ncol(pollution_data), nrow = k_markov-1)
  colnames(ww) <- colnames(pollution_data)
  rownames(ww) <- paste('t +', seq_len(nrow(ww)))
  ww
})

entropy_pns_p1 <- vapply(pns_p1_mat, function(par){-sum(par * log(par+.Machine$double.eps))}, 1.0)
entropy_pns_p2 <- vapply(pns_p2_mat, function(par){-sum(par * log(par+.Machine$double.eps))}, 1.0)
entropy_pns_p1_all <- vapply(pns_p1_mat_all, function(par){-sum(par * log(par+2*.Machine$double.eps))}, 1.0)
entropy_pns_p2_all <- vapply(pns_p2_mat_all, function(par){-sum(par * log(par+2*.Machine$double.eps))}, 1.0)
entropy_pn_p1 <- vapply(pn_p1_mat, function(par){-sum(par * log(par+.Machine$double.eps))}, 1.0)
entropy_pn_p2 <- vapply(pn_p2_mat, function(par){-sum(par * log(par+.Machine$double.eps))}, 1.0)
entropy_pn_p1_all <- vapply(pn_p1_mat_all, function(par){-sum(par * log(par+2*.Machine$double.eps))}, 1.0)
entropy_pn_p2_all <- vapply(pn_p2_mat_all, function(par){-sum(par * log(par+2*.Machine$double.eps))}, 1.0)
entropy_ps_p1 <- vapply(ps_p1_mat, function(par){-sum(par * log(par+.Machine$double.eps))}, 1.0)
entropy_ps_p2 <- vapply(ps_p2_mat, function(par){-sum(par * log(par+.Machine$double.eps))}, 1.0)
entropy_ps_p1_all <- vapply(ps_p1_mat_all, function(par){-sum(par * log(par+2*.Machine$double.eps))}, 1.0)
entropy_ps_p2_all <- vapply(ps_p2_mat_all, function(par){-sum(par * log(par+2*.Machine$double.eps))}, 1.0)

par(mfrow=c(2, 1), mar=c(4,5,3,3))
plot(lmbds, values_no_reg_pns_p1, log='x', type='o', ylim=c(0.02, .2), lwd=2, xlab='', cex.lab=1.3, ylab='PNS')
title(xlab=expression(lambda), cex.lab=1.8)
lines(lmbds, values_no_reg_pns_p2, type='o', lwd=2, col='blue')
lines(lmbds, values_no_reg_pns_all_p1, type='o', lwd=2, col='red')

plot(lmbds, values_no_reg_pns_p1/values_no_reg_pns_p1[1], 
     log='x', type='o', ylim=c(.0, 15), lwd=2, xlab='', cex.lab=1.3, ylab='PNS Factor of improvement')
title(xlab=expression(lambda), cex.lab=1.8)
lines(lmbds, values_no_reg_pns_p2/values_no_reg_pns_p2[1], type='o', lwd=2, col='blue')
lines(lmbds, values_no_reg_pns_all_p1/values_no_reg_pns_all_p1[1], type='o', lwd=2, col='red')

plot(lmbds, values_no_reg_pn_p1, log='x', type='o', lwd=2, ylim=c(.95, 1.05), xlab='', cex.lab=1.3, ylab='PNS')
title(xlab=expression(lambda), cex.lab=1.8)
lines(lmbds, values_no_reg_pn_p2, type='o', lwd=2, col='blue')


std_mat <- function(x, std=T){
  if(std){
    (x-min(x))/(max(x)-min(x))
  }else{
    x
  }
}

std_flag  <- T
tl.cex.val <- 1.5

# small examples
pdf('results/plots/causality/files/pns_pn_ps_all.pdf', height = 4, width = 12)
par(mfrow=c(1, 3), mar=c(3,5,3,3))
corrplot::corrplot(
  std_mat(pns_p1_mat_all[[1]], std=std_flag),
  cl.lim = c(-.Machine$double.eps,1), outline = "black", shade.lwd=20, cl.pos = "n", 
  method = "color", tl.cex=tl.cex.val, tl.col='black', tl.srt = 45
)
corrplot::corrplot(
  std_mat(pn_p1_mat_all[[1]], std=std_flag),
  cl.lim = c(-.Machine$double.eps,1), outline = "black", shade.lwd=20, cl.pos = "n", 
  method = "color", tl.cex=tl.cex.val, tl.col='black', tl.srt = 45
)
corrplot::corrplot(
  std_mat(ps_p1_mat_all[[1]], std=std_flag),
  cl.lim = c(-.Machine$double.eps,1), outline = "black", shade.lwd=20, cl.pos = "n", 
  method = "color", tl.cex=tl.cex.val, tl.col='black', tl.srt = 45
)
dev.off()
graphics.off()


pdf('results/plots/causality/files/with_without_target.pdf', height = 7, width = 12)
par(mfrow=c(2, 2), mar=c(3,5,3,3))
layout(
  matrix(c(1,2,3,4), nrow=2, ncol=2, byrow = T),
  width=c(1,1), heights = c(1, 1)
) 
corrplot::corrplot(
  std_mat(pns_p1_mat_all[[1]], std=std_flag),
  cl.lim = c(-.Machine$double.eps,1), outline = "black", shade.lwd=20, cl.pos = "n", 
  method = "color", tl.cex=tl.cex.val, tl.col='black', tl.srt = 45
)
corrplot::corrplot(
  std_mat(pns_p1_mat[[1]], std=std_flag),
  cl.lim = c(-.Machine$double.eps,1), outline = "black", shade.lwd=20, cl.pos = "n", 
  method = "color", tl.cex=tl.cex.val, tl.col='black', tl.srt = 45
)
barplot(apply(pns_p1_mat_all[[1]], 2, sum)*100, col=plotting_colours[2], cex.names = 1.4, cex.axis = 1.4)
barplot(apply(pns_p1_mat[[1]], 2, sum)*100, col=plotting_colours[2], cex.names = 1.4, cex.axis = 1.4)
dev.off()
graphics.off()



# PN Entropy and proba with reg
{
pdf('results/plots/causality/files/entropy_pn.pdf', height = 8, width = 11)
select_mats <- c(1, 3, 4, 6, 8)
par(mfrow=c(3, 1), mar=c(5,5,3,3))
layout(
  matrix(c(1,2,3), nrow=3, ncol=1, byrow = T),
  width=c(1,1,1), heights = c(3,3,1.5)
) 
lmbs_plot <- c(0.0001, unlist(lmbds[-1]))
plot(entropy_pn_p1_all / log(42), type = 'o', log='x', ylim=c(0.0,1), xaxt="n", ylab='Relative Entropy',
     cex.lab=cex_lab, xlab='', bty="n", lwd=3, col=plotting_colours[1], cex.axis=cex_axis, cex=2, pch=3)
text(x=seq_len(10), par("usr")[3], labels=c(0, lmbds[-1]), pos=1, xpd=T, 
     xlab=expression(lambda), cex=2, cex.axis=cex_axis)
title(xlab=expression(lambda), line=4, cex.lab=cex_lab*1.2)
lines( entropy_pn_p2_all / log(42), type = 'o', cex=2, pch=4,
       lwd=3, col=plotting_colours[2], lty=2)
lines( entropy_pn_p1 / log(36), type = 'o', cex=2, pch=5,
       lwd=3, col=plotting_colours[3], lty=3)
lines(entropy_pn_p2 / log(36), type = 'o', cex=2, pch=6,
      lwd=3, col=plotting_colours[4], lty=4)
# legend("bottomleft", legend=c('with v/hr L1', 'with v/hr L2', 'without v/hr L1', 'without v/hr L2'),
#        lty=1, lwd=4,  col=plotting_colours, cex=1.3, bty = "n", pt.cex = 2)

plot(values_no_reg_pn_p1_all, type = 'o', log='x', ylim=c(0.0,1), xaxt="n", ylab='PN',
     cex.lab=cex_lab, xlab='', cex=2, pch=3,
     bty="n", lwd=3, col=plotting_colours[1], cex.axis=cex_axis)
text(x=seq_len(10), par("usr")[3], labels=c(0, lmbds[-1]), pos=1, xpd=T, xlab=expression(lambda), cex=2)
title(xlab=expression(lambda), line=4, cex.lab=cex_lab*1.2)
lines(values_no_reg_pn_p2_all, type = 'o',
       lwd=3, col=plotting_colours[2], lty=2, cex=2, pch=4)
lines(values_no_reg_pn_p1, type = 'o',
       lwd=3, col=plotting_colours[3], lty=3, cex=2, pch=5)
lines(values_no_reg_pn_p2, type = 'o',
      lwd=3, col=plotting_colours[4], lty=4, cex=2, pch=6)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend=c('with v/hr L1', 'with v/hr L2', 'without v/hr L1', 'without v/hr L2'),
       lty=1:4, lwd=3, pch=3:6, col=plotting_colours, cex=1.8, bty = "n", pt.cex = 3, horiz=TRUE,  pt.bg = 'white')
dev.off()
graphics.off()
}

# PS Entropy and proba with reg
{
  pdf('results/plots/causality/files/entropy_ps.pdf', height = 8, width = 11)
  select_mats <- c(1, 3, 4, 6, 8)
  par(mfrow=c(3, 1), mar=c(5,5,3,3))
  layout(
    matrix(c(1,2,3), nrow=3, ncol=1, byrow = T),
    width=c(1,1,1), heights = c(3,3,1.5)
  ) 
  lmbs_plot <- c(0.0001, unlist(lmbds[-1]))
  plot(entropy_ps_p1_all / log(42), type = 'o', log='x', ylim=c(0.0,1), xaxt="n", ylab='Relative Entropy',
       cex.lab=cex_lab, xlab='', bty="n", lwd=3, col=plotting_colours[1], cex.axis=cex_axis, cex=2, pch=3)
  text(x=seq_len(10), par("usr")[3], labels=c(0, lmbds[-1]), pos=1, xpd=T, 
       xlab=expression(lambda), cex=2, cex.axis=cex_axis)
  title(xlab=expression(lambda), line=4, cex.lab=cex_lab*1.2)
  lines( entropy_ps_p2_all / log(42), type = 'o', cex=2, pch=4,
         lwd=3, col=plotting_colours[2], lty=2)
  lines( entropy_ps_p1 / log(36), type = 'o', cex=2, pch=5,
         lwd=3, col=plotting_colours[3], lty=3)
  lines(entropy_ps_p2 / log(36), type = 'o', cex=2, pch=6,
        lwd=3, col=plotting_colours[4], lty=4)
  # legend("bottomleft", legend=c('with v/hr L1', 'with v/hr L2', 'without v/hr L1', 'without v/hr L2'),
  #        lty=1, lwd=4,  col=plotting_colours, cex=1.3, bty = "n", pt.cex = 2)
  
  plot(values_no_reg_ps_p1_all, type = 'o', log='x', ylim=c(0.0,1), xaxt="n", ylab='PS',
       cex.lab=cex_lab, xlab='', cex=2, pch=3,
       bty="n", lwd=3, col=plotting_colours[1], cex.axis=cex_axis)
  text(x=seq_len(10), par("usr")[3], labels=c(0, lmbds[-1]), pos=1, xpd=T, xlab=expression(lambda), cex=2)
  title(xlab=expression(lambda), line=4, cex.lab=cex_lab*1.2)
  lines(values_no_reg_ps_p2_all, type = 'o',
        lwd=3, col=plotting_colours[2], lty=2, cex=2, pch=4)
  lines(values_no_reg_ps_p1, type = 'o',
        lwd=3, col=plotting_colours[3], lty=3, cex=2, pch=5)
  lines(values_no_reg_ps_p2, type = 'o',
        lwd=3, col=plotting_colours[4], lty=4, cex=2, pch=6)
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("center", legend=c('with v/hr L1', 'with v/hr L2', 'without v/hr L1', 'without v/hr L2'),
         lty=1:4, lwd=3, pch=3:6, col=plotting_colours, cex=1.8, bty = "n", pt.cex = 3, horiz=TRUE,  pt.bg = 'white')
  dev.off()
  graphics.off()
}

# PNS Entropy and proba with reg
{
pdf('results/plots/causality/files/entropy.pdf', height = 8, width = 11)
select_mats <- c(1, 3, 4, 6, 8)
par(mfrow=c(3, 1), mar=c(5,5,3,3))
layout(
  matrix(c(1,2,3), nrow=3, ncol=1, byrow = T),
  width=c(1,1,1), heights = c(3,3,1.5)
) 
lmbs_plot <- c(0.0001, unlist(lmbds[-1]))
plot(entropy_pns_p1_all / log(42), type = 'o', log='x', ylim=c(0.0,1), xaxt="n", ylab='Relative Entropy',
     cex.lab=cex_lab, xlab='', bty="n", lwd=3, col=plotting_colours[1], cex.axis=cex_axis, cex=2, pch=3)
text(x=seq_len(10), par("usr")[3], labels=c(0, lmbds[-1]), pos=1, xpd=T, 
     xlab=expression(lambda), cex=2, cex.axis=cex_axis)
title(xlab=expression(lambda), line=4, cex.lab=cex_lab*1.2)
lines( entropy_pns_p2_all / log(42), type = 'o', cex=2, pch=4,
       lwd=3, col=plotting_colours[2], lty=2)
lines( entropy_pns_p1 / log(36), type = 'o', cex=2, pch=5,
       lwd=3, col=plotting_colours[3], lty=3)
lines(entropy_pns_p2 / log(36), type = 'o', cex=2, pch=6,
      lwd=3, col=plotting_colours[4], lty=4)
# legend("bottomleft", legend=c('with v/hr L1', 'with v/hr L2', 'without v/hr L1', 'without v/hr L2'),
#        lty=1, lwd=4,  col=plotting_colours, cex=1.3, bty = "n", pt.cex = 2)

plot(values_no_reg_pns_p1_all, type = 'o', log='x', ylim=c(0.0,1), xaxt="n", ylab='PNS',
     cex.lab=cex_lab, xlab='', cex=2, pch=3,
     bty="n", lwd=3, col=plotting_colours[1], cex.axis=cex_axis)
text(x=seq_len(10), par("usr")[3], labels=c(0, lmbds[-1]), pos=1, xpd=T, xlab=expression(lambda), cex=2)
title(xlab=expression(lambda), line=4, cex.lab=cex_lab*1.2)
lines(values_no_reg_pns_p2_all, type = 'o',
      lwd=3, col=plotting_colours[2], lty=2, cex=2, pch=4)
lines(values_no_reg_pns_p1, type = 'o',
      lwd=3, col=plotting_colours[3], lty=3, cex=2, pch=5)
lines(values_no_reg_pns_p2, type = 'o',
      lwd=3, col=plotting_colours[4], lty=4, cex=2, pch=6)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend=c('with v/hr L1', 'with v/hr L2', 'without v/hr L1', 'without v/hr L2'),
       lty=1:4, lwd=3, pch=3:6, col=plotting_colours, cex=1.8, bty = "n", pt.cex = 3, horiz=TRUE,  pt.bg = 'white')
dev.off()
graphics.off()
}



# large examples
pdf('results/plots/causality/files/pns_matrices.pdf', height = 9, width = 20)
select_mats <- c(1, 3, 4, 6, 8, 9)
par(mfrow=c(4, length(select_mats) + 1), mai=c(0,0,0,0))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend=c('\tPNS \nwith v/hr - L1'), cex=2, bty = "n", pt.cex = 3, horiz=TRUE,  pt.bg = 'white')
lapply(pns_p1_mat_all[select_mats], function(mat){
  corrplot::corrplot(
    std_mat(mat, std=std_flag),
    cl.lim = c(-.Machine$double.eps,1), outline = "black", shade.lwd=20, cl.pos = "n",
    method = "color", tl.cex=tl.cex.val, tl.col='black', tl.srt = 45
  )
})
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend=c('\tPNS \nwith v/hr - L2'), cex=2, bty = "n", pt.cex = 3, horiz=TRUE,  pt.bg = 'white')
lapply(pns_p2_mat_all[select_mats], function(mat){
  corrplot::corrplot(
    std_mat(mat, std=std_flag),
    cl.lim = c(-.Machine$double.eps,1), outline = "black", shade.lwd=20, cl.pos = "n",
    method = "color", tl.cex=tl.cex.val, tl.col='black', tl.srt = 45, 
  )
})
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend=c('\tPNS \nwithout v/hr - L1'), cex=2, bty = "n", pt.cex = 3, horiz=TRUE,  pt.bg = 'white')
lapply(pns_p1_mat[select_mats], function(mat){
         corrplot::corrplot(
          std_mat(mat, std=std_flag),
           cl.lim = c(-.Machine$double.eps,1), outline = "black", shade.lwd=20, cl.pos = "n", 
           method = "color", tl.cex=tl.cex.val, tl.col='black', tl.srt = 45
         )
       })
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend=c('\tPNS \nwithout v/hr - L2'), cex=2, bty = "n", pt.cex = 3, horiz=TRUE,  pt.bg = 'white')
lapply(pns_p2_mat[select_mats], function(mat){
         corrplot::corrplot(
          std_mat(mat, std=std_flag),
           cl.lim = c(-.Machine$double.eps,1), outline = "black", shade.lwd=20, cl.pos = "n",
           method = "color", tl.cex=tl.cex.val, tl.col='black', tl.srt = 45
         )
       })
dev.off()
graphics.off()


lapply(pn_p1_mat[select_mats], function(mat){
         corrplot::corrplot(
          std_mat(mat, std=std_flag),
           cl.lim = c(-.Machine$double.eps,1), outline = "black", shade.lwd=20, cl.pos = "n",
           method = "shade", tl.cex=tl.cex.val, tl.col='black', tl.srt = 45
         )
       })
lapply(pn_p2_mat[select_mats], function(mat){
         corrplot::corrplot(
          std_mat(mat, std=std_flag),
           cl.lim = c(-.Machine$double.eps,1), outline = "black", shade.lwd=20, cl.pos = "n",
           method = "shade", tl.cex=tl.cex.val, tl.col='black', tl.srt = 45
         )
       })
lapply(ps_p1_mat[select_mats], function(mat){
         corrplot::corrplot(
          std_mat(mat, std=std_flag), 
           cl.lim = c(-.Machine$double.eps,1), outline = "black", shade.lwd=20, cl.pos = "n",
           method = "shade", tl.cex=tl.cex.val, tl.col='black', tl.srt = 45
         )
       })
lapply(ps_p2_mat[select_mats], function(mat){
         corrplot::corrplot(
          std_mat(mat, std=std_flag),
           cl.lim = c(-.Machine$double.eps,1), outline = "black", shade.lwd=20, cl.pos = "n",
           method = "shade", tl.cex=tl.cex.val, tl.col='black', tl.srt = 45
         )
       })
lapply(all_pns_p1_mat[select_mats], function(mat){
  corrplot::corrplot(
    std_mat(mat, std=std_flag),
    cl.lim = c(-.Machine$double.eps,1), outline = "black", shade.lwd=20, cl.pos = "n",
    method = "shade", tl.cex=tl.cex.val, tl.col='black', tl.srt = 45
  )
})

# aic_svines <- lapply(seq_len(12), function(k){
#   xvine::fit_markov_svines(
#     data = pollution_data_it[1:5000,],
#     k.markov = k,
#     family_set="archimedean",
#     selcrit="aic",
#     tree_crit="tau",
#     threshold=0.05
#   )
# })
# rlist::list.save(x = aic_svines, file = 'results/plots/causality/poc/aic_svines.RData')
aic_svines <- rlist::list.load('results/plots/causality/poc/aic_svines.RData')

par(mfrow=c(1,1), mar=c(3,3,3,3))
aic_values <- vapply(aic_svines, function(svine){2*svine$copula$npars - 2*svine$copula$loglik}, 1.0)
bic_values <- vapply(aic_svines, function(svine){2*svine$copula$npars * log(svine$copula$nobs) - 2*svine$copula$loglik}, 1.0)
mbicv_values <- vapply(aic_svines, function(svine){rvinecopulib::mBICV(svine$copula)}, 1.0)

fit_vals <- rbind(aic_values, bic_values, mbicv_values)
colnames(fit_vals) <- seq_len(12)
write.csv(x = round(fit_vals, 0), 'results/plots/causality/poc/aic_svines_fitting_criteria.csv')

plot(-log(abs(aic_values)), type='o', ylim = c(-10.08, -9.98))
lines(-log(abs(bic_values)), type='o')

plot(-log(abs(mbicv_values)), type='o')

plot(abs(aic_values)/max(abs(aic_values)), type='o')
lines(abs(bic_values)/max(abs(bic_values)))
lines(abs(mbicv_values)/max(abs(mbicv_values)))

aic_svine_argmin <- aic_svines[[which.min(aic_values)]]

w_t <- rep(.5, (ncol(pollution_data) ) * (k_markov - 1))
w_t <- w_t / sum(w_t)


