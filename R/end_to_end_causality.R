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

d_values <- apply(pollution_data, MARGIN = 2, FUN = function(x){qq <- fracdiff::fdGPH(x); qq$data <- x; qq})
all_data <- lapply(d_values, function(x){x$fd <- fracdiff::diffseries(x$data, d = x$d); x})
p_data <- do.call(cbind, lapply(all_data, function(x){x$fd}))

# pollution_data2 <- apply(
#   pollution_data, MARGIN = 2,
#   function(x){stats::stl(ts(x, frequency = 24), s.window = 24, t.window = 24*7)$time.series[1:max_length, 3]}
# )
# p_data <- pollution_data2


least_fracdiff <- function(data, level=0.05){
  low <- 0
  high <- 1
  vals <- aTSA::stationary.test(data, method='kpss', output = F) # tseries::adf.test(data)
  while (abs(high-low) > 1e-3) {
    mid <- (low+high)/2
    diff_data <- fracdiff::diffseries(data, mid)
    vals <- aTSA::stationary.test(data, method='kpss', output = F)
    p.val <- min(vals[,3]) # tseries::adf.test(diff_data)
    print(p.val)
    if(p.val < level){
        low <- mid
    }else{
        high <- mid
    }
  }
  
  return(list(d=low, p.value=min(vals[,3]), statistic=max(vals[,2])))
}

pollution_data <- t(t(p_data) / apply(p_data, MARGIN = 2, sd))
colnames(pollution_data) <- c(colnames(pollution_data)[1:6], 'v/hr')

# 
# pollution_data <- apply(
#   pollution_data, MARGIN = 2,
#   function(x){stats::stl(ts(x, frequency = 24), s.window = 24, t.window = 24*7)$time.series[1:max_length_val, 3]}
# )
# pollution_data_val <- pollution_data[(max_length+1):max_length_val,]
# pollution_data_val <- cbind(pollution_data_val, traffic_data_val)
# pollution_data <- pollution_data[1:max_length,]
# pollution_data <- cbind(pollution_data, traffic_data)
# pollution_data_val <- t(t(pollution_data_val) / apply(pollution_data, MARGIN = 2, sd))
# pollution_data <- t(t(pollution_data) / apply(pollution_data, MARGIN = 2, sd))
# 
# colnames(pollution_data) <- c(colnames(pollution_data)[1:6], 'v/hr')

set.seed(42)

# Hyperparamters
u0_single <- .95
k_markov <- 7
conditional_on_col <- 7

u0s_fixed <- rep(u0_single, ncol(pollution_data))
u0s_fixed[7] <- 0.85
pollution_data_it <- xvine::apply_integral_transform(pollution_data, u0s_fixed)
#unfold
params_gpd_ests <- pollution_data_it$par.ests
params_gpd_ses <- pollution_data_it$par.ses
pollution_data_it <- pollution_data_it$data

pollution_data_val_it <- xvine::apply_integral_transform(pollution_data_val, u0s_fixed)
#unfold
params_gpd_ests_val <- pollution_data_val_it$par.ests
params_gpd_ses_val <- pollution_data_val_it$par.ses
pollution_data_val_it <- pollution_data_val_it$data

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

rit_exp_stack_val <- lapply(
  data.frame(t(build_stack(xvine::apply_reverse_exponential(pollution_data_val_it)$data, k = k_markov))), 
  function(x) matrix(x, byrow = T, ncol=ncol(pollution_data_val_it))
) # returns a list
rit_exp_stack_val <- abind::abind(rit_exp_stack_val, along = 3)

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
    col_idx <- ncol(pollution_data) * 0:(k_markov-1)  + target
    
    cond_threshold <- quantile(stack_pollution_data[,conditional_on_col], u0s_fixed[target])
    cond_idx_above <- which(stack_pollution_data[,conditional_on_col] > cond_threshold)
    cond_target_threshold <- quantile(stack_pollution_data[,target], u0s_fixed[target])
    
    it_above <- xvine::integral_transform(
      as.vector(stack_pollution_data[cond_idx_above, col_idx]),
      u0=mean(as.vector(stack_pollution_data[cond_idx_above, col_idx]) < cond_target_threshold)
    )
    if(target == conditional_on_col){
      col_idx_under <- col_idx[-1]
    }else{col_idx_under <- col_idx}
    
    it_below <-  xvine::integral_transform(
      as.vector(stack_pollution_data[-cond_idx_above, col_idx_under]),
      u0=mean(as.vector(stack_pollution_data[-cond_idx_above, col_idx_under]) < cond_target_threshold)
    )
    return(list('above'=it_above, 'above_data'=as.vector(stack_pollution_data[cond_idx_above, col_idx]),
                'below'=it_below, 'below_data'=as.vector(stack_pollution_data[-cond_idx_above, col_idx_under]),
                'conditional_on_col'=conditional_on_col, 'target'=target))
  }
)
f_cf_params_matrix <- vapply(
  f_cf_params, 
  function(f_cf_target){
    prms <- as.vector(
      c(c(rbind(f_cf_target$above$par.ests, f_cf_target$above$par.ses)), f_cf_target$above$u0,
        c(rbind(f_cf_target$below$par.ests, f_cf_target$below$par.ses)), f_cf_target$below$u0)
    )
  },
  FUN.VALUE = rep(0, 10)
)
f_cf_params_matrix <- t(f_cf_params_matrix)

f_cf_params_matrix <- cbind(f_cf_params_matrix, vapply(d_values, function(x){x$d}, 1.0))
f_cf_params_matrix <- cbind(f_cf_params_matrix, vapply(d_values, function(x){x$sd.reg}, 1.0))
f_cf_params_matrix <- cbind(f_cf_params_matrix, apply(p_data, 2, function(x){unlist(tseries::adf.test(x, k = 12)$statistic)}))
colnames(f_cf_params_matrix) <- c('gamma_f', 'gamma_f_sd', 'sigma_f', 'sigma_f_sd', 'u0_f', 'gamma_cf', 'gamma_cf_sd', 'sigma_cf', 'sigma_cf_sd', 'u0_cf', 'd', 'd_sd', 'ADF statistic')
rownames(f_cf_params_matrix) <- colnames(pollution_data) 
write.csv(round(f_cf_params_matrix, 3), file = 'results/plots/causality/f_cf_params_matrix.csv')

par(mfrow=c(7,2), mar=c(3.9,2.2,0,0))
plotting_colours <- RColorBrewer::brewer.pal(n = ncol(pollution_data), name = 'Dark2')
max_quant <- 95
lapply(seq_len(ncol(pollution_data)), function(i){
  # above
  th_quantiles <- evir::qgpd(seq_len(max_quant)/100, mu = quantile(f_cf_params[[i]]$above_data, f_cf_params[[i]]$above$u0),
                             xi = f_cf_params[[i]]$above$par.ests[1], beta = f_cf_params[[i]]$above$par.ests[2])
  empi_quantiles <- quantile(
    f_cf_params[[i]]$above_data[f_cf_params[[i]]$above_data > quantile(f_cf_params[[i]]$above_data, f_cf_params[[i]]$above$u0)],
    seq_len(max_quant)/100)
  plot(th_quantiles, empi_quantiles,
    xlab='GPD quantiles', ylab='Emp. Quantiles', bty='n', cex.main=1.5, cex.lab=1.3, cex=1.4, col=plotting_colours[i], pch=16, cex.axis=1.5)
  lines(c(min(th_quantiles), 5*max(empi_quantiles)), c(min(th_quantiles), 5*max(empi_quantiles)), cex=1.5, lwd=2)
  
  # below
  th_quantiles <- evir::qgpd(seq_len(max_quant)/100, mu = quantile(f_cf_params[[i]]$below_data, f_cf_params[[i]]$below$u0),
                             xi = f_cf_params[[i]]$below$par.ests[1], beta = f_cf_params[[i]]$below$par.ests[2])
  empi_quantiles <- quantile(
    f_cf_params[[i]]$below_data[f_cf_params[[i]]$below_data > quantile(f_cf_params[[i]]$below_data, f_cf_params[[i]]$below$u0)],
    seq_len(max_quant)/100)
  plot(th_quantiles, empi_quantiles,
    xlab='GPD quantiles', ylab='Emp. Quantiles', bty='n', cex.main=1.5, cex.lab=1.3, cex=1.3, col=plotting_colours[i], pch=16, cex.axis=1.5)
  lines(c(min(th_quantiles), 5*max(empi_quantiles)), c(min(th_quantiles), 5*max(empi_quantiles)), cex=1.5, lwd=2)
})

# Exploration plots
layout(
  matrix(seq_len((1+ncol(pollution_data)) * 5), nrow=1+ncol(pollution_data), ncol=5, byrow = T),
  width=c(1,4,4,4,4,4), heights = c(1, rep(2, ncol(pollution_data)), mar=c(4,5,3,2))
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
  title(colnames(pollution_data)[i], cex.main=1.2, line=-4)
  
  plot(save_core_pollution_data[1:10000,i], xlab='Values', bty= 'n',
       main='', ylab='', cex.main=1.8, cex.lab=1.5, cex.axis=1.4, col=plotting_colours[i], type='l')
  hist(pollution_data[, i], probability = T, xlab='Values', breaks = 40, xlim = c(-6, 6),
       main='', ylab='', cex.main=1.8, cex.lab=1.5, cex.axis=1.4, col=plotting_colours[i])

  # above
  th_quantiles <- evir::qgpd(seq_len(max_quant)/100, mu = quantile(f_cf_params[[i]]$above_data, f_cf_params[[i]]$above$u0),
                             xi = f_cf_params[[i]]$above$par.ests[1], beta = f_cf_params[[i]]$above$par.ests[2])
  empi_quantiles <- quantile(
    f_cf_params[[i]]$above_data[f_cf_params[[i]]$above_data > quantile(f_cf_params[[i]]$above_data, f_cf_params[[i]]$above$u0)],
    seq_len(max_quant)/100)
  plot(th_quantiles, empi_quantiles, ylab='',
       xlab='GPD quantiles', bty='n', cex.main=1.5, cex.lab=1.3, cex=1.4, col=plotting_colours[i], pch=16, cex.axis=1.5)
  lines(c(min(th_quantiles), 5*max(empi_quantiles)), c(min(th_quantiles), 5*max(empi_quantiles)), cex=1.5, lwd=2)
  
  # below
  th_quantiles <- evir::qgpd(seq_len(max_quant)/100, mu = quantile(f_cf_params[[i]]$below_data, f_cf_params[[i]]$below$u0),
                             xi = f_cf_params[[i]]$below$par.ests[1], beta = f_cf_params[[i]]$below$par.ests[2])
  empi_quantiles <- quantile(
    f_cf_params[[i]]$below_data[f_cf_params[[i]]$below_data > quantile(f_cf_params[[i]]$below_data, f_cf_params[[i]]$below$u0)],
    seq_len(max_quant)/100)
  plot(th_quantiles, empi_quantiles, ylab='',
       xlab='GPD quantiles', bty='n', cex.main=1.5, cex.lab=1.3, cex=1.3, col=plotting_colours[i], pch=16, cex.axis=1.5)
  lines(c(min(th_quantiles), 5*max(empi_quantiles)), c(min(th_quantiles), 5*max(empi_quantiles)), cex=1.5, lwd=2)
}
}

par(mfrow=c(7,2), mar=c(3.9,2.2,0,0))
layout(
  matrix(seq_len((1+ncol(pollution_data)) * 5), nrow=1+ncol(pollution_data), ncol=5, byrow = T),
  width=c(1,4,4,4,4,4), heights = c(1, rep(2, ncol(pollution_data)), mar=c(4,5,3,2))
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
  title(colnames(pollution_data)[i], cex.main=1.2, line=-4)
  barplot(acf(pollution_data[,i], lag.max = 15, plot=F)$acf[,,1], names.arg=0:15,  bty='n', lwd=3, cex.lab=1.5, cex.axis=1.3, cex.names = 1.4,
          width = .3, space = .8, col=plotting_colours[i], xlab='Lag', border=1)
  abline(h=0, col='black')
  abline(h=qnorm((1 + .975)/2)/sqrt(nrow(pollution_data)), lty=2, col='blue')
  abline(h=-qnorm((1 + .975)/2)/sqrt(nrow(pollution_data)), lty=2, col='blue')
  
  barplot(pacf(pollution_data[,i], lag.max = 15, plot=F)$acf[,,1], names.arg=1:15,  bty='n', lwd=3, cex.lab=1.5, cex.axis=1.3, cex.names = 1.4,
          width = .3, space = .8, col=plotting_colours[i], xlab='Lag', border=1)
  abline(h=0, col='black')
  abline(h=qnorm((1 + .975)/2)/sqrt(nrow(pollution_data)), lty=2, col='blue')
  abline(h=-qnorm((1 + .975)/2)/sqrt(nrow(pollution_data)), lty=2, col='blue')
  
  barplot(acf_vals$acf[1:15, i, 7],names.arg=1:15,  bty='n', lwd=3, cex.lab=1.5, cex.axis=1.3, cex.names = 1.4,
          width = .3, space = .8, col=plotting_colours[i], xlab='Lag', border=1)
  abline(h=0, col='black')
  abline(h=qnorm((1 + .975)/2)/sqrt(nrow(pollution_data)), lty=2, col='blue')
  abline(h=-qnorm((1 + .975)/2)/sqrt(nrow(pollution_data)), lty=2, col='blue')
  
  barplot(pacf_vals$acf[1:15, i, 7],names.arg=1:15,  bty='n', lwd=3, cex.lab=1.5, cex.axis=1.3, cex.names = 1.4,
          width = .3, space = .8, col=plotting_colours[i], xlab='Lag', border=1)
  abline(h=0, col='black')
  abline(h=qnorm((1 + .975)/2)/sqrt(nrow(pollution_data)), lty=2, col='blue')
  abline(h=-qnorm((1 + .975)/2)/sqrt(nrow(pollution_data)), lty=2, col='blue')
}
}

# HANDMADE OPTIM PARALLEL
q_target <- qexp(u0_single)
q_source <- qexp(u0_single)

lmbds <- list(NULL, .001, .01, .1, .5, 1, 5, 10, 100, 200)
# PNS
{
# start values
penalty_weights <- 0
poc <- proba_necessary_sufficient_causation
wrap_pn_no_reg <-  xvine::wrapper_pn_all_but_source(
  data = rit_exp_stack,
  col_source = conditional_on_col, 
  u0_target = q_target,
  u0_source = q_source,
  poc = poc,
  lambda = NULL, p=1
)
w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
w_t <- w_t / sum(w_t)
cl <- parallel::makeCluster(parallel::detectCores()-1)
parallel::clusterExport(cl, c('wrap_pn_no_reg', 'wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
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
      data = rit_exp_stack,
      col_source = conditional_on_col, 
      poc = poc,
      u0_target = q_target,
      u0_source = q_source,
      lambda = l, p=1
    )
    w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
    w_t <- w_t / sum(w_t)
    fn_opt <- function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights}
    parallel::clusterExport(cl, c('wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
    system.time({
      all_p1_hot_start <- optimParallel::optimParallel(
        par=start_w_deoptim$optim$bestmem,
        fn = fn_opt,
        lower = rep(0, length(w_t)),
        upper = rep(1, length(w_t)),
        control = list(trace=0),
        parallel = list(cl=cl, forward=F, loginfo=F)  
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
      data = rit_exp_stack,
      col_source = conditional_on_col, 
      poc = poc,
      u0_target = q_target,
      u0_source = q_source,
      lambda = l, p=2
    )
    w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
    parallel::clusterExport(cl, c('wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
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

# PN
{
  # start values
  penalty_weights <- 0
  poc <- proba_necessary_causation
  wrap_pn_no_reg <-  xvine::wrapper_pn_all_but_source(
    data = rit_exp_stack,
    col_source = conditional_on_col, 
    u0_target = q_target,
    u0_source = q_source,
    poc = poc,
    lambda = NULL, p=1
  )
  w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
  w_t <- w_t / sum(w_t)
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c('wrap_pn_no_reg', 'wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
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
  aws_pn_p1 <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all_but_source(
        data = rit_exp_stack,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=1
      )
      w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
      w_t <- w_t / sum(w_t)
      parallel::clusterExport(cl, c('wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
      system.time({
        all_p1_hot_start <- optimParallel::optimParallel(
          par=start_w_deoptim$optim$bestmem,
          fn = function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights},
          lower = rep(0, length(w_t)),
          upper = rep(1, length(w_t)),
          control = list(trace=0),
          parallel = list(cl=cl, forward=F, loginfo=F)  
        )
      })
      all_p1_hot_start$value_no_reg <- wrap_pn_no_reg(all_p1_hot_start$par)
      return(all_p1_hot_start)
    }
  )
  rlist::list.save(x = aws_pn_p1, file = 'results/plots/causality/poc/aws_pn_p1.RData')
  tryCatch({parallel::stopCluster(cl)}) 
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  aws_pn_p2 <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all_but_source(
        data = rit_exp_stack,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=2
      )
      w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
      parallel::clusterExport(cl, c('wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
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
  rlist::list.save(x = aws_pn_p2, file = 'results/plots/causality/poc/aws_pn_p2.RData')
  tryCatch(parallel::stopCluster(cl)) 
}

# PS
{
  # start values
  penalty_weights <- 0
  poc <- proba_sufficient_causation
  wrap_pn_no_reg <-  xvine::wrapper_pn_all_but_source(
    data = rit_exp_stack,
    col_source = conditional_on_col, 
    u0_target = q_target,
    u0_source = q_source,
    poc = poc,
    lambda = NULL, p=1
  )
  w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
  w_t <- w_t / sum(w_t)
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c('wrap_pn_no_reg', 'wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
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
  
  tryCatch({parallel::stopCluster(cl)}) 
  wrap_pn_no_reg(start_w_deoptim$optim$bestmem)
  print(sum(start_w_deoptim$optim$bestmem))
  
  #p1 p2
  penalty_weights <- 0
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  aws_ps_p1 <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all_but_source(
        data = rit_exp_stack,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=1
      )
      w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
      w_t <- w_t / sum(w_t)
      parallel::clusterExport(cl, c('wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
      system.time({
        all_p1_hot_start <- optimParallel::optimParallel(
          par=start_w_deoptim$optim$bestmem,
          fn = function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights},
          lower = rep(0, length(w_t)),
          upper = rep(1, length(w_t)),
          control = list(trace=0),
          parallel = list(cl=cl, forward=F, loginfo=F)  
        )
      })
      all_p1_hot_start$value_no_reg <- wrap_pn_no_reg(all_p1_hot_start$par)
      return(all_p1_hot_start)
    }
  )
  rlist::list.save(x = aws_ps_p1, file = 'results/plots/causality/poc/aws_ps_p1.RData')
  tryCatch({parallel::stopCluster(cl)}) 
  
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  aws_ps_p2 <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all_but_source(
        data = rit_exp_stack,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=2
      )
      w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
      parallel::clusterExport(cl, c('wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
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
  rlist::list.save(x = aws_ps_p2, file = 'results/plots/causality/poc/aws_ps_p2.RData')
  tryCatch({parallel::stopCluster(cl)}) 
}


# PNS / PN / PS on all variables

# PNS
{
  # start values
  penalty_weights <- 0
  poc <- proba_necessary_sufficient_causation
  wrap_pn_no_reg <-  xvine::wrapper_pn_all(
    data = rit_exp_stack,
    col_source = conditional_on_col, 
    u0_target = q_target,
    u0_source = q_source,
    poc = poc,
    lambda = NULL, p=1
  )
  w_t <- rep(.5, (ncol(pollution_data) ) * (k_markov - 1))
  w_t <- w_t / sum(w_t)
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c('wrap_pn_no_reg', 'wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
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
  tryCatch({parallel::stopCluster(cl)})
  wrap_pn_no_reg(start_w_deoptim$optim$bestmem)
  print(sum(start_w_deoptim$optim$bestmem))
  
  penalty_weights <- 0
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  all_pns <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all(
        data = rit_exp_stack,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=1
      )
      w_t <- rep(.5, (ncol(pollution_data)) * (k_markov - 1))
      w_t <- w_t / sum(w_t)
      fn_opt <- function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights}
      parallel::clusterExport(cl, c('wrap_pn', 'wrap_pn_no_reg', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
      system.time({
        all_p1_hot_start <- optimParallel::optimParallel(
          par=start_w_deoptim$optim$bestmem,
          fn = fn_opt,
          lower = rep(0, length(w_t)),
          upper = rep(1, length(w_t)),
          control = list(trace=0),
          parallel = list(cl=cl, forward=F, loginfo=F)  
        )
      })
      all_p1_hot_start$value_no_reg <- wrap_pn_no_reg(all_p1_hot_start$par)
      return(all_p1_hot_start)
    }
  )
  rlist::list.save(x = all_pns, file = 'results/plots/causality/poc/all_pns.RData')
  tryCatch({parallel::stopCluster(cl)}) 
}

# PN
{
  # start values
  penalty_weights <- 0
  poc <- proba_necessary_causation
  wrap_pn_no_reg <-  xvine::wrapper_pn_all(
    data = rit_exp_stack,
    col_source = conditional_on_col, 
    u0_target = q_target,
    u0_source = q_source,
    poc = poc,
    lambda = NULL, p=1
  )
  w_t <- rep(.5, (ncol(pollution_data) ) * (k_markov - 1))
  w_t <- w_t / sum(w_t)
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c('wrap_pn_no_reg', 'wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
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
  tryCatch({parallel::stopCluster(cl)})
  wrap_pn_no_reg(start_w_deoptim$optim$bestmem)
  print(sum(start_w_deoptim$optim$bestmem))
  
  penalty_weights <- 0
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  all_pns <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all(
        data = rit_exp_stack,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=1
      )
      w_t <- rep(.5, (ncol(pollution_data)) * (k_markov - 1))
      w_t <- w_t / sum(w_t)
      fn_opt <- function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights}
      parallel::clusterExport(cl, c('wrap_pn', 'wrap_pn_no_reg', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
      system.time({
        all_p1_hot_start <- optimParallel::optimParallel(
          par=start_w_deoptim$optim$bestmem,
          fn = fn_opt,
          lower = rep(0, length(w_t)),
          upper = rep(1, length(w_t)),
          control = list(trace=0),
          parallel = list(cl=cl, forward=F, loginfo=F)  
        )
      })
      all_p1_hot_start$value_no_reg <- wrap_pn_no_reg(all_p1_hot_start$par)
      return(all_p1_hot_start)
    }
  )
  rlist::list.save(x = all_pns, file = 'results/plots/causality/poc/all_pn.RData')
  tryCatch({parallel::stopCluster(cl)}) 
}

# PS
{
  # start values
  penalty_weights <- 0
  poc <- proba_sufficient_causation
  wrap_pn_no_reg <-  xvine::wrapper_pn_all(
    data = rit_exp_stack,
    col_source = conditional_on_col, 
    u0_target = q_target,
    u0_source = q_source,
    poc = poc,
    lambda = NULL, p=1
  )
  w_t <- rep(.5, (ncol(pollution_data) ) * (k_markov - 1))
  w_t <- w_t / sum(w_t)
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c('wrap_pn_no_reg', 'wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
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
  tryCatch({parallel::stopCluster(cl)})
  wrap_pn_no_reg(start_w_deoptim$optim$bestmem)
  print(sum(start_w_deoptim$optim$bestmem))
  
  penalty_weights <- 0
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  all_pns <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all(
        data = rit_exp_stack,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=1
      )
      w_t <- rep(.5, (ncol(pollution_data)) * (k_markov - 1))
      w_t <- w_t / sum(w_t)
      fn_opt <- function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights}
      parallel::clusterExport(cl, c('wrap_pn', 'wrap_pn_no_reg', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
      system.time({
        all_p1_hot_start <- optimParallel::optimParallel(
          par=start_w_deoptim$optim$bestmem,
          fn = fn_opt,
          lower = rep(0, length(w_t)),
          upper = rep(1, length(w_t)),
          control = list(trace=0),
          parallel = list(cl=cl, forward=F, loginfo=F)  
        )
      })
      all_p1_hot_start$value_no_reg <- wrap_pn_no_reg(all_p1_hot_start$par)
      return(all_p1_hot_start)
    }
  )
  rlist::list.save(x = all_pns, file = 'results/plots/causality/poc/all_ps.RData')
  tryCatch({parallel::stopCluster(cl)}) 
}





aic_markov_level <- 10
svine_fit_aic <- xvine::fit_markov_svines(
  data = pollution_data_it,
  k.markov = aic_markov_level,
  family_set="archimedean",
  selcrit="aic",
  tree_crit="tau",
  threshold=0.05
)

# svine_fit_mbicv <- xvine::fit_markov_svines(
#   data = pollution_data_it,
#   k.markov = k_markov,
#   family_set="archimedean",
#   selcrit="mbicv",
#   tree_crit="tau",
#   threshold=0.05
# )

cond_sims_svines_aic <- lapply(
  1:dim(it_stack)[3],
  function(i){
    cond_sims_svines <- svines::svinecop_sim(
      n = k_markov-1, rep = 1, model = svine_fit_aic$copula, 
      past = rbind(matrix(runif((aic_markov_level-1)*ncol(pollution_data)), nrow=aic_markov_level-1, ncol=ncol(pollution_data)), it_stack[1,,i]),
      cores = 1)
    return(rbind(it_stack[1,,i], cond_sims_svines))
  }
)
cond_sims_svines_aic <- abind::abind(cond_sims_svines_aic, along = 3)
rlist::list.save('results/plots/causality/poc/cond_sims_svines_aic.RData')
cond_sims_svines_aic_exp <- xvine::apply_reverse_exponential(cond_sims_svines_aic)

# dual_reverse <- xvine::apply_dual_reverse_integral_transform(
#   data_unif = cond_sims_svines_aic,
#   data_source = pollution_data,
#   u0s = u0s_fixed,
#   conditional_on = 7,
#   shapes_f = f_cf_params_matrix[,1],
#   scales_f = f_cf_params_matrix[,3],
#   shapes_cf = f_cf_params_matrix[,6],
#   scales_cf = f_cf_params_matrix[,8]
# )
# 

cond_sims_svines_mbicv <- lapply(
  1:dim(it_stack)[3],
  function(i){
    cond_sims_svines <- svines::svinecop_sim(
      n = k_markov-1, rep = 1, model = svine_fit_mbicv$copula, 
      past = rbind(matrix(runif((k_markov-1)*ncol(pollution_data)), nrow=k_markov-1, ncol=ncol(pollution_data)), it_stack[1,,i]),
      cores = 1)
    return(rbind(it_stack[1,,i], cond_sims_svines))
  }
)
cond_sims_svines_mbicv <- abind::abind(cond_sims_svines_mbicv, along = 3)
cond_sims_svines_mbicv_exp <- xvine::apply_reverse_exponential(cond_sims_svines_mbicv)

tmp <- rlist::list.load('results/plots/causality/poc/cond_sims_svines_aic_exp.RData')
data_to_use <- tmp$data$data

q_target <- qexp(u0_single)
q_source <- qexp(u0_single)

# PNS
poc <- proba_necessary_sufficient_causation
{
  # start values
  penalty_weights <- 0
  poc <- proba_necessary_sufficient_causation
  wrap_pn_no_reg <-  xvine::wrapper_pn_all_but_source(
    data = data_to_use,
    col_source = conditional_on_col, 
    u0_target = q_target,
    u0_source = q_source,
    poc = poc,
    lambda = NULL, p=1
  )
  w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
  w_t <- w_t / sum(w_t)
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c('wrap_pn_no_reg', 'wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
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
  tryCatch(parallel::stopCluster(cl), error={}) 
  wrap_pn_no_reg(start_w_deoptim$optim$bestmem)
  print(sum(start_w_deoptim$optim$bestmem))
  
  #p1 p2
  penalty_weights <- 0
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  aws_pns_p1_svine  <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all_but_source(
        data = data_to_use,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=1
      )
      w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
      w_t <- w_t / sum(w_t)
      fn_opt <- function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights}
      parallel::clusterExport(cl, c('wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
      system.time({
        all_p1_hot_start <- optimParallel::optimParallel(
          par=start_w_deoptim$optim$bestmem,
          fn = fn_opt,
          lower = rep(0, length(w_t)),
          upper = rep(1, length(w_t)),
          control = list(trace=0),
          parallel = list(cl=cl, forward=F, loginfo=F)  
        )
      })
      all_p1_hot_start$value_no_reg <- wrap_pn_no_reg(all_p1_hot_start$par)
      return(all_p1_hot_start)
    }
  )
  rlist::list.save(x = aws_pns_p1_svine, file = 'results/plots/causality/poc/aws_pns_p1_svine_aic.RData')
  tryCatch(parallel::stopCluster(cl), error={}) 
  
  cl <- parallel::makeCluster(parallel::detectCores()-2)
  aws_pns_p2_svine <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all_but_source(
        data = data_to_use,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=2
      )
      w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
      parallel::clusterExport(cl, c('wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
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
  rlist::list.save(x = aws_pns_p2_svine, file = 'results/plots/causality/poc/aws_pns_p2_svine_aic.RData')
  
  parallel::stopCluster(cl)
}

# PN
poc <- proba_necessary_causation
{
  # start values
  penalty_weights <- 0
  poc <- proba_necessary_sufficient_causation
  wrap_pn_no_reg <-  xvine::wrapper_pn_all_but_source(
    data = data_to_use,
    col_source = conditional_on_col, 
    u0_target = q_target,
    u0_source = q_source,
    poc = poc,
    lambda = NULL, p=1
  )
  w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
  w_t <- w_t / sum(w_t)
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c('wrap_pn_no_reg', 'wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
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
  tryCatch(parallel::stopCluster(cl), error={}) 
  wrap_pn_no_reg(start_w_deoptim$optim$bestmem)
  print(sum(start_w_deoptim$optim$bestmem))
  
  #p1 p2
  penalty_weights <- 0
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  aws_pn_p1_svine  <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all_but_source(
        data = data_to_use,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=1
      )
      w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
      w_t <- w_t / sum(w_t)
      fn_opt <- function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights}
      parallel::clusterExport(cl, c('wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
      system.time({
        all_p1_hot_start <- optimParallel::optimParallel(
          par=start_w_deoptim$optim$bestmem,
          fn = fn_opt,
          lower = rep(0, length(w_t)),
          upper = rep(1, length(w_t)),
          control = list(trace=0),
          parallel = list(cl=cl, forward=F, loginfo=F)  
        )
      })
      all_p1_hot_start$value_no_reg <- wrap_pn_no_reg(all_p1_hot_start$par)
      return(all_p1_hot_start)
    }
  )
  rlist::list.save(x = aws_pn_p1_svine, file = 'results/plots/causality/poc/aws_pn_p1_svine_aic.RData')
  tryCatch(parallel::stopCluster(cl), error={}) 
  
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  aws_pn_p2_svine <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all_but_source(
        data = data_to_use,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=2
      )
      w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
      parallel::clusterExport(cl, c('wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
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
  rlist::list.save(x = aws_pn_p2_svine, file = 'results/plots/causality/poc/aws_pn_p2_svine_aic.RData')
  
  tryCatch(parallel::stopCluster(cl), error={}) 
}

# PS
poc <- proba_sufficient_causation
{
  # start values
  penalty_weights <- 0
  poc <- proba_necessary_sufficient_causation
  wrap_pn_no_reg <-  xvine::wrapper_pn_all_but_source(
    data = data_to_use,
    col_source = conditional_on_col, 
    u0_target = q_target,
    u0_source = q_source,
    poc = poc,
    lambda = NULL, p=1
  )
  w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
  w_t <- w_t / sum(w_t)
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c('wrap_pn_no_reg', 'wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
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
  tryCatch({parallel::stopCluster(cl)}) 
  wrap_pn_no_reg(start_w_deoptim$optim$bestmem)
  print(sum(start_w_deoptim$optim$bestmem))
  
  #p1 p2
  penalty_weights <- 0
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  aws_ps_p1_svine  <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all_but_source(
        data = data_to_use,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=1
      )
      w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
      w_t <- w_t / sum(w_t)
      fn_opt <- function(w){-wrap_pn(w) + (sum(w)-1)^2*penalty_weights}
      parallel::clusterExport(cl, c('wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
      system.time({
        all_p1_hot_start <- optimParallel::optimParallel(
          par=start_w_deoptim$optim$bestmem,
          fn = fn_opt,
          lower = rep(0, length(w_t)),
          upper = rep(1, length(w_t)),
          control = list(trace=0),
          parallel = list(cl=cl, forward=F, loginfo=F)  
        )
      })
      all_p1_hot_start$value_no_reg <- wrap_pn_no_reg(all_p1_hot_start$par)
      return(all_p1_hot_start)
    }
  )
  rlist::list.save(x = aws_ps_p1_svine, file = 'results/plots/causality/poc/aws_ps_p1_svine_aic.RData')
  tryCatch({parallel::stopCluster(cl)}) 
  
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  aws_ps_p2_svine <- lapply(
    lmbds,
    function(l){
      print(l)
      wrap_pn <-  xvine::wrapper_pn_all_but_source(
        data = data_to_use,
        col_source = conditional_on_col, 
        poc = poc,
        u0_target = q_target,
        u0_source = q_source,
        lambda = l, p=2
      )
      w_t <- rep(.5, (ncol(pollution_data) - 1) * (k_markov - 1))
      parallel::clusterExport(cl, c('wrap_pn', 'rit_exp_stack', 'conditional_on_col', 'q_target', 'q_source', 'penalty_weights'))
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
  rlist::list.save(x = aws_ps_p2_svine, file = 'results/plots/causality/poc/aws_ps_p2_svine_aic.RData')
  
  tryCatch({parallel::stopCluster(cl)}) 
}

# LOAD sparsity results

{
  a <- rlist::list.load('results/plots/causality/poc/aws_pns_p1.RData')
  b <- rlist::list.load('results/plots/causality/poc/aws_pns_p2.RData')
  a1 <- rlist::list.load('results/plots/causality/poc/aws_pns_p1_null.RData')
  b1 <- rlist::list.load('results/plots/causality/poc/aws_pns_p2_null.RData')
  
  aws_pns_p1 <- c(a, a1)
  aws_pns_p2 <- c(b, b1)
  
  a <- rlist::list.load('results/plots/causality/poc/aws_pn_p1.RData')
  b <- rlist::list.load('results/plots/causality/poc/aws_pn_p2.RData')
  a1 <- rlist::list.load('results/plots/causality/poc/aws_pn_p1_null.RData')
  b1 <- rlist::list.load('results/plots/causality/poc/aws_pn_p2_null.RData')
  
  aws_pn_p1 <- c(a, a1)
  aws_pn_p2 <- c(b, b1)
  
  a <- rlist::list.load('results/plots/causality/poc/aws_ps_p1.RData')
  b <- rlist::list.load('results/plots/causality/poc/aws_ps_p2.RData')
  a1 <- rlist::list.load('results/plots/causality/poc/aws_ps_p1_null.RData')
  b1 <- rlist::list.load('results/plots/causality/poc/aws_ps_p2_null.RData')
  
  aws_ps_p1 <- c(a, a1)
  aws_ps_p2 <- c(b, b1)
  
  a <- rlist::list.load('results/plots/causality/poc/aws_pns_p1_svine_aic.RData')
  b <- rlist::list.load('results/plots/causality/poc/aws_pns_p2_svine_aic.RData')
  a1 <- rlist::list.load('results/plots/causality/poc/aws_pns_p1_svine_aic_null.RData')
  b1 <- rlist::list.load('results/plots/causality/poc/aws_pns_p2_svine_aic_null.RData')
  
  aws_pns_p1_svine_aic <- c(a, a1)
  aws_pns_p2_svine_aic <- c(b, b1)
  
  a <- rlist::list.load('results/plots/causality/poc/aws_pn_p1_svine_aic.RData')
  b <- rlist::list.load('results/plots/causality/poc/aws_pn_p2_svine_aic.RData')
  a1 <- rlist::list.load('results/plots/causality/poc/aws_pn_p1_svine_aic_null.RData')
  b1 <- rlist::list.load('results/plots/causality/poc/aws_pn_p2_svine_aic_null.RData')
  
  aws_pn_p1_svine_aic <- c(a, a1)
  aws_pn_p2_svine_aic <- c(b, b1)
  
  aws_ps_p1_svine_aic <- rlist::list.load('results/plots/causality/poc/aws_ps_p1_svine_aic.RData')
  aws_ps_p2_svine_aic <- rlist::list.load('results/plots/causality/poc/aws_ps_p2_svine_aic.RData')
}



# all
# PNS all
q_target <- qexp(u0_single)
q_source <- qexp(u0_single)
system.time({
  r_time_w_empirical_pns_all <- xvine::maximise_pns(
    data = rit_exp_stack,
    col_source = conditional_on_col, 
    u0_target = q_target,
    u0_source = q_source,
    type = 'all', routine = 'optimParallel', lambda = 1
  )$weights
})

# all but source
# PNS all but source
q_target <- qexp(u0_single)
q_source <- qexp(u0_single)
system.time({
  r_time_w_empirical_pns_all_but_source <- xvine::maximise_pns(
    data = rit_exp_stack,
    col_source = conditional_on_col, 
    u0_target = q_target,
    u0_source = q_source,
    type = 'all-but-source', routine = 'optim', lambda = 1
  )$weights
})



colnames(r_time_w_empirical_pn) <- colnames(pollution_data)
rownames(r_time_w_empirical_pn) <- paste('t +', seq_len(nrow(r_time_w_empirical_pn)))
colnames(r_time_w_empirical_ps) <- colnames(pollution_data)
rownames(r_time_w_empirical_ps) <- paste('t +', seq_len(nrow(r_time_w_empirical_ps)))
colnames(r_time_w_empirical_pns) <- colnames(pollution_data)
rownames(r_time_w_empirical_pns) <- paste('t +', seq_len(nrow(r_time_w_empirical_pns)))

par(mfrow=c(2, 3))
corrplot::corrplot(
  r_time_w_empirical_pn,
  cl.lim = c(0,1), outline = "black", shade.lwd=20, 
  method = "shade", tl.cex=1.5, tl.col='black', tl.srt = 45
)
title(main='PN', cex.main=2)
corrplot::corrplot(
  r_time_w_empirical_ps,
  cl.lim = c(0,1), outline = "black", shade.lwd=20, 
  method = "shade", tl.cex=1.5, tl.col='black', tl.srt = 45
)
title(main='PS', cex.main=2)
corrplot::corrplot(
  r_time_w_empirical_pns,
  cl.lim = c(0,1), outline = "black", shade.lwd=20, 
  method = "shade", tl.cex=1.5, tl.col='black', tl.srt = 45
)
title(main='PNS', cex.main=2)

corrplot::corrplot(
  r_time_w_empirical_pn %*% diag(1/apply(r_time_w_empirical_pn, 2, max)),
  cl.lim = c(0,1), outline = "black", shade.lwd=20, 
  method = "shade", tl.cex=1.5, tl.col='black', tl.srt = 45
)
corrplot::corrplot(
  r_time_w_empirical_ps %*% diag(1/apply(r_time_w_empirical_ps, 2, max)),
  cl.lim = c(0,1), outline = "black", shade.lwd=20, 
  method = "shade", tl.cex=1.5, tl.col='black', tl.srt = 45
)
corrplot::corrplot(
  r_time_w_empirical_pns %*% diag(1/apply(r_time_w_empirical_pns, 2, max)),
  cl.lim = c(0,1), outline = "black", shade.lwd=20, 
  method = "shade", tl.cex=1.5, tl.col='black', tl.srt = 45
)






# r_time_w_empirical <- read.csv('results/r_time_w_empirical.csv')[,-1]
r_time_w_empirical <- as.matrix(r_time_w_empirical)
colnames(r_time_w_empirical) <- colnames(pollution_data)
rownames(r_time_w_empirical) <- paste('t +', seq_len(nrow(r_time_w_empirical)))
r_time_w_empirical

par(mfrow=c(1, 2))
corrplot::corrplot(
  r_time_w_empirical,
  cl.lim = c(0,1), outline = "black", shade.lwd=20, 
  method = "shade", tl.cex=1.5, tl.col='black', tl.srt = 45
)
corrplot::corrplot(
  t(t(r_time_w_empirical) / apply(r_time_w_empirical, 2, max)), 
  cl.lim = c(0,1),
  method = "shade", tl.cex=1.5, tl.col='black', tl.srt = 45
)

pairwise_idxs <- expand.grid(conditional_on_col, seq_len(ncol(pollution_data)))
system.time({
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c('pollution_data', 'u0_single', 'rit_stack'))
  r_time_w_empirical <- parallel::parApply(
    cl = cl, X=pairwise_idxs, MARGIN = 1, 
    FUN = function(x){
      q_target <- quantile(pollution_data[,x[2]], u0_single)
      q_source <- quantile(pollution_data[,x[1]], u0_single)
      w_t <- xvine::maximise_pn(
        data = rit_stack,
        target = x[2],
        col_source = x[1], 
        u0_target = q_target,
        u0_source = q_source,
        type = 'time', routine = 'deoptim'
      )
      return(w_t$weights)
    }
  )
  parallel::stopCluster(cl)
})
r_time_w_empirical

par(mar=c(4,4,2,1), mfrow=c(1,2))
matplot(r_time_w_empirical_o3, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
         main='Conditional on O3', ylab='Weights', xlab='Offset', cex.lab=1.3, cex=2, cex.axis=1.3)
matplot(r_time_w_empirical, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        main='Conditional on NO', ylab='Weights', xlab='Offset', cex.lab=1.3, cex=2, cex.axis=1.3)
legend(5, .5, colnames(pollution_data), plotting_colours)

colnames(r_time_w_empirical_o3) <- colnames(pollution_data)
write.csv(r_time_w_empirical_o3, file='results/plots/causality/r_time_w_empirical_o3.csv')

colnames(r_time_w_empirical) <- colnames(pollution_data)
write.csv(r_time_w_empirical, file='results/plots/causality/r_time_w_empirical_no.csv')

for(idx_target in 1:ncol(pollution_data)){
    q_target <- quantile(pollution_data[,idx_target], u0_single)
    q_source <- quantile(pollution_data[,idx_source], u0_single)
    w_t <- xvine::maximise_pn(
      data = rit_stack,
      target = idx_target,
      col_source = idx_source, 
      u0_target = q_target,
      u0_source = q_source,
      type = 'time', routine = 'deoptim'
    )
    cat(paste('\nSource:', colnames(pollution_data)[[idx_source]], '; Target', colnames(pollution_data)[[idx_target]], '\n'))
    cat(paste('Max weight:', round(max(w_t$weights) * 100, 2), '% (', which.max(w_t$weights) ,')\n'))
    cat(paste('Min weight:', round(min(w_t$weights) * 100, 2), '% (', which.min(w_t$weights) ,')\n'))
    cat(paste('Relative weight difference:', (max(w_t$weights) - min(w_t$weights))/min(w_t$weights) * 100, '%\n'))
    
    logits_w <- xvine::logits(w_t$weights)
    if(idx_target==1){
      plot(
        1:(k_markov-1), 
        logits_w,
        ylab = 'Weights (logits)',
        xlab = 'Timestep', main = paste('After an extreme in', colnames(pollution_data)[[idx_source]]),
        ylim=c(-2.5, -1),
        type='b', col=plotting_colours[idx_target],
        lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4, cex.main=1.6
      )
    }else{
      points(
        1:(k_markov-1), 
        logits_w,
        type='b',
        col=plotting_colours[idx_target],
        lwd=3
      )
    }
  legend(5, -1, colnames(pollution_data), plotting_colours)
}

# cross !
par(mfrow=c(2,3))
for(idx_source in 1:ncol(pollution_data)){
  for(time_target in 1:(k_markov-1)){
    w_t <- xvine::maximise_pn(
      data = it_stack,
      target = time_target,
      col_source = idx_source,
      u0_target = u0_single,
      u0_source = u0_single,
      type = 'cross'
    )
    print(w_t)
    logits_w <- xvine::logits(w_t$weights)
    if(time_target==1){
      barplot(
        logits_w,
        names.arg = colnames(pollution_data),
        ylab = 'Weights (logits)',
        xlab = 'Pollutant', main = paste('After an extreme in', colnames(pollution_data)[[idx_source]]),
        ylim=c(-6, .1),
        col=plotting_colours[idx_target],
        lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4, cex.main=1.6
      )
    }else{
      points(
        colnames(pollution_data), 
        logits_w,
        type='b',
        col=plotting_colours[idx_target],
        lwd=3
      )
    }
  }
  legend(1, -3.5, colnames(pollution_data), plotting_colours)
}


 # SIMULATION STUDY

# vine hpararms
hparams_vines <- list(
  family_set="archimedean",
  selcrit="aic",
  tree_crit="rho",
  trunc_lvl=NA,
  threshold=0.05
)
#  trunc_lvl=NA means trunc_level determined by mBICV


# fit R-Vines
rvine_fit <- xvine::fit_markov_rvines(
  data = pollution_data_it,
  k.markov = k_markov,
  hparams_vines
)

timewise_rvine_fit <- xvine::fit_markov_timewise_rvines(
  data = pollution_data_it,
  k.markov = k_markov,
  hparams_vines
)

cond_rvine_fit <- xvine::fit_markov_conditional_rvines(
  data = pollution_data_it,
  k.markov = k_markov,
  col_source = 7,
  u0 = u0_single,
  hparams_vines
)

cond_timewise_rvine_fit <- xvine::fit_markov_conditional_timewise_rvines(
  data = pollution_data_it,
  k.markov = k_markov,
  col_source = 7,
  u0 = u0_single,
  hparams_vines
)

svine_fit <- xvine::fit_markov_svines(
  data = pollution_data_it,
  k.markov = k_markov,
  family_set="archimedean",
  selcrit="aic",
  tree_crit="rho",
  threshold=0.05
)

plot(rvine_fit$copula, 1, edge_labels = "family")
plot(rvine_fit$copula, 26, edge_labels = "family")

# sims
n_sims <- 20000

rvine_sims <- model_simulation(n=n_sims, model=rvine_fit, qrng=F)
rvine_sims_reverse <- xvine::apply_reverse_integral_transform(
  data_unif = rvine_sims,
  data_source = pollution_data,
  u0s = u0s_fixed,
  shapes = params_gpd_ests[,1],
  scales = params_gpd_ests[,2]
)
rvine_sims_reverse_exp <- xvine::apply_reverse_exponential(
  data_unif = rvine_sims
)

timewise_rvine_sims <- model_simulation(n=n_sims, model=timewise_rvine_fit, qrng=F)
timewise_rvine_sims_reverse <- xvine::apply_reverse_integral_transform(
  data_unif = timewise_rvine_sims,
  data_source = pollution_data,
  u0s = u0s_fixed,
  shapes = params_gpd_ests[,1],
  scales = params_gpd_ests[,2]
)

cond_no_rvine_sims <- model_simulation(n=n_sims, model=cond_rvine_fit, qrng=F)
cond_no_rvine_sims_reversed <- xvine::apply_reverse_integral_transform(
  data_unif = cond_no_rvine_sims,
  data_source = pollution_data,
  u0s = u0s_fixed,
  shapes = params_gpd_ests[,1],
  scales = params_gpd_ests[,2]
)

cond_timewise_no_rvine_sims <- model_simulation(n=n_sims, model=cond_timewise_rvine_fit, qrng=F)
cond_timewise_no_rvine_sims_reversed <- xvine::apply_reverse_integral_transform(
  data_unif = cond_timewise_no_rvine_sims,
  data_source = pollution_data,
  u0s = u0s_fixed,
  shapes = params_gpd_ests[,1],
  scales = params_gpd_ests[,2]
)

svine_sims <- model_simulation(n=n_sims, model=svine_fit, qrng=F)
svine_sims <- svines::svine_sim(n = k_markov, rep = n_sims, model = svine_fit, cores = 1, qrng=F)
svine_sims_reverse <- xvine::apply_reverse_integral_transform(
  data_unif = svine_sims,
  data_source = pollution_data,
  u0s = u0s_fixed,
  shapes = params_gpd_ests[,1],
  scales = params_gpd_ests[,2]
)
svine_sims_reverse_exp <- xvine::apply_reverse_exponential(
  data_unif = svine_sims
)


# plot single-target PN
# plotting_colours <- viridis::(ncol(pollution_data), alpha = 1.)
par(mfrow=c(2,4))
plotting_colours <- RColorBrewer::brewer.pal(n = ncol(pollution_data), name = 'Dark2')
for(idx_source in 1:ncol(pollution_data)){
  for(idx_target in 1:ncol(pollution_data)){
    rvine_pst <- xvine::proba_single_target(
      data=rvine_sims,
      col_target = idx_target, col_source = idx_source,
      u0_target = u0_single, u0_source = u0_single
    )
    if(idx_target==1){
      plot(
        1:(k_markov-1), 
        xvine::proba_necessary_causation(
          rvine_pst$factual, rvine_pst$counterfactual
        ),
        ylim=c(0,1), ylab = 'Necessary causation',
        xlab = 'Offset', main = paste('After an extreme in', colnames(pollution_data)[[idx_source]]),
        type='b', col=plotting_colours[idx_target],
        lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4
      )
    }else{
      points(
        1:(k_markov-1), 
        xvine::proba_necessary_causation(
          rvine_pst$factual, rvine_pst$counterfactual
        ), ylim=c(0,1), type='b',
        col=plotting_colours[idx_target],
        lwd=3
      )
    }
    legend(4.2, .4, colnames(pollution_data), plotting_colours)
  }
}

# compute and plot maximum time-target PN
plotting_colours <- RColorBrewer::brewer.pal(n = ncol(pollution_data), name = 'Dark2')
par(mfrow=c(2,4))
for(idx_source in 1:ncol(pollution_data)){
  for(idx_target in 1:ncol(pollution_data)){
    w_t <- xvine::maximise_pn(
      data = rvine_sims,
      target = idx_target,
      col_source = idx_source,
      u0_target = u0_single,
      u0_source = u0_single,
      type = 'time'
    )
    cat(paste('\nSource:', colnames(pollution_data)[[idx_source]], '; Target', colnames(pollution_data)[[idx_target]], '\n'))
    cat(paste('Max weight:', max(w_t$weights) * 100, '%\n'))
    cat(paste('Min weight:', min(w_t$weights) * 100, '%\n'))
    cat(paste('Relative weight difference:', (max(w_t$weights) - min(w_t$weights))/min(w_t$weights) * 100, '%\n'))
    
    logits_w <- xvine::logits(w_t$weights)
    if(idx_target==1){
      plot(
        1:(k_markov-1), 
        logits_w,
        ylab = 'Weights (logits)',
        xlab = 'Timestep', main = paste('After an extreme in', colnames(pollution_data)[[idx_source]]),
        ylim=c(-3, -1),
        type='b', col=plotting_colours[idx_target],
        lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4, cex.main=1.6
      )
    }else{
      points(
        1:(k_markov-1), 
        logits_w,
        type='b',
        col=plotting_colours[idx_target],
        lwd=3
      )
    }
  }
  legend(1, -1, colnames(pollution_data), plotting_colours)
}

# compute and plot maximum time-target PN
# reversed!
plotting_colours <- RColorBrewer::brewer.pal(n = ncol(pollution_data), name = 'Dark2')
par(mfrow=c(2,4))
for(idx_source in 1:ncol(pollution_data)){
  for(idx_target in 1:ncol(pollution_data)){
    q_target <- quantile(pollution_data[,idx_target], u0s_fixed[idx_target])
    q_source <- quantile(pollution_data[,idx_source], u0s_fixed[idx_source])
    w_t <- xvine::maximise_pn(
      data = rvine_sims_reverse$data,
      target = idx_target,
      col_source = idx_source,
      u0_target = q_target,
      u0_source = q_source,
      type = 'time'
    )
    cat(paste('\nSource:', colnames(pollution_data)[[idx_source]], '; Target', colnames(pollution_data)[[idx_target]], '\n'))
    cat(paste('Max weight:', max(w_t$weights) * 100, '%\n'))
    cat(paste('Min weight:', min(w_t$weights) * 100, '%\n'))
    cat(paste('Relative weight difference:', (max(w_t$weights) - min(w_t$weights))/min(w_t$weights) * 100, '%\n'))
    
    logits_w <- xvine::logits(w_t$weights)
    if(idx_target==1){
      plot(
        1:(k_markov-1), 
        logits_w,
        ylab = 'Weights (logits)',
        xlab = 'Timestep', main = paste('(R) After an extreme in', colnames(pollution_data)[[idx_source]]),
        ylim=c(-3, -1),
        type='b', col=plotting_colours[idx_target],
        lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4, cex.main=1.6
      )
    }else{
      points(
        1:(k_markov-1), 
        logits_w,
        type='b',
        col=plotting_colours[idx_target],
        lwd=3
      )
    }
  }
  legend(1, -1, colnames(pollution_data), plotting_colours)
}


# TIMEWISE 
# compute and plot maximum time-target PN
plotting_colours <- RColorBrewer::brewer.pal(n = ncol(pollution_data), name = 'Dark2')
par(mfrow=c(2,4))
for(idx_source in 1:ncol(pollution_data)){
  for(idx_target in 1:ncol(pollution_data)){
    w_t <- xvine::maximise_pn(
      data = timewise_rvine_sims,
      target = idx_target,
      col_source = idx_source,
      u0_target = u0_single,
      u0_source = u0_single,
      type = 'time'
    )
    cat(paste('\nSource:', colnames(pollution_data)[[idx_source]], '; Target', colnames(pollution_data)[[idx_target]], '\n'))
    cat(paste('Max weight:', max(w_t$weights) * 100, '%\n'))
    cat(paste('Min weight:', min(w_t$weights) * 100, '%\n'))
    cat(paste('Relative weight difference:', (max(w_t$weights) - min(w_t$weights))/min(w_t$weights) * 100, '%\n'))
    
    logits_w <- xvine::logits(w_t$weights)
    if(idx_target==1){
      plot(
        1:(k_markov-1), 
        logits_w,
        ylab = 'Weights (logits)',
        xlab = 'Timestep', main = paste('After an extreme in', colnames(pollution_data)[[idx_source]]),
        ylim=c(-3, -1),
        type='b', col=plotting_colours[idx_target],
        lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4, cex.main=1.6
      )
    }else{
      points(
        1:(k_markov-1), 
        logits_w,
        type='b',
        col=plotting_colours[idx_target],
        lwd=3
      )
    }
  }
  legend(1, -1, colnames(pollution_data), plotting_colours)
}

# compute and plot maximum time-target PN
# reversed!
plotting_colours <- RColorBrewer::brewer.pal(n = ncol(pollution_data), name = 'Dark2')
par(mfrow=c(2,3))
for(idx_source in 1:ncol(pollution_data)){
  for(idx_target in 1:ncol(pollution_data)){
    q_target <- quantile(pollution_data[,idx_target], u0s_fixed[idx_target])
    q_source <- quantile(pollution_data[,idx_source], u0s_fixed[idx_source])
    w_t <- xvine::maximise_pn(
      data = timewise_rvine_sims_reverse$data,
      target = idx_target,
      col_source = idx_source,
      u0_target = q_target,
      u0_source = q_source,
      type = 'time'
    )
    cat(paste('\nSource:', colnames(pollution_data)[[idx_source]], '; Target', colnames(pollution_data)[[idx_target]], '\n'))
    cat(paste('Max weight:', max(w_t$weights) * 100, '%\n'))
    cat(paste('Min weight:', min(w_t$weights) * 100, '%\n'))
    cat(paste('Relative weight difference:', (max(w_t$weights) - min(w_t$weights))/min(w_t$weights) * 100, '%\n'))
    
    logits_w <- xvine::logits(w_t$weights)
    if(idx_target==1){
      plot(
        1:(k_markov-1), 
        logits_w,
        ylab = 'Weights (logits)',
        xlab = 'Timestep', main = paste('(R) After an extreme in', colnames(pollution_data)[[idx_source]]),
        ylim=c(-3, -1),
        type='b', col=plotting_colours[idx_target],
        lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4, cex.main=1.6
      )
    }else{
      points(
        1:(k_markov-1), 
        logits_w,
        type='b',
        col=plotting_colours[idx_target],
        lwd=3
      )
    }
  }
  legend(1, -1, colnames(pollution_data), plotting_colours)
}


# CONDITIONAL ON v/hr
plotting_colours <- RColorBrewer::brewer.pal(n = ncol(pollution_data), name = 'Dark2')
# cond single
# factual
par(mar=c(4,5,2,1), mfrow=c(2,3))
for(idx_source in 1:ncol(pollution_data)){
  for(idx_target in 1:ncol(pollution_data)){
    rvine_pst <- xvine::proba_single_target(
      data=cond_no_rvine_sims,
      col_target = idx_target, col_source = idx_source,
      u0_target = u0_single, u0_source = u0_single
    )
    if(idx_target==1){
      plot(
        1:(k_markov-1), 
        rvine_pst$factual,
        ylim=c(0,.8), ylab = 'Factual probability',
        xlab = 'Offset', main = paste('After an extreme in', colnames(pollution_data)[[idx_source]]),
        type='b', col=plotting_colours[idx_target],
        lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4
      )
    }else{
      points(
        1:(k_markov-1), 
        rvine_pst$factual, ylim=c(0,1), type='b',
        col=plotting_colours[idx_target],
        lwd=3
      )
    }
    abline(h=1-u0s_fixed[idx_target], lty=2, lwd=4)
    legend(8, .8, colnames(pollution_data), plotting_colours)
  }
}

# counterfactual
par(mar=c(4,5,2,1), mfrow=c(2,3))
for(idx_source in 1:ncol(pollution_data)){
  for(idx_target in 1:ncol(pollution_data)){
    rvine_pst <- xvine::proba_single_target(
      data=cond_no_rvine_sims,
      col_target = idx_target, col_source = idx_source,
      u0_target = u0_single, u0_source = u0_single
    )
    if(idx_target==1){
      plot(
        1:(k_markov-1), 
        rvine_pst$counterfactual,
        ylim=c(0,.4), ylab = 'Counterfactual probability',
        xlab = 'Offset', main = paste('After an extreme in', colnames(pollution_data)[[idx_source]]),
        type='b', col=plotting_colours[idx_target],
        lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4
      )
    }else{
      points(
        1:(k_markov-1), 
        rvine_pst$counterfactual, ylim=c(0,1), type='b',
        col=plotting_colours[idx_target],
        lwd=3
      )
    }
    abline(h=1-u0s_fixed[idx_target], lty=2, lwd=4)
    legend(8, 1., colnames(pollution_data), plotting_colours)
  }
}

# cond timewise
# factual
par(mfrow=c(2,3))
for(idx_source in 1:ncol(pollution_data)){
  for(idx_target in 1:ncol(pollution_data)){
    rvine_pst <- xvine::proba_single_target(
      data=cond_timewise_no_rvine_sims,
      col_target = idx_target, col_source = idx_source,
      u0_target = u0_single, u0_source = u0_single
    )
    if(idx_target==1){
      plot(
        1:(k_markov-1), 
        rvine_pst$factual,
        ylim=c(0,1), ylab = 'Factual probability',
        xlab = 'Offset', main = paste('After an extreme in', colnames(pollution_data)[[idx_source]]),
        type='b', col=plotting_colours[idx_target],
        lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4
      )
    }else{
      points(
        1:(k_markov-1), 
        rvine_pst$factual, ylim=c(0,1), type='b',
        col=plotting_colours[idx_target],
        lwd=3
      )
    }
    abline(h=1-u0s_fixed[idx_target], lty=2, lwd=4)
    legend(8, 1., colnames(pollution_data), plotting_colours)
  }
}

# counterfactual
par(mfrow=c(2,3))
for(idx_source in 1:ncol(pollution_data)){
  for(idx_target in 1:ncol(pollution_data)){
    rvine_pst <- xvine::proba_single_target(
      data=cond_timewise_no_rvine_sims,
      col_target = idx_target, col_source = idx_source,
      u0_target = u0_single, u0_source = u0_single
    )
    if(idx_target==1){
      plot(
        1:(k_markov-1), 
        rvine_pst$counterfactual,
        ylim=c(0,1), ylab = 'Counterfactual probability',
        xlab = 'Offset', main = paste('After an extreme in', colnames(pollution_data)[[idx_source]]),
        type='b', col=plotting_colours[idx_target],
        lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4
      )
    }else{
      points(
        1:(k_markov-1), 
        rvine_pst$counterfactual, ylim=c(0,1), type='b',
        col=plotting_colours[idx_target],
        lwd=3
      )
    }
    abline(h=1-u0s_fixed[idx_target], lty=2, lwd=4)
    legend(8, 1., colnames(pollution_data), plotting_colours)
  }
}


# Extremal path

#test
{
  par(mfrow=c(1,1))
cond_no_extremal_path <- cond_no_rvine_sims[,,which(cond_no_rvine_sims[1,3,] > u0_single)]
matplot(
  x = 0:(k_markov-1),
  apply(cond_no_extremal_path, c(1,2), mean) + apply(cond_no_extremal_path, c(1,2), function(x) 1.96 * sd(x) / sqrt(length(x))),
  col = plotting_colours, type = 'o', lwd = 2
)
matplot(
  x = 0:(k_markov-1),
  apply(cond_no_extremal_path, c(1,2), mean) - apply(cond_no_extremal_path, c(1,2), function(x) 1.96 * sd(x) / sqrt(length(x))),
  col = plotting_colours, type = 'o', lwd = 2, add = T
)
}

par(mfrow=c(2,2))
cond_no_extremal_path <- cond_no_rvine_sims[,,which(cond_no_rvine_sims[1,conditional_on_col,] > u0_single)]
matplot(
  x = 0:(k_markov-1),
  apply(cond_no_extremal_path, c(1,2), mean), main='Average response',
  col = plotting_colours, type = 'o', lty=1,
  ylab='Quantile', xlab='Time', lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4, cex.main=1.6
)
abline(h=u0_single, lwd=2)
matplot(
  x = 0:(k_markov-1),
  cond_no_extremal_path[,,1],
  col = plotting_colours, type = 'o', main='Example 1', lty=1,
  ylab='Quantile', xlab='Time', lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4, cex.main=1.6
)
abline(h=u0_single, lwd=2)
matplot(
  x = 0:(k_markov-1),
  cond_no_extremal_path[,,2],
  col = plotting_colours, type = 'o', main='Example 2', lty=1,
  ylab='Quantile', xlab='Time', lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4, cex.main=1.6
)
abline(h=u0_single, lwd=2)
matplot(
  x = 0:(k_markov-1),
  cond_no_extremal_path[,,3],
  col = plotting_colours, type = 'o', main='Example 3', lty=1,
  ylab='Quantile', xlab='Time', lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4, cex.main=1.6
)
abline(h=u0_single, lwd=2)


# extremal paths
it_stack_extremal_path <- it_stack[,,which(it_stack[1,conditional_on_col,] > u0_single)]
cond_no_extremal_path <- cond_no_rvine_sims[,,which(cond_no_rvine_sims[1,conditional_on_col,] > u0_single)]
cond_timewise_no_extremal_path <- cond_timewise_no_rvine_sims[,,which(cond_timewise_no_rvine_sims[1,conditional_on_col,] > u0_single)]
rvine_extremal_path <- rvine_sims[,,which(rvine_sims[1,conditional_on_col,] > u0_single)]
svine_extremal_path <- svine_sims[,,which(svine_sims[1,conditional_on_col,] > u0_single)]

par(mfrow=c(4,2))
lapply(seq_len(7), function(i) vioplot::vioplot(t(it_stack_extremal_path[,i,]), col=plotting_colours[i], main=colnames(pollution_data)[i]))
par(mfrow=c(4,2))
lapply(seq_len(7), function(i) vioplot::vioplot(t(rvine_extremal_path[,i,]), col=plotting_colours[i], main=colnames(pollution_data)[i]))
par(mfrow=c(4,2))
lapply(seq_len(7), function(i) vioplot::vioplot(t(cond_no_extremal_path[,i,]), col=plotting_colours[i], main=colnames(pollution_data)[i]))
par(mfrow=c(4,2))
lapply(seq_len(7), function(i) vioplot::vioplot(t(cond_timewise_no_extremal_path[,i,]), col=plotting_colours[i], main=colnames(pollution_data)[i]))


rit_stack_extremal_path <- xvine::apply_reverse_integral_transform(
  it_stack_extremal_path,
  data_source = pollution_data,
  u0s = u0s_fixed,
  shapes = params_gpd_ests[,1],
  scales = params_gpd_ests[,2]
)$data

r_rvine_extremal_path <- xvine::apply_reverse_integral_transform(
  rvine_extremal_path,
  data_source = pollution_data,
  u0s = u0s_fixed,
  shapes = params_gpd_ests[,1],
  scales = params_gpd_ests[,2]
)$data

r_cond_no_extremal_path <- xvine::apply_reverse_integral_transform(
  cond_no_extremal_path,
  data_source = pollution_data,
  u0s = u0s_fixed,
  shapes = params_gpd_ests[,1],
  scales = params_gpd_ests[,2]
)$data

r_cond_timewise_no_extremal_path <- xvine::apply_reverse_integral_transform(
  cond_timewise_no_extremal_path,
  data_source = pollution_data,
  u0s = u0s_fixed,
  shapes = params_gpd_ests[,1],
  scales = params_gpd_ests[,2]
)$data

# divergence
it_rvine_divergence <- xvine::apply_distr_distance(
  x = it_stack_extremal_path,
  y = rvine_extremal_path,
  method = 'divergence',
  per_bin=10
)

it_cond_divergence <- xvine::apply_distr_distance(
   x = it_stack_extremal_path,
   y = cond_no_extremal_path,
   method = 'divergence',
   per_bin=10
)

it_cond_timewise_divergence <- xvine::apply_distr_distance(
  x = it_stack_extremal_path,
  y = cond_timewise_no_extremal_path,
  method = 'divergence',
  per_bin=10
)

rit_rvine_divergence <- xvine::apply_distr_distance(
  x = rit_stack_extremal_path,
  y = r_rvine_extremal_path,
  method = 'divergence',
  per_bin=10
)

rit_cond_divergence <- xvine::apply_distr_distance(
  x = rit_stack_extremal_path,
  y = r_cond_no_extremal_path,
  method = 'divergence',
  per_bin=10
)

rit_cond_timewise_divergence <- xvine::apply_distr_distance(
  x = rit_stack_extremal_path,
  y = r_cond_timewise_no_extremal_path,
  method = 'divergence',
  per_bin=10
)


# KL
it_rvine_kl <- xvine::apply_distr_distance(
  x = it_stack_extremal_path,
  y = rvine_extremal_path,
  method = 'kullback-leibler',
  per_bin=10
)

it_cond_kl <- xvine::apply_distr_distance(
  x = it_stack_extremal_path,
  y = cond_no_extremal_path,
  method = 'kullback-leibler',
  per_bin=10
)

it_cond_timewise_kl <- xvine::apply_distr_distance(
  x = it_stack_extremal_path,
  y = cond_timewise_no_extremal_path,
  method = 'kullback-leibler',
  per_bin=10
)

rit_rvine_kl <- xvine::apply_distr_distance(
  x = rit_stack_extremal_path,
  y = r_rvine_extremal_path,
  method = 'kullback-leibler',
  per_bin=10
)

rit_cond_kl <- xvine::apply_distr_distance(
  x = rit_stack_extremal_path,
  y = r_cond_no_extremal_path,
  method = 'kullback-leibler',
  per_bin=10
)

rit_cond_timewise_kl <- xvine::apply_distr_distance(
  x = rit_stack_extremal_path,
  y = r_cond_timewise_no_extremal_path,
  method = 'kullback-leibler',
  per_bin=10
)

# DIVERGENCE / KL at uniform scale
layout(matrix(c(1:12, 13, 13, 13), nrow=5, ncol=3, byrow = T), width=c(1,1,1,1,2), heights = c(1,1,1,1,.4)) 
par(mar=c(2,4,2,1)) #No margin on the right side
matplot(it_rvine_divergence, type='o', pch=18, lty=1, lwd=2, col = plotting_colours, 
        log='y', main='R-Vine Divergence', ylab='Divergence (log)', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
par(mar=c(2,4,2,1)) #No margin on the right side
matplot(it_cond_divergence, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        log='y', main='Cond R-Vine Divergence', ylab='Divergence (log)', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
par(mar=c(2,4,2,1)) #No margin on the right side
matplot(it_cond_timewise_divergence, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        log='y', main='CT R-Vine Divergence', ylab='Divergence (log)', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)

par(mar=c(3,4,2,1)) #No margin on the right side
matplot(it_rvine_divergence/it_cond_divergence, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        main='R-Vine vs Cond R-Vine ratio KL', ylab='KL Ratio', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
par(mar=c(3,4,2,1)) #No margin on the right side
matplot(it_cond_divergence/it_cond_timewise_divergence, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        main='Cond vs CT R-Vine ratio KL', ylab='KL Ratio', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
par(mar=c(3,4,2,1)) #No margin on the right side
matplot(it_rvine_divergence/it_cond_timewise_divergence, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        main='R-Vine vs CT R-Vine ratio KL', ylab='KL Ratio', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)

par(mar=c(2,4,2,1)) #No margin on the right side
matplot(it_rvine_kl, type='o', pch=18, lty=1, lwd=2, col = plotting_colours, 
        log='y', main='R-Vine KL', ylab='KL (log)', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
par(mar=c(2,4,2,1)) #No margin on the right side
matplot(it_cond_kl, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        log='y', main='Cond R-Vine KL', ylab='KL (log)', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
par(mar=c(2,4,2,1)) #No margin on the right side
matplot(it_cond_timewise_kl, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        log='y', main='CT R-Vine KL', ylab='KL (log)', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
par(mar=c(3,4,2,1)) #No margin on the right side
matplot(it_rvine_kl/it_cond_kl, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        main='R-Vine vs Cond R-Vine ratio KL', ylab='KL Ratio', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
title(xlab = 'Horizon', line=2, cex.lab=1.3)
par(mar=c(3,4,2,1)) #No margin on the right side
matplot(it_cond_kl/it_cond_timewise_kl, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
         main='Cond vs CT R-Vine ratio KL', ylab='KL Ratio', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
title(xlab = 'Horizon', line=2, cex.lab=1.3)
par(mar=c(3,4,2,1)) #No margin on the right side
matplot(it_rvine_kl/it_cond_timewise_kl, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
         main='R-Vine vs CT R-Vine ratio KL', ylab='KL Ratio', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
title(xlab = 'Horizon', line=2, cex.lab=1.3)
par(mar=c(2,4,2,1)) #No margin on the right side
plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
legend("center", legend = colnames(pollution_data), col = plotting_colours, horiz = T, 
       cex = 1.5, pch=18, lty=1, lwd=2, bty = 'n', box.col = "white", pt.cex = 2)


# DIVERGENCE / KL at reversed scale
layout(matrix(c(1:12, 13, 13, 13), nrow=5, ncol=3, byrow = T), width=c(1,1,1,1,2), heights = c(1,1,1,1,.4)) 
par(mar=c(2,4,2,1)) #No margin on the right side
matplot(rit_rvine_divergence, type='o', pch=18, lty=1, lwd=2, col = plotting_colours, 
        log='y', main='R-Vine Divergence', ylab='Divergence (log)', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
par(mar=c(2,4,2,1)) #No margin on the right side
matplot(rit_cond_divergence, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        log='y', main='Cond R-Vine Divergence', ylab='Divergence (log)', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
par(mar=c(2,4,2,1)) #No margin on the right side
matplot(rit_cond_timewise_divergence, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        log='y', main='CT R-Vine Divergence', ylab='Divergence (log)', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)

par(mar=c(3,4,2,1)) #No margin on the right side
matplot(rit_rvine_divergence/rit_cond_divergence, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        main='R-Vine vs Cond R-Vine ratio Divergence', ylab='Divergence Ratio', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
par(mar=c(3,4,2,1)) #No margin on the right side
matplot(rit_cond_divergence/rit_cond_timewise_divergence, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        main='Cond vs CT R-Vine ratio Divergence', ylab='Divergence Ratio', xlab = 'Horizon', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
par(mar=c(3,4,2,1)) #No margin on the right side
matplot(rit_rvine_divergence/rit_cond_timewise_divergence, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        main='R-Vine vs CT R-Vine ratio Divergence', ylab='Divergence Ratio', xlab = 'Horizon', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)

par(mar=c(2,4,2,1)) #No margin on the right side
matplot(rit_rvine_kl, type='o', pch=18, lty=1, lwd=2, col = plotting_colours, 
        log='y', main='R-Vine KL', ylab='KL (log)', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
par(mar=c(2,4,2,1)) #No margin on the right side
matplot(rit_cond_kl, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        log='y', main='Cond R-Vine KL', ylab='KL (log)', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
par(mar=c(2,4,2,1)) #No margin on the right side
matplot(rit_cond_timewise_kl, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        log='y', main='CT R-Vine KL', ylab='KL (log)', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)

par(mar=c(2,4,2,1)) #No margin on the right side
matplot(rit_rvine_kl/rit_cond_kl, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        main='R-Vine vs Cond R-Vine ratio KL', ylab='KL Ratio', cex.lab=1.3, cex=2, cex.axis=1.3)
# title(xlab = 'Horizon', line=2, cex.lab=1.3)
abline(h=1, lwd=2)
par(mar=c(2,4,2,1)) #No margin on the right side
matplot(rit_cond_kl/rit_cond_timewise_kl, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        main='Cond vs CT R-Vine ratio KL', ylab='KL Ratio', xlab = 'Horizon', cex.lab=1.3, cex=2, cex.axis=1.3)
# title(xlab = 'Horizon', line=2, cex.lab=1.3)
abline(h=1, lwd=2)
par(mar=c(2,4,2,1)) #No margin on the right side
matplot(rit_rvine_kl/rit_cond_timewise_kl, type='o', pch=18, lty=1, lwd=2, col = plotting_colours,
        main='R-Vine vs CT R-Vine ratio KL', ylab='KL Ratio', xlab = 'Horizon', cex.lab=1.3, cex=2, cex.axis=1.3)
abline(h=1, lwd=2)
# title(xlab = 'Horizon', line=2, cex.lab=1.3)
par(mar=c(2,4,2,1)) #No margin on the right side
plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
legend("center", legend = colnames(pollution_data), col = plotting_colours, horiz = T, 
       cex = 1.5, pch=18, lty=1, lwd=2, bty = 'n', box.col = "white", pt.cex = 2)


# compute and plot maximum time-target PN
plotting_colours <- RColorBrewer::brewer.pal(n = ncol(pollution_data), name = 'Dark2')
par(mfrow=c(2,3))
for(idx_source in 1:ncol(pollution_data)){
  for(idx_target in 1:ncol(pollution_data)){
    w_t <- xvine::maximise_pn(
      data = cond_no_rvine_sims,
      target = idx_target,
      col_source = idx_source,
      u0_target = u0_single,
      u0_source = u0_single,
      type = 'time'
    )
    cat(paste('\nSource:', colnames(pollution_data)[[idx_source]], '; Target', colnames(pollution_data)[[idx_target]], '\n'))
    cat(paste('Max weight:', max(w_t$weights) * 100, '%\n'))
    cat(paste('Min weight:', min(w_t$weights) * 100, '%\n'))
    cat(paste('Relative weight difference:', (max(w_t$weights) - min(w_t$weights))/min(w_t$weights) * 100, '%\n'))
    
    logits_w <- xvine::logits(w_t$weights)
    if(idx_target==1){
      plot(
        1:(k_markov-1), 
        logits_w,
        ylab = 'Weights (logits)',
        xlab = 'Timestep', main = paste('After an extreme in', colnames(pollution_data)[[idx_source]]),
        ylim=c(-3, -1),
        type='b', col=plotting_colours[idx_target],
        lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4, cex.main=1.6
      )
    }else{
      points(
        1:(k_markov-1), 
        logits_w,
        type='b',
        col=plotting_colours[idx_target],
        lwd=3
      )
    }
  }
  legend(1, -1, colnames(pollution_data), plotting_colours)
}

# compute and plot maximum time-target PN
# reversed!
plotting_colours <- RColorBrewer::brewer.pal(n = ncol(pollution_data), name = 'Dark2')
par(mfrow=c(2,3))
for(idx_source in 1:ncol(pollution_data)){
  for(idx_target in 1:ncol(pollution_data)){
    q_target <- quantile(pollution_data[,idx_target], u0s_fixed[idx_target])
    q_source <- quantile(pollution_data[,idx_source], u0s_fixed[idx_source])
    w_t <- xvine::maximise_pn(
      data = cond_no_rvine_sims_reversed$data,
      target = idx_target,
      col_source = idx_source,
      u0_target = q_target,
      u0_source = q_source,
      type = 'time',
      routine = 'deoptim'
    )
    cat(paste('\nSource:', colnames(pollution_data)[[idx_source]], '; Target', colnames(pollution_data)[[idx_target]], '\n'))
    cat(paste('Max weight:', max(w_t$weights) * 100, '%\n'))
    cat(paste('Min weight:', min(w_t$weights) * 100, '%\n'))
    cat(paste('Relative weight difference:', (max(w_t$weights) - min(w_t$weights))/min(w_t$weights) * 100, '%\n'))
    
    logits_w <- xvine::logits(w_t$weights)
    if(idx_target==1){
      plot(
        1:(k_markov-1), 
        logits_w,
        ylab = 'Weights (logits)',
        xlab = 'Timestep', main = paste('(R) After an extreme in', colnames(pollution_data)[[idx_source]]),
        ylim=c(-3, -1),
        type='b', col=plotting_colours[idx_target],
        lwd=3, cex.lab=1.6, cex=1.3, cex.axis=1.4, cex.main=1.6
      )
    }else{
      points(
        1:(k_markov-1), 
        logits_w,
        type='b',
        col=plotting_colours[idx_target],
        lwd=3
      )
    }
  }
  legend(1, -1, colnames(pollution_data), plotting_colours)
}












