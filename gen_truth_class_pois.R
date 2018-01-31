library(tidyverse)

gen_arima <- function(n,a,b){

  while (TRUE){

    w_ar <- runif(1,.97,.99)
    w_ma <- runif(1,0,.99)
    i_order <- sample(c(0,1),1)
    timeseries <- try(arima.sim(n=n,
                                list(order=c(1,i_order,1),ar=w_ar,ma=w_ma),
                                sd=1),
                      silent=TRUE)
    if (class(timeseries) != 'try-error') break

  }

  timeseries <- scales::rescale(as.vector(timeseries),c(a,b))
  attr(timeseries,'w') <- list(ar=w_ar,i=i_order,ma=w_ma)

  return(timeseries)

}

# class
gen_table <- function(fl_sig=0,w_sig=6,
                      fl_bg=-6,w_bg=6,
                      bg_disp_mu=0,bg_disp_sigma=1,
                      sig_disp_mu1=0,sig_disp_sigma1=0,
                      sig_disp_mu2=0,sig_disp_sigma2=1,
                      n_clust=10,n_sig=10,n_tax_sig=1,n_bg=700,
                      len_arima=1000,len_ts=500,len_signal=300){

  params <- as.list(environment())

  idx_signal <- (len_ts - len_signal/2):(len_ts + len_signal/2 - 1)
  min_ts <- max(idx_signal) - len_ts + 1
  max_ts <- min(idx_signal)

  # generate background distribution
  mu_bg <- scales::rescale(rnorm(len_ts*n_bg,0,1),c(fl_bg,w_bg))
  background <- matrix(sapply(mu_bg,function(mu) rpois(1,exp(mu + rnorm(1,bg_disp_mu,bg_disp_sigma)))),
                       len_ts,n_bg)

  # generate features affacted by time series signal
  dat <- lapply(seq_len(n_clust), function(k){

    # generate n_sig features for cluster k

    timeseries_pure <- gen_arima(len_arima,fl_sig,w_sig) # gen signal for nb distriution means
    timeseries_noise1 <- rnorm(length(timeseries_pure),sig_disp_mu1,sig_disp_sigma1)
    timeseries <- sapply(seq_along(timeseries_pure),function(i) {
                           theta <- exp(timeseries_pure[i] + timeseries_noise1[i])
                           theta <- ifelse(theta > 20, 20, theta)
                           rpois(1,theta)
                           })
    attr(timeseries,'indexes') <- idx_signal

    # sample starting points within window to get indes for n_sig signals of length len_signal
    timesteps <- sapply(seq_len(n_sig),function(i){
      ts_start <- sample(min_ts:max_ts,1)
      ts_start:(ts_start+len_ts-1)
    })

    # generate n_sig time series via nb distribution within frame within window
    cluster <- sapply(seq_len(n_sig),function(i){
      y_tmp <- timeseries[timesteps[,i]]
      sapply(y_tmp,function(mu) {
        theta <- exp(rnorm(1,sig_disp_mu2,sig_disp_sigma2))
        theta <- ifelse(theta > 20, 20, theta)
        mu + rpois(n_tax_sig,theta)
        })
    })

    list(timeseries_pure=timeseries_pure,timeseries=timeseries,cluster=cluster,timesteps=timesteps)

  })
  abund <-  do.call(cbind,lapply(dat,function(x) x$cluster))
  final_table <- cbind(abund,background) # join signal features with background features

  colnames(final_table) <- c(paste0(rep('cl',n_clust*n_sig),rep(1:n_sig,each=n_clust),'_',rep('sig',n_clust*n_sig),rep(1:n_clust,n_sig)),
                             paste0('bg',1:ncol(background)))

  attr(final_table,'signals') <- dat
  attr(final_table,'background') <- background
  attr(final_table,'params') <- params

  class(final_table) <- 'mbts'

  return(final_table)

}


# methods
sig_cor <- function(object,...) UseMethod('sig_cor')
plot_sig <- function(object,...) UseMethod('plot_sig')
plot_sim <- function(object,...) UseMethod('plot_sim')
sparsity <- function(object,...) UseMethod('sparsity')
quantiles <- function(object,...) UseMethod('quantiles')


# print.mbts <- function(x,...){
#   cat('A mbts class object.')
# }

sig_cor.mbts <- function(x,i=1,method='spearman',round=3,...){

  z <- attr(x,'signals')[[i]]
  idx_signal <- attr(z$timeseries,'indexes')

  p <- cor(sapply(seq_len(ncol(z$cluster)),
                  function(j) z$cluster[z$timesteps[,j] %in% idx_signal,j]),
           method=method)

  p <- round(p,round)

  return(p)

}

plot_sig.mbts <- function(x,n=6,seed=sample.int(.Machine$integer.max,1)){

  params <- attr(x,'params')

  sims <- do.call(rbind,lapply(seq_len(n),function(x){
    w <- gen_arima(params$len_arima,params$fl_sig,params$w_sig)
    tibble(w=as.numeric(w),t=1:length(w),sim=x)
  }))

  # return(sims)

  par(mfrow=c(3,2))
  sapply(seq_len(n),function(x) plot(gen_arima(params$len_arima,params$fl_sig,params$w_sig)))
  par(mfrow=c(1,1))

}

plot_sim.mbts <- function(x,i=1){

  z <- attr(x,'signals')[[i]]
  idx_signal <- attr(z$timeseries,'indexes')

  df <- data.frame(ts='count',
                   t=matrix(z$timesteps),
                   count=matrix(z$cluster),
                   taxa=rep(1:ncol(z$cluster),each=nrow(z$cluster)),
                   stringsAsFactors=FALSE) %>%
    left_join(data.frame(t=min(.$t):max(.$t),signal=as.vector(z$timeseries)[min(.$t):max(.$t)]),by='t') %>%
    group_by(taxa) %>%
    mutate(signal=round(scales::rescale(signal,c(0,max(count)))),
           min_t=min(t),
           max_t=max(t))

  #return(df)

  p1 <- df %>%
    ggplot() +
    geom_rect(aes(xmin=min_t,xmax=max_t,ymin=-Inf,ymax=Inf),fill='gray') +
    annotate('rect',xmin=range(idx_signal)[1],xmax=range(idx_signal)[2],ymin=-Inf,ymax=Inf,alpha=0.2,fill='dodgerblue') +
    geom_line(aes(t,signal),linetype=3,color='red') +
    geom_line(aes(t,count),alpha=.5) +
    geom_point(aes(t,count),alpha=1) +
    facet_wrap(~taxa,ncol=2) +
    stat_smooth(aes(t,count),method='loess',color='green',se=FALSE,span=.1,size=.7,alpha=.5) +
    theme_classic() +
    xlim(range(df$t)) +
    labs(x='time',y='count')

  p1

}

quantiles.mbts <- function(x,min=.75,max=1,length=15){

  background <- attr(x,'background')
  signals <- attr(x,'signals')
  abund <-  do.call(cbind,lapply(signals,function(x) x$cluster))

  q <- cbind(quantile(abund,seq(min,max,length=length)),quantile(background,seq(min,max,length=length)))
  colnames(q) <- c('signal','background')

  print(q)

}

sparsity.mbts <- function(x,min=.75,max=1,length=15){

  background <- attr(x,'background')
  signals <- attr(x,'signals')
  abund <-  do.call(cbind,lapply(signals,function(x) x$cluster))

  s <- c(mean(abund==0),mean(background==0))
  names(s) <- c('signal','background')

  print(s)

}
