drm_emax2_bin <- function(d, ...){
  out <- rstan::sampling(stanmodels$drm01, data = d, ...)
  return(out)
}


