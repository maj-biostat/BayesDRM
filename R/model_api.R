drm_emax2_bin <- function(d, ...){
  out <- rstan::sampling(stanmodels$drm01, data = d, ...)
  return(out)
}

drm_emax3_bin <- function(d, ...){
  out <- rstan::sampling(stanmodels$drm02, data = d, ...)
  return(out)
}

drm_emax4_bin <- function(d, ...){
  out <- rstan::sampling(stanmodels$drm03, data = d, ...)
  return(out)
}
