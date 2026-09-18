select_jk_block_size <- function(N){
  size <- ceiling(N / 100)
  return(size)
}

jk_split <- function(data, jk_size){
  K <- ceiling(length(data) / jk_size)
  remain = length(data) %% jk_size
  block_size <- rep(jk_size, times = K)
  if(remain > 0){
    block_size[1] = remain
  }
  data_jk <- split(
    data,
    rep(1:K, times = block_size)
  )
  return(data_jk)
}
