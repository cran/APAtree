full_ca <- function(x){
  max_val <- apply(x, 2, max) + 1
  fac_val <- c(1, cumprod(max_val)[-ncol(x)])
  fac_mat <- matrix(rep(fac_val, nrow(x)), ncol = ncol(x), byrow = TRUE)
  count <- tabulate(rowSums(x * fac_mat) + 1)
  req_tab <- cbind(.feature_no = 1, expand.grid(lapply(max_val, seq_len)), area = 0)
  req_tab[1:length(count), "area"] <- count
  req_tab <- req_tab[do.call(what = order, as.data.frame(req_tab)[-ncol(x)-2]), ]
  req_tab
}