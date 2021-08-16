#include <Rcpp.h>
#include <cmath>        // std::abs
using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
IntegerMatrix neighborhood_ca(NumericMatrix x, NumericMatrix y, double radius) {
  int n_layer = y.ncol() - 2;
  
  NumericVector layer_size(n_layer);
  NumericVector layer_size_cumprod(n_layer + 1);
  layer_size_cumprod[layer_size_cumprod.size() - 1] = 1;
  int x_size = x.nrow();
  for (int i = n_layer-1; i != -1; i--) {
    NumericMatrix::Column layer_values = y( _, i + 2);
    layer_size[i] = max(layer_values) + 1;
    layer_size_cumprod[i] = layer_size_cumprod[i + 1] * layer_size[i] ;
  }
  int output_nrow = x_size * layer_size_cumprod[0];
  IntegerMatrix output(output_nrow, n_layer + 2);
  NumericVector level_vec;
  for (int i = n_layer - 1; i != -1; i--) {
    level_vec = rep_each(seq_len(layer_size[i]), layer_size_cumprod[i+1]);
    level_vec = rep_len(level_vec, output_nrow);
    IntegerMatrix::Column layer_col = output( _, i + 1);
    layer_col = level_vec;
  }
  IntegerMatrix::Column x_col = output( _, 0);
  x_col = rep_each(seq_len(x_size), layer_size_cumprod[0]);

  double x_diff, y_diff, xy_dist;
  int active_row, output_tick_idx = n_layer + 1;;
  for(int i = 0; i < x.nrow(); ++i){
    for(int j = 0; j < y.nrow(); ++j){
      x_diff = x(i, 0) - y(j, 0);
      if(std::abs(x_diff) > radius){
        continue;
      }
      y_diff = x(i, 1) - y(j, 1);
      if(std::abs(y_diff) > radius){
        continue;
      }
      xy_dist = sqrt(x_diff * x_diff + y_diff * y_diff);
      if(xy_dist > radius) continue;
      active_row = i * layer_size_cumprod[0];
      for(int k = 0; k < n_layer; ++k){
        active_row = active_row + y(j, 2 + k) * layer_size_cumprod(k + 1);
      }
      ++output(active_row, output_tick_idx);
    }
  }
  return(output);
}

