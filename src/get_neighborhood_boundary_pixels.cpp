#include <Rcpp.h>
#include <cmath>        // std::abs
using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
IntegerMatrix get_neighborhood_boundary_pixels(NumericMatrix x, NumericMatrix y, double radius) {
  std::vector<std::vector<int> > boundary_idx;
  std::vector<int> new_row(4);

  double x_diff_1, x_diff_2, y_diff_1, y_diff_2, xy_dist;
  for(int i = 0; i < x.nrow(); ++i){
    for(int j = 0; j < y.nrow(); j = j + 2){
      x_diff_1 = std::abs(x(i, 0) - y(j, 0));
      x_diff_2 = std::abs(x(i, 0) - y(j + 1, 0));
      if(x_diff_2 > x_diff_1) x_diff_1 = x_diff_2;
      if(x_diff_1 > radius) continue;
      y_diff_1 = std::abs(x(i, 1) - y(j, 1));
      y_diff_2 = std::abs(x(i, 1) - y(j + 1, 1));
      if(y_diff_2 > y_diff_1) y_diff_1 = y_diff_2;
      if(y_diff_1 > radius) continue;
      xy_dist = sqrt(x_diff_1 * x_diff_1 + y_diff_1 * y_diff_1);
      if(xy_dist > radius) continue;
      new_row[0] = i + 1;
      new_row[1] = y(j, 2);
      new_row[2] = y(j, 3);
      new_row[3] = y(j, 4);
      boundary_idx.push_back(new_row);
      new_row[1] = y(j + 1, 2);
      new_row[2] = y(j + 1, 3);
      new_row[3] = y(j + 1, 4);
      boundary_idx.push_back(new_row);
    }
  }
  IntegerMatrix output(boundary_idx.size(), 4);
  for(int i = 0; i < output.nrow(); i++){
    for(int j = 0; j < 4; j++){
      output(i, j) = boundary_idx[i][j];
    }
  }
  return(output);
}

