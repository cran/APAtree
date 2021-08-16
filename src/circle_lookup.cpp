#include <Rcpp.h>
#include <cmath>        // std::abs
using namespace Rcpp;

#include <vector>

//' @keywords internal
// [[Rcpp::export]]
NumericMatrix circle_lookup(NumericMatrix x, NumericMatrix y, double radius) {
  std::vector<std::vector<int> > vect;
  std::vector<int> new_row(2);
  double x_diff, y_diff, xy_dist;
  for(int i = 0; i < x.nrow(); ++i){
    for(int j = 0; j < y.nrow(); ++j){
      x_diff = x(i, 0) - y(j, 0);
      if(std::abs(x_diff) > radius) continue;
      y_diff = x(i, 1) - y(j, 1);
      if(std::abs(y_diff) > radius) continue;
      xy_dist = sqrt(x_diff * x_diff + y_diff * y_diff);
      if(xy_dist > radius) continue;
      new_row[0] = i;
      new_row[1] = j;
      vect.push_back(new_row);
    }
  }
  NumericMatrix joined(vect.size(), 2);
  for(int i = 0; i < joined.nrow(); i++){
      joined(i, 0) = vect[i][0] + 1;
      joined(i, 1) = vect[i][1] + 1;
  }
  return(joined);
}

