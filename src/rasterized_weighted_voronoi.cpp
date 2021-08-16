// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <Rinternals.h>
#include <cmath>        // std::abs

using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
List rasterized_weighted_voronoi(NumericMatrix rst, NumericMatrix points) {
  int rst_nrow = rst.nrow();
  int points_nrow = points.nrow();
  double squared_max_weight = pow(max(points( _, 2)), 2);

  // initialize variables needed within the loop
  double x_diff, y_diff, weight, weighted_dist_old, weighted_dist_new, critical_dist;
  int previous_tree = 0;
  NumericVector assigned_trees_vec(rst_nrow);
  LogicalVector critical(rst_nrow);

  for(int i = 0; i < rst_nrow; ++i){
    //initialize with previous tree as it is in many cases close to pixel i
    assigned_trees_vec(i) = previous_tree;
    x_diff = rst(i, 0) - points(previous_tree, 0);
    y_diff = rst(i, 1) - points(previous_tree, 1);
    weight = points(previous_tree, 2);
    weighted_dist_old = (x_diff * x_diff + y_diff * y_diff) / (weight * weight);
    critical_dist = sqrt(weighted_dist_old * squared_max_weight);
    for(int j = 0; j < points_nrow; ++j){
      x_diff = rst(i, 0) - points(j, 0);
      if(std::abs(x_diff) > critical_dist){
        continue;
      }
      y_diff = rst(i, 1) - points(j, 1);
      if(std::abs(y_diff) > critical_dist){
        continue;
      }
      weight = points(j, 2);
      weighted_dist_new = (x_diff * x_diff + y_diff * y_diff) / (weight * weight);
      if(weighted_dist_new < weighted_dist_old){
        weighted_dist_old = weighted_dist_new;
        assigned_trees_vec(i) = j;
        critical_dist = sqrt(weighted_dist_old * squared_max_weight);
      }
    }
    critical(i) = critical_dist >= rst(i, 2);
  }
  List result(2);
  result(0) = assigned_trees_vec + 1;
  result(1) = critical;
  return result;
}
