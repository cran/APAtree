// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/io/wkt/read.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/algorithms/append.hpp>
#include <boost/geometry/algorithms/within.hpp>
using namespace Rcpp;
namespace bg = boost::geometry;

//' @keywords internal
// [[Rcpp::export]]
IntegerMatrix subplot_ca(CharacterVector x_wkt, NumericMatrix y) {
  typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
  typedef bg::model::polygon<point_t> polygon_t;

  int polygon_number = x_wkt.size();
  std::vector<polygon_t> polygons(polygon_number);
  polygon_t new_poly;
  std::string polygon_wkt;
  for(unsigned int i = 0; i < polygon_number; ++i){
    polygon_wkt = x_wkt(i);
    bg::read_wkt(polygon_wkt, new_poly);
    polygons[i] = new_poly;
  }

  int n_layer = y.ncol() - 2;
  NumericVector layer_size(n_layer);
  NumericVector layer_size_cumprod(n_layer + 1);
  layer_size_cumprod[layer_size_cumprod.size() - 1] = 1;
  for (int i = n_layer-1; i != -1; i--) {
    NumericMatrix::Column layer_values = y( _, i + 2);
    layer_size[i] = max(layer_values) + 1;
    layer_size_cumprod[i] = layer_size_cumprod[i + 1] * layer_size[i] ;
  }
  int output_nrow = polygon_number * layer_size_cumprod[0];
  IntegerMatrix output(output_nrow, n_layer + 2);
  NumericVector level_vec;
  for (int i = n_layer - 1; i != -1; i--) {
    level_vec = rep_each(seq_len(layer_size[i]), layer_size_cumprod[i+1]);
    level_vec = rep_len(level_vec, output_nrow);
    IntegerMatrix::Column layer_col = output( _, i + 1);
    layer_col = level_vec;
  }
  IntegerMatrix::Column x_col = output( _, 0);
  x_col = rep_each(seq_len(polygon_number), layer_size_cumprod[0]);

  int active_row, output_tick_idx = n_layer + 1;
  point_t p;
  for(int i = 0; i < polygon_number; ++i){
    for(int j = 0; j < y.nrow(); ++j){
      p = point_t(y(j, 0), y(j, 1));
      if(!bg::within(p, polygons[i])) continue;
      active_row = i * layer_size_cumprod[0];
      for(int k = 0; k < n_layer; ++k){
        active_row = active_row + y(j, 2 + k) * layer_size_cumprod(k + 1);
      }
      ++output(active_row, output_tick_idx);
    }
  }
  return(output);
}

