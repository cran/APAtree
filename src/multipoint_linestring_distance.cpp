// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>
#include <Rinternals.h>
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/algorithms/append.hpp>
#include <boost/geometry/algorithms/distance.hpp>

using namespace Rcpp;

namespace bg = boost::geometry;

//' @keywords internal
// [[Rcpp::export]]

NumericVector multipoint_linestring_distance(NumericMatrix multipoint, NumericMatrix linestring){
  typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
  typedef bg::model::linestring<point_t> linestring_t;


  linestring_t ls;
  for(int i = 0; i < linestring.nrow(); ++i){
    bg::append(ls, point_t(linestring(i, 0), linestring(i, 1)));
  }

  point_t pt_i;
  NumericVector mp_ls_dist(multipoint.nrow()) ;
  for(int i = 0; i < multipoint.nrow(); ++i){
    pt_i = point_t(multipoint(i, 0), multipoint(i, 1));
    mp_ls_dist(i) = bg::distance(pt_i, ls);
  }
  //
  return mp_ls_dist;
}
