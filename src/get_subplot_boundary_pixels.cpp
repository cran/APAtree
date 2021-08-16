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
IntegerMatrix get_subplot_boundary_pixels(StringVector x_wkt, NumericMatrix y) {
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
  std::vector<std::vector<int> > boundary_idx;
  std::vector<int> new_row(4);
  point_t p1, p2;
  for(int i = 0; i < polygon_number; ++i){
    for(int j = 0; j < y.nrow(); j = j + 2){
      p1 = point_t(y(j, 0), y(j, 1));
      p2 = point_t(y(j+1, 0), y(j+1, 1));
      if(!bg::within(p1, polygons[i]) | !bg::within(p2, polygons[i])) continue;
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

