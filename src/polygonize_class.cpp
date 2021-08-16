#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <numeric>

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <Rinternals.h>
#include <boost/iterator/counting_iterator.hpp>

using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
List polygonize_class_cpp(NumericMatrix x, int direction) {
  //multi_poly: a Vector corresponding to the structure of what
  //sf::st_multipolygon() requires as input: A list of lists of nx2 matrices (in total
  //4 dimensions):
  std::vector<std::vector<std::vector<std::vector<int> > > > multi_poly;
  //An index which refers to each row number of `x`:
  std::vector<int> index( boost::counting_iterator<int>( 0 ),
                          boost::counting_iterator<int>( x.nrow()) );
  // `add_edge` willbe used to add one coordinate to a polygon element in
  // `multi_poly`. `new_edge is a (potential) new coordinate of the current
  // polygon while iterating.
  std::vector<int> add_edge(2), new_edge(2);
  // Preallocation : the angle of the last edge of the (potential) next edge of
  // the current polygon while iterating.
  int old_direction, next_direction;
  // Preallocation : the angular difference of the old and the next_direction:
  int direction_angle;
  // `direction_select` decides if raster points touching over corner (queens
  // case) should be connected. While iterating, the in `direction_select`
  // specified `direction_angle` indicates a correct edge. This is either 1 or
  // -1 as 2 does not exist and 0 always indicates a selection of this edge
  // (irrespective of queens or rooks case).
  int direction_select;
  // Preallocation: next_edge_i refers to the row index of `x` containing the
  // next (potential) edge.
  int next_edge_i;
  // Preallocation: Is the boundary outer (True) or inner (of a hole):
  bool is_outer_boundary;
  // Each time a coordinate (row) in `x` is assigned to a polygon, its row
  // number is erased from `index`. When `index` doesn't have any elements left,
  // the algorithm is done.
  while(index.size() > 0){
    // pick one of the most northern edges as starting coordinate (doesn't
    // matter which one):
    next_edge_i = 0;
    for(unsigned int i = 1; i < index.size(); ++i){
      if((x(index[i], 1) < x(index[next_edge_i], 1)) &
         (x(index[i], 3) == x(index[i], 1))) next_edge_i = i;
    }
    // copy first two starting coordinates
    std::vector<std::vector<int> > poly_boundary;
    add_edge[0] = x(index[next_edge_i], 1);
    add_edge[1] = x(index[next_edge_i], 2);
    poly_boundary.push_back(add_edge);
    add_edge[0] = x(index[next_edge_i], 3);
    add_edge[1] = x(index[next_edge_i], 4);
    poly_boundary.push_back(add_edge);
    old_direction = x(index[next_edge_i], 0);
    index.erase(index.begin() + next_edge_i);
    // assess if a new polygon is started or which polygon gets a new hole
    // punched in it
    is_outer_boundary = TRUE; // is only changed if the starting coordinate is within another polygon
    // to which outer border does an inner_border belong?
    int multi_poly_select = -1;
    // iterate towards north:
    for(int i = poly_boundary[0][0] - 2; i > 0; i = i - 2){
      // iterate over all already existing polygons:
      for(unsigned int j = 0; j < multi_poly.size(); ++j){
        // iterate through outer and inner boundary polygon coordinates:
        for(unsigned int k = 0; k < multi_poly[j].size(); ++k){
          // iterate through individual coordinates
          for(unsigned int l = 0; l < multi_poly[j][k].size() - 1; ++l){
            if((multi_poly[j][k][l][0] == i) & // same row
               (multi_poly[j][k][l + 1][0] == i) & //same row
               ((multi_poly[j][k][l][1] - poly_boundary[0][1]) * //same col
               (multi_poly[j][k][l + 1][1] - poly_boundary[0][1]) <= 0)){ //same col
              if(multi_poly[j][k][l][1] < multi_poly[j][k][l + 1][1]){
                is_outer_boundary = FALSE;
                multi_poly_select = j;
              }
              // as soon as the first boundary towards north is found, exit all loops
              goto endloop;
            }
          }
        }
      }
    }
    endloop: ;
    // inner boundaries need the opposite criteria for connectiviy of cornering
    // raster pixels.
    if(direction == 4) direction_select = 1; // When there are two potential ways to go, select clockwise
    else if(direction == 8) direction_select = -1;// When there are two potential ways to go, select counterclockwise
    if(is_outer_boundary == FALSE) direction_select = direction_select * -1;
    while((add_edge[0] != poly_boundary[0][0]) | (add_edge[1] != poly_boundary[0][1])){
      next_edge_i =  -1;
      for(unsigned int i = 0; i < index.size(); ++i){
        if((add_edge[0] == x(index[i], 1)) & (add_edge[1] == x(index[i], 2))){
          next_edge_i =  i;
          next_direction = x(index[i], 0);
          direction_angle = next_direction - old_direction;
          // transform a 270 degree direction angle (3 or -3) to 90 degree in
          // opposite direction (clockwise to counterclock wise and vice versa)
          if(direction_angle > 1) direction_angle = direction_angle/3 * -1;
          // Only continue to look for more edges at this point if another edge
          // would be the one to select
          if((direction_angle == direction_select) | (direction_angle == 0)) break;
        }
      }
      if(next_edge_i == -1) break;
      add_edge[0] = x(index[next_edge_i], 3);
      add_edge[1] = x(index[next_edge_i], 4);
      if(old_direction == next_direction){
        poly_boundary[poly_boundary.size() - 1] = add_edge;
      } else {
      poly_boundary.push_back(add_edge);
      }
      old_direction = next_direction;

      index.erase(index.begin() + next_edge_i);
    }
    if(multi_poly_select == -1){
      //add a new outer boundary polygon
      std::vector<std::vector<std::vector<int> > > add_poly;
      add_poly.push_back(poly_boundary);
      multi_poly.push_back(add_poly);
    }else{
      //add a new inner boundary to an existing outer boundary
      multi_poly[multi_poly_select].push_back(poly_boundary);
    }
  }
  List multi_poly_list(multi_poly.size());
  for(unsigned int i = 0; i < multi_poly.size(); ++i){
    List poly_list(multi_poly[i].size());
    for(unsigned int j = 0; j < multi_poly[i].size(); ++j){
      NumericMatrix poly_mat(multi_poly[i][j].size(), 2);
      for(unsigned int k = 0; k < multi_poly[i][j].size(); ++k){
        poly_mat(k, 0) = multi_poly[i][j][k][0];
        poly_mat(k, 1) = multi_poly[i][j][k][1];
      }
      poly_list(j) = poly_mat;
    }
    multi_poly_list(i) = poly_list;
  }
  return(multi_poly_list);
}
