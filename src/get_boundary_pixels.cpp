#include <Rcpp.h>
using namespace Rcpp;


//' @keywords internal
// [[Rcpp::export]]
IntegerMatrix get_boundary_pixels(IntegerMatrix x) {
  std::vector<std::vector<int> > boundary_idx;
  std::vector<int> new_row(3);
  IntegerVector active_cell_value(1), neighbor_cell_value(1);
  LogicalVector is_active_na(1);
  LogicalVector is_boundary_different(1);
  int neighbor_row, neighbor_col;
  int col_offset[4] = {-1, 0, 1, 0};
  int row_offset[4] = {0, 1, 0, -1};
  int direction[4] = {0, 1, 2, 3};
  for(int i = 1; i < x.nrow() - 1; ++i){
    for(int j = 1; j < x.ncol() - 1; ++j){
      active_cell_value = x(i, j);
      is_active_na = IntegerVector::is_na(active_cell_value[0]); // Na value for int
      if(is_active_na[0]) continue;
      for(int k = 0; k < 4; ++k){
        neighbor_row = i + col_offset[k];
        neighbor_col = j + row_offset[k];
        neighbor_cell_value = x(neighbor_row, neighbor_col);
        is_boundary_different = neighbor_cell_value != active_cell_value;
        if(is_boundary_different[0]){
          new_row[0] = 5; // active cell
          new_row[1] = i + 1;
          new_row[2] = j + 1;
          boundary_idx.push_back(new_row);
          new_row[0] = direction[k]; // boundary cell (5 - NA, 0 - North, 1 - East, ...)
          new_row[1] = neighbor_row + 1;
          new_row[2] = neighbor_col + 1;
          boundary_idx.push_back(new_row);
        }
      }
    }
  }
  IntegerMatrix output(boundary_idx.size(), 3);
  for(int i = 0; i < output.nrow(); i++){
    for(int j = 0; j < 3; j++){
      output(i, j) = boundary_idx[i][j];
    }
  }
  return(output);
}
