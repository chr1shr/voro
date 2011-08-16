#include "voro++_2d.hh"
#include <iostream>

int main() {
container_2d con(-1.1,1.1,-1.1,1.1,4,4,false,false,false,16);
con.import("/users/mac/voro/branches/2d_boundary/examples/bd_test");
con.setup();
con.debug_output();
con.draw_cells_gnuplot("bd_test.gnu");



}
