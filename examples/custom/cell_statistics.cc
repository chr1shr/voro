// Simple cell statistics demonstration code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 10th 2009

#include "voro++.cc"

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	double x,y,z;
	voronoicell v;

	// Initialize the Voronoi cell to be a cube of side length 2, centered
	// on the origin
	v.init(-1,1,-1,1,-1,1);
	
	// Remove one edge of the cell with a single plane cut
	v.plane(1,1,0,2);

	// Output the Voronoi cell to a file in gnuplot format
	v.draw_gnuplot("simple_cell.gnu",0,0,0);

	// Output vertex-based statistics
	cout << "Total vertices      : " << v.p << "\n";
	cout << "Vertex positions    : ";v.output_vertices(cout);cout << "\n";
	cout << "Vertex orders       : ";v.output_vertex_orders(cout);cout << "\n";
	cout << "Max rad. sq. vertex : " << 0.25*v.max_radius_squared() << "\n\n";
	
	// Output edge-based statistics
	cout << "Total edges         : " << v.number_of_edges() << "\n";
	cout << "Total edge distance : " << v.total_edge_distance() << "\n";
	cout << "Face perimeters     : ";v.output_face_perimeters(cout);cout << "\n\n";
	
	// Output face-based statistics
	cout << "Total faces         : " << v.number_of_faces() << "\n";
	cout << "Surface area        : " << v.surface_area() << "\n";
	cout << "Face freq. table    : ";v.output_face_freq_table(cout);cout << "\n";
	cout << "Face orders         : ";v.output_face_orders(cout);cout << "\n";
	cout << "Face areas          : ";v.output_face_areas(cout);cout << "\n";
	cout << "Face normals        : ";v.output_normals(cout);cout << "\n";
	cout << "Face vertices       : ";v.output_face_vertices(cout);cout << "\n\n";

	// Output volume-based statistics
	cout << "Volume              : " << v.volume() << "\n";
	v.centroid(x,y,z);
	cout << "Centroid vector     : (" << x << "," << y << "," << z << ")" << endl;

}
