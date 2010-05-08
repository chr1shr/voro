#include "cell_2d.cc"
#include "container_2d.cc"

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	double x,y,rsq,r;
	voronoicell_2d v;

	// Initialize the Voronoi cell to be a cube of side length 2, centered
	// on the origin
	v.init(-1,1,-1,1);

	// Cut the cell by 100 random planes which are all a distance 1 away
	// from the origin, to make an approximation to a sphere
	for(int i=0;i<100;i++) {
		x=2*rnd()-1;
		y=2*rnd()-1;
		rsq=x*x+y*y;
		if(rsq>0.01&&rsq<1) {
			r=1/sqrt(rsq);x*=r;y*=r;
			v.plane(x,y,1);
		}
	}

	// Print out several statistics about the computed cell
	cout << "Perimeter is " << v.perimeter() << endl;
	cout << "Area is " << v.area() << endl;
	v.centroid(x,y);
	cout << "Centroid is (" << x << "," << y << ")" << endl;

	// Output the Voronoi cell to a file, in the gnuplot format
	v.draw_gnuplot("single_cell.gnu",0,0);
}
