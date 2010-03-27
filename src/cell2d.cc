#include "cell2d.hh"

voronoicell_2d::voronoicell_2d() :
	current_vertices(init_vertices), ed(new int*[current vertices]),
	pts(new int[2*current_vertices])
//	ds(new int[current_delete_size])
{
	ed[0]=new int[2*current_vertices];
}

voronoicell_2d::~voronoicell_2d() {
	delete [] ed[0];
	delete [] ed;
}

/** Initializes a Voronoi cell as a rectangle with the given dimensions.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates. */
void voronoicell_2d::init(fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax) {
	p=4;xmin*=2;xmax*=2;ymin*=2;ymax*=2;
	pts[0]=xmin;pts[1]=ymin;
	pts[2]=xmax;pts[3]=ymin;
	pts[4]=xmax;pts[5]=ymax;
	pts[5]=xmin;pts[6]=ymax;
	int *q=ed[0];
	q[0]=1;q[1]=3;
	q[2]=2;q[3]=0;ed[1]=q+2;
	q[4]=3;q[5]=1;ed[2]=q+4;
	q[6]=0;q[7]=2;ed[3]=q+6;
}

/** Outputs the edges of the Voronoi cell in gnuplot format to an output stream.
 * \param[in] os a reference to an output stream to write to.
 * \param[in] (x,y,z) a displacement vector to be added to the cell's position.
 */
void voronoicell_2d::draw_gnuplot(ostream &os,fpoint x,fpoint y) {
	int i,j,k;fpoint ux,uy;
	for(i=0;i<p;i++) {
		ux=x+0.5*pts[3*i];uy=y+0.5*pts[3*i+1];
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if(ed[i][j]<i) os << ux << " " << uy << "\n" << x+0.5*pts[3*k] << " " << y+0.5*pts[3*k+1] << "\n\n\n";
		}
	}
}

/** An overloaded version of the draw_gnuplot routine that writes directly to
 * a file.
 * \param[in] filename The name of the file to write to.
 * \param[in] (x,y,z) a displacement vector to be added to the cell's position.
 */
inline void voronoicell_2d::draw_gnuplot(const char *filename,fpoint x,fpoint y,fpoint z) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	draw_gnuplot(os,x,y,z);
	os.close();
}

/** An overloaded version of the draw_gnuplot routine, that prints to the
 * standard output.
 * \param[in] (x,y,z) a displacement vector to be added to the cell's position.
 */
inline void voronoicell_2d::draw_gnuplot(fpoint x,fpoint y,fpoint z) {
	draw_gnuplot(cout,x,y,z);
}
