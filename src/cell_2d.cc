/** \file cell_2d.cc
 * \brief Function implementations for the voronoicell_2d class. */

#include "cell_2d.hh"

/** Constructs a 2D Voronoic cell and sets up the initial memory. */
voronoicell_2d::voronoicell_2d() :
	current_vertices(init_vertices), current_delete_size(init_delete_size),
	ed(new int*[current_vertices]), pts(new fpoint[2*current_vertices]),
	ds(new int[current_delete_size]) {
	ed[0]=new int[2*current_vertices];
	for(int i=1;i<current_vertices;i++) ed[i]=ed[i-1]+2;
}

/** The voronoicell_2d destructor deallocates all of the dynamic memory. */
voronoicell_2d::~voronoicell_2d() {
	delete [] ed[0];
	delete [] ed;
}

/** Doubles the storage for the vertices, by reallocating the pts and ed
 * arrays. If the allocation exceeds the absolute maximum set in max_vertices,
 * then the routine exits with a fatal error. */
void voronoicell_2d::add_memory_vertices() {
	int i=(current_vertices<<1);
	if(i>max_vertices) voropp_fatal_error("Vertex memory allocation exceeded absolute maximum",VOROPP_MEMORY_ERROR);
#if VOROPP_VERBOSE >=2
	cerr << "Vertex memory scaled up to " << i << endl;
#endif
	fpoint *ppts(new fpoint[2*i]);
	int j,*ped(new int[2*i]);
	for(j=0;j<2*current_vertices;j++) ped[j]=ed[0][j];
	delete [] ed[0];
	delete [] ed;
	ed=new int*[i];
	for(j=0;j<i;j++) ed[j]=ped+(2*j);
	for(j=0;j<2*current_vertices;j++) ppts[j]=pts[j];
	delete [] pts;
	pts=ppts;
	current_vertices=i;
}

/** Doubles the size allocation of the delete stack. If the allocation exceeds
 * the absolute maximum set in max_delete_size, then routine causes a fatal
 * error. */
void voronoicell_2d::add_memory_ds() {
	int i(current_delete_size<<1);
	if(i>max_delete_size) voropp_fatal_error("Delete stack memory allocation exceeded absolute maximum",VOROPP_MEMORY_ERROR);
#if VOROPP_VERBOSE >=2
	cerr << "Delete stack memory scaled up to " << i << endl;
#endif
	int j,*pds(new int[i]);
	for(j=0;j<current_delete_size;j++) pds[j]=ds[j];
	delete [] ds;ds=pds;
	current_delete_size=i;
}

/** Initializes a Voronoi cell as a rectangle with the given dimensions.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates. */
void voronoicell_2d::init(fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax) {
	p=4;xmin*=2;xmax*=2;ymin*=2;ymax*=2;
	pts[0]=xmin;pts[1]=ymin;
	pts[2]=xmax;pts[3]=ymin;
	pts[4]=xmax;pts[5]=ymax;
	pts[6]=xmin;pts[7]=ymax;
	int *q=ed[0];
	q[0]=1;q[1]=3;
	q[2]=2;q[3]=0;
	q[4]=3;q[5]=1;
	q[6]=0;q[7]=2;
}

/** Outputs the edges of the Voronoi cell in gnuplot format to an output
 * stream.
 * \param[in] os a reference to an output stream to write to.
 * \param[in] (x,y) a displacement vector to be added to the cell's position.
 */
void voronoicell_2d::draw_gnuplot(ostream &os,fpoint x,fpoint y) {
	if(p==0) return;
	int k=0;
	do {
		os << x+0.5*pts[2*k] << " " << y+0.5*pts[2*k+1] << "\n";
		k=ed[k][0];
	} while (k!=0);
	os << x+0.5*pts[0] << " " << y+0.5*pts[1] << "\n";
}

/** An overloaded version of the draw_gnuplot routine that writes directly to a
 * file.
 * \param[in] filename The name of the file to write to.
 * \param[in] (x,y) a displacement vector to be added to the cell's position.
 */
inline void voronoicell_2d::draw_gnuplot(const char *filename,fpoint x,fpoint y) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	draw_gnuplot(os,x,y);
	os.close();
}

/** An overloaded version of the draw_gnuplot routine, that prints to the
 * standard output.
 * \param[in] (x,y) a displacement vector to be added to the cell's position.
 */
inline void voronoicell_2d::draw_gnuplot(fpoint x,fpoint y) {
	draw_gnuplot(cout,x,y);
}

/** Outputs the edges of the Voronoi cell in POV-Ray format to an open file
 * stream, displacing the cell by given vector.
 * \param[in] os a output stream to write to.
 * \param[in] (x,y,z) a displacement vector to be added to the cell's position.
 */
void voronoicell_2d::draw_pov(ostream &os,fpoint x,fpoint y,fpoint z) {
	if(p==0) return;
	int k=0;
	do {
		os << "sphere{<" << x+0.5*pts[2*k] << "," << y+0.5*pts[2*k+1] << "," << z << ">,r}\n";
		os << "cylinder{<" << x+0.5*pts[2*k] << "," << y+0.5*pts[2*k+1] << "," << z << ">,<";
		k=ed[k][0];
		os << x+0.5*pts[2*k] << "," << y+0.5*pts[2*k+1] << "," << z << ">}\n";
	} while (k!=0);
}

/** An overloaded version of the draw_pov routine that writes directly to
 * a file.
 * \param[in] filename The name of the file to write to.
 * \param[in] (x,y,z) a displacement vector to be added to the cell's position.
 */
inline void voronoicell_2d::draw_pov(const char *filename,fpoint x,fpoint y,fpoint z) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	draw_pov(os,x,y,z);
	os.close();
}

/** An overloaded version of the draw_pov routine, that prints to the standard
 * output.
 * \param[in] (x,y,z) a displacement vector to be added to the cell's position.
 */
inline void voronoicell_2d::draw_pov(fpoint x,fpoint y,fpoint z) {
	draw_pov(cout,x,y,z);
}

/** Computes the maximum radius squared of a vertex from the center of the
 * cell. It can be used to determine when enough particles have been testing an
 * all planes that could cut the cell have been considered.
 * \return The maximum radius squared of a vertex.*/
fpoint voronoicell_2d::max_radius_squared() {
	int i;fpoint r,s;
	r=pts[0]*pts[0]+pts[1]*pts[1];
	for(i=2;i<2*p;i+=2) {
		s=pts[i]*pts[i]+pts[i+1]*pts[i+1];
		if(s>r) r=s;
	}
	return r;
}

/** Computes this distance of a Voronoi cell vertex to a plane.
 * \param[in] (x,y) the normal vector to the plane.
 * \param[in] rsq the distance along this vector of the plane.
 * \param[in] qp the index of the vertex to consider. */
inline fpoint voronoicell_2d::pos(fpoint x,fpoint y,fpoint rsq,int qp) {
	return x*pts[2*qp]+y*pts[2*qp+1]-rsq;
}

/** Cuts the Voronoi cell by a particle whose center is at a separation of
 * (x,y) from the cell center. The value of rsq should be initially set to
 * \f$x^2+y^2\f$.
 * \param[in] (x,y) the normal vector to the plane.
 * \param[in] rsq the distance along this vector of the plane.
 * \return False if the plane cut deleted the cell entirely, true otherwise. */
bool voronoicell_2d::plane(fpoint x,fpoint y,fpoint rsq) {
	int cp,lp,up=0,up2,up3,stack=0;fpoint fac,l,u,u2,u3;
	
	// Finish this section with an inside vertex if there is one
	u=pos(x,y,rsq,up);
	if(u<tolerance) {
		up2=ed[up][0];u2=pos(x,y,rsq,up2);
		up3=ed[up][1];u3=pos(x,y,rsq,up3);
		if(u2>u3) {
			while(u2<tolerance) {
				up2=ed[up2][0];
				u2=pos(x,y,rsq,up2);
				if(up2==up3) return true;
			}
			up=up2;u=u2;
		} else {
			while(u3<tolerance) {
				up3=ed[up3][1];
				u3=pos(x,y,rsq,up3);
				if(up2==up3) return true;
			}
			up=up3;u=u3;
		}
	}

	ds[stack++]=up;
	l=u;up2=ed[up][0];
	u2=pos(x,y,rsq,up2);
	while(u2>tolerance) {
		if(stack==current_delete_size) add_memory_ds();
		ds[stack++]=up2;
		up2=ed[up2][0];
		l=u2;
		u2=pos(x,y,rsq,up2);
		if(up2==up) return false;
	}
	
	if(u2>-tolerance) {
		// Adjust existing vertex
		cp=up2;
	} else {
		// Create new vertex
		if(p==current_vertices) add_memory_vertices();
		lp=ed[up2][1];
		fac=1/(u2-l);
		pts[2*p]=(pts[2*lp]*u2-pts[2*up2]*l)*fac;
		pts[2*p+1]=(pts[2*lp+1]*u2-pts[2*up2+1]*l)*fac;
		ed[p][0]=up2;
		ed[up2][1]=p;
		cp=p++;
	}

	l=u;up3=ed[up][1];u3=pos(x,y,rsq,up3);
	while(u3>tolerance) {
		if(stack==current_delete_size) add_memory_ds();
		ds[stack++]=up3;
		up3=ed[up3][1];
		l=u3;
		u3=pos(x,y,rsq,up3);
		if(up3==up2) break;
	}

	if(u3>tolerance) {
		// Adjust existing vertex
		ed[cp][1]=up3;
		ed[up3][0]=cp;
	} else {
		// Create new vertex
		if(p==current_vertices) add_memory_vertices();
		lp=ed[up3][0];
		fac=1/(u3-l);
		pts[2*p]=(pts[2*lp]*u3-pts[2*up3]*l)*fac;
		pts[2*p+1]=(pts[2*lp+1]*u3-pts[2*up3+1]*l)*fac;
		ed[p][0]=cp;
		ed[cp][1]=p;
		ed[p][1]=up3;
		ed[up3][0]=p++;
	}

	for(int i=0;i<stack;i++) ed[ds[i]][0]=-1;
	
	while(stack>0) {
		while(ed[--p][0]==-1);
		up=ds[--stack];
		if(up<p) {
			ed[ed[p][0]][1]=up;
			ed[ed[p][1]][0]=up;
			pts[2*up]=pts[2*p];
			pts[2*up+1]=pts[2*p+1];
			ed[up][0]=ed[p][0];
			ed[up][1]=ed[p][1];
		} else p++;
	}

	return true;
}

/** Calculates the perimeter of the Voronoi cell.
 * \return A floating point number holding the calculated distance. */
fpoint voronoicell_2d::perimeter() {
	if(p==0) return 0;
	int k=0,l;fpoint perim=0,dx,dy;
	do {
		l=ed[k][0];
		dx=pts[2*k]-pts[2*l];
		dy=pts[2*k+1]-pts[2*l+1];
		perim+=sqrt(dx*dx+dy*dy);
		k=l;
	} while (k!=0);
	return 0.5*perim;
}

/** Calculates the area of the Voronoi cell.
 * \return A floating point number holding the calculated distance. */
fpoint voronoicell_2d::area() {
	if(p==0) return 0;
	int k=ed[0][0];fpoint area=0,x=pts[0],y=pts[1],dx1,dy1,dx2,dy2;
	dx1=pts[2*k]-x;dy1=pts[2*k+1]-y;
	k=ed[k][0];
	while(k!=0) {
		dx2=pts[2*k]-x;dy2=pts[2*k+1]-y;
		area+=dx1*dy2-dx2*dy1;
		dx1=dx2;dy1=dy2;
		k=ed[k][0];
	}
	return 0.25*area;
}

/** Calculates the centroid of the Voronoi cell.
 * \param[out] (cx,cy) The coordinates of the centroid. */
void voronoicell_2d::centroid(fpoint &cx,fpoint &cy) {
	cx=cy=0;
	static const fpoint third=1/3.0;
	if(p==0) return;
	int k=ed[0][0];
	fpoint area,tarea=0,x=pts[0],y=pts[1],dx1,dy1,dx2,dy2;
	dx1=pts[2*k]-x;dy1=pts[2*k+1]-y;
	k=ed[k][0];
	while(k!=0) {
		dx2=pts[2*k]-x;dy2=pts[2*k+1]-y;
		area=dx1*dy2-dx2*dy1;
		tarea+=area;
		cx+=area*(dx1+dx2);
		cy+=area*(dy1+dy2);
		dx1=dx2;dy1=dy2;
		k=ed[k][0];
	}
	tarea=third/tarea;
	cx=0.5*(x+cx*tarea);
	cy=0.5*(y+cy*tarea);
}
