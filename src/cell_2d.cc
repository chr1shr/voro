/**d \file cell_2d.cc
 * \brief Function implementations for the voronoicell_2d class. */

#include "cell_2d.hh"
#include <iostream>
bool debugging=false;
/** Constructs a 2D Voronoi cell and sets up the initial memory. */
voronoicell_2d::voronoicell_2d() :
	current_vertices(init_vertices), current_delete_size(init_delete_size),
	ed(new int[2*current_vertices]), pts(new double[2*current_vertices]),
	ds(new int[current_delete_size]), stacke(ds+current_delete_size) {
}

/** The voronoicell_2d destructor deallocates all of the dynamic memory. */
voronoicell_2d::~voronoicell_2d() {
	delete [] ds;
	delete [] pts;
	delete [] ed;
}

/** Initializes a Voronoi cell as a rectangle with the given dimensions.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates. */
void voronoicell_2d::init(double xmin,double xmax,double ymin,double ymax) {
	p=4;xmin*=2;xmax*=2;ymin*=2;ymax*=2;
	*pts=xmin;pts[1]=ymin;
	pts[2]=xmax;pts[3]=ymin;
	pts[4]=xmax;pts[5]=ymax;
	pts[6]=xmin;pts[7]=ymax;
	int *q=ed;
	*q=1;q[1]=3;q[2]=2;q[3]=0;q[4]=3;q[5]=1;q[6]=0;q[7]=2;
}
/** Initialize a Voronoi cell as exactly fitting the non-convex domain specified in bnds
	* \param[in] bnds, a pointer to the beginning of the array specifing the intial vertices.
	The first point int bnds MUST be the origin, the second point in bnds MUST be the first edge of 
	the nonconvexity in the the counterclockwise direction. All subsequent points MUST define the cel
	l in a counterclockwise fashion.
	* \param[in] #bnds an int specifying the size of bnds. */
	
void voronoicell_2d::init_nonconvex(double bnds_loc[], int noofbnds){
	
	double x1, y1, x2, y2;
	p=noofbnds;  
	for(int i=0; i<(2*noofbnds); i+=2){
		pts[i]=2*bnds_loc[i];
		pts[i+1]=2*bnds_loc[i+1];
		if(i==0){
			ed[0]=noofbnds-1;
		}else{
			ed[i]=(i/2)-1;
		}if((i+2)==(noofbnds*2)){
			ed[i+1]=0;	
		}else{
			ed[i+1]=(i/2)+1;
		}
	}
	
	if (debugging){
		for(int j=0; j<(2*noofbnds); j+=2){
			cout << pts[j] << "  " << pts[j+1] << "  " << ed[2*j] << "   " << ed[(2*j)+1] << 
			"   ";
		}
	}


	x1=bnds_loc[2]; y1=bnds_loc[3]; x2=bnds_loc[(2*noofbnds)-2]; y2=bnds_loc[(2*noofbnds)-1];
	reg1[0]=x1; reg1[1]=y1; reg2[0]=x2; reg2[1]=y2;
	if(y1>0){
		reg1[2]=5;
		reg1[3]=-((5*x1)/y1);
	}if(y1<0){
		reg1[2]=-5;
		reg1[3]=((5*x1)/y1);
	}if(y1==0){
		if(x1<0){
			reg1[2]=0;
			reg1[3]=5;
		}else{
			reg1[2]=0;
			reg1[3]=-5;
		}
	}if(y2>0){
		reg2[2]=-5;
		reg2[3]=((5*x2)/y2);
	}if(y2<0){
		reg2[2]=5;
		reg2[3]=-((5*x2)/y2);
	}if(y2==0){
		if(x2<0){
			reg2[2]=0;
			reg2[3]=-5;
		}else{
			reg2[2]=0;
			reg2[3]=5;
		}
	}
}

/** Outputs the edges of the Voronoi cell in gnuplot format to an output
 * stream.
 * \param[in] (x,y) a displacement vector to be added to the cell's position.
 * \param[in] fp the file handle to write to. */
void voronoicell_2d::draw_gnuplot(double x,double y,FILE *fp) {
	if(p==0) return;
	int k=0;
	do {
		fprintf(fp,"%g %g\n",x+0.5*pts[2*k],y+0.5*pts[2*k+1]);
		k=ed[2*k];
	} while (k!=0);
	fprintf(fp,"%g %g\n\n",x+0.5*pts[0],y+0.5*pts[1]);
}

/** Outputs the edges of the Voronoi cell in POV-Ray format to an open file
 * stream, displacing the cell by given vector.
 * \param[in] (x,y,z) a displacement vector to be added to the cell's position.
 * \param[in] fp the file handle to write to. */
void voronoicell_2d::draw_pov(double x,double y,double z,FILE *fp) {
	if(p==0) return;
	int k=0;
	do {
		fprintf(fp,"sphere{<%g,%g,%g>,r}\ncylinder{<%g,%g,%g>,<"
			,x+0.5*pts[2*k],y+0.5*pts[2*k+1],z
			,x+0.5*pts[2*k],y+0.5*pts[2*k+1],z);
		k=ed[2*k];
		fprintf(fp,"%g,%g,%g>,r}\n",x+0.5*pts[2*k],y+0.5*pts[2*k+1],z);
	} while (k!=0);
}

/** Computes the maximum radius squared of a vertex from the center of the
 * cell. It can be used to determine when enough particles have been testing an
 * all planes that could cut the cell have been considered.
 * \return The maximum radius squared of a vertex.*/
double voronoicell_2d::max_radius_squared() {
	double r,s,*ptsp(pts+2),*ptse(pts+2*p);
	r=*pts*(*pts)+pts[1]*pts[1];
	while(ptsp<ptse) {
		s=*ptsp*(*ptsp);ptsp++;
		s+=*ptsp*(*ptsp);ptsp++;
		if(s>r) r=s;
	}
	return r;
}

/** Cuts the Voronoi cell by a particle whose center is at a separation of
 * (x,y) from the cell center. The value of rsq should be initially set to
 * \f$x^2+y^2\f$.
 * \param[in] (x,y) the normal vector to the plane.
 * \param[in] rsq the distance along this vector of the plane.
 * \return False if the plane cut deleted the cell entirely, true otherwise. */
bool voronoicell_2d::plane(double x,double y,double rsq) {
	int cp,lp,up=0,up2,up3,*stackp(ds);
	double fac,l,u,u2,u3;

	// First try and find a vertex that is within the cutting plane, if
	// there is one. If one can't be found, then the cell is not cut by
	// this plane and the routine immediately returns true.
	u=pos(x,y,rsq,up);
	if(u<tolerance) {
		up2=ed[2*up];u2=pos(x,y,rsq,up2);
		up3=ed[2*up+1];u3=pos(x,y,rsq,up3);
		if(u2>u3) {
			while(u2<tolerance) {
				up2=ed[2*up2];
				u2=pos(x,y,rsq,up2);
				if(up2==up3) return true;
			}
			up=up2;u=u2;
		} else {
			while(u3<tolerance) {
				up3=ed[2*up3+1];
				u3=pos(x,y,rsq,up3);
				if(up2==up3) return true;
			}
			up=up3;u=u3;
		}
	}

	// Add this point to the delete stack, and search clockwise
	// to find additional points that need to be deleted.
	*(stackp++)=up;
	l=u;up2=ed[2*up];
	u2=pos(x,y,rsq,up2);
	while(u2>tolerance) {
		if(stackp==stacke) add_memory_ds(stackp);
		*(stackp++)=up2;
		up2=ed[2*up2];
		l=u2;
		u2=pos(x,y,rsq,up2);
		if(up2==up) return false;
	}

	// Consider the first point that was found in the clockwise direction
	// that was not inside the cutting plane. If it lies on the cutting
	// plane then do nothing. Otherwise, introduce a new vertex.
	if(u2>-tolerance) {
		cp=up2;
	} else {
		if(p==current_vertices) add_memory_vertices();
		lp=ed[2*up2+1];
		fac=1/(u2-l);
		pts[2*p]=(pts[2*lp]*u2-pts[2*up2]*l)*fac;
		pts[2*p+1]=(pts[2*lp+1]*u2-pts[2*up2+1]*l)*fac;
		ed[2*p]=up2;
		ed[2*up2+1]=p;
		cp=p++;
	}

	// Search counter-clockwise for additional points that need to be
	// deleted
	l=u;up3=ed[2*up+1];u3=pos(x,y,rsq,up3);
	while(u3>tolerance) {
		if(stackp==stacke) add_memory_ds(stackp);
		*(stackp++)=up3;
		up3=ed[2*up3+1];
		l=u3;
		u3=pos(x,y,rsq,up3);
		if(up3==up2) break;
	}

	// Either adjust the existing vertex or create new one, and connect it
	// with the vertex found on the previous search in the clockwise
	// direction
	if(u3>tolerance) {
		ed[2*cp+1]=up3;
		ed[2*up3]=cp;
	} else {
		if(p==current_vertices) add_memory_vertices();
		lp=ed[2*up3];
		fac=1/(u3-l);
		pts[2*p]=(pts[2*lp]*u3-pts[2*up3]*l)*fac;
		pts[2*p+1]=(pts[2*lp+1]*u3-pts[2*up3+1]*l)*fac;
		ed[2*p]=cp;
		ed[2*cp+1]=p;
		ed[2*p+1]=up3;
		ed[2*up3]=p++;
	}

	// Mark points on the delete stack
	for(int *sp=ds;sp<stackp;sp++) ed[*sp*2]=-1;

	// Remove them from the memory structure
	while(stackp>ds) {
		while(ed[2*--p]==-1);
		up=*(--stackp);
		if(up<p) {
			ed[2*ed[2*p]+1]=up;
			ed[2*ed[2*p+1]]=up;
			pts[2*up]=pts[2*p];
			pts[2*up+1]=pts[2*p+1];
			ed[2*up]=ed[2*p];
			ed[2*up+1]=ed[2*p+1];
		} else p++;
	}
	return true;
}

/** The same as Plane excepts it handles nonconvex cells, for which the handling of the vertices is more complex. */

bool voronoicell_2d::plane_nonconvex(double x,double y,double rsq) {

//first determine which region we are in.
	if(((x*reg1[0]+y*reg1[1])>0)&&((x*reg1[2]+y*reg1[3])>0)){
		return halfplane(x, y, rsq, reg1[2], reg1[3]);
	}if(((x*reg2[0]+y*reg2[1])>0)&&((x*reg2[2]+y*reg2[3])>0)){
		return halfplane(x, y, rsq, reg2[2], reg2[3]);
	
	}else{
		return plane(x,y,rsq);
	}  
}
/** cuts the given cell by a half plane given
 * \param[in] (x1 ,y1) the normal vector to the plane.
 * \param[in] rsq the distance along this vector of the plane.
 * \param[in] (x2, y2) the normal vector to the plane from which we want the plane cut to start.
 * \return False if the plane cut deleted the cell entirely, true otherwise. */

bool voronoicell_2d::halfplane(double x1, double y1, double rsq, double x2, double y2) {
	int si=0, ci=0, ni, patch1, patch2, *stackp(ds);
	double cid=pos(x1,y1,rsq,ci), nid, fac;
	bool rightchunk=false;
	if(debugging) cout << "beginning";
//first find a vertex that is not being cut by the plane
	while(cid>tolerance){
		ci=ed[2*ci+1];
		cid=pos(x1, y1, rsq, ci);
		if(ci==si){
			return false;
		}
	}
	si=ci;



if(debugging)	cout << "part1     ";



//now circle around the vertices until we find the right chunk to cut.
//When we exit this loop, ci will not be cut, and ni will be the first index to be cut by the plane.
//if nothing is cut, return true.
	while(!rightchunk){
		ni=ed[ci*2];
		if(ni==si){
			cout << "error1";
			return true;
		}
		nid=pos(x1,y1,rsq,ni);
		if(nid<tolerance || (((pts[2*ni]*x2)+(pts[2*ni+1]*y2))<tolerance)) {
			ci=ni;
			cid=nid;
			continue;
			
		}
		else{
			rightchunk=true;
		if(debugging){
			cout << pts[ni*2];
			cout << "     ";
			cout << pts[ni*2+1];
			cout << "     ";
		}
		}
	}



if(debugging) cout << "part2";



if(stackp==stacke) add_memory_ds(stackp);
*(stackp++)=ni;

	
//see if ci lies on the cutting plane,  if it doesnt introduce a new vertex
	if(cid>-tolerance){
		patch1=ci;
	}else{
		if(p==current_vertices) add_memory_vertices();
		fac=1/(cid-nid);
		pts[2*p]=(pts[2*ni]*cid-pts[2*ci]*nid)*fac;
		pts[2*p+1]=(pts[2*ni+1]*cid-pts[2*ci+1]*nid)*fac;
		patch1=p++;
		ed[2*ci]=patch1;
		ed[2*patch1+1]=ci;
	}




if(debugging) cout << "part3";



//continue around the cell clockwise deleting points until ci is being cut and ni isn't.
ci=ni;
cid=nid;
ni=ed[2*ci];
nid=pos(x1,y1,rsq,ni);
	while(nid>tolerance){
		if(stackp==stacke) add_memory_ds(stackp);
		*(stackp++)=ni;
		ci=ni;
		cid=nid;
		ni=ed[2*ci];
		nid=pos(x1,y1,rsq,ni);
	}	



if(debugging) cout << "part4";



//now ci is being cut and ni is not. If ni lies on the cutting plane, then do nothing, if it does not, introduce
//a new vertex.
	if(nid>-tolerance){
		patch2=ni;
		ed[2*patch1]=patch2;
		ed[2*patch2+1]=patch1;
	}else{
		if(p==current_vertices) add_memory_vertices();
		fac=1/(nid-cid);
		pts[2*p]=(pts[2*ci]*nid-pts[2*ni]*cid)*fac;
		pts[2*p+1]=(pts[2*ci+1]*nid-pts[2*ni+1]*cid)*fac;
		ed[2*p]=ni;
		ed[2*ni+1]=p;
		patch2=p++;
		ed[patch2*2+1]=patch1;
		ed[2*patch1]=patch2;
	}
	// Mark points on the delete stack
	for(int *sp=ds;sp<stackp;sp++) ed[*sp*2]=-1;

	// Remove them from the memory structure
	while(stackp>ds) {
		while(ed[2*--p]==-1);
		int up=*(--stackp);
		if(up<p) {
			ed[2*ed[2*p]+1]=up;
			ed[2*ed[2*p+1]]=up;
			pts[2*up]=pts[2*p];
			pts[2*up+1]=pts[2*p+1];
			ed[2*up]=ed[2*p];
			ed[2*up+1]=ed[2*p+1];
		} else p++;
	}

if(debugging) cout << "part5";

	return true;

}

bool voronoicell_2d::wallcut(double wx1,double wy1,double wx2,double wy2){
	double wox, woy, wpx, wpy, wpl,nl, pcx, pcy, rs;
	if((wx1==0 && wy1==0) || (wx2==0 && wy2==0)) return true; 
	wox=wx2-wx1; woy=wy2-wy1;
	wpx=-woy; wpy=wox;
	wpl=pow((pow(wpx,2.0)+pow(wpy,2.0)),0.5);
	wpx=wpx/wpl;
	wpy=wpy/wpl;
	nl=wx1*wpx+wy1*wpy;
	pcx=wpx*nl;//MULTIPLY BY 2?
	pcy=wpy*nl;//MULTIPLY BY 2?
	rs=pcx*pcx+pcy*pcy;
	this.plane(pcx,pcy,rs);	
		
	

	return true;
}


/** Calculates the perimeter of the Voronoi cell.
 * \return A floating point number holding the calculated distance. */
double voronoicell_2d::perimeter() {
	if(p==0) return 0;
	int k=0,l;double perim=0,dx,dy;
	do {
		l=ed[2*k];
		dx=pts[2*k]-pts[2*l];
		dy=pts[2*k+1]-pts[2*l+1];
		perim+=sqrt(dx*dx+dy*dy);
		k=l;
	} while (k!=0);
	return 0.5*perim;
}

/** Calculates the area of the Voronoi cell.
 * \return A floating point number holding the calculated distance. */
double voronoicell_2d::area() {
	if(p==0) return 0;
	int k(*ed);double area=0,x=*pts,y=pts[1],dx1,dy1,dx2,dy2;
	dx1=pts[2*k]-x;dy1=pts[2*k+1]-y;
	k=ed[2*k];
	while(k!=0) {
		dx2=pts[2*k]-x;dy2=pts[2*k+1]-y;
		area+=dx1*dy2-dx2*dy1;
		dx1=dx2;dy1=dy2;
		k=ed[2*k];
	}
	return 0.125*area;
}

/** Calculates the centroid of the Voronoi cell.
 * \param[out] (cx,cy) The coordinates of the centroid. */
void voronoicell_2d::centroid(double &cx,double &cy) {
	cx=cy=0;
	static const double third=1/3.0;
	if(p==0) return;
	int k(*ed);
	double area,tarea=0,x=*pts,y=pts[1],dx1,dy1,dx2,dy2;
	dx1=pts[2*k]-x;dy1=pts[2*k+1]-y;
	k=ed[2*k];
	while(k!=0) {
		dx2=pts[2*k]-x;dy2=pts[2*k+1]-y;
		area=dx1*dy2-dx2*dy1;
		tarea+=area;
		cx+=area*(dx1+dx2);
		cy+=area*(dy1+dy2);
		dx1=dx2;dy1=dy2;
		k=ed[2*k];
	}
	tarea=third/tarea;
	cx=0.5*(x+cx*tarea);
	cy=0.5*(y+cy*tarea);
}

/** Computes the Voronoi cells for all particles in the container, and for each
 * cell, outputs a line containing custom information about the cell structure.
 * The output format is specified using an input string with control sequences
 * similar to the standard C printf() routine.
 * \param[in] format the format of the output lines, using control sequences to
 *                   denote the different cell statistics.
 * \param[in] i the ID of the particle associated with this Voronoi cell.
 * \param[in] (x,y) the position of the particle associated with this Voronoi
 *                    cell.
 * \param[in] r a radius associated with the particle.
 * \param[in] fp the file handle to write to. */
void voronoicell_2d::output_custom(const char *format,int i,double x,double y,double r,FILE *fp) {
	char *fmp(const_cast<char*>(format));
	while(*fmp!=0) {
		if(*fmp=='%') {
			fmp++;
			switch(*fmp) {

				// Particle-related output
				case 'i': fprintf(fp,"%d",i);break;
				case 'x': fprintf(fp,"%g",x);break;
				case 'y': fprintf(fp,"%g",y);break;
				case 'q': fprintf(fp,"%g %g",x,y);break;
				case 'r': fprintf(fp,"%g",r);break;

				// Vertex-related output
				case 'w': fprintf(fp,"%d",p);break;
				case 'm': fprintf(fp,"%g",0.25*max_radius_squared());break;

				// Edge-related output
				case 'p': fprintf(fp,"%g",perimeter());break;

				// Area-related output
				case 'a': fprintf(fp,"%g",area());break;
				case 'c': {
						  double cx,cy;
						  centroid(cx,cy);
						  fprintf(fp,"%g %g",cx,cy);
					  } break;
				case 'C': {
						  double cx,cy;
						  centroid(cx,cy);
						  fprintf(fp,"%g %g",x+cx,y+cy);
					  } break;

				// End-of-string reached
				case 0: fmp--;break;

				// The percent sign is not part of a
				// control sequence
				default: putc('%',fp);putc(*fmp,fp);
			}
		} else putc(*fmp,fp);
		fmp++;
	}
	fputs("\n",fp);
}

/** Doubles the storage for the vertices, by reallocating the pts and ed
 * arrays. If the allocation exceeds the absolute maximum set in max_vertices,
 * then the routine exits with a fatal error. */
void voronoicell_2d::add_memory_vertices() {
	double *ppe(pts+2*current_vertices);
	int *ede(ed+2*current_vertices);

	// Double the memory allocation and check it is within range
	current_vertices<<=1;
	if(current_vertices>max_vertices) voropp_fatal_error("Vertex memory allocation exceeded absolute maximum",VOROPP_MEMORY_ERROR);
#if VOROPP_VERBOSE >=2
	fprintf(stderr,"Vertex memory scaled up to %d\n",current_vertices);
#endif

	// Copy the vertex positions
	double *npts(new double[2*current_vertices]),*npp(npts),*pp(pts);
	while(pp<ppe) *(npp++)=*(pp++);
	delete [] pts;pts=npts;

	// Copy the edge table
	int *ned(new int[2*current_vertices]),*nep(ned),*edp(ed);
	while(edp<ede) *(nep++)=*(edp++);
	delete [] ed;ed=ned;
}

/** Doubles the size allocation of the delete stack. If the allocation exceeds
 * the absolute maximum set in max_delete_size, then routine causes a fatal
 * error.
 * \param[in] stackp a reference to the current stack pointer. */
void voronoicell_2d::add_memory_ds(int *&stackp) {
	current_delete_size<<=1;
	if(current_delete_size>max_delete_size) voropp_fatal_error("Delete stack 1 memory allocation exceeded absolute maximum",VOROPP_MEMORY_ERROR);
#if VOROPP_VERBOSE >=2
	fprintf(stderr,"Delete stack 1 memory scaled up to %d\n",current_delete_size);
#endif
	int *dsn(new int[current_delete_size]),*dsnp(dsn),*dsp(ds);
	while(dsp<stackp) *(dsnp++)=*(dsp++);
	delete [] ds;ds=dsn;stackp=dsnp;
	stacke=ds+current_delete_size;
}
