// Voronoi calculation code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : February 27th 2008

#include "cell.hh"
#include "container.hh"

/** Container constructor. The first six arguments set the corners of the box to
 * be (xa,ya,za) and (xb,yb,zb). The box is then divided into an nx by ny by nz
 * grid of blocks, set by the following three arguments. The next three
 * arguments are booleans, which set the periodicity in each direction. The
 * final argument sets the amount of memory allocated to each block. */
container::container(fpoint xa,fpoint xb,fpoint ya,fpoint yb,fpoint za,fpoint zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi,int isz)
	: length_scale(1),sz(isz),ax(xa),bx(xb),ay(ya),by(yb),az(za),bz(zb),
	xsp(xn/(xb-xa)),ysp(yn/(yb-ya)),zsp(zn/(zb-za)),
	nx(xn),ny(yn),nz(zn),nxy(xn*yn),nxyz(xn*yn*zn),
	hx(xper?2*xn+1:xn),hy(yper?2*yn+1:yn),hz(zper?2*zn+1:zn),hxy(hx*hy),hxyz(hx*hy*hz),
	xperiodic(xper),yperiodic(yper),zperiodic(zper),
	mv(0),max_radius(0),wall_number(0) {
	int l;
	co=new int[nxyz];
	for(l=0;l<nxyz;l++) co[l]=0;
	mem=new int[nxyz];
	for(l=0;l<nxyz;l++) mem[l]=memi;
	id=new int*[nxyz];
	for(l=0;l<nxyz;l++) id[l]=new int[memi];
	p=new fpoint*[nxyz];
	for(l=0;l<nxyz;l++) p[l]=new fpoint[sz*memi];
	mask=new unsigned int[hxyz];
	for(l=0;l<hxyz;l++) mask[l]=0;
	s_size=30*(hxy+hz*(hx+hy));if (s_size<18) s_size=18;
	sl=new int[s_size];
	walls=new wall*[current_wall_size];
}

/** Container constructor. The first six arguments set the corners of the box to
 * be (xa,ya,za) and (xb,yb,zb). The box is then divided into an nx by ny by nz
 * grid of blocks, set by the following three arguments. The next three
 * arguments are booleans, which set the periodicity in each direction. The
 * final argument sets the amount of memory allocated to each block. */
container::container(fpoint xa,fpoint xb,fpoint ya,fpoint yb,fpoint za,fpoint zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi)
	: length_scale(1),sz(3),ax(xa),bx(xb),ay(ya),by(yb),az(za),bz(zb),
	xsp(xn/(xb-xa)),ysp(yn/(yb-ya)),zsp(zn/(zb-za)),
	nx(xn),ny(yn),nz(zn),nxy(xn*yn),nxyz(xn*yn*zn),
	hx(xper?2*xn+1:xn),hy(yper?2*yn+1:yn),hz(zper?2*zn+1:zn),hxy(hx*hy),hxyz(hx*hy*hz),
	xperiodic(xper),yperiodic(yper),zperiodic(zper),
	mv(0),max_radius(0),wall_number(0),current_wall_size(init_wall_size) {
	int l;
	co=new int[nxyz];
	for(l=0;l<nxyz;l++) co[l]=0;
	mem=new int[nxyz];
	for(l=0;l<nxyz;l++) mem[l]=memi;
	id=new int*[nxyz];
	for(l=0;l<nxyz;l++) id[l]=new int[memi];
	p=new fpoint*[nxyz];
	for(l=0;l<nxyz;l++) p[l]=new fpoint[sz*memi];
	mask=new unsigned int[hxyz];
	for(l=0;l<hxyz;l++) mask[l]=0;
	s_size=30*(hxy+hz*(hx+hy));if (s_size<18) s_size=18;
	sl=new int[s_size];
	walls=new wall*[current_wall_size];
}

/** Container destructor - free memory. */
container::~container() {
	int l;
	for(l=0;l<nxyz;l++) delete [] p[l];
	for(l=0;l<nxyz;l++) delete [] id[l];
	delete [] p;
	delete [] id;
	delete [] mem;
	delete [] co;	
}

/** Dumps all the particle positions and identifies to a file. */
void container::dump(char *filename) {
	int c,l,i;
	ofstream file;
	file.open(filename,ofstream::out|ofstream::trunc);
	for(l=0;l<nxyz;l++) {
		for (c=0;c<co[l];c++)
			file << id[l][c];
			for(i=sz*c;i<sz*(c+1);i++) file << " " << p[l][i];
			file << endl;
	}
	file.close();
}

/** Put a particle into the correct region of the container. */
void container::put(int n,fpoint x,fpoint y,fpoint z) {
	if(x>ax&&y>ay&&z>az) {
		int i,j,k;
		i=int((x-ax)*xsp);j=int((y-ay)*ysp);k=int((z-az)*zsp);
		if(i<nx&&j<ny&&k<nz) {
			i+=nx*j+nxy*k;
			if(co[i]==mem[i]) add_particle_memory(i);
			p[i][sz*co[i]]=x;p[i][sz*co[i]+1]=y;p[i][sz*co[i]+2]=z;
			id[i][co[i]++]=n;
		}
	}
}

/** Increase memory for a particular region. */
void container::add_particle_memory(int i) {
	int *idp;fpoint *pp;
	int l,nmem=2*mem[i];
	if (nmem>max_particle_memory) throw fatal_error("Absolute maximum memory allocation exceeded");
	idp=new int[nmem];
	for(l=0;l<co[i];l++) idp[l]=id[i][l];
	pp=new fpoint[sz*nmem];
	for(l=0;l<sz*co[i];l++) pp[l]=p[i][l];
	mem[i]=nmem;
	delete [] id[i];id[i]=idp;
	delete [] p[i];p[i]=pp;
}

/** Add list memory. */
inline void container::add_list_memory() {
	cout << "addmem\n";
	int i,j=0,*ps;
	s_size*=2;
	ps=new int[s_size];
	if(s_start<=s_end) {
		for(i=s_start;i<s_end;i++) ps[j++]=sl[i];
	} else {
		cout << s_start << " " << s_end << " " << s_size << endl;
		for(i=s_start;i<s_size;i++) ps[j++]=sl[i];
		for(i=0;i<s_end;i++) ps[j++]=sl[i];
	}
	s_start=0;s_end=j;
	delete [] sl;sl=ps;
}

/** Import a list of particles from standard input. */
void container::import(istream &is) {
	int n;fpoint x,y,z;
	is >> n >> x >> y >> z;
	while(!is.eof()) {
		put(n,x,y,z);
		is >> n >> x >> y >> z;
	}
}

/** An overloaded version of the import routine, that reads the standard input.
 */
inline void container::import() {
	import(cin);
}

/** An overloaded version of the import routine, that reads in particles from
 * a particular file.
 * \param[in] filename The name of the file to read from. */
inline void container::import(char *filename) {
	ifstream is;
	is.open(filename,ifstream::in);
	import(is);
	is.close();
}

/** Outputs the number of particles within each region. */
void container::region_count() {
	int i,j,k,ijk=0;
	for(k=0;k<nz;k++) {
		for(j=0;j<ny;j++) {
			for(i=0;i<nx;i++) cout << "Region (" << i << "," << j << "," << k << "): " << co[ijk++] << " particles" << endl;
		}
	}
}

/** Clears a container of particles. */
void container::clear() {
	for(int ijk=0;ijk<nxyz;ijk++) co[ijk]=0;
	max_radius=0;
}

/** Guess length scale by dividing the total volume of the container by the
 * total number of particles, and then taking the cube root. */
void container::guess_length_scale() {
	const fpoint third=1/3.0;
	int sp=0;
	for(int ijk=0;ijk<nxyz;ijk++) sp+=co[ijk];
	if (sp>0) {
		length_scale=(fpoint) sp;
		length_scale=pow(abs((bx-ax)*(by-ay)*(bz-az)/length_scale),third);
	}
}

/** Computes the Voronoi cells for all particles within a box with corners
 * (xmin,ymin,zmin) and (xmax,ymax,zmax), and saves the output in a format
 * that can be read by gnuplot. */
void container::draw_gnuplot(char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax) {
	fpoint x,y,z,px,py,pz;
	facets_loop l1(this);
	int q,s;
	voronoicell c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	s=l1.init(xmin,xmax,ymin,ymax,zmin,zmax,px,py,pz);
	do {
		for(q=0;q<co[s];q++) {
			x=p[s][sz*q]+px;y=p[s][sz*q+1]+py;z=p[s][sz*q+2]+pz;
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax&&z>zmin&&z<zmax) {
				compute_cell(c,l1.ip,l1.jp,l1.kp,s,q,x,y,z);
				c.dump_gnuplot(os,x,y,z);
			}
		}
	} while ((s=l1.inc(px,py,pz))!=-1);
	os.close();
}

/** If only a filename is supplied to draw_gnuplot(), then assume that we are
 * calculating the entire simulation region. */
void container::draw_gnuplot(char *filename) {
	draw_gnuplot(filename,ax,bx,ay,by,az,bz);
}

/** Computes the Voronoi cells for all particles within a box with corners
 * (xmin,ymin,zmin) and (xmax,ymax,zmax), and saves the output in a format
 * that can be read by gnuplot.*/
void container::draw_pov(char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax) {
	fpoint x,y,z,px,py,pz;
	facets_loop l1(this);
	int q,s;
	voronoicell c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	os << "#declare voronoi=union{\n";
	s=l1.init(xmin,xmax,ymin,ymax,zmin,zmax,px,py,pz);
	do {
		for(q=0;q<co[s];q++) {
			x=p[s][sz*q]+px;y=p[s][sz*q+1]+py;z=p[s][sz*q+2]+pz;
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax&&z>zmin&&z<zmax) {
				compute_cell(c,l1.ip,l1.jp,l1.kp,s,q,x,y,z);
				c.dump_pov(os,x,y,z);
			}
		}
	} while ((s=l1.inc(px,py,pz))!=-1);
	os << "}\n";
	os.close();
}

/** If only a filename is supplied to draw_pov(), then assume that we are
 * calculating the entire simulation region.*/
void container::draw_pov(char *filename) {
	draw_pov(filename,ax,bx,ay,by,az,bz);
}


/** Computes the Voronoi volumes for all the particles, and stores the
 * results according to the particle label in the fpoint array bb.*/
void container::store_cell_volumes(fpoint *bb) {
	voronoicell c;
	int i,j,k,ijk=0,q;
	for (k=0;k<nz;k++) {
		for(j=0;j<ny;j++) {
			for(i=0;i<nx;i++) {
				for(q=0;q<co[ijk];q++) {
					compute_cell(c,i,j,k,ijk,q);
					bb[id[ijk][q]]=c.volume();
				}
				ijk++;
			}
		}
	}
}

/** Prints a list of all particle labels, positions, and Voronoi volumes to the
 * standard output. */
template<class n_option>
inline void container::print_all(ostream &os,voronoicell_base<n_option> &c) {
	fpoint x,y,z;
	int i,j,k,ijk=0,q;
	for (k=0;k<nz;k++) {
		for(j=0;j<ny;j++) {
			for(i=0;i<nx;i++) {
				for(q=0;q<co[ijk];q++) {
					x=p[ijk][sz*q];y=p[ijk][sz*q+1];z=p[ijk][sz*q+2];
					compute_cell(c,i,j,k,ijk,q,x,y,z);
					os << id[ijk][q] << " " << x << " " << y << " " << z;
					if (sz==4) os << " " << p[ijk][4*q+3];
					os << " " << c.volume();
					c.neighbors(os);		
					os << endl;
				}
				ijk++;
			}
		}
	}
}

/** Prints a list of all particle labels, positions, and Voronoi volumes to the
 * standard output. */
void container::print_all(ostream &os) {
	voronoicell c;
	print_all(os,c);
}

/** An overloaded version of print_all(), which just prints to standard output. */
void container::print_all() {
	voronoicell c;
	print_all(cout);
}

/** An overloaded version of print_all(), which outputs the result to a particular
 * file.
 * \param[in] filename The name fo the file to write to. */
inline void container::print_all(char* filename) {
	voronoicell c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	print_all(os,c);
	os.close();
}

/** Prints a list of all particle labels, positions, Voronoi volumes, and a list
 * of neighboring particles to an output stream.
 * \param[in] os The output stream to print to.*/
void container::print_all_neighbor(ostream &os) {
	voronoicell_neighbor c;
	print_all(os,c);
}

/** An overloaded version of print_all_neighbor(), which just prints to
 * standard output. */
void container::print_all_neighbor() {
	voronoicell_neighbor c;
	print_all(cout,c);
}

/** An overloaded version of print_all_neighbor(), which outputs the result to a
 * particular file
 * \param[in] filename The name of the file to write to. */
inline void container::print_all_neighbor(char* filename) {
	voronoicell_neighbor c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	print_all(os,c);
	os.close();
}

/** Initialize the Voronoi cell to be the entire container. For non-periodic
 * coordinates, this is set by the position of the walls. For periodic
 * coordinates, the space is equally divided in either direction from the
 * particle's initial position. That makes sense since those boundaries would
 * be made by the neighboring periodic images of this particle. */
template<class n_option>
inline void container::initialize_voronoicell(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z) {
	float x1,x2,y1,y2,z1,z2;
	if (xperiodic) x1=-(x2=0.5*(bx-ax));else {x1=ax-x;x2=bx-x;}
	if (yperiodic) y1=-(y2=0.5*(by-ay));else {y1=ay-y;y2=by-y;}
	if (zperiodic) z1=-(z2=0.5*(bz-az));else {z1=az-z;z2=bz-z;}
	c.init(x1,x2,y1,y2,z1,z2);
	for(int j=0;j<wall_number;j++) walls[j]->cut_cell(c,x,y,z);
}

bool container::point_inside(fpoint x,fpoint y,fpoint z) {
	if(x<ax||x>bx||y<ay||y>by||z<az||z>bz) return false;
	return point_inside_walls(x,y,z);
}

bool container::point_inside_walls(fpoint x,fpoint y,fpoint z) {
	for(int j=0;j<wall_number;j++) if (!walls[j]->point_inside(x,y,z)) return false;
	return true;
}

/** Computes a single Voronoi cell in the container. This routine can be run by
 * the user, and it is also called multiple times by the functions vprintall,
 * store_cell_volumes() and draw(). */ 
template<class n_option>
void container::compute_cell_slow(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s,fpoint x,fpoint y,fpoint z) {
	fpoint x1,y1,z1,qx,qy,qz,lr=0,lrs=0,ur,urs,rs;
	int q,t;
	facets_loop l(this);
	initialize_voronoicell(c,x,y,z);
	// Now the cell is cut by testing neighboring particles in concentric
	// shells. Once the test shell becomes twice as large as the Voronoi
	// cell we can stop testing.
	//
	// TODO: This initial test, to figure out if this is a polydisperse
	// case or not, is ugly. It's probably best to recompile using
	// templates.
	if (sz==3) {
		while(lrs<c.maxradsq()) {
			ur=lr+0.5*length_scale;urs=ur*ur;
			t=l.init(x,y,z,ur,qx,qy,qz);
			do {
				for(q=0;q<co[t];q++) {
					x1=p[t][sz*q]+qx-x;y1=p[t][sz*q+1]+qy-y;z1=p[t][sz*q+2]+qz-z;
					rs=x1*x1+y1*y1+z1*z1;
					if (lrs-tolerance<rs&&rs<urs&&(q!=s||ijk!=t))
						c.nplane(x1,y1,z1,rs,id[t][q]);
				}
			} while ((t=l.inc(qx,qy,qz))!=-1);
			lr=ur;lrs=urs;
		}
	} else {
		fpoint crad=p[s][sz*i+3];
		const fpoint mul=1+(crad*crad-max_radius*max_radius)/((max_radius+crad)*(max_radius+crad));
		crad*=crad;
		while(lrs*mul<c.maxradsq()) {
			ur=lr+0.5*length_scale;urs=ur*ur;
			t=l.init(x,y,z,ur,qx,qy,qz);
			do {
				for(q=0;q<co[t];q++) {
					x1=p[t][sz*q]+qx-x;y1=p[t][sz*q+1]+qy-y;z1=p[t][sz*q+2]+qz-z;
					rs=x1*x1+y1*y1+z1*z1;
					if (lrs-tolerance<rs&&rs<urs&&(q!=s||ijk!=t))
						c.nplane(x1,y1,z1,rs+crad-p[t][sz*q+3]*p[t][sz*q+3],id[t][q]);
				}
			} while ((t=l.inc(qx,qy,qz))!=-1);
			lr=ur;lrs=urs;
		}
	}
}

/** A overloaded version of compute_cell_slow, that sets up the x, y, and z variables. */
template<class n_option>
inline void container::compute_cell_slow(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s) {
	fpoint x=p[ijk][sz*s],y=p[ijk][sz*s+1],z=p[ijk][sz*s+2];
	compute_cell_slow(c,i,j,k,ijk,s,x,y,z);
}

template<class n_option>
void container::compute_cell(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s,fpoint x,fpoint y,fpoint z) {
	const fpoint boxx=(bx-ax)/nx,boxy=(by-ay)/ny,boxz=(bz-az)/nz;
	fpoint x1,y1,z1,qx=0,qy=0,qz=0;
	fpoint xlo,ylo,zlo,xhi,yhi,zhi,rs;
	int ci,cj,ck,cijk,di,dj,dk,dijk,ei,ej,ek,eijk,q;
	// Initialize the Voronoi cell to fill the entire container
	initialize_voronoicell(c,x,y,z);
	int fuc=0;fpoint mrs;

	// Test all particles in the particle's local region first
	for(q=0;q<s;q++) {
		x1=p[ijk][sz*q]-x;
		y1=p[ijk][sz*q+1]-y;
		z1=p[ijk][sz*q+2]-z;
		rs=x1*x1+y1*y1+z1*z1;
		c.nplane(x1,y1,z1,rs,id[i][j]);
	}
	while(q<co[ijk]) {
		x1=p[ijk][sz*q]-x;
		y1=p[ijk][sz*q+1]-y;
		z1=p[ijk][sz*q+2]-z;
		rs=x1*x1+y1*y1+z1*z1;
		c.nplane(x1,y1,z1,rs,id[i][j]);
		q++;
	}

	mrs=c.maxradsq();
	// Update the mask counter, and if it has wrapped around, then
	// reset the mask
	mv++;
	if (mv==0) {
		for(q=0;q<hxyz;q++) mask[q]=0;
		mv=1;
	}
	s_start=s_end=0;

	ci=xperiodic?nx:i;
	cj=yperiodic?ny:j;
	ck=zperiodic?nz:k;
	fpoint fx=x-ax-boxx*(i-ci),fy=y-ay-boxy*(j-cj),fz=z-az-boxz*(k-ck);
	cijk=ci+hx*(cj+hy*ck);
	mask[cijk]=mv;

	if (ci>0) {mask[cijk-1]=mv;sl[s_end++]=ci-1;sl[s_end++]=cj;sl[s_end++]=ck;};
	if (cj>0) {mask[cijk-hx]=mv;sl[s_end++]=ci;sl[s_end++]=cj-1;sl[s_end++]=ck;};
	if (ck>0) {mask[cijk-hxy]=mv;sl[s_end++]=ci;sl[s_end++]=cj;sl[s_end++]=ck-1;};
	if (ci<hx-1) {mask[cijk+1]=mv;sl[s_end++]=ci+1;sl[s_end++]=cj;sl[s_end++]=ck;};
	if (cj<hy-1) {mask[cijk+hx]=mv;sl[s_end++]=ci;sl[s_end++]=cj+1;sl[s_end++]=ck;};
	if (ck<hz-1) {mask[cijk+hxy]=mv;sl[s_end++]=ci;sl[s_end++]=cj;sl[s_end++]=ck+1;};
	while(s_start!=s_end) {
		if (++fuc==4) {
			fuc=0;mrs=c.maxradsq();
		}
		if(s_start==s_size) s_start=0;
		di=sl[s_start++];dj=sl[s_start++];dk=sl[s_start++];
		xlo=di*boxx-fx;xhi=xlo+boxx;
		ylo=dj*boxy-fy;yhi=ylo+boxy;
		zlo=dk*boxz-fz;zhi=zlo+boxz;
		if(di>ci) {
			if(dj>cj) {
				if(dk>ck) {
					if (corner_test(c,xlo,ylo,zlo,xhi,yhi,zhi)) continue;
				} else if(dk<ck) {
					if (corner_test(c,xlo,ylo,zhi,xhi,yhi,zlo)) continue;
				} else {
					if (edge_z_test(c,xlo,ylo,zlo,xhi,yhi,zhi)) continue;
				}
			} else if(dj<cj) {
				if(dk>ck) {
					if (corner_test(c,xlo,yhi,zlo,xhi,ylo,zhi)) continue;
				} else if(dk<ck) {
					if (corner_test(c,xlo,yhi,zhi,xhi,ylo,zlo)) continue;
				} else {
					if (edge_z_test(c,xlo,yhi,zlo,xhi,ylo,zhi)) continue;
				}
			} else {
				if(dk>ck) {
					if (edge_y_test(c,xlo,ylo,zlo,xhi,yhi,zhi)) continue;
				} else if(dk<ck) {
					if (edge_y_test(c,xlo,ylo,zhi,xhi,yhi,zlo)) continue;
				} else {
					if (face_x_test(c,xlo,ylo,zlo,yhi,zhi)) continue;
				}
			}
		} else if(di<ci) {
			if(dj>cj) {
				if(dk>ck) {
					if (corner_test(c,xhi,ylo,zlo,xlo,yhi,zhi)) continue;
				} else if(dk<ck) {
					if (corner_test(c,xhi,ylo,zhi,xlo,yhi,zlo)) continue;
				} else {
					if (edge_z_test(c,xhi,ylo,zlo,xlo,yhi,zhi)) continue;
				}
			} else if(dj<cj) {
				if(dk>ck) {
					if (corner_test(c,xhi,yhi,zlo,xlo,ylo,zhi)) continue;
				} else if(dk<ck) {
					if (corner_test(c,xhi,yhi,zhi,xlo,ylo,zlo)) continue;
				} else {
					if (edge_z_test(c,xhi,yhi,zlo,xlo,ylo,zhi)) continue;
				}
			} else {
				if(dk>ck) {
					if (edge_y_test(c,xhi,ylo,zlo,xlo,yhi,zhi)) continue;
				} else if(dk<ck) {
					if (edge_y_test(c,xhi,ylo,zhi,xlo,yhi,zlo)) continue;
				} else {
					if (face_x_test(c,xhi,ylo,zlo,yhi,zhi)) continue;
				}
			}
		} else {
			if(dj>cj) {
				if(dk>ck) {
					if (edge_x_test(c,xlo,ylo,zlo,xhi,yhi,zhi)) continue;
				} else if(dk<ck) {
					if (edge_x_test(c,xlo,ylo,zhi,xhi,yhi,zlo)) continue;
				} else {
					if (face_y_test(c,xlo,ylo,zlo,xhi,zhi)) continue;
				}
			} else if(dj<cj) {
				if(dk>ck) {
					if (edge_x_test(c,xlo,yhi,zlo,xhi,ylo,zhi)) continue;
				} else if(dk<ck) {
					if (edge_x_test(c,xlo,yhi,zhi,xhi,ylo,zlo)) continue;
				} else {
					if (face_y_test(c,xlo,yhi,zlo,xhi,zhi)) continue;
				}
			} else {
				if(dk>ck) {
					if (face_z_test(c,xlo,ylo,zlo,xhi,yhi)) continue;
				} else if(dk<ck) {
					if (face_z_test(c,xlo,ylo,zhi,xhi,yhi)) continue;
				} else {
					cout << "error\n";
				}
			}
		}

		if(xperiodic) {ei=i+di-nx;if (ei<0) {qx=ax-bx;ei+=nx;} else if (ei>=nx) {qx=bx-ax;ei-=nx;} else qx=0;} else ei=di;
		if(yperiodic) {ej=j+dj-ny;if (ej<0) {qy=ay-by;ej+=ny;} else if (ej>=ny) {qy=by-ay;ej-=ny;} else qy=0;} else ej=dj;
		if(zperiodic) {ek=k+dk-nz;if (ek<0) {qz=az-bz;ek+=nz;} else if (ek>=nz) {qz=bz-az;ek-=nz;} else qz=0;} else ek=dk;

		cout << "newblock\n";
		eijk=ei+nx*(ej+ny*ek);
		for(q=0;q<co[eijk];q++) {
			x1=p[eijk][sz*q]+qx-x;
			y1=p[eijk][sz*q+1]+qy-y;
			z1=p[eijk][sz*q+2]+qz-z;
			rs=x1*x1+y1*y1+z1*z1;
			if (rs<mrs) {
				cout << x1 << y1 << z1;
			c.nplane(x1,y1,z1,rs,id[eijk][q]);}
		}

		if((s_start<=s_end?s_size-s_end+s_start:s_end-s_start)<18) add_list_memory();

		dijk=di+hx*(dj+hy*dk);
		if(di>0) if(mask[dijk-1]!=mv) {if(s_end==s_size) s_end=0;mask[dijk-1]=mv;sl[s_end++]=di-1;sl[s_end++]=dj;sl[s_end++]=dk;}
		if(dj>0) if(mask[dijk-hx]!=mv) {if(s_end==s_size) s_end=0;mask[dijk-hx]=mv;sl[s_end++]=di;sl[s_end++]=dj-1;sl[s_end++]=dk;}
		if(dk>0) if(mask[dijk-hxy]!=mv) {if(s_end==s_size) s_end=0;mask[dijk-hxy]=mv;sl[s_end++]=di;sl[s_end++]=dj;sl[s_end++]=dk-1;}
		if(di<hx-1) if(mask[dijk+1]!=mv) {if(s_end==s_size) s_end=0;mask[dijk+1]=mv;sl[s_end++]=di+1;sl[s_end++]=dj;sl[s_end++]=dk;}
		if(dj<hy-1) if(mask[dijk+hx]!=mv) {if(s_end==s_size) s_end=0;mask[dijk+hx]=mv;sl[s_end++]=di;sl[s_end++]=dj+1;sl[s_end++]=dk;}
		if(dk<hz-1) if(mask[dijk+hxy]!=mv) {if(s_end==s_size) s_end=0;mask[dijk+hxy]=mv;sl[s_end++]=di;sl[s_end++]=dj;sl[s_end++]=dk+1;}
	}
}

template<class n_option>
inline bool container::corner_test(voronoicell_base<n_option> &c,fpoint xl,fpoint yl,fpoint zl,fpoint xh,fpoint yh,fpoint zh) {
	if (c.plane_intersects_guess(xh,yl,zl,xl*xh+yl*yl+zl*zl)) return false;
	if (c.plane_intersects(xh,yh,zl,xl*xh+yl*yh+zl*zl)) return false;
	if (c.plane_intersects(xl,yh,zl,xl*xl+yl*yh+zl*zl)) return false;
	if (c.plane_intersects(xl,yh,zh,xl*xl+yl*yh+zl*zh)) return false;
	if (c.plane_intersects(xl,yl,zh,xl*xl+yl*yl+zl*zh)) return false;
	if (c.plane_intersects(xh,yl,zh,xl*xh+yl*yl+zl*zh)) return false;
	return true;
}

template<class n_option>
inline bool container::edge_x_test(voronoicell_base<n_option> &c,fpoint x0,fpoint yl,fpoint zl,fpoint x1,fpoint yh,fpoint zh) {
	if (c.plane_intersects_guess(x0,yl,zh,yl*yl+zl*zh)) return false;
	if (c.plane_intersects(x1,yl,zh,yl*yl+zl*zh)) return false;
	if (c.plane_intersects(x1,yl,zl,yl*yl+zl*zl)) return false;
	if (c.plane_intersects(x0,yl,zl,yl*yl+zl*zl)) return false;
	if (c.plane_intersects(x0,yh,zl,yl*yh+zl*zl)) return false;
	if (c.plane_intersects(x1,yh,zl,yl*yh+zl*zl)) return false;
	return true;
}

template<class n_option>
inline bool container::edge_y_test(voronoicell_base<n_option> &c,fpoint xl,fpoint y0,fpoint zl,fpoint xh,fpoint y1,fpoint zh) {
	if (c.plane_intersects_guess(xl,y0,zh,xl*xl+zl*zh)) return false;
	if (c.plane_intersects(xl,y1,zh,xl*xl+zl*zh)) return false;
	if (c.plane_intersects(xl,y1,zl,xl*xl+zl*zl)) return false;
	if (c.plane_intersects(xl,y0,zl,xl*xl+zl*zl)) return false;
	if (c.plane_intersects(xh,y0,zl,xl*xh+zl*zl)) return false;
	if (c.plane_intersects(xh,y1,zl,xl*xh+zl*zl)) return false;
	return true;
}

template<class n_option>
inline bool container::edge_z_test(voronoicell_base<n_option> &c,fpoint xl,fpoint yl,fpoint z0,fpoint xh,fpoint yh,fpoint z1) {
	if (c.plane_intersects_guess(xl,yh,z0,xl*xl+yl*yh)) return false;
	if (c.plane_intersects(xl,yh,z1,xl*xl+yl*yh)) return false;
	if (c.plane_intersects(xl,yl,z1,xl*xl+yl*yl)) return false;
	if (c.plane_intersects(xl,yl,z0,xl*xl+yl*yl)) return false;
	if (c.plane_intersects(xh,yl,z0,xl*xh+yl*yl)) return false;
	if (c.plane_intersects(xh,yl,z1,xl*xh+yl*yl)) return false;
	return true;
}

template<class n_option>
inline bool container::face_x_test(voronoicell_base<n_option> &c,fpoint xl,fpoint y0,fpoint z0,fpoint y1,fpoint z1) {
	if (c.plane_intersects_guess(xl,y0,z0,xl*xl)) return false;
	if (c.plane_intersects(xl,y0,z1,xl*xl)) return false;
	if (c.plane_intersects(xl,y1,z1,xl*xl)) return false;
	if (c.plane_intersects(xl,y1,z0,xl*xl)) return false;
	return true;
}

template<class n_option>
inline bool container::face_y_test(voronoicell_base<n_option> &c,fpoint x0,fpoint yl,fpoint z0,fpoint x1,fpoint z1) {
	if (c.plane_intersects_guess(x0,yl,z0,yl*yl)) return false;
	if (c.plane_intersects(x0,yl,z1,yl*yl)) return false;
	if (c.plane_intersects(x1,yl,z1,yl*yl)) return false;
	if (c.plane_intersects(x1,yl,z0,yl*yl)) return false;
	return true;
}

template<class n_option>
inline bool container::face_z_test(voronoicell_base<n_option> &c,fpoint x0,fpoint y0,fpoint zl,fpoint x1,fpoint y1) {
	if (c.plane_intersects_guess(x0,y0,zl,zl*zl)) return false;
	if (c.plane_intersects(x0,y1,zl,zl*zl)) return false;
	if (c.plane_intersects(x1,y1,zl,zl*zl)) return false;
	if (c.plane_intersects(x1,y0,zl,zl*zl)) return false;
	return true;
}

/** A overloaded version of compute_cell, that sets up the x, y, and z variables. */
template<class n_option>
inline void container::compute_cell(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s) {
	fpoint x=p[ijk][sz*s],y=p[ijk][sz*s+1],z=p[ijk][sz*s+2];
	compute_cell(c,i,j,k,ijk,s,x,y,z);
}

/** Creates a facets_loop object, by pulling the necesssary constants about the container
 * geometry from a pointer to the current container class. */
facets_loop::facets_loop(container *q) : sx(q->bx-q->ax), sy(q->by-q->ay), sz(q->bz-q->az),
	xsp(q->xsp),ysp(q->ysp),zsp(q->zsp),
	ax(q->ax),ay(q->ay),az(q->az),
	nx(q->nx),ny(q->ny),nz(q->nz),nxy(q->nxy),nxyz(q->nxyz),
	xperiodic(q->xperiodic),yperiodic(q->yperiodic),zperiodic(q->zperiodic) {}

/** Initializes a facets_loop object, by finding all blocks which are within a distance
 * r of the vector (vx,vy,vz). It returns the first block which is to be
 * tested, and sets the periodic displacement vector (px,py,pz) accordingly. */
inline int facets_loop::init(fpoint vx,fpoint vy,fpoint vz,fpoint r,fpoint &px,fpoint &py,fpoint &pz) {
	ai=step_int((vx-ax-r)*xsp);
	bi=step_int((vx-ax+r)*xsp);
	if (!xperiodic) {
		if(ai<0) {ai=0;if (bi<0) bi=0;}
		if(bi>=nx) {bi=nx-1;if (ai>=nx) ai=nx-1;} 
	}
	aj=step_int((vy-ay-r)*ysp);
	bj=step_int((vy-ay+r)*ysp);
	if (!yperiodic) {
		if(aj<0) {aj=0;if (bj<0) bj=0;}
		if(bj>=ny) {bj=ny-1;if (aj>=ny) aj=ny-1;} 
	}
	ak=step_int((vz-az-r)*zsp);
	bk=step_int((vz-az+r)*zsp);
	if (!zperiodic) {
		if(ak<0) {ak=0;if (bk<0) bk=0;}
		if(bk>=nz) {bk=nz-1;if (ak>=nz) ak=nz-1;} 
	}
	i=ai;j=aj;k=ak;
	aip=ip=step_mod(i,nx);apx=px=step_div(i,nx)*sx;
	ajp=jp=step_mod(j,ny);apy=py=step_div(j,ny)*sy;
	akp=kp=step_mod(k,nz);apz=pz=step_div(k,nz)*sz;
	inc1=aip-step_mod(bi,nx);
	inc2=nx*(ny+ajp-step_mod(bj,ny))+inc1;
	inc1+=nx;
	s=aip+nx*(ajp+ny*akp);
	return s;
}

/** Initializes a facets_loop object, by finding all blocks which overlap the box with
 * corners (xmin,ymin,zmin) and (xmax,ymax,zmax). It returns the first block
 * which is to be tested, and sets the periodic displacement vector (px,py,pz)
 * accordingly. */
inline int facets_loop::init(fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax,fpoint &px,fpoint &py,fpoint &pz) {
	ai=step_int((xmin-ax)*xsp);
	bi=step_int((xmax-ax)*xsp);
	if (!xperiodic) {
		if(ai<0) {ai=0;if (bi<0) bi=0;}
		if(bi>=nx) {bi=nx-1;if (ai>=nx) ai=nx-1;} 
	}
	aj=step_int((ymin-ay)*ysp);
	bj=step_int((ymax-ay)*ysp);
	if (!yperiodic) {
		if(aj<0) {aj=0;if (bj<0) bj=0;}
		if(bj>=ny) {bj=ny-1;if (aj>=ny) aj=ny-1;} 
	}
	ak=step_int((zmin-az)*zsp);
	bk=step_int((zmax-az)*zsp);
	if (!zperiodic) {
		if(ak<0) {ak=0;if (bk<0) bk=0;}
		if(bk>=nz) {bk=nz-1;if (ak>=nz) ak=nz-1;} 
	}
	i=ai;j=aj;k=ak;
	aip=ip=step_mod(i,nx);apx=px=step_div(i,nx)*sx;
	ajp=jp=step_mod(j,ny);apy=py=step_div(j,ny)*sy;
	akp=kp=step_mod(k,nz);apz=pz=step_div(k,nz)*sz;
	inc1=aip-step_mod(bi,nx);
	inc2=nx*(ny+ajp-step_mod(bj,ny))+inc1;
	inc1+=nx;
	s=aip+nx*(ajp+ny*akp);
	return s;
}

/** Returns the next block to be tested in a loop, and updates the periodicity
 * vector if necessary. */
inline int facets_loop::inc(fpoint &px,fpoint &py,fpoint &pz) {
	if (i<bi) {
		i++;
		if (ip<nx-1) {ip++;s++;} else {ip=0;s+=1-nx;px+=sx;}
		return s;
	} else if (j<bj) {
		i=ai;ip=aip;px=apx;j++;
		if (jp<ny-1) {jp++;s+=inc1;} else {jp=0;s+=inc1-nxy;py+=sy;}
		return s;
	} else if (k<bk) {
		i=ai;ip=aip;j=aj;jp=ajp;px=apx;py=apy;k++;
		if (kp<nz-1) {kp++;s+=inc2;} else {kp=0;s+=inc2-nxyz;pz+=sz;}
		return s;
	} else return -1;
}

/** Custom int function, that gives consistent stepping for negative numbers.
 * With normal int, we have (-1.5,-0.5,0.5,1.5) -> (-1,0,0,1).
 * With this routine, we have (-1.5,-0.5,0.5,1.5) -> (-2,-1,0,1).*/
inline int facets_loop::step_int(fpoint a) {
	return a<0?int(a)-1:int(a);
}

/** Custom mod function, that gives consistent stepping for negative numbers. */
inline int facets_loop::step_mod(int a,int b) {
	return a>=0?a%b:b-1-(b-1-a)%b;
}

/** Custom div function, that gives consistent stepping for negative numbers. */
inline int facets_loop::step_div(int a,int b) {
	return a>=0?a/b:-1+(a+1)/b;
}

/** Put a particle into the correct region of the container.
 * \param[in] n The numerical ID of the inserted particle.
 * \param[in] (x,y,z) The position vector of the inserted particle.
 * \param[in] r The radius of the particle.*/
void container_poly::put(int n,fpoint x,fpoint y,fpoint z,fpoint r) {
	if(x>ax&&y>ay&&z>az) {
		int i,j,k;
		i=int((x-ax)*xsp);j=int((y-ay)*ysp);k=int((z-az)*zsp);
		if(i<nx&&j<ny&&k<nz) {
			i+=nx*j+nxy*k;
			if(co[i]==mem[i]) add_particle_memory(i);
			p[i][sz*co[i]]=x;p[i][sz*co[i]+1]=y;p[i][sz*co[i]+2]=z;p[i][sz*co[i]+3]=r;
			if (r>max_radius) max_radius=r;
			id[i][co[i]++]=n;
		}
	}
}

/** If the radius argument is not supplied to the polydisperse put() routine
 * then assume that we're using a particle with a unit diameter.
 * \param[in] n The numerical ID of the inserted particle.
 * \param[in] (x,y,z) The position vector of the inserted particle. */
void container_poly::put(int n,fpoint x,fpoint y,fpoint z) {
	put(n,x,y,z,0.5);
}

/** Import a list of particles.
 * \param[in] is An open input stream to read from. */
void container_poly::import(istream &is) {
	int n;fpoint x,y,z;
	fpoint r;
	is >> n >> x >> y >> z >> r;
	while(!is.eof()) {
		put(n,x,y,z,r);
		is >> n >> x >> y >> z >> r;
	}
}

/** An overloaded version of the import routine, that reads the standard input.
 */
inline void container_poly::import() {
	import(cin);
}

/** An overloaded version of the import routine, that reads in particles from
 * a particular file.
 * \param[in] filename The name of the file to open and import. */
inline void container_poly::import(char *filename) {
	ifstream is;
	is.open(filename,ifstream::in);
	import(is);
	is.close();
}

/** Adds a wall to the container. */
void container::add_wall(wall& w) {
	if (wall_number==current_wall_size) {
		current_wall_size*=2;
		if (current_wall_size>max_wall_size) throw fatal_error("Wall memory allocation exceeded absolute maximum");
		wall **pwall;
		pwall=new wall*[current_wall_size];
		for(int i=0;i<wall_number;i++) pwall[i]=walls[i];
		delete [] walls;walls=pwall;		
	}
	walls[wall_number++]=&w;
}
