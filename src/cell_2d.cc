// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file cell_2d.cc
 * \brief Function implementations for the voronoicell_2d class. */

#include "cell_2d.hh"

#include <cmath>
#include <cstring>

namespace voro {

/** Constructs a 2D Voronoi cell and sets up the initial memory.
 * \param[in] max_len_sq a maximum square length scale to set the internal
 *                   tolerance for floating point vertex calculations. */
voronoicell_base_2d::voronoicell_base_2d(double max_len_sq) :
    current_vertices(init_vertices), current_delete_size(init_delete_size),
    tol(tolerance*max_len_sq), ed(new int[2*current_vertices]),
    pts(new double[2*current_vertices]), ds(new int[current_delete_size]),
    stacke(ds+current_delete_size) {}

/** Constructs a 2D Voronoi cell by copying the constants from another.
 * \param[in] c a pointer to the Voronoi cell information to copy. */
voronoicell_base_2d::voronoicell_base_2d(voronoicell_base_2d *cp) :
    current_vertices(cp->current_vertices),
    current_delete_size(init_delete_size), p(cp->p), tol(cp->tol),
    ed(new int[2*current_vertices]), pts(new double[2*current_vertices]),
    ds(new int[current_delete_size]), stacke(ds+current_delete_size) {

    // Copy the vertex information across
    memcpy(ed,cp->ed,2*sizeof(int)*p);
    memcpy(pts,cp->pts,2*sizeof(double)*p);
}

/** The voronoicell_2d destructor deallocates all of the dynamic memory. */
voronoicell_base_2d::~voronoicell_base_2d() {
    delete [] ds;
    delete [] pts;
    delete [] ed;
}

/** Initializes a Voronoi cell as a rectangle with the given dimensions.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates. */
void voronoicell_base_2d::init_base(double xmin,double xmax,double ymin,double ymax) {
    p=4;xmin*=2;xmax*=2;ymin*=2;ymax*=2;
    *pts=xmin;pts[1]=ymin;
    pts[2]=xmax;pts[3]=ymin;
    pts[4]=xmax;pts[5]=ymax;
    pts[6]=xmin;pts[7]=ymax;
    int *q=ed;
    *q=1;q[1]=3;q[2]=2;q[3]=0;q[4]=3;q[5]=1;q[6]=0;q[7]=2;
}

/** Outputs the edges of the Voronoi cell in Gnuplot format to an output
 * stream.
 * \param[in] (x,y) a displacement vector to be added to the cell's position.
 * \param[in] fp the file handle to write to. */
void voronoicell_base_2d::draw_gnuplot(double x,double y,FILE *fp) {
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
 * \param[in] (x,y) a displacement vector to be added to the cell's position.
 * \param[in] fp the file handle to write to. */
void voronoicell_base_2d::draw_pov(double x,double y,FILE *fp) {
    if(p==0) return;
    int k=0;
    do {
        fprintf(fp,"sphere{<%g,%g,0>,r}\ncylinder{<%g,%g,0>,<",
                x+0.5*pts[2*k],y+0.5*pts[2*k+1],
                x+0.5*pts[2*k],y+0.5*pts[2*k+1]);
        k=ed[2*k];
        fprintf(fp,"%g,%g,0>,r}\n",x+0.5*pts[2*k],y+0.5*pts[2*k+1]);
    } while (k!=0);
}

/** Computes the maximum radius squared of a vertex from the center of the
 * cell. It can be used to determine when enough particles have been testing an
 * all planes that could cut the cell have been considered.
 * \return The maximum radius squared of a vertex.*/
double voronoicell_base_2d::max_radius_squared() {
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
bool voronoicell_base_2d::plane_intersects(double x,double y,double rsq) {
    int up=0,up2,up3;
    double u,u2,u3;

    // First try and find a vertex that is within the cutting plane, if there
    // is one. If one can't be found, then the cell is not cut by this plane
    // and the routine immediately returns true.
    u=pos(x,y,rsq,up);
    if(u<tol) {
        up2=ed[2*up];u2=pos(x,y,rsq,up2);
        up3=ed[2*up+1];u3=pos(x,y,rsq,up3);
        if(u2>u3) {
            while(u2<tol) {
                up2=ed[2*up2];
                u2=pos(x,y,rsq,up2);
                if(up2==up3) return false;
            }
            up=up2;u=u2;
        } else {
            while(u3<tol) {
                up3=ed[2*up3+1];
                u3=pos(x,y,rsq,up3);
                if(up2==up3) return false;
            }
            up=up3;u=u3;
        }
    }
    return true;
}

/** Cuts the Voronoi cell by a particle whose center is at a separation of
 * (x,y) from the cell center. The value of rsq should be initially set to
 * \f$x^2+y^2\f$.
 * \param[in] (x,y) the normal vector to the plane.
 * \param[in] rsq the distance along this vector of the plane.
 * \return False if the plane cut deleted the cell entirely, true otherwise. */
template<class vc_class>
bool voronoicell_base_2d::nplane(vc_class &vc,double x,double y,double rsq,int p_id) {
    int up=0,up2,up3;
    double u,u2,u3;

    // First try and find a vertex that is within the cutting plane, if there
    // is one. If one can't be found, then the cell is not cut by this plane
    // and the routine immediately returns true.
    u=pos(x,y,rsq,up);
    if(u<tol) {
        up2=ed[2*up];u2=pos(x,y,rsq,up2);
        up3=ed[2*up+1];u3=pos(x,y,rsq,up3);
        if(u2>u3) {
            while(u2<tol) {
                up2=ed[2*up2];
                u2=pos(x,y,rsq,up2);
                if(up2==up3) return true;
            }
            up=up2;u=u2;
        } else {
            while(u3<tol) {
                up3=ed[2*up3+1];
                u3=pos(x,y,rsq,up3);
                if(up2==up3) return true;
            }
            up=up3;u=u3;
        }
    }

    return nplane_cut(vc,x,y,rsq,p_id,u,up);
}

template<class vc_class>
bool voronoicell_base_2d::nplane_cut(vc_class &vc,double x,double y,double rsq,int p_id,double u,int up) {
    int cp,lp,up2,up3,*stackp=ds;
    double fac,l,u2,u3;

    // Add this point to the delete stack, and search counter-clockwise to find
    // additional points that need to be deleted.
    *(stackp++)=up;
    l=u;up2=ed[2*up];
    u2=pos(x,y,rsq,up2);
    while(u2>tol) {
        if(stackp==stacke) add_memory_ds(stackp);
        *(stackp++)=up2;
        up2=ed[2*up2];
        l=u2;
        u2=pos(x,y,rsq,up2);
        if(up2==up) return false;
    }

    // Consider the first point that was found in the counter-clockwise
    // direction that was not inside the cutting plane. If it lies on the
    // cutting plane then do nothing. Otherwise, introduce a new vertex.
    if(u2>-tol) cp=up2;
    else {
        if(p==current_vertices) add_memory_vertices(vc);
        lp=ed[2*up2+1];
        fac=1/(u2-l);
        pts[2*p]=(pts[2*lp]*u2-pts[2*up2]*l)*fac;
        pts[2*p+1]=(pts[2*lp+1]*u2-pts[2*up2+1]*l)*fac;
        vc.n_copy(p,lp);
        ed[2*p]=up2;
        ed[2*up2+1]=p;
        cp=p++;
    }

    // Search clockwise for additional points that need to be deleted
    l=u;up3=ed[2*up+1];u3=pos(x,y,rsq,up3);
    while(u3>tol) {
        if(stackp==stacke) add_memory_ds(stackp);
        *(stackp++)=up3;
        up3=ed[2*up3+1];
        l=u3;
        u3=pos(x,y,rsq,up3);
        if(up3==up2) break;
    }

    // Either adjust the existing vertex or create new one, and connect it with
    // the vertex found on the previous search in the counter-clockwise
    // direction
    if(u3>-tol) {
        ed[2*cp+1]=up3;
        ed[2*up3]=cp;
        vc.n_set(up3,p_id);
    } else {
        if(p==current_vertices) add_memory_vertices(vc);
        lp=ed[2*up3];
        fac=1/(u3-l);
        pts[2*p]=(pts[2*lp]*u3-pts[2*up3]*l)*fac;
        pts[2*p+1]=(pts[2*lp+1]*u3-pts[2*up3+1]*l)*fac;
        ed[2*p]=cp;
        ed[2*cp+1]=p;
        vc.n_set(p,p_id);
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
            vc.n_copy(up,p);
            ed[2*up]=ed[2*p];
            ed[2*up+1]=ed[2*p+1];
        } else p++;
    }
    return p>=3;
}

/** Returns a vector of the vertex vectors using the local coordinate system.
 * \param[out] v the vector to store the results in. */
void voronoicell_base_2d::vertices(std::vector<double> &v) {
    v.resize(2*p);
    double *ptsp=pts;
    for(int i=0;i<2*p;i+=2) {
        v[i]=*(ptsp++)*0.5;
        v[i+1]=*(ptsp++)*0.5;
    }
}

/** Outputs the vertex vectors using the local coordinate system.
 * \param[out] fp the file handle to write to. */
void voronoicell_base_2d::output_vertices(FILE *fp) {
    if(p>0) {
        fprintf(fp,"(%g,%g)",*pts*0.5,pts[1]*0.5);
        for(double *ptsp=pts+2;ptsp<pts+2*p;ptsp+=2) fprintf(fp," (%g,%g)",*ptsp*0.5,ptsp[1]*0.5);
    }
}

/** Outputs the vertex vectors using the local coordinate system.
 * \param[in] pr the precision to output the numbers to.
 * \param[out] fp the file handle to write to. */
void voronoicell_base_2d::output_vertices(int pr,FILE *fp) {
    if(p>0) {
        fprintf(fp,"(%.*g,%.*g)",pr,*pts*0.5,pr,pts[1]*0.5);
        for(double *ptsp=pts+2;ptsp<pts+2*p;ptsp+=2) fprintf(fp," (%.*g,%.*g)",pr,*ptsp*0.5,pr,ptsp[1]*0.5);
    }
}

/** Returns a vector of the vertex vectors in the global coordinate system.
 * \param[out] v the vector to store the results in.
 * \param[in] (x,y,z) the position vector of the particle in the global
 *                    coordinate system. */
void voronoicell_base_2d::vertices(double x,double y,std::vector<double> &v) {
    v.resize(2*p);
    double *ptsp=pts;
    for(int i=0;i<2*p;i+=2) {
        v[i]=x+*(ptsp++)*0.5;
        v[i+1]=y+*(ptsp++)*0.5;
    }
}

/** Outputs the vertex vectors using the global coordinate system.
 * \param[out] fp the file handle to write to.
 * \param[in] (x,y,z) the position vector of the particle in the global
 *                    coordinate system. */
void voronoicell_base_2d::output_vertices(double x,double y,FILE *fp) {
    if(p>0) {
        fprintf(fp,"(%g,%g)",x+*pts*0.5,y+pts[1]*0.5);
        for(double *ptsp=pts+2;ptsp<pts+2*p;ptsp+=2) fprintf(fp," (%g,%g)",x+*ptsp*0.5,y+ptsp[1]*0.5);
    }
}

/** Outputs the vertex vectors using the global coordinate system.
 * \param[in] pr the precision to output the numbers to.
 * \param[out] fp the file handle to write to.
 * \param[in] (x,y,z) the position vector of the particle in the global
 *                    coordinate system. */
void voronoicell_base_2d::output_vertices(int pr,double x,double y,FILE *fp) {
    if(p>0) {
        fprintf(fp,"(%.*g,%.*g)",pr,x+*pts*0.5,pr,y+pts[1]*0.5);
        for(double *ptsp=pts+2;ptsp<pts+2*p;ptsp+=2) fprintf(fp," (%.*g,%.*g)",pr,x+*ptsp*0.5,pr,y+ptsp[1]*0.5);
    }
}

/** Calculates the perimeter of the Voronoi cell.
 * \return A floating point number holding the calculated distance. */
double voronoicell_base_2d::perimeter() {
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

/** Calculates the perimeter of the Voronoi cell.
 * \return A floating point number holding the calculated distance. */
void voronoicell_base_2d::edge_lengths(std::vector<double> &vd) {
    if(p==0) {vd.clear();return;}
    vd.resize(p);
    int i=0,k=0,l;double dx,dy;
    do {
        l=ed[2*k];
        dx=pts[2*k]-pts[2*l];
        dy=pts[2*k+1]-pts[2*l+1];
        vd[i++]=sqrt(dx*dx+dy*dy);
        k=l;
    } while (k!=0);
}

/** Calculates the perimeter of the Voronoi cell.
 * \return A floating point number holding the calculated distance. */
void voronoicell_base_2d::normals(std::vector<double> &vd) {
    if(p==0) {vd.clear();return;}
    vd.resize(2*p);
    int i=0,k=0,l;double dx,dy,nor;
    do {
        l=ed[2*k];
        dx=pts[2*k]-pts[2*l];
        dy=pts[2*k+1]-pts[2*l+1];
        nor=dx*dx+dy*dy;
        if(nor>tol) {
            nor=1/sqrt(nor);
            vd[i++]=dy*nor;
            vd[i++]=-dx*nor;
        } else {
            vd[i++]=0;
            vd[i++]=0;
        }
        k=l;
    } while (k!=0);
}

/** Calculates the area of the Voronoi cell.
 * \return A floating point number holding the calculated distance. */
double voronoicell_base_2d::area() {
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
void voronoicell_base_2d::centroid(double &cx,double &cy) {
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
    if(fabs(tarea)>tol) {
        tarea=third/tarea;
        cx=0.5*(x+cx*tarea);
        cy=0.5*(y+cy*tarea);
    } else cx=cy=0;
}

/** Computes the Voronoi cells for all particles in the container, and for each
 * cell, outputs a line containing custom information about the cell structure.
 * The output format is specified using an input string with control sequences
 * similar to the standard C printf() routine.
 * \param[in] format the format of the output lines, using control sequences to
 *                   denote the different cell statistics.
 * \param[in] i the ID of the particle associated with this Voronoi cell.
 * \param[in] (x,y) the position of the particle associated with this Voronoi
 *                  cell.
 * \param[in] r a radius associated with the particle.
 * \param[in] fp the file handle to write to. */
void voronoicell_base_2d::output_custom(const char *format,int i,double x,double y,double r,FILE *fp) {
    char *fmp(const_cast<char*>(format));
    std::vector<int> vi;
    std::vector<double> vd;
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
                case 'p': output_vertices(fp);break;
                case 'P': output_vertices(x,y,fp);break;
                case 'm': fprintf(fp,"%g",0.25*max_radius_squared());break;

                // Edge-related output
                case 'g': fprintf(fp,"%d",p);break;
                case 'E': fprintf(fp,"%g",perimeter());break;
                case 'e': edge_lengths(vd);voro_print_vector(vd,fp);break;
                case 'l': normals(vd);
                      voro_print_positions_2d(vd,fp);
                      break;
                case 'n': neighbors(vi);
                      voro_print_vector(vi,fp);
                      break;

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

                // Precision is specified
                case '.': {
                        int pr;
                        if(voro_read_precision(fp,fmp,pr)) switch(*fmp) {

                            // Particle-related output
                            case 'x': fprintf(fp,"%.*g",pr,x);break;
                            case 'y': fprintf(fp,"%.*g",pr,y);break;
                            case 'q': fprintf(fp,"%.*g %.*g",pr,x,pr,y);break;
                            case 'r': fprintf(fp,"%g",r);break;

                            // Vertex-related output
                            case 'p': output_vertices(pr,fp);break;
                            case 'P': output_vertices(pr,x,y,fp);break;
                            case 'm': fprintf(fp,"%.*g",pr,0.25*max_radius_squared());break;

                            // Edge-related output
                            case 'E': fprintf(fp,"%.*g",pr,perimeter());break;
                            case 'e': edge_lengths(vd);voro_print_vector(pr,vd,fp);break;
                            case 'l': normals(vd);
                                  voro_print_positions_2d(pr,vd,fp);
                                  break;

                            // Area-related output
                            case 'a': fprintf(fp,"%.*g",pr,area());break;
                            case 'c': {
                                      double cx,cy;
                                      centroid(cx,cy);
                                      fprintf(fp,"%.*g %.*g",pr,cx,pr,cy);
                                  } break;
                            case 'C': {
                                      double cx,cy;
                                      centroid(cx,cy);
                                      fprintf(fp,"%.*g %.*g",pr,x+cx,pr,y+cy);
                                  } break;

                            // End-of-string reached
                            case 0: fprintf(fp,"%%.%d",pr);fmp--;break;

                            // This is not a valid control sequence
                            default: fprintf(fp,"%%.%d%c",pr,*fmp);
                          }
                      } break;

                // End-of-string reached
                case 0: putc('%',fp);fmp--;break;

                // The percent sign is not part of a
                // control sequence
                default: putc('%',fp);putc(*fmp,fp);
            }
        } else putc(*fmp,fp);
        fmp++;
    }
    putc('\n',fp);
}

/** Doubles the storage for the vertices, by reallocating the pts and ed
 * arrays. If the allocation exceeds the absolute maximum set in max_vertices,
 * then the routine exits with a fatal error. */
template<class vc_class>
void voronoicell_base_2d::add_memory_vertices(vc_class &vc) {
    int ocv=current_vertices;

    // Double the memory allocation and check it is within range
    current_vertices<<=1;
    if(current_vertices>max_vertices) voro_fatal_error("Vertex memory allocation exceeded absolute maximum",VOROPP_MEMORY_ERROR);
#if VOROPP_VERBOSE >=2
    fprintf(stderr,"Vertex memory scaled up to %d\n",current_vertices);
#endif

    // Copy the vertex positions
    double *npts=new double[2*current_vertices];
    memcpy(npts,pts,sizeof(double)*2*ocv);
    delete [] pts;pts=npts;

    // Copy the edge table
    int *ned=new int[2*current_vertices];
    memcpy(ned,ed,sizeof(int)*2*ocv);
    delete [] ed;ed=ned;

    // Double the neighbor information if necessary
    vc.n_add_memory_vertices(ocv);
}

/** Constructs the 2D Voronoi cell by copying the information from an existing one.
 * \param[in] c_ a reference to the class to copy from. */
voronoicell_neighbor_2d::voronoicell_neighbor_2d(voronoicell_neighbor_2d &c_) :
	voronoicell_base_2d((voronoicell_base_2d*) &c_), ne(new int[current_vertices]) {
    memcpy(ne,c_.ne,sizeof(int)*p);
}

inline void voronoicell_neighbor_2d::n_add_memory_vertices(int ocv) {
    int *nne=new int[current_vertices];
    memcpy(nne,ne,sizeof(int)*ocv);
    delete [] ne;ne=nne;
}

void voronoicell_neighbor_2d::neighbors(std::vector<int> &v) {
    v.resize(p);
    for(int i=0;i<p;i++) v[i]=ne[i];
}

void voronoicell_neighbor_2d::neighbors_sorted(std::vector<int> &v) {
    v.resize(p);
    if(p==0) return;
    int k=0,i=0;
    do {
        v[i++]=ne[k];
        k=ed[2*k];
    } while (k!=0);
}

void voronoicell_neighbor_2d::init(double xmin,double xmax,double ymin,double ymax) {
    init_base(xmin,xmax,ymin,ymax);
    *ne=-3;ne[1]=-2;ne[2]=-4;ne[3]=-1;
}

/** Doubles the size allocation of the delete stack. If the allocation exceeds
 * the absolute maximum set in max_delete_size, then routine causes a fatal
 * error.
 * \param[in] stackp a reference to the current stack pointer. */
void voronoicell_base_2d::add_memory_ds(int *&stackp) {
    current_delete_size<<=1;
    if(current_delete_size>max_delete_size) voro_fatal_error("Delete stack memory allocation exceeded absolute maximum",VOROPP_MEMORY_ERROR);
#if VOROPP_VERBOSE >=2
    fprintf(stderr,"Delete stack memory scaled up to %d\n",current_delete_size);
#endif
    int *dsn=new int[current_delete_size],*dsnp=dsn,*dsp=ds;
    while(dsp<stackp) *(dsnp++)=*(dsp++);
    delete [] ds;ds=dsn;stackp=dsnp;
    stacke=ds+current_delete_size;
}

// Explicit instantiation
template bool voronoicell_base_2d::nplane(voronoicell_2d&,double,double,double,int);
template bool voronoicell_base_2d::nplane(voronoicell_neighbor_2d&,double,double,double,int);
template bool voronoicell_base_2d::nplane_cut(voronoicell_2d&,double,double,double,int,double,int);
template bool voronoicell_base_2d::nplane_cut(voronoicell_neighbor_2d&,double,double,double,int,double,int);
template void voronoicell_base_2d::add_memory_vertices(voronoicell_2d&);
template void voronoicell_base_2d::add_memory_vertices(voronoicell_neighbor_2d&);

}
