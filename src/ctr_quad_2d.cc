#include "ctr_quad_2d.hh"

namespace voro {

container_quad_2d::container_quad_2d(double ax_,double bx_,double ay_,double by_) :
	quadtree((ax_+bx_)*0.5,(ay_+by_)*0.5,(bx_-ax_)*0.5,(by_-ay_)*0.5),
	ax(ax_), bx(bx_), ay(ay_), by(by_) {

}

quadtree::quadtree(double cx_,double cy_,double lx_,double ly_) :
	cx(cx_), cy(cy_), lx(lx_), ly(ly_), ps(2), id(new int[qt_max]),
	p(new double[ps*qt_max]), co(0) {

}

quadtree::~quadtree() {
	if(id==NULL) {
		delete qne;delete qnw;
		delete qse;delete qsw;
	} else {
		delete [] p;
		delete [] id;
	}
}

void quadtree::split() {
	double hx=0.5*lx,hy=0.5*ly;
	qsw=new quadtree(cx-hx,cy-hy,hx,hy);
	qse=new quadtree(cx+hx,cy-hy,hx,hy);
	qnw=new quadtree(cx-hx,cy+hy,hx,hy);
	qne=new quadtree(cx+hx,cy+hy,hx,hy);
	for(int i=0;i<co;i++)
		(p[ps*i]<cx?(p[ps*i+1]<cy?qsw:qnw)
			   :(p[ps*i+1]<cy?qse:qne))->quick_put(id[i],p[ps*i],p[ps*i+1]);
	delete [] id;id=NULL;
	delete [] p;
}

void quadtree::put(int i,double x,double y) {
	if(id!=NULL) {
		if(co==qt_max) split();
		else {
			quick_put(i,x,y);
			return;
		}
	}
	(x<cx?(y<cy?qsw:qnw):(y<cy?qse:qne))->put(i,x,y);
}

void quadtree::draw_cross(FILE *fp) {
	if(id==NULL) {
		fprintf(fp,"%g %g\n%g %g\n\n\n%g %g\n%g %g\n\n\n",
			cx-lx,cy,cx+ly,cy,cx,cy-ly,cx,cy+ly);
		qsw->draw_cross(fp);
		qse->draw_cross(fp);
		qnw->draw_cross(fp);
		qne->draw_cross(fp);
	}
}

void container_quad_2d::draw_quadtree(FILE *fp) {
	fprintf(fp,"%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n",ax,ay,bx,ay,bx,by,ax,by,ax,ay);
	draw_cross(fp);
}

void quadtree::draw_particles(FILE *fp) {
	if(id==NULL) {
		qsw->draw_particles(fp);
		qse->draw_particles(fp);
		qnw->draw_particles(fp);
		qne->draw_particles(fp);
	} else for(int i=0;i<co;i++)
		fprintf(fp,"%d %g %g\n",id[i],p[ps*i],p[ps*i+1]);
}

}
