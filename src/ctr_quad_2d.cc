#include "ctr_quad_2d.hh"
#include "quad_march.hh"

#include <vector>

namespace voro {

container_quad_2d::container_quad_2d(double ax_,double bx_,double ay_,double by_) :
	quadtree((ax_+bx_)*0.5,(ay_+by_)*0.5,(bx_-ax_)*0.5,(by_-ay_)*0.5),
	ax(ax_), bx(bx_), ay(ay_), by(by_) {

}

quadtree::quadtree(double cx_,double cy_,double lx_,double ly_) :
	cx(cx_), cy(cy_), lx(lx_), ly(ly_), ps(2), id(new int[qt_max]),
	p(new double[ps*qt_max]), co(0), nco(0), nmax(0) {

}

quadtree::~quadtree() {
	if(id==NULL) {
		delete qne;delete qnw;
		delete qse;delete qsw;
	} else {
		delete [] p;
		delete [] id;
	}
	if(nmax>0) delete [] nei;
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

void quadtree::draw_neighbors(FILE *fp) {
	for(int i=0;i<nco;i++)
		fprintf(fp,"%g %g %g %g\n",cx,cy,nei[i]->cx-cx,nei[i]->cy-cy);
	if(id==NULL) {
		qsw->draw_neighbors(fp);
		qse->draw_neighbors(fp);
		qnw->draw_neighbors(fp);
		qne->draw_neighbors(fp);
	}
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

void quadtree::setup_neighbors() {
	if(id==NULL) {
		qsw->setup_neighbors();
		qse->setup_neighbors();
		qnw->setup_neighbors();
		qne->setup_neighbors();
		we_neighbors(qsw,qse);
		we_neighbors(qnw,qne);
		ns_neighbors(qsw,qnw);
		ns_neighbors(qse,qne);
	}
}

void quadtree::we_neighbors(quadtree *qw,quadtree *qe) {
	const int ma=1<<30;
	quad_march<0> mw(qw);
	quad_march<1> me(qe);
	while(mw.s<ma||me.s<ma) {
		mw.cu()->add_neighbor(me.cu());
		me.cu()->add_neighbor(mw.cu());
		if(mw.ns>me.ns) me.step();
		else {
			if(mw.ns==me.ns) me.step();
			mw.step();
		}
	}
}

void quadtree::ns_neighbors(quadtree *qs,quadtree *qn) {
	const int ma=1<<30;
	quad_march<2> ms(qs);
	quad_march<3> mn(qn);
	while(ms.s<ma||mn.s<ma) {
		ms.cu()->add_neighbor(mn.cu());
		mn.cu()->add_neighbor(ms.cu());
		if(ms.ns>mn.ns) mn.step();
		else {
			if(ms.ns==mn.ns) mn.step();
			ms.step();
		}
	}
}

void quadtree::add_neighbor_memory() {
	if(nmax==0) {
		nmax=4;
		nei=new quadtree*[nmax];
	}
	if(nmax>16777216) {
		fputs("Maximum quadtree neighbor memory exceeded\n",stderr);
		exit(1);
	}
	nmax<<=1;
	quadtree** pp=new quadtree*[nmax];
	for(int i=0;i<nco;i++) pp[i]=nei[i];
	delete [] nei;
	nei=pp;
}

}
