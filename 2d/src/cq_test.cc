#include <cstdio>
#include <cmath>

#include "voro++_2d.hh"
using namespace voro;

// Set up timing routine. If code was compiled with OpenMP, then use the
// accurate wtime function. Otherwise use the clock function in the ctime
// library.
#ifdef _OPENMP
#include "omp.h"
inline double wtime() {return omp_get_wtime();}
#else
#include <ctime>
inline double wtime() {return double(clock())/CLOCKS_PER_SEC;}
#endif

double rnd(double a,double b) {return a+(b-a)/RAND_MAX*static_cast<double>(rand());}

int main() {
    int i,l,n=10;double x,y,r,t1,t2;

    while(n<10000000) {
        container_quad_2d con1(-1,1,-1,1);
        l=int(sqrt(double(n))/3.46)+1;
        container_2d con2(-1,1,-1,1,l,l,false,false,8);

        for(i=0;i<n;i++) {
            x=rnd(-1,1);
            y=rnd(-1,1);
            r=1;//(x*x+y*y)*0.5;
            con1.put(i,x*r,y*r);
            con2.put(i,x*r,y*r);
        }

        con1.setup_neighbors();
        t2=wtime();
        con2.compute_all_cells();
        t2=wtime()-t2;
        t1=wtime();
        con1.compute_all_cells();
        t1=wtime()-t1;

        printf("%d %g %g\n",n,t1,t2);
        n+=n>>2;
    }
}
