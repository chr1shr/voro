// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

#include <cstdlib>
#include <cmath>

#include "voro++.cc"
using namespace voro;

// Set up timing routine. If code was compiled with OpenMP, then use the
// accurate wtime function. Otherwise use the clock function in the ctime
// library.
#ifdef _OPENMP
#include "omp.h"
inline double wtime_() {return omp_get_wtime();}
#else
inline double wtime_() {return 1./CLOCKS_PER_SEC*static_cast<double>(clock());}
#endif

// Returns a random double that is uniformly distributed between 0 and 1
inline double rnd() {return 1./RAND_MAX*static_cast<double>(rand());}

// Prints out a message about the syntax of the command-line utility
void syntax_message() {
    puts("Syntax: ./timing_test <num> <block_lo> <block_hi> <reps> <prd>\n"
         "Arguments:\n"
         "<num>         The number of particles       [100000]\n"
         "<block_lo>    The lower limit of blocks     [10]\n"
         "<block_hi>    The upper limit of blocks     [40]\n"
         "<reps>        The number of repeat trails   [1]\n"
         "<prd>         Whether to use periodicity    [0]\n\n"
         "If any argument is missing, the default value in the square brackets is used\n");
    exit(1);
}

int main(int argc,char **argv) {

    // Check for a valid number of command-line arguments
    if(argc>6) syntax_message();

    // Read the command-line arguments, and check that they are valid
    int blo=10,bhi=40,num=100000,reps=1;
    bool prd=false;
    if(argc>1) {
        num=atoi(argv[1]);
        if(num<=0) syntax_message();
        if(argc>2) {
            blo=atoi(argv[2]);
            if(blo<=0) syntax_message();
            if(argc>3) {
                bhi=atoi(argv[3]);
                if(bhi<blo) syntax_message();
                if(argc>4) {
                    reps=atoi(argv[4]);
                    if(reps<=0) syntax_message();
                    if(argc>5) {
                        int k=atoi(argv[5]);
                        if(k!=0&&k!=1) syntax_message();
                        prd=k==1;
                    }
                }
            }
        }
    }

    // Loop over the block sizes to test
    for(int b=blo;b<=bhi;b++) {
        container_3d con(0,1,0,1,0,1,b,b,b,prd,prd,prd,8);

        // Perform the repeat trials
        double st=0,stt=0,sc=0,scc=0,t0,t1;
        for(int l=0;l<reps;l++) {

            // Randomly insert the particles into the container
            if(l!=0) con.clear();
            t0=wtime_();
            for(int i=0;i<num;i++) con.put(i,rnd(),rnd(),rnd());

            // Perform the dummy computation to compute all the Voronoi cells
            t1=wtime_();
            con.compute_all_cells();

            // Store the timing statistics
            t0=t1-t0;t1=wtime_()-t1;
            printf("%g %g\n",t0,t1);
            st+=t0;stt+=t0*t0;
            sc+=t1;scc+=t1*t1;
        }

        // Output the timing information
        st/=reps;stt=stt/reps-st*st;
        sc/=reps;scc=scc/reps-sc*sc;
        printf("%d %g %g %g %g\n",b,st,stt<0?0:sqrt(stt),sc,scc<0?0:sqrt(scc));
        con.clear();
    }
}
