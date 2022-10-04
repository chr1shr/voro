#include <cstdio>
#include <cstdlib>

#include "voro++.hh"
using namespace voro;

// Returns floating point number uniformly distributed over [0,1]
inline double rnd() {return (1./RAND_MAX)*static_cast<double>(rand());}

int main() {

    // Number of parallel threads
    int num_t=4;

    // Number of particles to use
    int N=100000;

    // Construct a 2D container as a periodic unit square divided into a
    // 160x160 grid of blocks. Each block initially holds up to 8 particles. The
    // final arguments set the number of voro_compute objects for use by the
    // threads.
    container_2d con(0.0,1.0,0.0,1.0,160,160,true,true,8,num_t);

    // Add particles to the container
    for(int i=0;i<N;i++) con.put(i,rnd(),rnd());

    // Declare iterator
    container_2d::iterator cli;

    // Parallel Voronoi computation to compute the average Voronoi cell
    // perimeter
    double tperim=0.;
#pragma omp parallel num_threads(num_t)
    {
        // Thread-private Voronoi cell object and perimeter counter
        voronoicell_2d c(con);
        double perim=0.;

        // Iterate through the particles
#pragma omp for
        for(cli=con.begin();cli<con.end();cli++)
            if(con.compute_cell(c,cli))
                perim+=c.perimeter();

        // Add local perimeter counter to global perimeter counter using atomic
        // operation to prevent race condition
#pragma omp atomic
        tperim+=perim;
    }

    // Print average Voronoi cell perimeter
    printf("Average Voronoi cell perimeter is %.12g\n",tperim/N);
}
