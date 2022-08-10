// Voro++ example code demonstrating find_voronoi_cell function
// By Chris H. Rycroft and the Rycroft Group

#include "voro++.hh"
using namespace voro;

// The sampling distance for the grids of find_voronoi_cell calls
const double h=0.05;

// The cube of the sampling distance, corresponding the amount of volume
// associated with a sample point
const double hcube=h*h*h;

// Set the number of particles that are going to be randomly introduced
const int particles=20;

// This function returns a random double between 0 and 1
double rnd() {return 1./RAND_MAX*static_cast<double>(rand());}

int main() {
    int i;
    double x,y,z,rx,ry,rz;

    // Create a unit cube and make it non-periodic in each of the three
    // coordinates. Allocate space for eight particles within each
    // computational block
    container_3d con(0,1,0,1,0,1,5,5,5,false,false,false,8);

    // Randomly add particles into the container
    for(i=0;i<particles;i++) {
        x=rnd();
        y=rnd();
        z=rnd();
        con.put(i,x,y,z);
    }

    // Output the particle positions in gnuplot format
    con.draw_particles("find_voro_cell_p.gnu");

    // Scan a 2D slice in the container, and for each point in the slice, find
    // the Voronoi cell that the point is in. Store a vector
    FILE *f1=safe_fopen("find_voro_cell.vec","w");
    for(x=0.5*h;x<1;x+=h) for(y=0.5*h;y<1;y+=h) {
        if(con.find_voronoi_cell(x,y,0.5,rx,ry,rz,i))
            fprintf(f1,"%g %g %g %g %g %g %g\n",x,y,0.5,rx-x,ry-y,rz-0.5,
                sqrt((rx-x)*(rx-x)+(ry-y)*(ry-y)+(rz-0.5)*(rz-0.5)));
        else fprintf(stderr,"# find_voronoi_cell error for %g %g 0.5\n",x,y);
    }
    fclose(f1);

    // Create a blank array for storing the sampled Voronoi volumes
    int samp_v[particles];
    for(i=0;i<particles;i++) samp_v[i]=0;

    // Scan over a grid covering the entire container, finding which Voronoi
    // cell each point is in, and tallying the result as a method of sampling
    // the volume of each Voronoi cell
    for(z=0.5*h;z<1;z+=h) for(y=0.5*h;y<1;y+=h) for(x=0.5*h;x<1;x+=h) {
        if(con.find_voronoi_cell(x,y,z,rx,ry,rz,i)) samp_v[i]++;
        else fprintf(stderr,"# find_voronoi_cell error for %g %g %g\n",x,y,z);
    }

    // Output the Voronoi cells in gnuplot format and a file with the
    // comparisons between the Voronoi cell volumes and the sampled volumes
    f1=safe_fopen("find_voro_cell.vol","w");
    FILE *f2=safe_fopen("find_voro_cell_v.gnu","w");
    container_3d::iterator ci;
    voronoicell_3d c(con);
    for(ci=con.begin();ci<con.end();ci++) if(con.compute_cell(c,ci)) {

        // Get the position and ID information for the particle currently being
        // considered by the loop. Ignore the radius information.
        con.pos(ci,x,y,z);

        // Save and entry to the .vol file, storing both the computed
        // Voronoi cell volume, and the sampled volume based on the
        // number of grid points that were inside the cell
        fprintf(f1,"%d %g %g %g %g %g\n",i,x,y,z,c.volume(),samp_v[i]*hcube);

        // Draw the Voronoi cell
        c.draw_gnuplot(x,y,z,f2);
    }
    fclose(f1);
    fclose(f2);
}
