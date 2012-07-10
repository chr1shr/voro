#include "voro++_2d.hh"


using namespace voro;
#include <stdio.h>
#include <iostream>
#include <vector>
using namespace std;
/*
void add_memory_array(double* &old,int size){
	double *newa= new double[2*size];
	for(int i=0;i<(size);i++){
		newa[9-i]=old[i];
	}
	delete [] old;
	old=newa;
}*/

int main(int argc,char **argv){

	v_connect connect;
	char *outfn;

	FILE *fp=fopen(argv[1],"r");
	connect.import(fp);
	connect.assemble_vertex();
	connect.assemble_gen_ed();
	connect.assemble_boundary();
//	connect.lloyds(.0001);
	
	outfn=new char[1000];
	sprintf(outfn,"%s.gv",argv[1]);
	connect.draw_vtg_gnu(outfn);
	delete [] outfn;
	
	outfn=new char[1000];
	sprintf(outfn,"%s.vertlabel",argv[1]);
	connect.label_vertices(outfn);
	delete [] outfn;

	outfn=new char[1000];
	sprintf(outfn,"%s.genlabel",argv[1]);
	connect.label_generators(outfn);
	delete [] outfn;

	outfn=new char[1000];
	sprintf(outfn,"%s.edlabel",argv[1]);
	connect.label_edges(outfn);
	delete [] outfn;

	outfn=new char[1000];
	sprintf(outfn,"%s.neighbors",argv[1]);
	connect.draw_gen_gen(outfn);
	delete [] outfn;

	outfn=new char[1000];
	sprintf(outfn,"%s.gen_to_ed_table",argv[1]);
	connect.print_gen_to_ed_table(outfn);
	delete outfn;

	outfn=new char[1000];
	sprintf(outfn,"%s.gen_to_vert_table",argv[1]);
	connect.print_gen_to_vert_table(outfn);
	delete [] outfn;

	outfn=new char[1000];
	sprintf(outfn,"%s.vert_to_gen_table",argv[1]);
	connect.print_vert_to_gen_table(outfn);
	delete [] outfn;

	outfn=new char[1000];
	sprintf(outfn,"%s.ed_to_gen_table",argv[1]);
	connect.print_ed_to_gen_table(outfn);
	delete []outfn;

	outfn=new char[1000];
	sprintf(outfn,"%s.vert_to_ed_table",argv[1]);
	connect.print_vert_to_ed_table(outfn);
	delete [] outfn;

	outfn=new char[1000];
	sprintf(outfn,"%s.vert_boundary_table",argv[1]);
	connect.print_vert_boundary(outfn);
	delete [] outfn;

	outfn=new char[1000];
	sprintf(outfn,"%s.ed_boundary_table",argv[1]);
	connect.print_ed_boundary(outfn);
	delete [] outfn;

	outfn=new char[1000];
	sprintf(outfn,"%s.gnu",argv[1]);
	connect.draw_gnu(outfn);
	delete [] outfn;

	outfn=new char[1000];
	sprintf(outfn,"%s.centroidlabel",argv[1]);
	connect.label_centroids(outfn);
	delete [] outfn;

	outfn=new char[1000];
	sprintf(outfn,"%s.dualmesh",argv[1]);
	connect.draw_median_mesh(outfn);
	delete [] outfn;

	outfn=new char[1000];
	sprintf(outfn,"%s.ascii",argv[1]);
	connect.ascii_output(outfn);
	delete [] outfn;

	if(argc==4){
		double x=atof(argv[2]); 
		double y=atof(argv[3]);
		outfn=new char[1000];
		sprintf(outfn,"%s.closest",argv[1]);
		connect.draw_closest_generator(outfn,x,y);
		delete [] outfn;
	}
/*
	double *test=new double[10];
	for(int i=0;i<10;i++){
		test[i]=i;
	}
	for(int i=0;i<10;i++){
		cout << test[i] << endl;
	}
	cout << "\n\n\n" << endl;
	add_memory_array(test,10);
	for(int i=0;i<20;i++){
		cout << test[i] << endl;
	}
*/	
	
}























		
	










