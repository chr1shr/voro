#include "voro++_2d.hh"


using namespace voro;
#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

void test(int &x,int &y){
	x=1;y=2;
}

int main(int argc,char **argv){
	char *outfn(new char[1000]);
	sprintf(outfn,"circle");
	FILE *fp=safe_fopen(outfn,"w");
	for(int i=0;i<200;i++){
		 fprintf(fp,"%i %f %f\n",i+4,cos((double)2 * 3.14159 * ((double)i)/200),sin((double)2*3.14159 * ((double)i)/200));
	}
	cout << "the sin is" << sin (30.0 * 3.14159/180) << endl;
	return 0;





}

