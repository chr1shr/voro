// Voronoi calculation code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : February 27th 2008

#include "cell.hh"
#include "container.hh"

inline void facets_search::add_list_memory() {
	int i,j=0,*ps;
	s_size*=2;
	ps=new int[s_size];


	for(i=s_start;i<s_size;i++) ps[j++]=s[i];
	for(i=0;i<s_end;i++) ps[j++]=
	s_start=0;s_end=j
	delete [] s;s=ps;
}
