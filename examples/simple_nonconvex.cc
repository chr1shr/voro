#include "voro++_2d.hh"
#include <iostream>
int main() {
	double up[14];
	int noofbnds=7;
	voronoicell_2d up1;
	voronoicell_2d up2;
	up[0]=0; up[1]=0; up[2]=-2; up[3]=5; up[4]=-5; up[5]=5;
	up[6]=-5; up[7]=-5; up[8]=5; up[9]=-5; up[10]=5; up[11]=5;
	up[12]=2; up[13]=5;
	up1.init_nonconvex(up,noofbnds);
	up1.draw_gnuplot(0,0,"up.gnu");
	up1.nonconvexity=true;
	up1.plane_nonconvex(0,-3,9);
	up1.draw_gnuplot(0,0,"up_cut1.gnu");
	up1.plane_nonconvex(-2,5,29);
	up1.draw_gnuplot(0,0,"up_cut2.gnu");
	up2.init_nonconvex(up, noofbnds);
	up2.plane_nonconvex(2,5,29);
	up2.draw_gnuplot(0,0,"up_cut3.gnu");

	double right[14];
	voronoicell_2d right1;
	voronoicell_2d right2;
	right[0]=0; right[1]=0; right[2]=5; right[3]=2; right[4]=5; right[5]=5;
	right[6]=-5; right[7]=5; right[8]=-5; right[9]=-5; right[10]=5; right[11]=-5;
	right[12]=5; right[13]=-2;
	right1.init_nonconvex(right, noofbnds);
	right2.init_nonconvex(right, noofbnds);
	right1.draw_gnuplot(0,0,"right.gnu");
	right1.plane_nonconvex(-3,0,9);
	right1.draw_gnuplot(0,0,"right_cut1.gnu");
	right1.plane_nonconvex(2,2,4);
	right1.draw_gnuplot(0,0,"right_cut2.gnu");
	right2.plane_nonconvex(2,-2,4);
	right2.draw_gnuplot(0,0,"right_cut3.gnu");

	double down[14];
	voronoicell_2d down1;
	voronoicell_2d down2;
	down[0]=0; down[1]=0; down[2]=2; down[3]=-5; down[4]=5; down[5]=-5;
	down[6]=5; down[7]=5; down[8]=-5; down[9]=5; down[10]=-5; down[11]=-5;
	down[12]=-2; down[13]=-5;
	down1.init_nonconvex(down, noofbnds);
	down2.init_nonconvex(down, noofbnds);

	down1.draw_gnuplot(0,0,"down.gnu");
	down1.plane_nonconvex(0,3,9);
	down1.draw_gnuplot(0,0,"down_cut1.gnu");
	down1.plane_nonconvex(2,-2,4);
	down1.draw_gnuplot(0,0,"down_cut2.gnu");
	down2.plane_nonconvex(-2,-2,4);
	down2.draw_gnuplot(0,0,"down_cut3.gnu");
	
	double left[14];
	voronoicell_2d left1;
	voronoicell_2d left2;
	left[0]=0; left[1]=0; left[2]=-5; left[3]=-2; left[4]=-5; left[5]=-5;
	left[6]=5; left[7]=-5; left[8]=5; left[9]=5; left[10]=-5; left[11]=5;
	left[12]=-5; left[13]=2;
	left1.init_nonconvex(left, noofbnds);
	left2.init_nonconvex(left, noofbnds);
	left1.draw_gnuplot(0,0,"left.gnu");
	left1.init_nonconvex(left, noofbnds);
	left1.plane_nonconvex(3,0,9);
	left1.draw_gnuplot(0,0,"left_cut1.gnu");
	left1.plane_nonconvex(-2,2,4);
	left1.draw_gnuplot(0,0,"left_cut2.gnu");
	left2.plane_nonconvex(-2,-2,4);
	left2.draw_gnuplot(0,0,"left_cut3.gnu");

	cout << up1.reg1[0] << "   " << up1.reg1[1] << "  " << up1.reg1[2] << "   " <<
	up1.reg1[3]<< "    ";
	cout << up1.reg2[0] << "   " << up1.reg2[1] << "  " << up1.reg2[2] << "   " <<
	up1.reg2[3]<< "    ";

	cout << right1.reg1[0] << "   " << right1.reg1[1] << "  " << right1.reg1[2] << "   " <<
	right1.reg1[3]<< "    ";
	cout << right1.reg2[0] << "   " << right1.reg2[1] << "  " << right1.reg2[2] << "   " <<
	right1.reg2[3]<< "    ";

	cout << down1.reg1[0] << "   " << down1.reg1[1] << "  " << down1.reg1[2] << "   " <<
	down1.reg1[3]<< "    ";
	cout << down1.reg2[0] << "   " << down1.reg2[1] << "  " << down1.reg2[2] << "   " <<
	down1.reg2[3]<< "    ";

	cout << left1.reg1[0] << "   " << left1.reg1[1] << "  " << left1.reg1[2] << "   " <<
	left1.reg1[3]<< "    ";
	cout << left1.reg2[0] << "   " << left1.reg2[1] << "  " << left1.reg2[2] << "   " <<
	left1.reg2[3]<< "    ";
	

	

	
	return 0;
}











