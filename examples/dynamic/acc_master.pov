#version 3.6;

#include "colors.inc"
#include "metals.inc"
#include "textures.inc"
#include "glass.inc"

global_settings{
	max_trace_level 64
}

camera {
	location <13,17,-34>
	right 0.205*x*image_width/image_height
	up 0.205*y
	look_at <0,0,0>
}

background{rgb 1}

light_source{<-8,30,-20> color rgb <0.73,0.71,0.71>}
light_source{<20,5,-15> color rgb <0.38,0.40,0.40>}

#declare r=0.016;
#declare s=0.5;

union{
#include "acc_sphere_p.pov"
	rotate <270,0,0>
//	pigment{rgb <0.75,0.88,1>}
//	finish{F_MetalB}
	pigment{rgb <1,0.93,0.5>} finish{reflection 0.1 specular 0.3 ambient 0.45}
}

/*union{
#include "acc_sphere_vi.pov"
	rotate <270,0,0>
	interior{I_Glass1}
	finish{F_Glass5}
	pigment{color rgbf <0.97, 0.99, 0.98, 0.92>}
	pigment{rgbft <0.91,0.36,0.53,0,0>} finish{specular 0.36 ambient 0.44}
}*/

#declare r=0.01;
union{
//#include "acc_sphere_v.pov"
	rotate <270,0,0>
//	texture{T_Silver_3C}
//	pigment{rgb <0.91,0.36,0.53>} finish{specular 0.56 ambient 0.24 metallic reflection 0.28}
	pigment{rgb <0.41,0.76,0.93>} finish{specular 0.56 ambient 0.24 metallic reflection 0.28}
}
