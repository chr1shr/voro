#version 3.6;

#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

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

#declare r=0.02;
#declare s=0.5;

union{
#include "acc_sphere_p.pov"
	rotate <270,0,0>
//	texture{T_Silver_2C}
	pigment{rgb <1,0.95,0.5>} finish{reflection 0.1 specular 0.3 ambient 0.42}
}


union{
#include "acc_sphere_v.pov"
	rotate <270,0,0>
	pigment{rgb <0.91,0.36,0.53>} finish{specular 0.56 ambient 0.24 metallic reflection 0.28}
}
