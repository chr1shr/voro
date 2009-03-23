#version 3.6;
#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

camera {
	location <30,25,-50>
	right 0.15*x*image_width/image_height
	up 0.15*y
	look_at <0,2.8,0>
}

background{rgb 1}

light_source{<-8,30,-20> color rgb <0.77,0.75,0.75>}
light_source{<25,12,-12> color rgb <0.43,0.45,0.45>}

#declare r=0.08;
#declare s=0.5;

union{
#include "pack_six_cube_p.pov"
	rotate <270,0,0>
	pigment{rgb <0.92,0.65,1>} finish{reflection 0.17 specular 0.3 ambient 0.42}
}

union{
#include "pack_six_cube_v.pov"
	rotate <270,0,0>
	pigment{rgb <0.7,0.95,1>} finish{specular 0.5 ambient 0.42}
}
