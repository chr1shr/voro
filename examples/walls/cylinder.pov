#version 3.6;

#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

camera {
	location <20,30,-50>
	right 0.4*x*image_width/image_height
	up 0.4*y
	look_at <0,9.25,0>
}

background{rgb 1}

light_source{<-8,30,-20> color rgb <0.77,0.75,0.5>}
light_source{<20,5,-15> color rgb <0.38,0.40,0.40>}

#declare r=0.08;
INCLUDE
//#include "cylinder_p.pov"
//#include "cylinder_v.pov"

union{
	particles
	rotate <270,0,0>
	rotate <0,ROT,0>
	pigment{rgb <1,0.95,0.5>} finish{reflection 0.1 specular 0.3 ambient 0.42}
}

union{
	voronoi
	rotate <270,0,0>
	rotate <0,ROT,0>
	pigment{rgb <0.5,0.8,1>} finish{specular 0.5 ambient 0.42}
}
