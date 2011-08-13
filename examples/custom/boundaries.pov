#version 3.6;
#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

camera {
	location <30,25,-50>
	right 0.15*x*image_width/image_height
	up 0.15*y
	look_at <0,2.6,0>
}

background{rgb 1}

light_source{<-8,30,-20> color rgb <0.77,0.75,0.75>}
light_source{<25,12,-12> color rgb <0.43,0.45,0.45>}

#declare r=0.06;
#declare s=0.5;

union{
#include "boundaries_p.pov"
	rotate <270,0,0>
	pigment{rgb <1,0.82,0.42>} finish{reflection 0.10 specular 0.4 ambient 0.36}
}

union{
#include "boundaries_v.pov"
	rotate <270,0,0>
	pigment{rgb <1,0.3,0.4>} finish{specular 0.5 ambient 0.42}
}
