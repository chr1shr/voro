#version 3.6;

#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

camera {
	location <20,30,-50>
	right 0.3*x*image_width/image_height
	up 0.3*y
	look_at <0,3,0>
}

background{rgb 1}

light_source{<-8,30,-20> color rgb <0.77,0.75,0.75>}
light_source{<20,5,-15> color rgb <0.38,0.40,0.40>}

#declare r=0.025;
#declare s=0.1;

union{
#include "frustum_p.pov"
	rotate <270,0,0>
	scale 10
	pigment{rgb <0.4,0.85,0.95>} finish{reflection 0.1 specular 0.3 ambient 0.42 metallic}
}

union{
#include "frustum_v.pov"
	rotate <270,0,0>
	scale 10
	pigment{rgb <0.5,0.5,0.51>} finish{specular 0.3 ambient 0.42 reflection 0.4 metallic}
}
