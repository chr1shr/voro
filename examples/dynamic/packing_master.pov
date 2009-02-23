#version 3.6;

#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

camera {
	orthographic
	location <0,58,-100>
	right 120*x*image_width/image_height
	up 120*y
	look_at <0,58,0>
}

background{rgb 1}

light_source{<-56,100,-140> color rgb <0.77,0.75,0.75>}
light_source{<80,-20,-80> color rgb <0.38,0.40,0.40>}

#declare r=0.07;
#declare s=0.5;

/*union{
	torus{2.5,r}
	torus{3.5,r}
	rotate <90,0,0>
	translate <0,0,-10>
	pigment{rgb <0.4,0.6,0.8>} finish{specular 0.5 ambient 0.42}
}*/

union{
#include "temp.pov"
	rotate <270,0,0>
//	texture{T_Chrome_3C}
	pigment{rgb <1,0.95,0.5>} finish{reflection 0.1 specular 0.3 ambient 0.42}
}


union{
#include "temp2.pov"
	pigment{rgb <0.4,0.6,0.8>} finish{specular 0.5 ambient 0.42}
}
