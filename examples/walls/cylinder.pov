#version 3.6;

#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

camera {
	location <20,-50,30>
	sky <0,0,1>
	right -0.4*x*image_width/image_height
	up 0.4*z
	look_at <0,0,9.25>
}

background{rgb 1}

light_source{<-8,-20,30> color rgb <0.77,0.75,0.75>}
light_source{<20,-15,5> color rgb <0.38,0.40,0.40>}

#declare r=0.08;
#declare s=0.5;

union{
#include "cylinder_p.pov"
	pigment{rgb <1,0.95,0.5>} finish{reflection 0.1 specular 0.3 ambient 0.42}
}

union{
#include "cylinder_v.pov"
	pigment{rgb <0.5,0.8,1>} finish{specular 0.5 ambient 0.42}
}
