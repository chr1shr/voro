#version 3.6;

#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

global_settings {
	max_trace_level 64
}

camera {
	location <0,34,-40>
	right 0.4*x*image_width/image_height
	up 0.4*y
	look_at <0,-1.31,0>
}

background{rgb 1}

light_source{<-16,43,-20> color rgb <0.72,0.69,0.69>}
light_source{<30,12,-15> color rgb <0.34,0.37,0.37>}

#declare s=0.5;
#declare r=0.06;

union{
	#include "torus_p.pov"
	rotate <90,0,0>
	pigment{rgb 0.97}
	finish{reflection 0.1 ambient 0.30 specular 0.3}
}

union{
	#include "torus_v.pov"
	rotate <90,0,0>
	pigment{radial pigment_map {
			[0 rgb <0.5,0.7,1>]
			[0.25 rgb <0.38,0.82,0.92>]
			[0.5 rgb <0.5,0.7,1>]
			[0.75 rgb <0.65,0.4,1>]
			[1 rgb <0.5,0.7,1>]}}
	finish{specular 0.3 ambient 0.3 reflection 0.1}
}
