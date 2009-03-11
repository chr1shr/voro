#version 3.6;

#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

global_settings {
	max_trace_level 256
}

camera {
	location <40/1.3,35*1.3,-80/1.3>
	right 0.77*x*image_width/image_height
	up 0.77*y
	look_at <0,-5,0>
}

background{rgb 1}

light_source{<0,100,-150> color rgb <0.72,0.70,0.70>}
light_source{<80,-30,-90> color rgb <0.38,0.40,0.40>}

#declare s=0.5;
#declare r=0.07;

union{
#include "helix_p.pov"
	rotate <270,0,0>
	rotate <0,0,28>
	texture{T_Chrome_3C}
//	pigment{rgb <0.9,0.85,0.35>} finish{reflection 0.1 specular 0.3 ambient 0.42}
}

union{
#include "helix_v.pov"
	rotate <270,0,0>
	rotate <0,0,28>
//	texture{T_Chrome_3C}
	pigment{rgb <0.4,0.7,1>} finish{specular 0.4 ambient 0.42}
}
