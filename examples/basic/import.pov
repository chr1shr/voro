#version 3.6;
#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

global_settings {
	max_trace_level 64
}

camera {
	location <30,25,-50>
	right 0.24*x*image_width/image_height
	up 0.24*y
	look_at <0,4.5,0>
}

background{rgb 1}

light_source{<-8,30,-20> color rgb <0.77,0.75,0.75>}
light_source{<25,12,-12> color rgb <0.38,0.40,0.40>}

#declare r=0.08;
#declare s=0.5;

union{
#include "pack_ten_cube.pov"
	rotate <270,0,0>
	pigment{rgb 0.95} finish{reflection 0.1 specular 0.3 ambient 0.42}
}

union{
#include "cells_ten_cube.pov"
	rotate <270,0,0>
	pigment{rgb <1,0.4,0.45>} finish{specular 0.5 ambient 0.42}
}
