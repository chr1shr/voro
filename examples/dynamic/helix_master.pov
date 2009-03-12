#version 3.6;

#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

global_settings {
	max_trace_level 256
}

camera {
	location <30,35,-65>
	right 0.60*x*image_width/image_height
	up 0.60*y
	look_at <9,0,0>
}

background{rgb 1}

light_source{<-20,100,-150> color rgb <0.72,0.70,0.70>}
light_source{<80,-10,-90> color rgb <0.38,0.40,0.40>}

#declare s=0.5;
#declare r=0.07;

union{
#include "helix_p.pov"
	rotate <270,0,0>
	rotate <0,-32,95>
	rotate <0,-8,0>
	pigment{rgb 0.90}
finish {
    ambient 0.35
    brilliance 3
    diffuse 0.5
    metallic
    specular 0.35
    roughness 1/60
   // reflection 0.35
    reflection 0.29
}
	//	pigment{rgb <0.9,0.85,0.35>} finish{reflection 0.1 specular 0.3 ambient 0.42}
}

union{
#include "helix_v.pov"
	rotate <270,0,0>
	pigment{radial pigment_map {
			[0 rgb <0.3,0.6,0.9>]
			[0.33333333 rgb <0.3,0.3,0.9>]
			[0.66666666 rgb <0.6,0.8,0.9>]
			[1 rgb <0.3,0.6,0.9>]}}
	finish{specular 0.4 ambient 0.42 reflection 0.1}
	rotate <0,-32,95>
	rotate <0,-8,0>
//	texture{T_Chrome_3C}
}
