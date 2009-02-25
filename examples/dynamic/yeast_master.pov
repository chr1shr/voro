#version 3.6;

#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

camera {
	location <40,35,-80>
	right 0.65*x*image_width/image_height
	up 0.65*y
	look_at <0,-5,0>
}

background{rgb 1}

light_source{<0,100,-150> color rgb <0.72,0.70,0.70>}
light_source{<80,-30,-90> color rgb <0.38,0.40,0.40>}

#declare r=0.2;
#declare s=0.5;
#declare b=20;

union {
	cylinder{<-b,-b,-b>,<-b,-b,b>,r}
	cylinder{<b,-b,-b>,<b,-b,b>,r}
	cylinder{<-b,b,-b>,<-b,b,b>,r}
	cylinder{<b,b,-b>,<b,b,b>,r}
	cylinder{<-b,-b,-b>,<-b,b,-b>,r}
	cylinder{<b,-b,-b>,<b,b,-b>,r}
	cylinder{<-b,-b,b>,<-b,b,b>,r}
	cylinder{<b,-b,b>,<b,b,b>,r}
	cylinder{<-b,-b,-b>,<b,-b,-b>,r}
	cylinder{<-b,b,-b>,<b,b,-b>,r}
	cylinder{<-b,-b,b>,<b,-b,b>,r}
	cylinder{<-b,b,b>,<b,b,b>,r}
	sphere{<-b,-b,-b>,r}
	sphere{<-b,-b,b>,r}
	sphere{<-b,b,-b>,r}
	sphere{<-b,b,b>,r}
	sphere{<b,-b,-b>,r}
	sphere{<b,-b,b>,r}
	sphere{<b,b,-b>,r}
	sphere{<b,b,b>,r}
	pigment{rgb <0.4,0.6,0.8>} finish{specular 0.5 ambient 0.42}
}

union{
#include "temp.pov"
	rotate <270,0,0>
//	texture{T_Chrome_3C}
	pigment{rgb <0.9,0.85,0.35>} finish{reflection 0.1 specular 0.3 ambient 0.42}
}
