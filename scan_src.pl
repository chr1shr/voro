#!/usr/bin/perl
open A,"$ARGV[0]";
open B,">temp";

$c=0;$vr=0;
while(<A>) {
	$c++ if s/[ \t]+$//;
	s/voropp_order/particle_order/g;
	s/voropp_/voro_/g;
	s/voro_safe_fopen/safe_fopen/g;
	s/v_loop/c_loop/g;
	print B;
}
close A;
close B;
`mv temp $ARGV[0]`;

print "$c tabs/spaces removed from $ARGV[0] ".($vr==1?"(check)\n":"(fail)\n");
