#!/usr/bin/perl 

#$verb=">/dev/null 2>/dev/null"
$nodes=1;

$h=0;

$dir="output";

$g=0;$ng="0000";

while (-e "$dir/${ng}_p.pov.gz") {
	open T,">rtemp$h.pov";

	open A,"yeast_master.pov";
	while (<A>) {
		last if m/\#include "temp.pov"/;
		print T;
	}
	
	$fn="$dir/${ng}_p.pov";
	system "gunzip -c $fn.gz >gztemp.pov";

	open B,"gztemp.pov";
	print T while <B>;
	close B;system "rm gztemp.pov";
	
#	while (<A>) {
#		last if m/\#include "temp2.pov"/;
#		print T;
#	}
#
#	$fn="$dir/${ng}_v.pov";
#	system "gunzip $fn.gz";
#	
#	$a=0;
#	open B,"$fn";
#	while(<B>) {
#		m/^cylinder{<([-\.\d]+),([-\.\d]+),([-\.\d]+)>,<\1,\2,\3>,r}/?$a++:(print T);
#	}
#	close B;system "rm $fn";
#	print "$a degenerate cylinders stripped\n" if $a>0;

	print T while (<A>);
	close A;close T;

	print "Rendering movie frames\n";
	print "Frame $g, timestep $ts to $h\n";
	exec "povray +SU +Ofr_$ng.png +W800 +H600 +A0.3 +R3 -J rtemp$h.pov $verb; mv fr_$ng.png output" if (($pid[$h]=fork)==0);
	if ($queue) {
		print "Waiting...\n";
		$piddone=wait;$h=0;
		$h++ while $piddone!=$pid[$h]&&$h<=$nodes;
		die "PID return error!\n" if $h>$nodes;
	} else {
		$h++;$queue=1 if $h==$nodes;
	}
	$g++;
	$ng=sprintf("%04d",$g);
}

wait foreach 1..$nodes;

if ($g>0) {
	system "qt_export --sequencerate=15 output/fr_0000.png --loadsettings=/Users/chr/misc/qtprefs/qt --replacefile $ARGV[0]";
	system "rm output/fr_????.png";
} else {
	print "No files to process\n";
}
