#!/usr/bin/perl 

@nodes=("yuba.dhcp.lbl.gov","yuba.dhcp.lbl.gov","jordan.lbl.gov","tigris.lbl.gov","jordan.lbl.gov","ganges.lbl.gov","volga.lbl.gov","jordan.lbl.gov","volga.lbl.gov","volga.lbl.gov","tigris.lbl.gov","yuba.dhcp.lbl.gov","yuba.dhcp.lbl.gov","mekong.lbl.gov","mekong.lbl.gov");

$h=0;

$dir="output";

foreach $g (0..200) {
	$ng=sprintf("%04d",$g);

	open T,">rtemp.pov";

	open A,"dyn_master.pov";
	while (<A>) {
		last if m/\#include "temp.pov"/;
		print T;
	}
	
	$fn="$dir/${ng}_p.pov";
	system "gunzip $fn.gz";

	open B,"$fn";
	print T while <B>;
	close B;system "rm $fn";
	
	while (<A>) {
		last if m/\#include "temp2.pov"/;
		print T;
	}

	$fn="$dir/${ng}_v.pov";
	system "gunzip $fn.gz";

	open B,"$fn";
	print T while <B>;
	close B;system "rm $fn";
	
	print T while (<A>);
	close A;close T;

	print "Rendering movie frames\n";
	print "Frame $f, timestep $ts to $nodes[$h]\n";
	`rsync -z rtemp.pov $nodes[$h]:cgran/rtemp$h.pov`;
	$nice=(@nodes[$h]=~m/yuba/)?"nice -n 19":"nice +19";
	exec "ssh @nodes[$h] \"cd cgran;$nice povray +SU +Ofr_$ng.png +W512 +H512 +A0.3 +R3 -J rtemp$h.pov\" >/dev/null 2>/dev/null; rsync -rz @nodes[$h]:cgran/fr_$ng.png output ; ssh @nodes[$h] \"rm cgran/fr_$ng.png\"\n" if (($pid[$h]=fork)==0);
	if ($queue) {
		print "Waiting...\n";
		$piddone=wait;$h=0;
		$h++ while $piddone!=$pid[$h]&&$h<=$#nodes;
		die "PID return error!\n" if $h>$#nodes;
	} else {
		$h++;$queue=1 if $h==$#nodes;
	}
	$f++;
}

wait foreach 1..$#nodes;

system "qt_export --sequencerate=15 output/fr_0000.png --loadsettings=/Users/chr/misc/qtprefs/qt --replacefile $ARGV[0]";
system "rm output/fr_????.png";
