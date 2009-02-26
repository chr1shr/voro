#!/usr/bin/perl 

@nodes=("yuba.dhcp.lbl.gov","yuba.dhcp.lbl.gov","jordan.lbl.gov","tigris.lbl.gov","jordan.lbl.gov","ganges.lbl.gov","volga.lbl.gov","jordan.lbl.gov","volga.lbl.gov","volga.lbl.gov","tigris.lbl.gov","yuba.dhcp.lbl.gov","yuba.dhcp.lbl.gov","mekong.lbl.gov","mekong.lbl.gov","ganges.lbl.gov","kaczawa.lbl.gov","kaczawa.lbl.gov","kaczawa.lbl.gov","rio-grande.lbl.gov","rio-grande.lbl.gov");

$dir=$#ARGV==0?"output":$ARGV[1];

$h=0;
$voronoi=0;
$verb=" >/dev/null 2>/dev/null";$verb="";

$g=0;$ng="0000";

while (-e "$dir/${ng}_p.pov.gz") {
	if (-e "$dir/fr_$ng.png") {
		$g++;
		$ng=sprintf("%04d",$g);
		next;
	}

	open T,">rtemp.pov";

	open A,"packing_master.pov";
	while (<A>) {
		last if m/\#include "temp.pov"/;
		print T;
	}
	
	$fn="$dir/${ng}_p.pov";
	system "gunzip -c $fn.gz >gztemp.pov";
	open B,"gztemp.pov";
	print T while <B>;
	close B;system "rm gztemp.pov";

	if ($voronoi==1) {
		
		while (<A>) {
			last if m/\#include "temp2.pov"/;
			print T;
		}	
	
		$fn="$dir/${ng}_v.pov";
		system "gunzip $fn.gz";
		
		$a=0;
		open B,"$fn";
		while(<B>) {
			m/^cylinder{<([-\.\d]+),([-\.\d]+),([-\.\d]+)>,<\1,\2,\3>,r}/?$a++:(print T);
		}
		close B;system "rm $fn";
		print "$a degenerate cylinders stripped\n" if $a>0;
	}
	
	print T while (<A>);
	close A;close T;

	print "Frame $g to $nodes[$h]\n";
	`rsync -z rtemp.pov $nodes[$h]:cgran/rtemp$h.pov`;
	$nice=(@nodes[$h]=~m/yuba/)?"nice -n 19":"nice +19";
	exec "ssh @nodes[$h] \"cd cgran;$nice povray +SU +Ofr_$ng.png +W300 +H800 +A0.001 +R3 -J rtemp$h.pov\"$verb; rsync -rz @nodes[$h]:cgran/fr_$ng.png $dir ; ssh @nodes[$h] \"rm cgran/fr_$ng.png\"\n" if (($pid[$h]=fork)==0);
#	exec "ssh @nodes[$h] \"cd cgran;$nice povray +SU +Ofr_$ng.png +W800 +H700 +A0.001 +R3 -J rtemp$h.pov\" >/dev/null 2>/dev/null; rsync -rz @nodes[$h]:cgran/fr_$ng.png $dir ; ssh @nodes[$h] \"rm cgran/fr_$ng.png\"\n" if (($pid[$h]=fork)==0);
	if ($queue) {
		print "Waiting...\n";
		$piddone=wait;$h=0;
		$h++ while $piddone!=$pid[$h]&&$h<=$#nodes;
		die "PID return error!\n" if $h>$#nodes;
	} else {
		$h++;$queue=1 if $h==$#nodes;
	}

	$g++;
	$ng=sprintf("%04d",$g);
}

wait foreach 1..$#nodes;

if ($g>0) {
	system "qt_export --sequencerate=15 $dir/fr_0000.png --loadsettings=/Users/chr/misc/qtprefs/qt --replacefile $ARGV[0]";
#	system "rm $dir/fr_????.png";
} else {
	print "No files to process\n";
}
