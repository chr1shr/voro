#!/usr/bin/perl 

@nodes=("yuba.dhcp.lbl.gov","yuba.dhcp.lbl.gov","jordan.lbl.gov","tigris.lbl.gov","jordan.lbl.gov","ganges.lbl.gov","volga.lbl.gov","jordan.lbl.gov","volga.lbl.gov","volga.lbl.gov","tigris.lbl.gov","yuba.dhcp.lbl.gov","yuba.dhcp.lbl.gov","mekong.lbl.gov","mekong.lbl.gov");

$h=0;

$dir="output";

foreach $g (0..200) {
	$p=$g*4;
	open T,">rtemp.pov";
	open A,"cylinder.pov";
	while (<A>) {
		last if m/INCLUDE/;
		print T;
	}
	$o=1;
	open B,"cylinder_p.pov";
	while (<B>) {
		if (m/id (\d*)/) {
			$o=$1<$p?1:0;next;
		}
		print T if $o;
	}
	print T "}\n" unless $o;
	close B;
	$o=1;
	open B,"cylinder_v.pov";
	while (<B>) {
		if (m/cell (\d*)/) {
			$o=$1<$p?1:0;next;
		}
		print T if $o;
	}
	print T "}\n" unless $o;
	close B;
	while (<A>) {
		$ty=$g/2;
		s/ROT/$ty/g;
		print T;
	}
	close T;
	print "Rendering movie frames\n";
	print "Frame $f, timestep $ts to $nodes[$h]\n";
	`rsync -z rtemp.pov $nodes[$h]:cgran/rtemp$h.pov`;
	$nice=(@nodes[$h]=~m/yuba/)?"nice -n 19":"nice +19";
	exec "ssh @nodes[$h] \"cd cgran;$nice povray +SU +Ofr_$g.png +W200 +H320 +A0.0001 +R9 -J rtemp$h.pov\" >/dev/null 2>/dev/null; rsync -rz @nodes[$h]:cgran/fr_$g.png . ; ssh @nodes[$h] \"rm cgran/fr_$g.png\"\n" if (($pid[$h]=fork)==0);
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
