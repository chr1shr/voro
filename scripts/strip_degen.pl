#!/usr/bin/perl
# This perl script can remove degenerate cylinders from the POV files.
# Ideally, the code just wouldn't include these in the first place.
use File::Copy;

$a=0;
open A,@ARGV[0];
open B,">@ARGV[0].temp";
while (<A>) {
	m/^cylinder{<([-\.\d]+),([-\.\d]+),([-\.\d]+)>,<\1,\2,\3>,r}/?$a++:(print B);
}
print "$a cylinders stripped\n";
move("@ARGV[0].temp",@ARGV[0]);
