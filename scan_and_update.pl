#!/usr/bin/perl
open A,"$ARGV[0]";
open B,">temp";

$c=0;
while(<A>) {
	$c++ if s/[ \t]+$//;
	if(/^\/\/ Date/) {
		print B "// Date     : August 30th 2011\n";
		next;
	}
	print B;
}
close A;
close B;
`mv temp $ARGV[0]`;

print "$c tabs/spaces removed from $ARGV[0]\n";
