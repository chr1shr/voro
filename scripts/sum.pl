#!/usr/bin/perl
open A,"@ARGV[0]" or die "Can't open input file";
while(<A>) {
	@A=split;
	$c+=@A[5];
}
print "$c\n";
