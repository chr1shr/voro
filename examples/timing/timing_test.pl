@range=(10..40);
$tries=3;

foreach $r (@range) {
	system "g++ -I../../src -DNNN=$r -o timing_test -O3 -fast -pedantic -Winline --param large-function-growth=1000 --param max-inline-insns-single=2000 -Wall timing_test.cc";
	$st=$stt=0;
	foreach $t (1..$tries) {
		system "./timing_test >time_temp";
		open F,"time_temp" or die "Can't open timing file: $!";
		($t)=split ' ',<F>;
		$st+=$t;$stt+=$t*$t;
		close F;
	}
	$st/=$tries;
	$stt=$stt/$tries-$st*$st;$stt=$stt>0?sqrt($stt):0;
	print "$r $st $stt\n";
}
