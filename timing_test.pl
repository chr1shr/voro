@source=("worklist.cc","worklist.hh","wall_test.cc","wall.cc","cell.cc","container.cc","container.hh","config.hh","cell.hh","wall.hh");
@range=(10..40);
$tries=1;
$compile='g++ -o wall_test -O3 wall_test.cc -pedantic -Winline --param large-function-growth=600 --param max-inline-insns-single=1000 -fast -Wall';
$executable="./wall_test";

foreach $r (@range) {
	foreach $s (@source) {
		open S,$s or die "Can't open source file: $!";
		open T,">.timingtemp/$s" or die "Can't open temporary file: $!";
		while(<S>) {
			s/NNN/$r/g;
			print T;
		}
		close T;
		close S;
	}
	
	chdir ".timingtemp";
	system $compile;
	$tr=$tu=$ts=0;
	$tr2=$tu2=$ts2=0;
	foreach $t (1..$tries) {
		system "eval time $executable 2>ti";
		open F,"ti" or die "Can't open timing file: $!";<F>;
		<F>=~m/(\d*)m([\d\.]*)s/;$rr=$1*60+$2;$tr+=$rr;$tr2+=$rr*$rr;
		<F>=~m/(\d*)m([\d\.]*)s/;$ru=$1*60+$2;$tu+=$ru;$tu2+=$ru*$ru;
		<F>=~m/(\d*)m([\d\.]*)s/;$rs=$1*60+$2;$ts+=$rs;$ts2+=$rs*$rs;
	}
	$tr/=$tries;$tu/=$tries;$ts/=$tries;
	$tr2=$tr2/$tries-$tr*$tr;$tr2=$tr2>0?sqrt($tr2):0;
	$tu2=$tu2/$tries-$tu*$tu;$tu2=$tu2>0?sqrt($tu2):0;
	$ts2=$ts2/$tries-$ts*$ts;$ts2=$ts2>0?sqrt($ts2):0;
	chdir "..";
	print "$r $tr $tu $ts $tr2 $tu2 $ts2\n";
}
