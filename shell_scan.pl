#!/usr/bin/perl
$r=6;
$ls=92;$pad=$ls*4;
$d=8;
use Math::Trig;

$maxlen=0;
for($kk=0;$kk<6;$kk++) {
	for($jj=0;$jj<6;$jj++) {
		for($ii=0;$ii<6;$ii++) {
			$in=$ii+$r*($jj+$r*$kk);
			$len[$in]=worklist($ii,$jj,$kk);
			$maxlen=$len[$in] if $len[$in]>$maxlen;
		}
	}
}

open W,">worklist.cc";
print W "const int sgrid=$r;\n";
print W "const int seq_length=$maxlen;\n";
print W "static const wl[]=(\n";
foreach (0..($r*$r*$r-1)) {
	print W "\t";
	for($q=0;$q<$maxlen;$q++) {
		print W ($q<$len[$_]?"$c[$_][$q]":"0").($q==$maxlen-1?"\n":",");
	}
}
print W ");\n";
close W;

sub worklist {
	print "@_[0] @_[1] @_[2]\n";
	$ind=@_[0]+$r*(@_[1]+$r*@_[2]);
	$ac=0;$v++;
	$xp=$yp=$zp=0;
	$x=(@_[0]+0.5)/$r;
	$y=(@_[1]+0.5)/$r;
	$z=(@_[2]+0.5)/$r;
	$m[$d][$d][$d]=$v;
	add(1,0,0);add(0,1,0);add(0,0,1);
	add(-1,0,0);add(0,-1,0);add(0,0,-1);
	foreach $l (1..$ls) {
		$minwei=1e9;
		foreach (0..$ac-1) {
			$xt=@a[3*$_];$yt=@a[3*$_+1];$zt=@a[3*$_+2];
#			$wei=dis($x,$y,$z,$xt,$yt,$zt)+1*acos(($xt*$xp+$yt*$yp+$zt*$zp)/($xt*$xt+$yt*$yt+$zt*$zt)*($xp*$xp+$yp*$yp+$zp*$zp));
			$wei=adis($x,$y,$z,$xt,$yt,$zt)+0.1*sqrt(($xt-$xp)**2+($yt-$yp)**2+($zt-$zp)**2);
			$nx=$_,$minwei=$wei if $wei<$minwei;
		}
		$xp=@a[3*$nx];$yp=@a[3*$nx+1];$zp=@a[3*$nx+2];
		add($xp+1,$yp,$zp);add($xp,$yp+1,$zp);add($xp,$yp,$zp+1);
		add($xp-1,$yp,$zp);add($xp,$yp-1,$zp);add($xp,$yp,$zp-1);
		push @b,(splice @a,3*$nx,3);$ac--;
	}
	$v++;
	for($i=0;$i<$#b;$i+=3) {
		$xt=@b[$i];$yt=@b[$i+1];$zt=@b[$i+2];
		$m[$d+$xt][$d+$yt][$d+$zt]=$v;
	}
	$m[$d][$d][$d]=$v;
	for($i=0;$i<$#b;$i+=3) {
		$xt=@b[$i];$yt=@b[$i+1];$zt=@b[$i+2];
		last if $m[$d+$xt+1][$d+$yt][$d+$zt]!=$v;
		last if $m[$d+$xt][$d+$yt+1][$d+$zt]!=$v;
		last if $m[$d+$xt][$d+$yt][$d+$zt+1]!=$v;
		last if $m[$d+$xt-1][$d+$yt][$d+$zt]!=$v;
		last if $m[$d+$xt][$d+$yt-1][$d+$zt]!=$v;
		last if $m[$d+$xt][$d+$yt][$d+$zt-1]!=$v;	
	}
	$j=0;
	$c[$ind][$j++]=$i;
	while ($#b>0) {
		$i-=3;
		$xt=shift @b;$yt=shift @b;$zt=shift @b;
		$o=0;$oo=0;
		$oo++,$o|=1 if $m[$d+$xt+1][$d+$yt][$d+$zt]!=$v;
		$oo++,$o|=2 if $m[$d+$xt][$d+$yt+1][$d+$zt]!=$v;
		$oo++,$o|=4 if $m[$d+$xt][$d+$yt][$d+$zt+1]!=$v;
		$oo++,$o|=8 if $m[$d+$xt-1][$d+$yt][$d+$zt]!=$v;
		$oo++,$o|=16 if $m[$d+$xt][$d+$yt-1][$d+$zt]!=$v;
		$oo++,$o|=32 if $m[$d+$xt][$d+$yt][$d+$zt-1]!=$v;
		$c[$ind][$j++]=$xt;
		$c[$ind][$j++]=$yt;
		$c[$ind][$j++]=$zt;
		$c[$ind][$j++]=$o if $i<0;
	}
	undef @a;
	undef @b;
	return $j;
}



sub add {
	if ($m[$d+@_[0]][$d+@_[1]][$d+@_[2]]!=$v) {
		$ac++;
		push @a,@_[0],@_[1],@_[2];
		$m[$d+@_[0]][$d+@_[1]][$d+@_[2]]=$v;
	}
}

sub dis {
	$xl=@_[3]+0.25-@_[0];$xh=@_[3]+0.75-@_[0];
	$yl=@_[4]+0.25-@_[1];$yh=@_[4]+0.75-@_[1];
	$zl=@_[5]+0.25-@_[2];$zh=@_[5]+0.75-@_[2];
	$dis=(abs($xl)<abs($xh)?$xl:$xh)**2
		+(abs($yl)<abs($yh)?$yl:$yh)**2
		+(abs($zl)<abs($zh)?$zl:$zh)**2;
	return sqrt $dis;
}

sub adis {
	$xl=@_[3]-@_[0];$xh=@_[3]+1-@_[0];
	$yl=@_[4]-@_[1];$yh=@_[4]+1-@_[1];
	$zl=@_[5]-@_[2];$zh=@_[5]+1-@_[2];
	$dis=(abs($xl)<abs($xh)?$xl:$xh)**2
		+(abs($yl)<abs($yh)?$yl:$yh)**2
		+(abs($zl)<abs($zh)?$zl:$zh)**2;
	return sqrt $dis;
}
