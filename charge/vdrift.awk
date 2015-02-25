#!/usr/bin/awk -f
#
# Usage: vdrift.awk <epositions.txt>
#
# simple analysis of epositions.txt to compute drift speed using average
# of drift distances and times.  Assumes particles start at (0,0,0).
#
# 20150109  Michael Kelsey
# 20150219  Deal with new, more complicated format for epositions.txt

BEGIN { FS=","; Ne=0; De=0; Ze=0; Te=0; Nh=0; Dh=0; Zh=0; Th=0; }

$2==1 {
    Nh++; x=$9-$6; y=$10-$7; z=$11-$8; Th += $4;
    Dh += sqrt(x*x+y*y+z*z); Zh += abs(z);
}
$2==-1 {
    Ne++; x=$9-$6; y=$10-$7; z=$11-$8; Te += $4;
    De += sqrt(x*x+y*y+z*z); Ze += abs(z);
}

function abs(x) { return x<0?-x:x; }

END {
  if (Nh>0) {
    Dh /= Nh; Th /= Nh; Zh /= Nh;
    Vh = 1e6 * Dh/Th;
    Vzh = 1e6 * abs(Zh)/Th;
    print "    holes:",Dh*100,"cm @",Th,"ns:",Vh,"km/s (along z",Vzh,")";
  }

  if (Ne>0) {
    De /= Ne; Te /= Ne; Ze /= Ne;
    Ve = 1e6 * De/Te;
    Vze = 1e6 * abs(Ze)/Te;
    print "electrons:",De*100,"cm @",Te,"ns:",Ve,"km/s (along z",Vze,")";
  }
}
