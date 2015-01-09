#!/usr/bin/awk -f
#
# Usage: vdrift.awk <epositions.txt>
#
# simple analysis of epositions.txt to compute drift speed using average
# of drift distances and times.  Assumes particles start at (0,0,0).
#
# 20150109  Michael Kelsey

BEGIN { Ne=0; De=0; Te=0; Nh=0; Dh=0; Th=0; }

/^1 /  { Nh++; Dh += sqrt($2*$2+$3*$3+$4*$4); Th += $5; }
/^-1 / { Ne++; De += sqrt($2*$2+$3*$3+$4*$4); Te += $5; }

END {
  if (Nh>0) {
    Dh /= Nh; Th /= Nh;
    Vh = 1e6 * Dh/Th;
    print "    holes:",Dh*100,"cm @",Th,"ns:",Vh,"km/s";
  }

  if (Ne>0) {
    De /= Ne; Te /= Ne;
    Ve = 1e6 * De/Te;
    print "electrons:",De*100,"cm @",Te,"ns:",Ve,"km/s";
  }
}
