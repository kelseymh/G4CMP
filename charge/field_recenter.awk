#!/usr/bin/awk -f
#
# Usage: centerfield <Epot-file> > <new-file>
#
# Takes a CDMS DMC-format electric field file, and shifts the Z coordinate
# by -1.27 cm.  This puts the coordinate origin at the center of the field,
# as expected for GEANT4 use.
#
# $Id$
# 20140522  Michael Kelsey

### Every line is <x> <y> <z> <V>, with dimensions in meters
{
    x = $1; y = $2; z = $3 - 0.0127; V = $4;
    printf "%23.16e %23.16e %23.16e %23.16e\n",x,y,z,V;
}
