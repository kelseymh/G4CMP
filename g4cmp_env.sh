# $Id: 035eab37646641ac52d1a4c35a8cd0540867c993 $
#
# Configure's users Bourne shell (/bin/sh) or BASH environment for G4CMP
#
# Usage: . g4cmp_env.sh
#
# 20161006  Add G4WORKDIR to (DY)LD_LIBRARY_PATH
# 20170509  Define G4CMPLIB and G4CMPINCLUDE relative to G4CMPINSTALL
# 20200719  Set undefined *LD_LIBRARY_PATH; use $() in place of ``
# 20220331  G4CMP-293: Remove G4CMPORDPARAMTABLE; not using RegisterProcess()

# Identify location of script from user command (c.f. geant4make.sh)

if [ -z "$G4CMPINSTALL" ]; then
  if [ -n "$BASH_VERSION" ]; then
    g4cmp_dir=$(dirname ${BASH_ARGV[0]})
    export G4CMPINSTALL=$(cd $g4cmp_dir > /dev/null; pwd)
  elif [ -f g4cmp_env.sh ]; then
    export G4CMPINSTALL=$(pwd)
  fi
fi

# Ensure that G4CMP installation is known

if [ -z "$G4CMPINSTALL" ]; then
  echo "ERROR: g4cmp_env.sh could not self-locate G4CMP installation."
  echo "Please cd to the installation area and source script again."
  return 1
fi

# If running script from source directory, assume GMake build

if [ -r $G4CMPINSTALL/README.md ]; then
  export G4CMPLIB=$G4WORKDIR/lib/$G4SYSTEM
  export G4CMPINCLUDE=$G4CMPINSTALL/library/include
elif [ $(basename $(dirname $G4CMPINSTALL)) = "share" ]; then
  topdir=$(dirname $(dirname $G4CMPINSTALL))
  export G4CMPLIB=$topdir/lib
  export G4CMPINCLUDE=$topdir/include/G4CMP
fi

# Extend library path to include G4CMP library location

export LD_LIBRARY_PATH=${G4CMPLIB}${LD_LIBRARY_PATH:+:}$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=${G4CMPLIB}${DYLD_LIBRARY_PATH:+:}$DYLD_LIBRARY_PATH

# Assign environment variables for runtime configuraiton

export G4LATTICEDATA=$G4CMPINSTALL/CrystalMaps
