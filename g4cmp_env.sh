# $Id: 035eab37646641ac52d1a4c35a8cd0540867c993 $
#
# Configure's users Bourne shell (/bin/sh) or BASH environment for G4CMP
#
# Usage: . g4cmp_env.sh

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
  echo "ERROR: g4cmp_env.sh could self-locate G4CMP installation."
  echo "Please cd to the installation area and source script again."
  return 1
fi

# Assign environment variables for runtime configuraiton

export G4LATTICEDATA=$G4CMPINSTALL/CrystalMaps
export G4ORDPARAMTABLE=$G4CMPINSTALL/G4CMPOrdParamTable.txt
