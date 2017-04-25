#! /bin/sh

if [ -z ${G4CMPINSTALL+x} ]; then
  #G4CMPINSTALL is not set so let's set it
  if [ -z ${BASH_VERSION+x} ]; then
    #The shell is not bash
    if [ ! -f test.sh ]; then #test if in the same directory as file
      printf "ERROR: %s could not self-locate G4CMP\n" "test.sh";
      printf "installation. This is most likely because you \n";
      printf "are using ksh, zsh or similar. To fix this issue,\n";
      printf "cd to the directory containing this script.\n";
      printf "Then source %s from there.\n" "test.sh";
      return 1
    else
      #Not using bash but in the directory of the script
      G4CMPINSTALL=$(pwd);
    fi
  else
    #The shell is bash so use BASH_SOURCE to locate
    SCRIPT_DIR="$(dirname "$BASH_SOURCE")";
    #the two step process ensures the G4CMP_BIN is not set as
    #relative path.
    G4CMPINSTALL=$(cd $SCRIPT_DIR > /dev/null ; pwd);
  fi
  
  if test "x$PATH" = "x" ; then
    export PATH="$G4CMPINSTALL";
  else
    export PATH="$G4CMPINSTALL":${PATH};
  fi
  
  if test "x$LD_LIBRARY_PATH" = "x" ; then
    export LD_LIBRARY_PATH="`cd $G4CMPINSTALL/../lib > /dev/null ; pwd`";
  else
    export LD_LIBRARY_PATH="`cd $G4CMPINSTALL/../lib > /dev/null ; pwd`":${LD_LIBRARY_PATH};
  fi
  
  if test "x$DYLD_LIBRARY_PATH" = "x" ; then
    export DYLD_LIBRARY_PATH="`cd $G4CMPINSTALL/../lib > /dev/null ; pwd`";
  else
    export DYLD_LIBRARY_PATH="`cd $G4CMPINSTALL/../lib > /dev/null ; pwd`":${DYLD_LIBRARY_PATH};
  fi

  #whatever method was used to set it. Now export it.
  export G4CMPINSTALL;
  export G4LATTICEDATA="`cd $G4CMPINSTALL/../share/G4CMP/CrystalMaps/ > /dev/null ; pwd`";
  export G4ORDPARAMTABLE="$G4CMPINSTALL/G4CMPOrdParamTable.txt";

else
  printf "G4CMPINSTALL is already set to %s.\n" $G4CMPINSTALL;
  printf "This might indicate a problem.\n";
  printf "This script will stop without settting anything.\n";
fi


