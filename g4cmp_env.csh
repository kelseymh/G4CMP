# $Id: 3ddff4078d93e1dac811f0aeac166bf1ac512cd0 $
#
# Configure's users C shell (/bin/csh) or TCSH environment for G4CMP
#
# Usage: source g4cmp_env.csh

# Identify location of script from user command (c.f. geant4make.csh)

if (! $?G4CMPINSTALL) then
  unset g4cmp_dir

  set ARGS = ($_)		# Pull directory from terminal command
  if ("$ARGS" != "") then
    if ("$ARGS[2]" =~ */g4cmp_env.csh) then
      setenv G4CMPINSTALL `dirname $ARGS[2]`
    endif
  endif

  if (! $?G4CMPINSTALL) then	# Sourced within script, command not avail
    if (-e g4cmp_env.csh) then
      setenv G4CMPINSTALL `pwd`
    else if ("$1" != "") then	# TCSH allows arguments to source
      if (-e $1/g4cmp_env.csh) then
        setenv G4CMPINSTALL $1
      endif
    endif
  endif
endif

# Ensure that G4CMP installation is known

if (! $?G4CMPINSTALL) then
  echo "ERROR: g4cmp_env.sh could self-locate G4CMP installation."
  echo "Please cd to the installation area and source script again."
  return 1
endif

# Extend library path to include G4CMP library location

set g4cmplib = $G4WORKDIR/lib/$G4SYSTEM
if ($?LD_LIBRARY_PATH) then
  setenv LD_LIBRARY_PATH ${g4cmplib}:$LD_LIBRARY_PATH
endif
if ($?DYLD_LIBRARY_PATH) then
  setenv DYLD_LIBRARY_PATH ${g4cmplib}:$DYLD_LIBRARY_PATH
endif

# Assign environment variables for runtime configuraiton

setenv G4LATTICEDATA $G4CMPINSTALL/CrystalMaps
setenv G4ORDPARAMTABLE $G4CMPINSTALL/G4CMPOrdParamTable.txt
