# $Id: 3ddff4078d93e1dac811f0aeac166bf1ac512cd0 $
#
# Configure's users C shell (/bin/csh) or TCSH environment for G4CMP
#
# Usage: source g4cmp_env.csh
#
# 20170509  Define G4CMPLIB and G4CMPINCLUDE relative to G4CMPINSTALL
# 20180917  Fix initialization of G4CMPINSTALL to use absolute path
# 20200719  Set undefined *LD_LIBRARY_PATH
# 20220331  G4CMP-293: Remove G4CMPORDPARAMTABLE; not using RegisterProcess()

# Identify location of script from user command (c.f. geant4make.csh)

if (! $?G4CMPINSTALL) then
  set ARGS = ($_)		# Pull directory from terminal command
  if ("$ARGS" != "") then
    if ("$ARGS[2]" =~ */g4cmp_env.csh) then
      set g4cmp_dir = `dirname $ARGS[2]`
      setenv G4CMPINSTALL `cd $g4cmp_dir; pwd`
    endif
  endif

  if (! $?G4CMPINSTALL) then	# Sourced within script, command not avail
    if (-e g4cmp_env.csh) then
      setenv G4CMPINSTALL `pwd`
    else if ("$1" != "") then	# TCSH allows arguments to source
      if (-e $1/g4cmp_env.csh) then
        set g4cmp_dir = `cd $1; pwd`
        setenv G4CMPINSTALL $1
      endif
    endif
  endif

  unset g4cmp_dir
endif

# Ensure that G4CMP installation is known

if (! $?G4CMPINSTALL) then
  echo "ERROR: g4cmp_env.sh could not self-locate G4CMP installation."
  echo "Please cd to the installation area and source script again."
  return 1
endif

# If running script from source directory, assume GMake build

if (-r $G4CMPINSTALL/README.md) then
  setenv G4CMPLIB     $G4WORKDIR/lib/$G4SYSTEM
  setenv G4CMPINCLUDE $G4CMPINSTALL/library/include
else if (`dirname $G4CMPINSTALL|xargs basename` == "share") then
  set topdir = `dirname $G4CMPINSTALL|xargs dirname`
  setenv G4CMPLIB     $topdir/lib
  setenv G4CMPINCLUDE $topdir/include/G4CMP
  unset topdir
endif

# Extend library path to include G4CMP library location

if ($?LD_LIBRARY_PATH) then
  setenv LD_LIBRARY_PATH ${G4CMPLIB}:$LD_LIBRARY_PATH
else
  setenv LD_LIBRARY_PATH ${G4CMPLIB}
endif

if ($?DYLD_LIBRARY_PATH) then
  setenv DYLD_LIBRARY_PATH ${G4CMPLIB}:$DYLD_LIBRARY_PATH
else
  setenv DYLD_LIBRARY_PATH ${G4CMPLIB}
endif

# Assign environment variables for runtime configuraiton

setenv G4LATTICEDATA   $G4CMPINSTALL/CrystalMaps
