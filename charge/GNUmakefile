# G4CMP/charge/GNUmakefile -- for building charge-carrier propagation demo
# $Id$
#
# 20140509  Support different executables by looking at target name
# 20160518  Use G4CMPINSTALL instead of ".." to find includes

ifndef G4CMP_NAME		# Two-pass method
  .PHONY : g4cmpCharge luke_dist electron_Epv ionize latticeVecs

  g4cmpCharge luke_dist electron_Epv ionize latticeVecs :
	@$(MAKE) G4CMP_NAME=$@

  G4CMP_NAME := g4cmpCharge
endif

include $(G4CMPINSTALL)/g4cmp.gmk
