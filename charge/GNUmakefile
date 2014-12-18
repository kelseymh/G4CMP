# G4CMP/charge/GNUmakefile -- for building charge-carrier propagation demo
# $Id$
#
# 20140509  Support different executables by looking at target name

ifndef G4CMP_NAME		# Two-pass method
  .PHONY : g4cmpCharge luke_dist electron_Epv

  g4cmpCharge luke_dist electron_Epv :
	@$(MAKE) G4CMP_NAME=$@ bin

  G4CMP_NAME := g4cmpCharge
endif

include ../g4cmp.gmk
