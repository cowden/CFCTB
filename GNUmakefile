# $Id: GNUmakefile,v 1.2 2003-01-23 15:31:39 maire Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := bigModule
G4TARGET := $(name)
G4EXLIB := true

CPPFLAGS += $(shell root-config --cflags) -I$(EXTERN_BASE)/include -g
EXTRALIBS = $(shell root-config --glibs) -L$(EXTERN_BASE)/lib -lHepMC -lCLHEP -g

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
