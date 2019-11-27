#!/bin/env make
# 
# created by Sebastian Ehlert in March 2018 (v 2.1)
# 2.0: recreated from scratch
# 2.1: refined the generation of dependencies, since it crashed if not
#      done in the main directory
#
###------------------------------------------------------------------------###
# a nice trick to print all 'important' targets this Makefile can produce
.PHONY: help

# would be made by default if there is no .DEFAULT_GOAL set
help:
	@echo "This Makefile kindly provides to you:"                     && \
	$(MAKE) --print-data-base --question |                               \
	awk '/^[^.%][-A-Za-z0-9]*:/ {print substr($$1, 1, length($$1)-1)}' | \
	sort | pr --omit-pagination --width=80 --columns=4

##############################################################################
## this file is part of the xtb build system                                ##
## for more information on this have a look in the AK-Wiki or               ##
## consult the README in MAKE/                                              ##
##############################################################################

###------------------------------------------------------------------------###
# specifiy where you are in the source tree
MAINDIR := .
# here resides the build system
MAKEDIR := $(MAINDIR)/MAKE

include $(MAKEDIR)/Makeconfig

##############################################################################
## put all programs here, we use this silly syntax for you, to comment
## out parts of it more easily
PROGRAMS := xtb4stda
###------------------------------------------------------------------------###
# initialized empty, then incremental build up for easy commenting
SCRIPTS := 
###------------------------------------------------------------------------###
# initialized empty, then incremental build up for easy commenting
PARAMETER :=
PARAMETER += .param_gbsa_acetone
PARAMETER += .param_gbsa_acetonitrile
PARAMETER += .param_gbsa_benzene
PARAMETER += .param_gbsa_ch2cl2
PARAMETER += .param_gbsa_chcl3
PARAMETER += .param_gbsa_cs2
PARAMETER += .param_gbsa_dmso
PARAMETER += .param_gbsa_ether
PARAMETER += .param_gbsa_h2o
PARAMETER += .param_gbsa_methanol
PARAMETER += .param_gbsa_thf
PARAMETER += .param_gbsa_toluene
PARAMETER += .param_stda1.xtb
PARAMETER += .param_stda2.xtb
##############################################################################
## only edit below if you know what you are doing                           ##
##############################################################################

###------------------------------------------------------------------------###
# set the targets
.PHONY: all main libs progs pack setup depend clean distclean veryclean
.PHONY: $(PROGRAMS)

# this is the first, so it would be made by default
all: setup libs progs

# it may not be desireable to make all as default,
# therefore, we set xtb as default
.DEFAULT_GOAL = main
main: setup $(EXEDIR)/xtb4stda

###------------------------------------------------------------------------###
# generate the name of the program and the library
LIBNAME  := $(patsubst %,$(LIBDIR)/lib%.a,$(PROGRAMS))
PROGNAME := $(patsubst %,$(EXEDIR)/%,$(PROGRAMS))
# this might be buggy, so don't use
#ifeq ($(BUILD_IN_BIN),yes)
#   PROGNAME := $(patsubst %, $(HOME)/bin/$(NAME), $(PROGRAMS))
#endif # BIN
# everything should also dependent on the local Makefile
ifeq ($(MAKEFILE_DEP),yes)
   MAKEDEP := Makefile
endif # DEP

###------------------------------------------------------------------------###
# there is also the possibility to pack the xtb via the Makefile
pack: all
	tar czvf xtb4stda-$$(date +'%y%m%d').tgz \
		$(PROGNAME)\
		$(patsubst %,$(MAINDIR)/scripts/%,$(SCRIPTS))\
		$(PARAMETER)\
		README.md\
		COPYING\
		COPYING.LESSER\
		.xtb4stdarc

###------------------------------------------------------------------------###
# set up all directories
setup: $(EXEDIR) $(LIBDIR) $(MODDIR) $(OUTDIR) $(DEPDIR)

$(EXEDIR):
	@echo "creating directory for executables"  && \
	mkdir -p $@

$(LIBDIR):
	@echo "creating directory for libraries"    && \
	mkdir -p $@

$(MODDIR):
	@echo "creating directory for modules"      && \
	mkdir -p $@

$(OUTDIR):
	@echo "creating directory for object files" && \
	mkdir -p $@

$(DEPDIR):
	@echo "creating directory for dependencies" && \
	mkdir -p $@

###------------------------------------------------------------------------###
# set up the targets for all
program: $(PROGNAME)

library: $(LIBNAME)

###------------------------------------------------------------------------###
# now the rules to actually start the submakefiles
$(EXEDIR)/xtb4stda:
	$(MAKE) -C src program

$(LIBDIR)/libxtb4stda.a:
	$(MAKE) -C src library

###------------------------------------------------------------------------###
# clean up, at the right places
clean:
	$(RM) $(OUTDIR)/*.o

distclean: veryclean
	$(RM) $(MODDIR)/*.mod $(PROGNAME) $(LIBNAME)

# in case you use the dependency build you can get rid of those files to
veryclean: clean
	$(RM) $(DEPDIR)/*.d

###------------------------------------------------------------------------###
