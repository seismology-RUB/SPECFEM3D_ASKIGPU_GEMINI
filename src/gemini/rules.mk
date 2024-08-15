#-------------------------------------------------------------
#   Makefile for geminiUnified
#----------------------------------------------------------------------------
#   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
#
#   This file is part of GEMINI_UNIFIED version 1.0.
#
#   GEMINI_UNIFIED version 1.0 is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   GEMINI_UNIFIED version 1.0 is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with GEMINI_UNIFIED version 1.0.  If not, see <http://www.gnu.org/licenses/>.
#---------------------------------------------------------------------------------------
## compilation directories
S := ${S_TOP}/src/gemini
$(gemini_OBJECTS): S = ${S_TOP}/src/gemini

#######################################

gemini_TARGETS = \
	$(gemini_OBJECTS) \
	$(gemini_SHARED_OBJECTS) \
	$(EMPTY_MACRO)


gemini_OBJECTS = \
	$O/anyRankArray.gem.o \
	$O/anyRankRealArray.gem.o \
	$O/anyRankIntegerArray.gem.o \
	$O/axesRotation.gem.o \
	$O/hdfWrapper.gemmpi.o \
	$O/errorMessage.gem.o \
	$O/mathConstants.gem.o \
	$O/string.gem.o \
	$O/realloc.gem.o \
	$O/travelTimes.gem.o \
	$O/locatePoint.gem.o \
	$(EMPTY_MACRO)


gemini_MODULES = \
	$(FC_MODDIR)/anyrankarray.$(FC_MODEXT) \
	$(FC_MODDIR)/anyrankrealarray.$(FC_MODEXT) \
	$(FC_MODDIR)/anyrankintegerarray.$(FC_MODEXT) \
	$(FC_MODDIR)/axesrotation.$(FC_MODEXT) \
	$(FC_MODDIR)/hdfwrapper.$(FC_MODEXT) \
	$(FC_MODDIR)/errormessage.$(FC_MODEXT) \
	$(FC_MODDIR)/mathconstants.$(FC_MODEXT) \
	$(FC_MODDIR)/string.$(FC_MODEXT) \
	$(FC_MODDIR)/realloc.$(FC_MODEXT) \
	$(FC_MODDIR)/traveltimes.$(FC_MODEXT) \
	$(FC_MODDIR)/locatepoint.$(FC_MODEXT) \
	$(EMPTY_MACRO)


###
### MPI
###
gemini_OBJECTS += $(COND_MPI_OBJECTS)
ifeq ($(MPI),yes)
gemini_MODULES += $(FC_MODDIR)/my_mpi.$(FC_MODEXT)
endif

####
#### rule for each .o file below
####

##
## gemini dependencies
##
$O/anyRankRealArray.gem.o: $O/anyRankArray.gem.o
$O/anyRankIntegerArray.gem.o: $O/anyRankArray.gem.o
$O/hdfWrapper.gemmpi.o: $O/anyRankRealArray.gem.o $O/anyRankIntegerArray.gem.o $O/errorMessage.gem.o $O/string.gem.o 
$O/string.gem.o: $O/errorMessage.gem.o
$O/errorMessage.gem.o: $O/realloc.gem.o
$O/axesRotation.gem.o: $O/mathConstants.gem.o
$O/travelTimes.gem.o: $O/anyRankRealArray.gem.o $O/hdfWrapper.gemmpi.o $O/errorMessage.gem.o $O/locatePoint.gem.o

##
##  rules for compilation
##

$O/%.gem.o: $S/%.f90 
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.gemmpi.o: $S/%.f90 
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -I$(HDF_INCLUDES) -c -o $@ $<
