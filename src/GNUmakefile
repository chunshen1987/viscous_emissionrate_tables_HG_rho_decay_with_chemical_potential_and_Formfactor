# ===========================================================================
#  Makefile iSS                                    Chun Shen Mar. 19, 2013
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install	make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := g++
CFLAGS= -O3 $(shell gsl-config --cflags)

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS) $(shell gsl-config --libs)
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	cal_photon_emissionrate_tables.e
endif

SRC		=	main.cpp Arsenal.cpp gauss_quadrature.cpp \
                  Table2D.cpp chemical_potential.cpp Formfactor.cpp \
                  HG_1to3_decay.cpp ParameterReader.cpp Table.cpp

INC		= 	Arsenal.h Physicalconstants.h Stopwatch.h \
                  gauss_quadrature.h Table2D.h chemical_potential.h \
                  Formfactor.h HG_1to3_decay.h ParameterReader.h \
                  Table.h

# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	../

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(LDFLAGS) $(OBJECTS) -o $(TARGET)
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
main.cpp : Arsenal.h Stopwatch.h gauss_quadrature.h Physicalconstants.h HG_1to3_decay.h ParameterReader.h
Arsenal.cpp : gauss_quadrature.h
Table2D.cpp : Arsenal.h
chemical_potential.cpp : Arsenal.h Table2D.h Physicalconstants.h
Formfactor.cpp : Physicalconstants.h
ParameterReader.cpp : Arsenal.h
HG_1to3_decay.cpp : ParameterReader.h Arsenal.h Physicalconstants.h gauss_quadrature.h chemical_potential.h
