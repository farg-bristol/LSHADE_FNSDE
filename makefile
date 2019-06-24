##########################################################
#VARIABLES

OBJS=	MAIN.o\
	de.o\
        de_config.o\
        niching_func_cons.o\
        perfanal.o

EXEC =  de.exe

FC = gfortran

FCFLAGS = -O2
#FCFLAGS = -O0 -Wall -Wextra -Wconversion -pedantic -fbacktrace -fimplicit-none -fimplicit-none

##########################################################
#BUILD INSTRUCTIONS

EXEC: $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXEC) $(OBJS) $(COBJS) $(CPOBJS)

%.o:%.f90
	$(FC) -c $(FCFLAGS) $<




