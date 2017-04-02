CC=gcc
FC=gfortran
AR=ar
FFLAGS=-fPIC -Ofast
ARFLAGS         =      ru
INCLUDES = 

F2PY_FLAGS= -I./ 

SRCS = cs_stat.f90 ARWpost.f90

F2PY = f2py
OBJS = $(SRCS:.f90=.so)
all:    $(OBJS)
	@echo  ${MODULE}.so has been compiled

%.so: %.f90 
	$(F2PY) $(F2PY_FLAGS) --f90flags=-fPIC  -m $* -c $< 



#$(OBJS) #${MODULE}.F90 


# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
#.F90.o:
#	$(FC) $(FFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.a *.o *.mod  *.so

#depend: $(SRCS)
#	makedepend $(INCLUDES) $^


# DO NOT DELETE THIS LINE -- make depend needs it
