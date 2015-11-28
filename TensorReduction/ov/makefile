# Makefile routine.
VPATH		= $(PWD)/obj
INCPATH  	= $(PWD)/../Include
OUTPUT_OPTION	= -o $(PWD)/obj/$@

FC = gfortran
FFLAGS 	= -O2 -g -I$(INCPATH) -Iobj -Jobj

LIBRARY = libov.a

# -----------------------------------------------------------------------------
# Specify the object files. 

OBJS = \
TRsettensorcontrol.o \
del2.o \
del3.o \
del4.o \
doBtensor.o \
doCtensor.o \
doDtensor.o \
ludcmp.o \
mprove.o \
ovArraysetup.o \
ovBtensor.o \
ovCcheck.o \
ovCtensor.o \
ovDcheck.o \
ovDtensor.o \
ovE0scalar.o \
ovEtensor.o \
ovIpiksrt.o \
ovdot.o \
ovdsvdcmp.o \
ovindex.o \
ovdpythag.o \
ovswitch.o \
ovw2.o \
ovw3.o \
pvfndd.o \
zlubksb.o \
zsvbksb.o
      
      
$(LIBRARY):	$(OBJS) 
		ar cru $(LIBRARY) $(patsubst %,obj/%,$(OBJS))
		ranlib $(LIBRARY)
 
# -----------------------------------------------------------------------------
# Specify other options.


clean:
	- rm -f *.o obj/*.o *.s *.prj *~ core

# -----------------------------------------------------------------------------

# DO NOT DELETE


# LocalWords:  Makefile
