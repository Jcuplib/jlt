TOPDIR = ../..
include ${TOPDIR}/Mkinclude
include files.mk

JLTMODS = $(JLTOBJS:.o=.mod)

test:
	$(MAKE) $(JLTLIB)
	$(MAKE) modules
	$(FC) prg_remapping_table.f90 -I$(MODDIR) $(JLTLIB)

all:
	$(MAKE) $(JLTLIB)
	$(MAKE) modules
	@echo "Complete making JcupLT."


$(JLTLIB): $(JLTOBJS)
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@

modules: $(JLTOBJS)
	$(INSTALL) $(JLTMODS) $(MODDIR)

allclean: clean
	$(RM) -f $(LIBDIR)/$(LIBJLT)

clean:
	$(RM) -f *.o *.mod *.lst *.L *.L


.SUFFIXES:
.SUFFIXES: .o .f90 .mod

.f90.o:
	$(FC) $(FFLAGS) -c $< -I$(MODDIR)

%.mod: %.f90
	make $(patsubst %.f90,%.o,$<)
