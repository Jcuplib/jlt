TOPDIR = ..
include $(TOPDIR)/Mkinclude

TESTSRCDIR = $(TOPDIR)/src

all: setup libjlt tests

setup:
	test -d $(LIBDIR) || mkdir -p $(LIBDIR)
	test -d $(MODDIR) || mkdir -p $(MODDIR)

libjlt: setup
	$(MAKE) -C $(JLTSRCDIR) TOPDIR=$(abspath $(TOPDIR)) all

tests: libjlt
	$(MAKE) -C $(TESTSRCDIR)
	cp $(TESTSRCDIR)/comp1 .
	cp $(TESTSRCDIR)/comp2 .

clean:
	$(MAKE) -C $(TESTSRCDIR) clean
	-$(RM) -f comp1 comp2

distclean: clean
	$(MAKE) -C $(JLTSRCDIR) TOPDIR=$(abspath $(TOPDIR)) clean
	$(RM) -rf $(LIBDIR) $(MODDIR)
