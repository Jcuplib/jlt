TOPDIR = .
include Mkinclude

all: setup libjlt tests

setup:
	test -d $(BINDIR) || mkdir -p $(BINDIR)
	test -d $(LIBDIR) || mkdir -p $(LIBDIR)
	test -d $(MODDIR) || mkdir -p $(MODDIR)

libjlt: setup
	$(MAKE) -C $(JLTSRCDIR) TOPDIR=$(abspath $(TOPDIR)) all

tests: libjlt
	$(MAKE) -C src

allclean: clean
	$(RM) -rf $(BINDIR)
	$(RM) -rf $(LIBDIR)
	$(RM) -rf $(MODDIR)

clean:
	$(MAKE) -C src clean
	$(MAKE) -C $(JLTSRCDIR) TOPDIR=$(abspath $(TOPDIR)) clean
