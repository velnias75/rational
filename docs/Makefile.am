if HAVE_DOXYGEN
doxyfile.stamp: Doxyfile
	$(DOXYGEN) Doxyfile
	echo Timestamp > doxyfile.stamp

CLEANFILES = doxyfile.stamp

MAINTAINER_CLEANFILES = Doxyfile

all-local: doxyfile.stamp

install-data-local: doxyfile.stamp
	mkdir -p $(DESTDIR)$(docdir)/html
	find $(top_builddir)/docs/html -type f -exec $(INSTALL) -m 644 "{}" $(DESTDIR)$(docdir)/html \;

uninstall-local:
	rm -rf $(DESTDIR)$(docdir)/html

clean-local:
	rm -rf $(top_builddir)/docs/html
endif
