SUBDIRS = src

ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = README.md rational.pc.in

DISTCLEANFILES = rational.pc

doc_DATA = README.md

pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = rational.pc

AM_DISTCHECK_CONFIGURE_FLAGS = --enable-cln --enable-gmp

if WITH_CLN
PKGCFG_WITH_CLN = cln >= 1.3.2
endif

if WITH_GMPXX
PKGCFG_GMP_CFLAGS = ${GMP_CFLAGS}
PKGCFG_GMP_LIBS = ${GMP_LIBS}
endif

if WITH_INFINT
PKGCFG_GMP_CFLAGS = ${INFINT_CFLAGS}
endif

rational.pc: $(top_srcdir)/rational.pc.in
	$(AM_V_GEN)$(SED) \
	-e 's|@requires_cln[@]|$(PKGCFG_WITH_CLN)|g' \
	-e 's|@INFINT_CFLAGS[@]|$(PKGCFG_INFINT_CFLAGS)|g' \
	-e 's|@GMP_CFLAGS[@]|$(PKGCFG_GMP_CFLAGS)|g' \
	-e 's|@GMP_LIBS[@]|$(PKGCFG_GMP_LIBS)|g' \
	-e 's|@prefix[@]|$(prefix)|g' \
	-e 's|@exec_prefix[@]|$(exec_prefix)|g' \
	-e 's|@PACKAGE_VERSION[@]|$(PACKAGE_VERSION)|g' \
	-e 's|@includedir[@]|$(includedir)|g' \
	-e 's|@libdir[@]|$(libdir)|g' < $< > $@
