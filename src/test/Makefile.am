check_PROGRAMS = test_rational

TESTS = $(check_PROGRAMS)

noinst_HEADERS = rationaltest.h gmptest.h clntest.h infinttest.h exprtest.h exprtest_cln.h \
	exprtest_infint.h unlogicaltest.h

test_rational_CXXFLAGS = -fexceptions -finline-functions -I$(top_srcdir)/src/rational \
	--param large-function-growth=400 $(CPPUNIT_CFLAGS)
test_rational_CPPFLAGS = -Wno-inline
test_rational_SOURCES = rationaltest.cpp exprtest.cpp test_rational.cpp

if WITH_CLN
test_rational_CXXFLAGS += $(CLN_CFLAGS)
test_rational_SOURCES += clntest.cpp exprtest_cln.cpp
endif

if WITH_UNLOGICAL
test_rational_CXXFLAGS += --param inline-unit-growth=180
else
if WITH_GMPXX
test_rational_CXXFLAGS += --param inline-unit-growth=60
endif
endif

if WITH_GMPXX
test_rational_CXXFLAGS += $(GMP_CFLAGS)
test_rational_SOURCES += gmptest.cpp
endif

if WITH_INFINT
test_rational_CXXFLAGS += $(INFINT_CFLAGS)
test_rational_SOURCES += infinttest.cpp exprtest_infint.cpp
endif

if WITH_UNLOGICAL
test_rational_CXXFLAGS += $(UNLOGICAL_CFLAGS)
test_rational_SOURCES += unlogicaltest.cpp
endif

test_rational_LDFLAGS = -no-install
test_rational_LDADD = $(CPPUNIT_LIBS)

if WITH_CLN
test_rational_LDADD += $(CLN_LIBS)
endif

if WITH_GMPXX
test_rational_LDADD += $(GMP_LIBS)
endif

if WITH_UNLOGICAL
test_rational_LDADD += $(UNLOGICAL_LIBS)
endif

AUTOMAKE_OPTIONS = serial-tests
