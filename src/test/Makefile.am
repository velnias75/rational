check_PROGRAMS = test_rational

TESTS = $(check_PROGRAMS)

noinst_HEADERS = rationaltest.h gmptest.h clntest.h exprtest.h exprtest_cln.h

test_rational_CXXFLAGS = -std=c++0x -fexceptions -I$(top_srcdir)/src/librational $(CPPUNIT_CFLAGS)
test_rational_SOURCES = rationaltest.cpp exprtest.cpp test_rational.cpp

if WITH_CLN
test_rational_CXXFLAGS += $(CLN_CFLAGS)
test_rational_SOURCES += clntest.cpp exprtest_cln.cpp
endif

if WITH_GMPXX
test_rational_CXXFLAGS += --param inline-unit-growth=60
test_rational_SOURCES += gmptest.cpp
endif

test_rational_LDFLAGS = -no-install
test_rational_LDADD = $(CPPUNIT_LIBS)

if WITH_CLN
test_rational_LDADD += $(CLN_LIBS)
endif

if WITH_GMPXX
test_rational_LDADD += -lgmpxx -lgmp
endif
