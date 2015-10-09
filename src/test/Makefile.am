check_PROGRAMS = test_rational

TESTS = $(check_PROGRAMS)

noinst_HEADERS = rationaltest.h gmptest.h exprtest.h

test_rational_CXXFLAGS = -std=c++98 -fexceptions -I$(top_srcdir)/src/librational $(CPPUNIT_CFLAGS)
test_rational_SOURCES = rationaltest.cpp exprtest.cpp test_rational.cpp

if WITH_GMPXX
test_rational_CXXFLAGS += --param inline-unit-growth=60
test_rational_SOURCES += gmptest.cpp
endif

test_rational_LDFLAGS = -no-install
test_rational_LDADD = $(CPPUNIT_LIBS)

if WITH_GMPXX
test_rational_LDADD += -lgmpxx -lgmp
endif
