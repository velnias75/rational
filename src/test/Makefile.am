check_PROGRAMS = test_rational

TESTS = $(check_PROGRAMS)

noinst_HEADERS = rationaltest.h

test_rational_CXXFLAGS = -std=c++98 -I$(top_srcdir)/src/librational $(CPPUNIT_CFLAGS)
test_rational_SOURCES = rationaltest.cpp test_rational.cpp
test_rational_LDFLAGS = -no-install
test_rational_LDADD = $(CPPUNIT_LIBS)
