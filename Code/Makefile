# Makefile for programs for Bennett-Gherga-Rechnitzer algorithm

# eclib is a requirement.  If installed in the usual place /usr/local
# the following line might not be necessary.  If installed anywhere
# else set the install directory here:

ECLIB_BASE=/usr/local
#ECLIB_BASE=$(HOME)/eclib
INCDIR = $(ECLIB_BASE)/include
LIBDIR = $(ECLIB_BASE)/lib

GCC=g++ -std=c++11 -fmax-errors=1
CC = $(GCC)


# The type of integers used for components of Quad, Qideal, RatQuad
# can be either long or bigint (=NTL's ZZ), and is typedef'd to QUINT
# in the code.  By default the type is long and QUINT_IS_long is
# defined; this can result in overflow for large levels over larger
# fields.  Change this to 1 (and make clean and rebuild) to compile
# using ZZ as base integer type for Quads instead.  That results in
# slower code, but it does not overflow!
BASE_TYPE_ZZ=0
ifeq ($(BASE_TYPE_ZZ), 1)
 BASE_TYPE_FLAG = -D QUINT_IS_ZZ
endif

# NB If used with a multithreaded build of eclib then you MUST define
# USE_BOOST=1 below so that the correct compiler and linker stuff is
# appended below.  Otherwise set USE_BOOST=0 (or do not set it).

USE_BOOST=1
ifeq ($(USE_BOOST), 1)
 BOOST_ASIO_LIB = -lboost_system-mt
 BOOST_CPPFLAGS =   -DECLIB_MULTITHREAD -DHAVE_STDCXX_0X=/\*\*/ -DHAVE_TR1_UNORDERED_MAP=/\*\*/ -DHAVE_STDCXX_0X=/\*\*/ -DHAVE_UNORDERED_MAP=/\*\*/# -pthread -I/usr/include
 BOOST_SYSTEM_LIB = -lboost_system
 BOOST_THREAD_LIB = -lboost_thread
 BOOST_LDFLAGS = -L/usr/lib -pthread $(BOOST_SYSTEM_LIB) $(BOOST_THREAD_LIB)
endif

# to disable checking of assert() use the following:
#OPTFLAG = -DNDEBUG -O3 -fPIC
# to enable checking of assert() use the following:
OPTFLAG = -O3 -Wall -fPIC

# for profiling:
#CFLAGS = -c -g -pg $(OPTFLAG) $(BOOST_CPPFLAGS) $(BASE_TYPE_FLAG) -I$(INCDIR)
#LFLAGS = -pg -lec -lntl -lstdc++  -L$(LIBDIR) -Wl,-rpath -Wl,$(LIBDIR) $(BOOST_LDFLAGS)

#for coverage:
#CFLAGS = -c -g --coverage $(BOOST_CPPFLAGS) $(BASE_TYPE_FLAG) -I$(INCDIR)
#LFLAGS = --coverage -fprofile-arcs -lec -lntl -lstdc++  -L$(LIBDIR) -Wl,-rpath -Wl,$(LIBDIR) $(BOOST_LDFLAGS)

#for normal use:
CFLAGS = -c -g $(OPTFLAG) $(BOOST_CPPFLAGS) $(BASE_TYPE_FLAG) -I$(INCDIR)
LFLAGS = -lec -lntl -lstdc++  -L$(LIBDIR) -Wl,-rpath -Wl,$(LIBDIR) $(BOOST_LDFLAGS)

all: tests

sources: ccs headers
	chmod a+r *.h *.cc

ccs: cubic_utils.cc TME.cc cubic_forms_support.cc cubic_forms_support.cc

headers: cubic_utils.h TME.h

%.o:   %.cc
	$(CC) $(CFLAGS) $<

TESTS = cubic_forms_conductor cubic_forms_support
tests: $(TESTS)

clean:
	rm -f $(TESTS)
	rm -f *.o *~

cubic_forms_conductor: cubic_forms_conductor.o cubic_utils.o TME.o
	$(CC) -o cubic_forms_conductor cubic_forms_conductor.o cubic_utils.o TME.o $(LFLAGS)

cubic_forms_support: cubic_forms_support.o cubic_utils.o
	$(CC) -o cubic_forms_support cubic_forms_support.o cubic_utils.o $(LFLAGS)

cubic_forms_conductor.o: cubic_forms_conductor.cc cubic_utils.h TME.h
cubic_forms_support.o: cubic_forms_support.cc cubic_utils.h
cubic_utils.o: cubic_utils.h cubic_utils.cc
TME.o: cubic_utils.h TME.h

