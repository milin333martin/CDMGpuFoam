# Morton-key adaptive mesh engine
#
# One source, four back ends, exactly as src/gpu does it:
#
#   make            serial, any C++20 compiler
#   make parallel   host threads via g++ and TBB
#   make multicore  host threads via nvc++
#   make gpu        the device, via nvc++ -stdpar=gpu
#
# Only "make gpu" uses the GPU.  The others are for validating that the
# same code gives the same answer wherever it runs, which is how the flux
# back end's bugs were caught.

SRC      := cdmAmr.cpp
TEST     := test/testAmr.cpp
BIN      := test/testAmr

GPUFLAGS ?= -gpu=ccall

.PHONY: all parallel multicore gpu foam test carry clean

all:
	g++ -std=c++20 -O3 -march=native -Wall -Wextra -I. $(TEST) $(SRC) -o $(BIN)
	@echo "built serial"

parallel:
	g++ -std=c++20 -O3 -march=native -Wall -Wextra -DCDM_STDPAR -I. \
	    $(TEST) $(SRC) -o $(BIN) -ltbb
	@echo "built with host threads (TBB)"

multicore:
	nvc++ -std=c++20 -O3 -stdpar=multicore -DCDM_STDPAR -I. \
	    $(TEST) $(SRC) -o $(BIN)
	@echo "built with host threads (nvc++)"

gpu:
	nvc++ -std=c++20 -O3 -stdpar=gpu $(GPUFLAGS) -DCDM_STDPAR -I. \
	    $(TEST) $(SRC) -o $(BIN)
	@echo "built for the device"

test: $(BIN)
	./$(BIN)

# Refinement and coarsening in the same step.  Separate from testAmr
# because it needs no back end: it is a statement about indices, and it
# is the same on the host and on the device.
carry:
	g++ -std=c++20 -O3 -Wall -Wextra -I. test/testCarry.cpp $(SRC) -o test/testCarry
	./test/testCarry

clean:
	rm -f $(BIN) test/testCarry

# A shared library the OpenFOAM side can link against.  Built with the same
# compiler OpenFOAM was built with, deliberately: the engine's interface is
# C++ classes and std::vector, so it must not cross a compiler boundary the
# way src/gpu's plain-C interface safely can.  The device build of the
# engine is exercised by "make gpu" and testAmr, not from inside OpenFOAM.
# The standard must match the one OpenFOAM's own wmake uses, because
# cdmAmr::Mesh is defined in a header that both this library and
# libCDMDynamicMesh2D include.  The OpenFOAM-facing target deliberately
# keeps OpenFOAM's C++ standard (default c++14).  Concurrent flag writes
# use OpenMP atomics in this target; std::atomic_ref is compiled only
# for the standalone C++20 stdpar targets.
FOAM_STD ?= c++14

foam:
	@test -n "$(FOAM_USER_LIBBIN)" || { echo "source the OpenFOAM environment first"; exit 1; }
	g++ -std=$(FOAM_STD) -O3 -march=native -fPIC -Wall -Wextra -fopenmp -DCDM_OPENMP -shared $(SRC) \
	    -o $(FOAM_USER_LIBBIN)/libCDMAmr.so
	@echo "built $(FOAM_USER_LIBBIN)/libCDMAmr.so"
