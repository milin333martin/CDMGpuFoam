# ---------------------------------------------------------------------------
# libCDMGpu - the GPU back end for CDMGpuFoam
#
#   make              serial, system compiler.  Builds and tests anywhere.
#   make parallel     every CPU core, g++ + Intel TBB (libtbb-dev)
#   make multicore    every CPU core, nvc++ -stdpar=multicore
#   make gpu          one GPU,        nvc++ -stdpar=gpu
#   make test         run the ZND validation on whichever build is present
#
# This library is built separately from the rest of CDMGpuFoam and on
# purpose: it is the one piece that may be compiled by nvc++ rather than the
# compiler OpenFOAM was built with, and its interface is plain C so that the
# two can never disagree about an ABI.  Nothing here includes an OpenFOAM
# header.
#
# nvc++ comes with the NVIDIA HPC SDK, which is a separate free download and
# is NOT part of the CUDA toolkit:  https://developer.nvidia.com/hpc-sdk
# Without it, 'make parallel' still gives you every core, and plain 'make'
# still builds and validates.
# ---------------------------------------------------------------------------

CXX    ?= g++
NVCXX  ?= nvc++
STD     = -std=c++17
WARN    = -Wall -Wextra
LIB     = $(FOAM_USER_LIBBIN)/libCDMGpu.so
SRC     = cdmGpuKernels.cpp

# Device code generation.  ccall emits every compute capability the SDK knows
# plus PTX, so one .so runs on a workstation card and on whatever the cluster
# has - at the cost of a slow compile.  Override for a faster local build, or
# to pin the CUDA version when the driver is older than the SDK's default:
#
#   make gpu GPUFLAGS="-gpu=cc86"              # an RTX A2000 only, quick
#   make gpu GPUFLAGS="-gpu=ccall,cuda12.9"    # portable, CUDA 12.9 runtime
#
# A driver that reports "CUDA Version: 12.x" in nvidia-smi cannot run code
# built against CUDA 13; pin cuda12.x in that case and use a cuda_multi
# install of the SDK.
GPUFLAGS ?= -gpu=ccall

.PHONY: all parallel multicore gpu test clean

all:
	$(CXX) $(STD) -O3 -march=native -fPIC $(WARN) -shared $(SRC) -o libCDMGpu.so

parallel:
	$(CXX) $(STD) -O3 -march=native -fPIC $(WARN) -DCDM_STDPAR -shared $(SRC) \
	    -o libCDMGpu.so -ltbb

multicore:
	$(NVCXX) $(STD) -O3 -fPIC -stdpar=multicore -DCDM_STDPAR -shared $(SRC) \
	    -o libCDMGpu.so

gpu:
	$(NVCXX) $(STD) -O3 -fPIC -stdpar=gpu $(GPUFLAGS) -DCDM_STDPAR -shared $(SRC) \
	    -o libCDMGpu.so

install:
	@test -n "$(FOAM_USER_LIBBIN)" || { echo "source the OpenFOAM environment first"; exit 1; }
	install -D -m 0755 libCDMGpu.so $(LIB)
	@echo "installed $(LIB)"

test:
	$(CXX) $(STD) -O3 -march=native $(WARN) -I. test/testZND.cpp $(SRC) -o test/testZND
	./test/testZND

clean:
	rm -f libCDMGpu.so test/testZND
