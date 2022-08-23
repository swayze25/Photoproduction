# PYTHIA configuration file.
# Generated on Tuesday 16 August 2022 11:30:46 AM IST with the user supplied options:
# --with-fastjet3=/home/siddharth/HEP/Fastjet/fastjet-install
# --with-root=/home/siddharth/ROOT/root-6.24.06-install

# Install directory prefixes.
PREFIX_BIN=/home/siddharth/Pythia8/pythia8306/bin
PREFIX_INCLUDE=/home/siddharth/Pythia8/pythia8306/include
PREFIX_LIB=/home/siddharth/Pythia8/pythia8306/lib
PREFIX_SHARE=/home/siddharth/Pythia8/pythia8306/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=true
FASTJET3_CONFIG=fastjet-config
FASTJET3_BIN=/home/siddharth/HEP/Fastjet/fastjet-install/bin/
FASTJET3_INCLUDE=-I/home/siddharth/HEP/Fastjet/fastjet-install/include
FASTJET3_LIB=-L/home/siddharth/HEP/Fastjet/fastjet-install/lib -Wl,-rpath,/home/siddharth/HEP/Fastjet/fastjet-install/lib -lfastjet

HEPMC2_USE=false
HEPMC2_CONFIG=
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

HEPMC3_USE=false
HEPMC3_CONFIG=
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=false
LHAPDF6_CONFIG=
LHAPDF6_BIN=
LHAPDF6_INCLUDE=
LHAPDF6_LIB=

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

ROOT_USE=true
ROOT_CONFIG=root-config
ROOT_BIN=/home/siddharth/ROOT/root-6.24.06-install/bin/
ROOT_INCLUDE=-I/home/siddharth/ROOT/root-6.24.06-install/include
ROOT_LIB=-L/home/siddharth/ROOT/root-6.24.06-install/lib -Wl,-rpath,/home/siddharth/ROOT/root-6.24.06-install/lib -lCore

GZIP_USE=false
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=false
PYTHON_CONFIG=
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=

MPICH_USE=false
MPICH_CONFIG=
MPICH_BIN=
MPICH_INCLUDE=
MPICH_LIB=

HDF5_USE=false
HDF5_CONFIG=
HDF5_BIN=
HDF5_INCLUDE=
HDF5_LIB=

HIGHFIVE_USE=false
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=
HIGHFIVE_INCLUDE=
HIGHFIVE_LIB=
