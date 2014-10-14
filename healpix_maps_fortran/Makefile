FC=gfortran

HEALPIX_INC_DIR=/arxiv/libraries/source/Healpix_3.11/include
HEALPIX_LIB_DIR=/arxiv/libraries/source/Healpix_3.11/lib
CFITSIO_LIB_DIR=/arxiv/libraries/ubuntu/gcc/cfitsio3350/lib

LIBS=-lgif -lhealpix -lhpxgif -lsharp_healpix_f -lcfitsio -lgomp

FFLAGS=-I$(HEALPIX_INC_DIR) -Wall -Wextra -fopenmp
LDFLAGS=-L$(HEALPIX_LIB_DIR) -L$(CFITSIO_LIB_DIR)

default:write_healpix_maps.exe

%.exe:%.o
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)
%.o:%.f95
	$(FC) $(FFLAGS) -c $<

clean:
	rm *.exe
