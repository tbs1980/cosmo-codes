CXX=g++

HEALPIX_INC_DIR=/share/splinter/cosmos/modules/nov_2014/install_dir/Healpix_3.20_cfitsio3370/include
HEALPIX_LIB_DIR=/share/splinter/cosmos/modules/nov_2014/install_dir/Healpix_3.20_cfitsio3370/lib
CFITSIO_LIB_DIR=/share/splinter/cosmos/modules/nov_2014/install_dir/cfitsio3370/lib

LIBS=-lhealpix_cxx -lcxxsupport -lsharp -lc_utils -lfftpack -lcfitsio -lstdc++

CFLAGS=-I$(HEALPIX_INC_DIR) -std=c++98 -Wall -Wextra -fopenmp
LDFLAGS=-L$(HEALPIX_LIB_DIR) -L$(CFITSIO_LIB_DIR)

default:read_hp_maps.exe

%.exe:%.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
%.o:%.cpp
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.exe
