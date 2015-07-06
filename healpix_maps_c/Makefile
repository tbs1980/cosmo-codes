CC=gcc

HEALPIX_INC_DIR=/share/splinter/cosmos/modules/nov_2014/install_dir/Healpix_3.20_cfitsio3370/include
HEALPIX_LIB_DIR=/share/splinter/cosmos/modules/nov_2014/install_dir/Healpix_3.20_cfitsio3370/lib
CFITSIO_LIB_DIR=/share/splinter/cosmos/modules/nov_2014/install_dir/cfitsio3370/lib

LIBS=-lchealpix -lcfitsio

CFLAGS=-I$(HEALPIX_INC_DIR) -std=c99 -Wall -Wextra -fopenmp
LDFLAGS=-L$(HEALPIX_LIB_DIR) -L$(CFITSIO_LIB_DIR)

default:read_hp_maps.exe

%.exe:%.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
%.o:%.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.exe
