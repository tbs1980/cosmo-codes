# Write `healpix` maps with minimal header using frotran

An example showing how to write a healpix map with minimal header using fortran.
The example will write the following header to a file called `test_map.fts`.


    XTENSION= 'BINTABLE'           / binary table extension
    BITPIX  =                    8 / 8-bit bytes
    NAXIS   =                    2 / 2-dimensional binary table
    NAXIS1  =                 8192 / width of table in bytes
    NAXIS2  =                  768 / number of rows in table
    PCOUNT  =                    0 / size of special data area
    GCOUNT  =                    1 / one data group (required keyword)
    TFIELDS =                    1 / number of fields in each row
    COMMENT
    COMMENT  -----------------------------------------------
    COMMENT  Sky Map Pixelisation Specific Keywords
    COMMENT  -----------------------------------------------
    PIXTYPE = 'HEALPIX '           / HEALPIX Pixelisation
    ORDERING= 'RING    '           / Pixel ordering scheme, either RING or NESTED
    NSIDE   =                  256 / Resolution parameter for HEALPIX
    FIRSTPIX=                    0 / First pixel # (0 based)
    LASTPIX =               786431 / Last pixel # (0 based)
    COORDSYS= 'unknown '           / Pixelisation coordinate system
    COMMENT  G = Galactic, E = ecliptic, C = celestial = equatorial
    BAD_DATA=  -1.637500000000E+30 / Sentinel value given to bad pixels
    COMMENT
    COMMENT  -----------------------------------------------
    COMMENT  Planck Simulation Specific Keywords
    COMMENT  -----------------------------------------------
    EXTNAME = 'FULL SKY MAP'
    COMMENT
    POLCCONV= 'COSMO   '           / Coord. convention for polarisation (COSMO/IAU)
    COMMENT  -----------------------------------------------
    COMMENT  Data Description Specific Keywords
    COMMENT  -----------------------------------------------
    COMMENT
    COMMENT  Full sky data
    OBJECT  = 'FULLSKY '
    INDXSCHM= 'IMPLICIT'           / Indexing : IMPLICIT or EXPLICIT
    GRAIN   =                    0 / Grain of pixel indexing
    COMMENT  GRAIN=0 : no indexing of pixel data                         (IMPLICIT)
    COMMENT  GRAIN=1 : 1 pixel index -> 1 pixel data                     (EXPLICIT)
    COMMENT  GRAIN>1 : 1 pixel index -> data of GRAIN consecutive pixels (EXPLICIT)
    COMMENT
    POLAR   =                    F / Polarisation included (True/False)
    DERIV   =                    0 / Derivative included (0, 1 or 2)
    COMMENT
    TTYPE1  = 'TEMPERATURE'        / Temperature map
    TFORM1  = '1024D   '           / data format of field: 8-byte REAL
    TUNIT1  = 'unknown '           / map unit
    COMMENT
    MYCARD  = 'Test a card'
    END
