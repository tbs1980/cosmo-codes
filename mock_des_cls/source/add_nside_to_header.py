import pyfits as pf

#open base

types = ["counts","noise"]
models = ["fix","reg","std"]


for typ in types:
    filename = "filemap_base_"+typ+".fits"
    print "opening ",filename
    #open file
    hdulist = pf.open(filename,mode='update')
    #add stuff
    hdulist[1].header['PIXTYPE'] = ('HEALPIX','HEALPIX pixelisation')
    hdulist[1].header['ORDERING'] = ('RING','Pixel ordering scheme, either RING or NESTED')
    hdulist[1].header['EXTNAME'] = ('xtension','name of this binary table extension')
    hdulist[1].header['NSIDE'] = (128,'Resolution parameter of HEALPIX')
    hdulist[1].header['FIRSTPIX'] = (0,'First pixel # (0 based)')
    hdulist[1].header['LASTPIX'] = (196607,'Last pixel # (0 based)')
    hdulist[1].header['INDXSCHM'] = ('IMPLICIT','Indexing: IMPLICIT or EXPLICIT')
    hdulist.flush()
    hdulist.close()

for model in models:
    for typ in types:
        for ind in range(101,126):
            filename = "filemap_"+model+"_"+typ+"_"+str(ind)+".fits"
            print "opening ",filename
            hdulist = pf.open(filename,mode='update')
            #add stuff
            hdulist[1].header['PIXTYPE'] = ('HEALPIX','HEALPIX pixelisation')
            hdulist[1].header['ORDERING'] = ('RING','Pixel ordering scheme, either RING or NESTED')
            hdulist[1].header['EXTNAME'] = ('xtension','name of this binary table extension')
            hdulist[1].header['NSIDE'] = (128,'Resolution parameter of HEALPIX')
            hdulist[1].header['FIRSTPIX'] = (0,'First pixel # (0 based)')
            hdulist[1].header['LASTPIX'] = (196607,'Last pixel # (0 based)')
            hdulist[1].header['INDXSCHM'] = ('IMPLICIT','Indexing: IMPLICIT or EXPLICIT')
            hdulist.flush()
            hdulist.close()
