#include <string.h>
#include <stdio.h>
#include <fitsio.h>

static void CheckErrors(int status)
{
    if(status)
    {
        char err_text[32];
        fits_get_errstatus(status, err_text);
        printf("Error: %s\n",err_text);
    }
}

int main(void)
{
    fitsfile *fptr;
    int status = 0;
    char* file_name = "!mybintab.fits";

    //create a fits file
    CheckErrors( fits_create_file(&fptr,file_name, &status) );

    // specify what is in each column, see the link for more details
    //http://heasarc.gsfc.nasa.gov/docs/software/fitsio/quick/node10.html
    const long nrows=1025;
    const int tfields=2;
    char *ttype[tfields], *tform[tfields], *tunit[tfields];

    for (int ii = 0; ii < tfields; ii++)
    {
        ttype[ii] = (char *) malloc(20);
        tform[ii] = (char *) malloc(20);
        tunit[ii] = (char *) malloc(20);
    }

    // name for the column
    strcpy(ttype[0], "l");
    strcpy(ttype[1], "C_l");

    // data format
    strcpy(tform[0], "I11");
    strcpy(tform[1], "E13.5");

    // units
    strcpy(tunit[0], "");
    strcpy(tunit[1], "mK^2");

    // name of the HDU extension
    char *extname ="POWER SPECTRUM";

    //insert an empy table first
    CheckErrors(fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype,tform, tunit,extname,&status) );

    //note that everything here starts with 1 rather than 0 (unlike the C arryas)
    long firstrow =1;
    long firstelem = 1;
    long nelements = nrows;
    int* ell = malloc(nrows*sizeof(int));
    double* c_ell=malloc(nrows*sizeof(double));

    // create data
    for(long i=0;i<nrows;++i)
    {
        ell[i] = (int) i;
        c_ell[i] = (double) i*100;
    }

    //first column datatype TINT
    int colnum = 1;
    CheckErrors( fits_write_col(fptr, TINT, colnum, firstrow, firstelem, nelements, ell,&status) );

    //second column datatype is TDOUBLE
    colnum = 2 ;
    CheckErrors( fits_write_col(fptr, TDOUBLE, colnum, firstrow, firstelem, nelements, c_ell,&status) );

    //free memory
    for (int ii = 0; ii < tfields; ii++)
    {
        free(ttype[ii]);
        free(tform[ii]);
        free(tunit[ii]);
    }
    free(ell);
    free(c_ell);

    //close fits file
    CheckErrors( fits_close_file(fptr, &status) );
    return 0;
}
