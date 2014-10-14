program write_healpix_maps
    use healpix_types
    use pix_tools, only : nside2npix, npix2nside
    use head_fits, only : add_card, write_minimal_header
    use fitstools, only : output_map

    integer(kind=i4b) nside
    integer(kind=i8b) npix
    integer(kind=i8b) i
    real(kind=DP), dimension(:,:), allocatable :: map
    character(len=80), dimension(1:120) :: header
    integer(kind=i4b) nmap
    integer(i4b) :: nlheader
    character(len=filenamelen) :: outfile

    nside = 256
    npix = nside2npix(nside)
    nmap = 1


    allocate(map(0:npix-1,1:nmap))

    do i=0,npix-1
        map(i,1) = real(i)
    enddo

    ! create a minimal header
    call write_minimal_header(header, 'MAP', nside=nside, ordering='RING')
    ! write something of my own information
    call add_card(header, 'MYCARD', 'Test a card')

    nlheader = size(header)
    outfile = "!test_map.fits"

    ! write ouput
    call output_map(map,header,outfile)

    deallocate(map)
end program write_healpix_maps
