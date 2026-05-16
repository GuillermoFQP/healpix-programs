! Intrinsic Mode Internal Linear Combination for full-sky HEALPix maps

program mapstacker

use healpix_modules
use hp_globals
use hp_utils

implicit none

!======================================================================================
real(DP), allocatable :: map_in(:,:), freqmaps(:,:)
real(DP)              :: fwhm, bl_min, fwhm_rad
integer               :: nside, nside_in, npix_in, nmaps, ord_in, n
integer               :: i, lcut
integer               :: nfreqmaps
character(len=80)     :: fout
character(len=16)     :: key
character(len=32)     :: label
!======================================================================================

! Read number of frequency maps
write(*,'(/,X,A)',advance='no') "Enter number of frequency maps: "
read(*,*) nfreqmaps

! Read common resolution parameter for the frequency maps
write(*,'(/,X,A)',advance='no') "Enter resulution parameter (Nside): "
read(*,*) nside

n = 12 * nside**2 - 1

! Provide storage for frequency maps array
allocate(freqmaps(0:n,nfreqmaps))

do i = 1, nfreqmaps
	! Read input file name of frequency map
	write(*,'(/,X,"File name of map ", I0,": ")',advance='no') i
	read(*,*) fin
	
	! Store frequency map in array "map_in"
	npix_in = getsize_fits(fin, nmaps=nmaps, nside=nside_in, ordering=ord_in)
	allocate(map_in(0:12*nside_in**2-1,1))
	call input_map(fin, map_in, npix_in, 1)
	
	! Change ordering to RING
	if (ord_in == NEST) call convert_nest2ring(nside_in, map_in(:,1))

	! Change frequency maps resolution to the given resolution parameter
	if (nside_in == nside) freqmaps(:,i) = map_in(:,1)
	if (nside_in /= nside) call udgrade_ring(map_in(:,1), nside_in, freqmaps(:,i), nside)

	deallocate(map_in)
end do

! Read smoothing FWHM parameter
write(*,'(/,X,A)',advance='no') "Enter smoothing FWHM parameter (arcmin): "
read(*,*) fwhm

! Smoothing per channel
do i = 1, nfreqmaps
	call smoothing(nside, RING, 3*nside-1, freqmaps(:,i), fwhm)
end do

! Output IMFs
call write_minimal_header(header, "map", nside=nside, order=RING, creator="MAPSTACKER", version="1.0", fwhm_degree=fwhm)

! Generate and assign column names (TTYPEi) to each map
do i = 1, nfreqmaps
	! Create FITS keyword name: TTYPE1, TTYPE2, ...
	write(key, '(A,I0)') 'TTYPE', i
	
	! Create a label for each map (customize as you like)
	write(label, '("CHANNEL_", I0)') i
	
	! Add to header
	call add_card(header, trim(key), trim(label), 'Frequency map label')
end do

! Read smoothing FWHM parameter
write(*,'(/,X,A)',advance='no') "Enter output file name: "
read(*,*) fout

call output_map(freqmaps, header, fout)
!write_map(fout, M, nside, ord, pol, vec, creator, version)

write (*,*) "Output file generated successfully."

deallocate(freqmaps)

end program
