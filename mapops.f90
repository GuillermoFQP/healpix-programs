! Intrinsic Mode Internal Linear Combination for full-sky HEALPix maps
program mapops

use healpix_modules
use hp_globals
use hp_utils

implicit none

!======================================================================================
real(DP), allocatable :: map_in(:,:), map_out(:,:)
integer               :: nside, nside_in, npix_in, nmaps, ord_in, n
integer               :: i
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
allocate(map_out(0:n,1), source=0.0)

do i = 1, nfreqmaps
	! Read input file name of frequency map
	write(*,'(/,X,"File name of map ", I0,": ")',advance='no') i
	read(*,*) fin
	
	! Store frequency map in array "map_in"
	npix_in = getsize_fits(fin, nmaps=nmaps, nside=nside_in, ordering=ord_in)
	allocate(map_in(0:12*nside_in**2-1,1))
	call input_map(fin, map_in, npix_in, 1)
	map_out = map_out + map_in
	deallocate(map_in)
end do

! Output IMFs
call write_minimal_header(header, "map", nside=nside, order=ord_in, creator="MAPSTACKER", version="1.0")

! Read smoothing FWHM parameter
write(*,'(/,X,A)',advance='no') "Enter output file name: "
read(*,*) fout

call output_map(map_out, header, fout)
!write_map(fout, M, nside, ord, pol, vec, creator, version)

write (*,*) "Output file generated successfully."

deallocate(map_out)

end program
