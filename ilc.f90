! Intrinsic Mode Internal Linear Combination for full-sky HEALPix maps
program ilc

use healpix_modules
use hp_globals
use hp_utils

implicit none

!======================================================================================
integer               :: n, nfreqmaps, nside, npix, ord_in, nmaps
integer               :: i, j, l
real(DP), allocatable :: map_in(:,:), freqmaps(:,:), map_out(:,:)
real(DP), allocatable :: cov(:,:), cov_inv(:,:), w(:)
character(len=80)     :: filename
!======================================================================================

! Read number of frequency maps
write(*,'(/,X,A)',advance='no') "Enter number of frequency maps: "
read(*,*) nfreqmaps

! Read resolution parameter
write(*,'(/,X,A)',advance='no') "Enter resolution parameter (Nside): "
read(*,*) nside

n = 12*nside**2-1

! Allocate full input map (all frequencies at once)
allocate(map_in(0:n,nfreqmaps), freqmaps(0:n,nfreqmaps), source=0.0)

! Read single input file containing all frequency maps
write(*,'(/,X,A)',advance='no') "Enter input file name: "
read(*,*) filename

npix = getsize_fits(filename, nmaps=nmaps, nside=nside, ordering=ord_in)

call input_map(filename, map_in, n+1, nfreqmaps)

! Change ordering to nested if necessary
if (ord_in == RING) then
	do i = 1, nfreqmaps
		call convert_ring2nest(nside, map_in(:,i))
	end do
end if

freqmaps = map_in

! Allocate covariance structures
allocate(cov(nfreqmaps,nfreqmaps), cov_inv(nfreqmaps,nfreqmaps), w(nfreqmaps), source=0.0)

! Compute weighted covariance matrix
call cov_mtx(nfreqmaps, n+1, freqmaps, cov)

! Invert covariance
call mtx_inv(nfreqmaps, cov, cov_inv)

! Compute weights
w = sum(cov_inv, dim=2) / sum(cov_inv)

allocate(map_out(0:n,1), source=0.0)

! Compute CMB component
do i = 1, nfreqmaps
	freqmaps(:,i) = freqmaps(:,i) * w(i)
end do

map_out(:,1) = sum(freqmaps, dim=2)

! Output result
write(*,'(/,X,A)',advance='no') "Enter output file name: "
read(*,*) filename

call write_minimal_header(header, 'map', nside=nside, order=NEST, creator="MEMD", version="1.0")
call output_map(map_out, header, filename)

! Free memory
deallocate(map_out, freqmaps, map_in, cov, cov_inv, w)

end program
