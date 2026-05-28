! Intrinsic Mode Internal Linear Combination for full-sky HEALPix maps
program imilc

use healpix_modules
use hp_globals
use hp_utils

implicit none

!======================================================================================
integer                        :: n, nfreqmaps, nside, npix, ord_in, nmaps
integer                        :: i, j, p, q, l, nneigh
integer, allocatable           :: neigh(:)
real(DP), allocatable          :: map_in(:,:), freqmaps(:,:), map_out(:,:)
real(DP), allocatable          :: cov(:,:), cov_inv(:,:), w(:)
real(DP), allocatable          :: local_data(:,:), weights_kernel(:)
real(DP)                       :: radius, radius_rad, sigma, dist, vp(3), vq(3)
character(len=80)              :: filename
!======================================================================================

! Read number of frequency maps
write(*,'(/,X,A)',advance='no') "Enter number of frequency maps: "
read(*,*) nfreqmaps

! Read resolution parameter
write(*,'(/,X,A)',advance='no') "Enter resolution parameter (Nside): "
read(*,*) nside

n = 12*nside**2-1

! Read Gaussian kernel radius (in radians)
write(*,'(/,X,A)',advance='no') "Enter Gaussian kernel radius (arcmin): "
!write(*,'(/,X,A)',advance='no') "Enter top-hat kernel radius (arcmin): "
read(*,*) radius

radius_rad = (radius/60.0) * (pi/180.0) ! FWHM parameter in radians
sigma      = radius_rad / 2.0           ! Standard deviation for Gaussian kernel

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

! Allocate output
allocate(map_out(0:n,1), source=0.0)

! Allocate covariance structures
allocate(cov(nfreqmaps,nfreqmaps), cov_inv(nfreqmaps,nfreqmaps), w(nfreqmaps), source=0.0)

! Loop over pixels
!$OMP PARALLEL PRIVATE(p, vp, l, neigh, nneigh, local_data, weights_kernel, cov, cov_inv, w, i, j, q, vq, dist) SHARED(n, nside, radius_rad, sigma, freqmaps, map_out, nfreqmaps)
!$OMP DO SCHEDULE(DYNAMIC)
do p = 0, n
	! Find vector pointing to the center of pixel "p"
	call pix2vec_nest(nside, p, vp)

	! Estimate the size of the array containing all the neighbour pixels
	l = 12 * nside**2 * sin(radius_rad/2.0)**2
	allocate(neigh(0:3*l/2))
	
	! Find neighbors within radius
	call query_disc(nside, vp, radius_rad, neigh, nneigh, nest=1)
	
	! Allocate local data
	allocate(local_data(0:nneigh-1,nfreqmaps), weights_kernel(0:nneigh-1))
	
	! Build local weighted dataset
	do i = 0, nneigh-1
		! Select a pixel in neighborhood
		q = neigh(i)
		
		! Find vector pointing to the center of pixel "q"
		call pix2vec_nest(nside, q, vq)
		
		! Calculate angular distance between it and the central pixel
		call angdist(vp, vq, dist)

		! Calculate the weight for that pixel based on its distance from the center
		weights_kernel(i) = exp( - (dist**2) / (2.0 * sigma**2) )
!		weights_kernel(i) = 1.0
		
		do j = 1, nfreqmaps
			local_data(i,j) = freqmaps(q,j)
		end do
	end do
	
	! Compute weighted covariance matrix
	call weighted_cov_mtx(nfreqmaps, nneigh, local_data, weights_kernel, cov)

	! Invert covariance
	call mtx_inv(nfreqmaps, cov, cov_inv)

	! Compute weights
	w = sum(cov_inv, dim=2) / sum(cov_inv)

	! Apply weights at central pixel only
	do j = 1, nfreqmaps
		map_out(p,1) = map_out(p,1) + w(j) * freqmaps(p,j)
	end do

	deallocate(local_data, weights_kernel, neigh)
end do
!$OMP END DO
!$OMP END PARALLEL

! Change ordering to ring if necessary
!do i = 1, nfreqmaps
!	call convert_nest2ring(nside, map_out(:,1))
!end do

! Output result
write(*,'(/,X,A)',advance='no') "Enter output file name: "
read(*,*) filename

call write_minimal_header(header, 'map', nside=nside, order=NEST, creator="MEMD", version="1.0")
call output_map(map_out, header, filename)

! Free memory
deallocate(map_out, freqmaps, map_in, cov, cov_inv, w)

end program
