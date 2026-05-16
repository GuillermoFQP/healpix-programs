! Empirical Mode Decomposition for HEALPix maps

program empirical_mode_decomposition

use healpix_modules
use hp_globals
use hp_utils

implicit none

!======================================================================================
real(DP), allocatable          :: map_in(:,:), map_out(:,:,:)
real(DP)                       :: fwhm, bl_min, fwhm_rad
integer                        :: nside, npix, nmaps, ord, n, imf, ch, i, stff, tens, lmax, lcut
character(len=80), allocatable :: fout(:)
!======================================================================================

! Hilbert-Huang transform
call getArgument(1, fin) ! Input file name
call getArgument(2, arg) ! Number of IMFs
read(arg,*) imf
call getArgument(3, arg) ! Stiffness
read(arg,*) stff
call getArgument(4, arg) ! Tension parameter
read(arg,*) tens
call getArgument(5, arg) ! FWHM for Gaussian smoothing (in arcmin)
read(arg,*) fwhm

! Output file names
ch = imf
allocate(fout(2*ch))

do i = 1, ch
	write (fout(i),'("imf",I0,".fits")') i
	write (fout(ch+i),'("res",I0,".fits")') i
end do

! Parameters of the FITS file containing the input map
npix = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord)

n        = nside2npix(nside) - 1                                         ! Total number of pixels minus one
fwhm_rad = (fwhm/60.0) * (pi/180.0)                                      ! FWHM parameter in radians
bl_min   = 1.0d-7                                                        ! Cutoff value for Gaussian beam
lcut     = int(sqrt(0.25 - 16.0*log(bl_min)*log(2.0)/fwhm_rad**2) - 0.5) ! Cutoff value for "l" due to Gaussian beam
lmax     = min(3*nside - 1, lcut, 2*2048)                                ! Maximum "l" for interpolation

write (*,'(/,X, "Input map Nside = ", I0)') nside

! Allocating arrays
allocate(map_in(0:n,nmaps), map_out(0:n,nmaps,ch), source=0.0)

! Reading input map
call input_map(fin, map_in, npix, nmaps)
write (*,'(/,X,A)') "Map read successfully."

! Gaussian smoothing the input map
call smoothing(nside, ord, 3*nside-1, map_in(:,1), fwhm)

! All the following subroutines are designed for maps in NESTED ordering
if (ord == 1) call convert_ring2nest(nside, map_in)

! Find extrema or perform empirical mode decomposition
write (*,'(/,X,"Interpolation with stiffness parameter ",I0," and tension parameter ",I0,".")') stff, tens
call emd(nside, map_in(:,1), imf, stff, tens, lmax, map_out(:,1,:))
write (*,'(/, X, A)') "Empirical Mode Decomposition completed successfully."

! Go back to RING ordering if necessary
if (ord == 1) then
	call convert_nest2ring(nside, map_in)
	do i = 1, ch
		call convert_nest2ring(nside, map_out(:,:,i))
	end do
end if

! Generating output files
call write_minimal_header(header, 'map', nside=nside, order=ord)
do i = 1, ch
	call output_map(map_out(:,:,i), header, fout(i))
end do

! Generating files containing the residues and the EMD
do i = 1, ch
	call output_map(map_in - sum(map_out(:,:,1:i),dim=3), header, fout(ch+i))
end do

write (*,*) "Output files generated successfully."

deallocate(map_in, map_out, fout)

end program
