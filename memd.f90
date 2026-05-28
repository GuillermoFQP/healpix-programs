! Multivariate Empirical Mode Decomposition for HEALPix maps
program multivariate_empirical_mode_decomposition

use healpix_modules
use hp_globals
use hp_utils

implicit none

!======================================================================================
real(DP), allocatable          :: map_in(:,:), map_out(:,:,:), freqmaps(:,:), mask_in(:,:), mask(:,:)
real(DP)                       :: fwhm, bl_min, fwhm_rad, strip_bound
integer                        :: nside, nside_in, npix, npix_in, nmaps, ord, ord_in, n
integer                        :: nside_mask, npix_mask, nmaps_mask, ord_mask
integer                        :: nimf, i, stff, tens, lmax, lcut
integer                        :: nfreqmaps, ndir
character(len=80)              :: fmask
character(len=80), allocatable :: fout(:)
!======================================================================================

! Read number of frequency maps
write(*,'(/,X,A)',advance='no') "Enter number of frequency maps: "
read(*,*) nfreqmaps

! Read common resolution parameter for the frequency maps
write(*,'(/,X,A)',advance='no') "Enter resulution parameter (Nside): "
read(*,*) nside

! Provide storage for frequency maps array
allocate(freqmaps(0:12*nside**2-1,nfreqmaps), mask(0:12*nside**2-1,1), source=0.0)

! Store mask map in array "mask"
!fmask = "COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits"
fmask = "COM_Mask_CMB-Inpainting-Mask-Int_2048_R3.00.fits"
npix_mask = getsize_fits(fmask, nmaps=nmaps_mask, nside=nside_mask, ordering=ord_mask)
allocate(mask_in(0:12*nside_mask**2-1,1))
call input_map(fmask, mask_in, npix_mask, 1)
if (ord_mask == RING) call convert_ring2nest(nside_mask, mask_in(:,1))
if (nside_mask == nside) mask = mask_in
if (nside_mask /= nside) call udgrade_nest(mask_in(:,1), nside_mask, mask(:,1), nside)

do i = 1, nfreqmaps
	! Read input file name of frequency map
	write(*,'(/,X,"File name of map ", I0,": ")',advance='no') i
	read(*,*) fin
	
	! Store frequency map in array "map_in"
	npix_in = getsize_fits(fin, nmaps=nmaps, nside=nside_in, ordering=ord_in)
	allocate(map_in(0:12*nside_in**2-1,1))
	call input_map(fin, map_in, npix_in, 1)
	
	! Change ordering to NEST
	if (ord_in == RING) call convert_ring2nest(nside_in, map_in(:,1))
	
!	strip_bound = sin(2.5 * DEG2RAD)
!	call apply_mask(map_in, NEST, zbounds=[strip_bound, -strip_bound])

	! Change frequency maps resolution to the given resolution parameter
	if (nside_in == nside) freqmaps(:,i) = map_in(:,1)
	if (nside_in /= nside) call udgrade_nest(map_in(:,1), nside_in, freqmaps(:,i), nside)
	
	deallocate(map_in)
end do

! Mask all the frequency maos
call apply_mask(freqmaps, NEST, mask=mask)

! Read number of intrinsic modes
write(*,'(/,X,A)',advance='no') "Enter number of intrinsic modes: "
read(*,*) nimf

! Read interpolation stiffness parameter
write(*,'(/,X,A)',advance='no') "Enter interpolation stiffness parameter: "
read(*,*) stff

! Read interpolation tension parameter
write(*,'(/,X,A)',advance='no') "Enter interpolation tension parameter: "
read(*,*) tens

! Read smoothing FWHM parameter
write(*,'(/,X,A)',advance='no') "Enter smoothing FWHM parameter: "
read(*,*) fwhm

! Read smoothing FWHM parameter
write(*,'(/,X,A)',advance='no') "Enter number of projections for the computation of the MEMD mean envelopes: "
read(*,*) ndir

! Output file names
allocate(fout(nimf+1))

do i = 1, nimf+1
	if (i <= nimf  ) write (fout(i), '("imf_",I0,".fits")') i
	if (i == nimf+1) write (fout(i), '("res.fits")')
end do

n        = nside2npix(nside) - 1
fwhm_rad = (fwhm/60.0) * (pi/180.0)
bl_min   = 1.0d-7
lcut     = int(sqrt(0.25 - 16.0*log(bl_min)*log(2.0) / fwhm_rad**2) - 0.5)
lmax     = min(3*nside-1, lcut, 2*2048)

! Smoothing per channel
do i = 1, nfreqmaps
	call smoothing(nside, NEST, 3*nside-1, freqmaps(:,i), fwhm)
end do

write (*,'(/,X,"Multivariate Empirical Mode Decomposition with ", I0, " directions.")') ndir

! Allocate arrays
allocate(map_out(0:n,nfreqmaps,nimf), source=0.0)

! Perform MEMD
call memd(nside, freqmaps, nfreqmaps, nimf, stff, tens, lmax, ndir, map_out, mask=mask)

write (*,'(/,X,A)') "Multivariate EMD completed successfully."

! Convert to RING
call convert_nest2ring(nside, freqmaps)
do i = 1, nimf
	call convert_nest2ring(nside, map_out(:,:,i))
end do

! Output IMFs
call write_minimal_header(header, "map", nside=nside, order=RING, creator="MEMD", version="1.0", fwhm_degree=fwhm)

do i = 1, nimf+1
	if (i <= nimf  ) call output_map(map_out(:,:,i), header, fout(i))
	if (i == nimf+1) call output_map(freqmaps - sum(map_out(:,:,1:nimf), dim=3), header, fout(i))
end do

write (*,*) "Output files generated successfully."

deallocate(freqmaps, fout, map_out)

end program
