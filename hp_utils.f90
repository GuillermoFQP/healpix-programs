! Empirical Mode Decomposition for HEALPix maps

module hp_utils

use healpix_modules
use hp_globals

implicit none

!======================================================================================
! Parameters for the calculations of spherical harmonics
integer, parameter  :: LOG2LG = 100, RSMAX = 20, RSMIN = -20
real(DP), parameter :: FL_LARGE = 2.0**LOG2LG, FL_SMALL = 2.0**(-LOG2LG)
real(DP), parameter :: ALN2_INV = 1.4426950408889634073599246810 ! 1/log(2)
real(DP)            :: rescale_tab(RSMIN:RSMAX)
!======================================================================================

contains

! Gaussian smoothing process
subroutine smoothing(nside, ord, lmax, map_in, fwhm)
	integer, intent(in)       :: nside, ord
	real(DP), intent(in)      :: fwhm
	real(DP), intent(inout)   :: map_in(0:12*nside**2-1)
	complex(DPC), allocatable :: alm(:,:,:)
	integer                   :: lmax

	write (*,'(/,X,"Smoothing input map with Gaussian beam.")')

	allocate(alm(1,0:lmax,0:lmax))

	! The subroutines used here are designed for maps in RING ordering
	if (ord == NEST) call convert_nest2ring(nside, map_in)

	call map2alm(nside, lmax, lmax, map_in, alm)
	call alter_alm(nside, lmax, lmax, fwhm, alm)
	call alm2map(nside, lmax, lmax, alm, map_in)

	deallocate(alm)

	! Go back to NESTED ordering if necessary
	if (ord == NEST) call convert_ring2nest(nside, map_in)

end subroutine smoothing

! Positions and values of the local extrema of an input map
subroutine local_extrema(nside, map_in, nmax, nmin, imax, imin, mask)
	integer, intent(in)   :: nside
	real(DP), intent(in)  :: map_in(0:12*nside**2-1)
	real(DP), intent(in), optional :: mask(0:12*nside**2-1,1)
	integer, intent(out)  :: nmax, nmin, imax(12*nside**2/9), imin(12*nside**2/9)
	integer               :: n, i, j, l, nlist
	integer, allocatable  :: list(:)
	real(DP)              :: rpix, radius, vi(3)
	real(DP)              :: theta, phi, lat_deg, strip_bound
	logical               :: use_disc
	
	n = nside2npix(nside) - 1
	
	nmin = 0 ! Local minima counter
	nmax = 0 ! Local maxima counter
	
	! Use "query_disc" subroutine to find local extrema, otherwise use "neighbours_nest"
	use_disc = .false.
	
	if (use_disc) then
		!======================================================================================
		! Finding local extrema using disk with given radius
		!======================================================================================
		rpix = sqrt(4.0*pi/(n+1)) / sqrt(2.0) ! Approximate radius of one pixel
		radius = 5.0*rpix                      ! Disk radius
		l = (2.0*radius)**2 / (4.0*pi/(n+1))  ! Approximate number of pixels inside the disk
		
		allocate(list(0:l))
		
		do i = 0, n
			! Skip pixels within strip around equator
			call pix2ang_nest(nside, i, theta, phi)
!			lat_deg = 90.0 - (theta * RAD2DEG)
!			strip_bound = 2.5
!			if (lat_deg > -strip_bound .and. lat_deg < strip_bound) cycle
			
			if (present(mask) .and. mask(i, 1) == 0.0) cycle

			! Pixel indices of a neighborhood around pixel "i"
			call pix2vec_nest(nside, i, vi)
			call query_disc(nside, vi, radius, list, nlist, nest=1)
			
			! Find local extrema
			if (map_in(i) >= maxval(map_in(list(1:nlist)))) then
				nmax = nmax + 1
				imax(nmax) = i
			end if
			if (map_in(i) <= minval(map_in(list(1:nlist)))) then
				nmin = nmin + 1
				imin(nmin) = i
			end if
		end do
		!======================================================================================
	else
		!======================================================================================
		! Finding local extrema using nearest neighbours
		!======================================================================================
		allocate(list(8))
		
		do i = 0, n
			! Skip pixels within ±15° latitude of the equator
			call pix2ang_nest(nside, i, theta, phi)
!			lat_deg = 90.0 - (theta * RAD2DEG)
!			strip_bound = 2.5
!			if (lat_deg > -strip_bound .and. lat_deg < strip_bound) cycle
			
			if (present(mask) .and. mask(i, 1) == 0.0) cycle
			
			! Pixel indices of a neighborhood around pixel "i"
			call neighbours_nest(nside, i, list, nlist)
			
			! Find local extrema
			if (map_in(i) >= maxval(map_in(list(1:nlist)))) then
				nmax = nmax + 1
				imax(nmax) = i
			end if
			if (map_in(i) <= minval(map_in(list(1:nlist)))) then
				nmin = nmin + 1
				imin(nmin) = i
			end if
		end do
		!======================================================================================
	end if
	
	deallocate(list)
	
end subroutine local_extrema

! Solve the real system of "n" symmetric linear equations in "n" unknowns in the form "A*X=B"
subroutine lsolve(n, A, B)
	integer, intent(in)     :: n
	real(DP), intent(in)    :: A(n,n)
	real(DP), intent(inout) :: B(n)
	real(DP), allocatable   :: work(:)
	real(DP)                :: get_lwork(1)
	integer                 :: pivot(n), stat, lwork
	
	stat = 0 ! Status indicator ("stat/=0" indicates an error)
	
	! Returns the optimal size of the WORK array as the first entry of the GET_LWORK array
	call dsysv('U', n, 1, A, n, pivot, B, n, get_lwork, -1, stat)
	
	lwork = get_lwork(1) ! Parameter to calculate the optimal size of the WORK array
	
	allocate(work(lwork))
	
	! LAPACK subroutine to get "X=A^(-1)*B" for a symmetric matrix "A"
	call dsysv('U', n, 1, A, n, pivot, B, n, work, lwork, stat)
	
	! Stop the program if necessary
	if (stat /= 0) call fatal_error('Singular matrix in subroutine LSOLVE()')
	
	deallocate(work)
	
end subroutine lsolve

! First factor for interpolation in harmonic space
subroutine interp_alms(c, ang, lmax, mmax, alm)
	integer, intent(in)       :: lmax, mmax
	real(DP), intent(in)      :: c(:), ang(size(c),2)
	complex(DPC), intent(out) :: alm(1,0:lmax,0:mmax)
	real(DP)                  :: mfac(0:mmax), recfac(0:1,0:lmax), lam_lm(0:lmax)
	integer                   :: l, m, i
	
	alm = (0.0, 0.0)
	
	! Recursion factors used in "lambda_mm" calculation for all "m" in "0<=m<=m_max"
	call gen_mfactor(mmax, mfac)
	
	!$OMP PARALLEL PRIVATE(m, i, recfac, lam_lm, l) SHARED(mmax, lmax, c, ang, mfac, alm)
	!$OMP DO SCHEDULE(DYNAMIC)
	do m = 0, mmax
		! Generate recursion factors useful for "lambda_lm" for a given "m"
		call gen_recfactor(lmax, m, recfac)
		
		do i = 1, size(c)
			! Compute "lam_lm(theta_i)" for all "l>=m" for a given "m"
			call do_lambda_lm(lmax, m, abs(cos(ang(i,1))), sin(ang(i,1)), mfac(m), recfac, lam_lm)
			if (cos(ang(i,1)) < 0.0) forall (l=m:lmax) lam_lm(l) = (-1.0)**(l+m) * lam_lm(l)
			
			! Compute numerators for all "l>=m" for a given "m"
			alm(1,m:lmax,m) = alm(1,m:lmax,m) + cmplx(c(i) * lam_lm(m:lmax), kind=DPC) * cdexp((0.0, -1.0) * m * ang(i,2))
		end do
	end do
	!$OMP END DO
	!$OMP END PARALLEL
	
end subroutine interp_alms

! Global interpolation by spherical spline (Perrin et al, 1988)
subroutine ss_interp(nside, lmax, stff, tens, next, LUT, iext, map_out)
	integer, intent(in)       :: nside, lmax, stff, tens, next, iext(next)
	real(DP), intent(in)      :: LUT(next)
	real(DP), intent(out)     :: map_out(0:12*nside**2-1)
	real(DP)                  :: vec(next,3), ang(next,2), vi(3), fwhm
	integer                   :: i, j, n, l, lmin, omega_pix
	real(DP), allocatable     :: A(:,:), B(:), bl1(:,:), bl2(:,:), wl(:,:)
	complex(DPC), allocatable :: alm(:,:,:)
	logical                   :: harmonic_space
	
	n = nside2npix(nside) - 1
	fwhm = 20.0 ! 2.0 * acos(1.0 - 2.0/next) * 180.0 * 60.0 / pi
	
	allocate(A(next,next), B(next), source=0.0)                      ! System of linear equations
	allocate(bl1(0:lmax,1), bl2(0:lmax,1), wl(0:lmax,1), source=1.0) ! Beams
	
	! Positions of the local extrema
	do i = 1, next
		call pix2vec_nest(nside, iext(i), vec(i,:))
		call pix2ang_nest(nside, iext(i), ang(i,1), ang(i,2))
	end do
	
	! Initialize RESCALE_TAB useful for calculation of spherical harmonics
	call init_rescale_tab()
	
	! Inverse of differential operator in harmonic space
	if (stff /= 0 .or. tens == 0) then
		lmin = 1
		bl1(0,1) = 1.0
	else
		lmin = 0
		bl1(0,1) = -1.0 / dble(tens)
	end if
	
	do l = 1, lmax
		bl1(l,1) = (-1.0)**(stff+1) / (dble(l)*dble(l+1) + dble(tens)) / (dble(l)*dble(l+1))**stff
	end do
	
	! Pixel window function
	call pixel_window(wl, nside)
	
	! Gaussian beam for "stff=0"
!	if (stff == 0) call gaussbeam(fwhm, lmax, bl2)
	
	write (*,'(/, X, "--- Computing interpolation matrix of dimension ", I0, ".")') next
	
	! Computing "A"
	!$OMP PARALLEL PRIVATE(i, j) SHARED(A, vec, lmin, lmax)
	!$OMP DO SCHEDULE(DYNAMIC)
	do j = 1, next
		do i = 1, j
			A(i,j) = Aij(dot_product(vec(i,:),vec(j,:)), wl*bl1*bl2, lmin, lmax)
!			A(i,j) = G(dot_product(vec(i,:),vec(j,:)), wl*bl1*bl2, lmin, lmax)
			if (i /= j) A(j,i) = A(i,j)
		end do
	end do
	!$OMP END DO
	!$OMP END PARALLEL
	
	! Computing "B"
	B = LUT - lmin * sum(LUT)/next

	write (*,'(X, A)') "--- Solving system of equations."

	! Calculating and storing the interpolation coefficients "X=A^(-1)*B" in array "B"
	call lsolve(next, A, B)
	
	deallocate(A)
	
	harmonic_space = .true.
	
	if (harmonic_space) then
		! Interpolation in spherical harmonic space
		write (*,'(X, A)') "--- Calculating spherical harmonic coefficients."
		
		allocate(alm(1,0:lmax,0:lmax), source=(0.0, 0.0))
		
		! Calculating the sum of "C_i*Y^(*)_lm(theta_i,phi_i)" over all pixels "i" with values in LUT
		call interp_alms(B, ang, lmax, lmax, alm)
		
		! Convolution in harmonic space and generation of output map
		do l = lmin, lmax
			alm(1,l,0:l) = alm(1,l,0:l) * wl(l,1) * bl1(l,1) * bl2(l,1)
		end do
		
		! "a_00" takes a different value when "stff>0" or "tens=0"
		if (lmin == 1) alm(1,0,0) = sqrt(4.0*pi) * sum(LUT) / next
		
		write (*,'(X, A)') "--- Generating map."
		
		call alm2map(nside, lmax, lmax, alm, map_out) ! Generates map in RING ordering
		call convert_ring2nest(nside, map_out)        ! Goes back to NESTED ordering

		write (*,'(X, "--- Max. and min. interpolation error: ", E10.4, X, E10.4)') maxval(abs(LUT-map_out(iext))), minval(abs(LUT-map_out(iext)))
		
		deallocate(B, wl, bl1, bl2, alm)
	else
		! Interpolation in real space
		write (*,'(X, A)') "--- Interpolation in real space started."
		
		map_out = lmin * sum(LUT)/next
		
		!$OMP PARALLEL PRIVATE (i, vi) SHARED(n, nside, map_out, lmin, lmax)
		!$OMP DO SCHEDULE(DYNAMIC)
		do i = 0, n
			! Position of the pixel on the map
			call pix2vec_nest(nside, i, vi)
			
			! Calculation of the interpolated value
			do j = 1, next
				map_out(i) = map_out(i) + B(j) * Aij(dot_product(vi,vec(j,:)), wl*bl1*bl2, lmin, lmax)
				!map_out(i) = map_out(i) + B(j) * G(dot_product(vi,vec(j,:)), wl*bl1*bl2, lmin, lmax)
			end do
		end do
		!$OMP END DO
		!$OMP END PARALLEL
		
		deallocate(B, wl, bl1, bl2)
		
		write (*,'(X, "--- Max. and min. interpolation error: ", E10.4, X, E10.4)') maxval(abs(LUT-map_out(iext))), minval(abs(LUT-map_out(iext)))
	end if
	
end subroutine ss_interp

! Global interpolation by spherical spline (Perrin et al, 1988) [OPTIMIZED]
subroutine ss_interp_precomp(nside, lmax, stff, tens, next, LUT, iext, map_out)
	integer, intent(in)       :: nside, lmax, stff, tens, next, iext(next)
	real(DP), intent(in)      :: LUT(next)
	real(DP), intent(out)     :: map_out(0:12*nside**2-1)
	real(DP)                  :: vec(next,3), ang(next,2), vi(3), fwhm, mfac0, mfac_arr(0:0)
	integer                   :: i, j, n, l, lmin, omega_pix
	real(DP), allocatable     :: A(:,:), B(:), recfac(:,:), bl1(:,:), bl2(:,:), wl(:,:), wgt(:), bl_eff(:,:)
	complex(DPC), allocatable :: alm(:,:,:)

	n = nside2npix(nside) - 1
	fwhm = 20.0 ! 2.0 * acos(1.0 - 2.0/next) * 180.0 * 60.0 / pi

	allocate(A(next,next), B(next), source=0.0)                      ! System of linear equations
	allocate(recfac(0:1,0:lmax), source=0.0)                         ! Precomputed factors
	allocate(bl1(0:lmax,1), bl2(0:lmax,1), wl(0:lmax,1), source=1.0) ! Beams
	allocate(wgt(0:lmax), bl_eff(0:lmax,1), source=1.0)              ! Beams

	! Positions of the local extrema
	do i = 1, next
		call pix2vec_nest(nside, iext(i), vec(i,:))
		call pix2ang_nest(nside, iext(i), ang(i,1), ang(i,2))
	end do

	! Initialize RESCALE_TAB useful for calculation of spherical harmonics
	call init_rescale_tab()

	! Inverse of differential operator in harmonic space
	if (stff /= 0 .or. tens == 0) then
		lmin = 1
		bl1(0,1) = 1.0
	else
		lmin = 0
		bl1(0,1) = -1.0 / dble(tens)
	end if

	do l = 1, lmax
		bl1(l,1) = (-1.0)**(stff+1) / (dble(l)*dble(l+1) + dble(tens)) / (dble(l)*dble(l+1))**stff
	end do
	
	! Pixel window function
	call pixel_window(wl, nside)

	! Gaussian beam for "stff=0"
!	if (stff == 0) call gaussbeam(fwhm, lmax, bl2)

	! Precompute effective beam
	bl_eff(:,1) = wl(:,1) * bl1(:,1) * bl2(:,1)
	
	! Precompute weight array
	do l = 0, lmax
		wgt(l) = sqrt(dble(2*l+1) / (4.0*pi)) * bl_eff(l,1)
	end do
	
	! Precompute recursion quantities for m=0  (called ONCE instead of
	! once per Aij invocation, i.e. next*(next+1)/2 + (n+1)*next fewer calls)
	call gen_mfactor(0, mfac_arr)
	mfac0 = mfac_arr(0)
	call gen_recfactor(lmax, 0, recfac)

!	write (*,'(/, X, "--- Computing interpolation matrix of dimension ", I0, ".")') next

	! Computing "A"
	!$OMP PARALLEL PRIVATE(i, j) SHARED(A, vec, lmin, lmax, mfac0, recfac, wgt)
	!$OMP DO SCHEDULE(DYNAMIC)
	do j = 1, next
		do i = 1, j
			A(i,j) = Aij_precomp(dot_product(vec(i,:), vec(j,:)), lmin, lmax, mfac0, recfac, wgt)
			if (i /= j) A(j,i) = A(i,j)
		end do
	end do
	!$OMP END DO
	!$OMP END PARALLEL

	! Computing "B"
	B = LUT - lmin * sum(LUT)/next

!	write (*,'(X, A)') "--- Solving system of equations."

	! Calculating and storing the interpolation coefficients "X=A^(-1)*B" in array "B"
	call lsolve(next, A, B)

	deallocate(A)
	
!	write (*,'(X, A)') "--- Calculating spherical harmonic coefficients."
	
	allocate(alm(1,0:lmax,0:lmax), source=(0.0, 0.0))
	
	! Calculating the sum of "C_i*Y^(*)_lm(theta_i,phi_i)" over all pixels "i" with values in LUT
	call interp_alms(B, ang, lmax, lmax, alm)
	
	! Convolution in harmonic space and generation of output map
	do l = lmin, lmax
		alm(1,l,0:l) = alm(1,l,0:l) * wl(l,1) * bl1(l,1) * bl2(l,1)
	end do
	
	! "a_00" takes a different value when "stff>0" or "tens=0"
	if (lmin == 1) alm(1,0,0) = sqrt(4.0*pi) * sum(LUT) / next
	
!	write (*,'(X, A)') "--- Generating map."
	
	call alm2map(nside, lmax, lmax, alm, map_out) ! Generates map in RING ordering
	call convert_ring2nest(nside, map_out)        ! Goes back to NESTED ordering

!	write (*,'(X, "--- Max. and min. interpolation error: ", E10.4, X, E10.4)') maxval(abs(LUT-map_out(iext))), minval(abs(LUT-map_out(iext)))
	
	deallocate(B, wl, bl1, bl2, alm)
	
end subroutine ss_interp_precomp

! Computes the upper and lower envelopes of an input map
subroutine extrema(nside, map_in, stff, tens, lmax, map_out)
	integer, intent(in)   :: nside, stff, tens, lmax
	real(DP), intent(in)  :: map_in(0:12*nside**2-1)
	real(DP), intent(out) :: map_out(0:12*nside**2-1,4)
	integer               :: nmax, nmin, imax(12*nside**2/9), imin(12*nside**2/9)
	
	! Finding the positions and values of the local extrema of the input map
	write (*,'(/, X, A)') "- Finding local extrema."
	call local_extrema(nside, map_in, nmax, nmin, imax, imin)
	
	map_out(:,1:2) = HPX_DBADVAL
	
	map_out(imax(1:nmax),1) = map_in(imax(1:nmax))
	map_out(imin(1:nmin),2) = map_in(imin(1:nmin))
	
	write (*,'(/, X, A)') "- Computing upper envelope."
	call ss_interp(nside, lmax, stff, tens, nmax, map_in(imax(1:nmax)), imax(1:nmax), map_out(:,3))
	write (*,'(/, X, A)') "- Computing lower envelope."
	call ss_interp(nside, lmax, stff, tens, nmin, map_in(imin(1:nmin)), imin(1:nmin), map_out(:,4))
	
end subroutine extrema

! Empirical Mode Decomposition (EMD) process
subroutine emd(nside, map_in, nimf, stff, tens, lmax, imf)
	integer, intent(in)   :: nside, nimf, stff, tens, lmax
	real(DP), intent(in)  :: map_in(0:12*nside**2-1)
	real(DP), intent(out) :: imf(0:12*nside**2-1,nimf)
	integer               :: i, j, k, n, nmax, nmin, imax(12*nside**2/9), imin(12*nside**2/9)
	real(DP), allocatable :: inp(:), Emax(:), Emin(:)
	real(DP)              :: mean_stdv, last_mean_stdv, stop_crit
	
	n = nside2npix(nside) - 1
	
	allocate(inp(0:n), Emax(0:n), Emin(0:n), source=0.0)
	
	! IMF number
	do i = 1, nimf
		write (*,'(/, X, "- Computing Intrinsic Mode Function ", I0, "/", I0, ".")') i, nimf
		
		j = 1           ! Iteration number
		mean_stdv = 1.0 ! Standard deviation of the mean envelope
		stop_crit = 0.0 ! Stoppage criterion (sifting process will stop if MEAN_STDV is less than STOP_CRIT)
		
		! Sifting process
		do while (mean_stdv >= stop_crit)
			write (*,'(/, X, "-- Iteration in progress: " I0, ".")') j
			
			if (j == 1) then
				inp = map_in
				if (i /= 1) then
					do k = 1, i - 1
						inp = inp - imf(:,k)
					end do
				end if
			else
				inp = imf(:,i)
			end if
			
			! Finding the positions and values of the local extrema of the input map
			call local_extrema(nside, inp, nmax, nmin, imax, imin)

			! Computing smooth upper and lower envelopes
			call ss_interp(nside, lmax, stff, tens, nmax, inp(imax(1:nmax)), imax(1:nmax), Emax)
			call ss_interp(nside, lmax, stff, tens, nmin, inp(imin(1:nmin)), imin(1:nmin), Emin)
			
			! Update standard deviation of mean envelope
			mean_stdv = sqrt(sum(abs((Emax + Emin) / 2.0)**2) / dble(n+1))
			
			! Stoppage criterion
			if (j == 2) stop_crit = mean_stdv / 2.0
			
			write (*,'(/, X,          "-- Mean envelope SD = ", E10.4)') mean_stdv
			if (j >= 2) write (*,'(X, "-- Stoppage criteria SD < ", E10.4)') stop_crit
			
			! Force the standard deviation to decrease, otherwise exit loop
			if (j >= 2 .and. mean_stdv >= last_mean_stdv) exit
		 
			! Update last mean standard deviation of mean envelope
			last_mean_stdv = mean_stdv
			
			! Update IMF
			imf(:,i) = inp - (Emax + Emin) / 2.0
			
			! Update iteration number
			j = j + 1
		end do
	end do
	
	deallocate(inp, Emax, Emin)
	
end subroutine emd

! Halton sequence
pure function halton(index, base) result(halton_seq)
	integer, intent(in) :: index, base
	integer             :: i
	real(DP)            :: f, halton_seq

	halton_seq = 0.0
	f          = 1.0 / base
	i          = index

	do while (i > 0)
		halton_seq = halton_seq + f * mod(i, base)
		i = i / base
		f = f / base
	end do

end function halton

! First primes
subroutine first_primes(n, primes)
	integer, intent(in)  :: n
	integer, intent(out) :: primes(n)
	integer              :: count, num, i
	logical              :: is_prime

	count = 0
	num   = 2

	do while (count < n)
		is_prime = .true.
		
		do i = 2, int(sqrt(real(num)))
			if (mod(num, i) == 0) then
				is_prime = .false.
				exit
			end if
		end do

		if (is_prime) then
			count         = count + 1
			primes(count) = num
		end if

		num = num + 1

	end do

end subroutine first_primes

! Generator of unit vectors uniformly distributed on the (N_nu-1)-dimensional unit sphere
subroutine generate_directions(N_nu, ndir, dir)
	integer, intent(in)   :: N_nu, ndir
	real(DP), intent(out) :: dir(N_nu, ndir)
	integer               :: d, k
	real(DP)              :: u1, u2, r, theta
	real(DP)              :: vec(N_nu), norm
	integer, allocatable  :: primes(:)
	
	! Generate first N_nu prime numbers (bases for Halton)
	allocate(primes(N_nu))
	call first_primes(N_nu, primes)
	
	! Generate quasi-random directions
	do d = 1, ndir

		! Build Gaussian vector using Halton + Box-Muller
		do k = 1, N_nu, 2

			u1 = halton(d, primes(k))
			if (k+1 <= N_nu) then
				u2 = halton(d, primes(k+1))
			else
				u2 = halton(d+ndir, primes(1)) ! fallback
			end if

			! Avoid log(0)
			if (u1 <= 1d-12) u1 = 1d-12

			r      = sqrt(-2.0 * log(u1))
			theta  = 2.0 * pi * u2
			vec(k) = r * cos(theta)
			
			if (k+1 <= N_nu) vec(k+1) = r * sin(theta)

		end do

		! Normalize to unit sphere
		norm = sqrt(sum(vec**2))
		if (norm == 0.0) then
			vec    = 0.0
			vec(1) = 1.0
			norm   = 1.0
		end if

		dir(:,d) = vec / norm

	end do

	deallocate(primes)

end subroutine generate_directions

! Multivariate Empirical Mode Decomposition (MEMD) process
subroutine memd(nside, map_in, N_nu, nimf, stff, tens, lmax, ndir, imf, mask)
	integer, intent(in)   :: nside, N_nu, nimf, stff, tens, lmax, ndir
	real(DP), intent(in)  :: map_in(0:12*nside**2-1, N_nu)
	real(DP), intent(in), optional :: mask(0:12*nside**2-1,1)
	real(DP), intent(out) :: imf(0:12*nside**2-1, N_nu, nimf)
	integer               :: i, j, k, d, n, nmax, nmin
	integer               :: imax(12*nside**2/9), imin(12*nside**2/9)
	real(DP), allocatable :: inp(:,:), proj(:), Emax(:), Emin(:)
	real(DP), allocatable :: mean_env(:,:), dir(:,:)
	real(DP)              :: mean_stdv, last_mean_stdv, stop_crit

	n = nside2npix(nside) - 1

	allocate(inp(0:n, N_nu), proj(0:n), Emax(0:n), Emin(0:n))
	allocate(mean_env(0:n, N_nu), dir(N_nu, ndir))
	
	! Generate projection directions (uniform on hypersphere)
	call generate_directions(N_nu, ndir, dir)

	! IMF loop
	do i = 1, nimf
		write (*,'(/, X, "- Computing Multivariate IMF ", I0, "/", I0, ".")') i, nimf

		j = 1
		mean_stdv = 1.0
		stop_crit = 0.0

		do while (mean_stdv >= stop_crit)
			write (*,'(/, X, "-- Iteration: ", I0, ".")') j

			! Residual input
			if (j == 1) then
				inp = map_in
				if (i /= 1) then
					do k = 1, i-1
						inp = inp - imf(:,:,k)
					end do
				end if
			else
				inp = imf(:,:,i)
			end if

			mean_env = 0.0

			! Loop over projection directions
			do d = 1, ndir
				write(*,'(/, X, "-- Projection ", I0, "/", I0, " with direction vector ", A, ".")') d, ndir, trim(vec_str(dir(:,d)))

				! Calculate projection to get scalar signal
				proj = 0.0
				do k = 1, N_nu
					proj = proj + dir(k,d) * inp(:,k)
				end do

				! Find extrema of projected signal
				call local_extrema(nside, proj, nmax, nmin, imax, imin, mask=mask)
				
				write(*,'(4X, "Number of local maxima = ", I0)') nmax
				write(*,'(4X, "Number of local minima = ", I0)') nmin

				! Envelopes in projected space
				call ss_interp_precomp(nside, lmax, stff, tens, nmax, proj(imax(1:nmax)), imax(1:nmax), Emax)
				call ss_interp_precomp(nside, lmax, stff, tens, nmin, proj(imin(1:nmin)), imin(1:nmin), Emin)

				! Mean envelope of projection
				proj = (Emax + Emin) / 2.0

				! Lift back to multivariate space
				do k = 1, N_nu
					mean_env(:,k) = mean_env(:,k) + proj * dir(k,d)
				end do

			end do

			! Average over directions
			mean_env = mean_env / dble(ndir)
			
			! Stopping criterion (vector version)
			mean_stdv = sqrt(sum(mean_env**2) / dble((n+1)*N_nu))

			if (j == 2) stop_crit = mean_stdv / 2.0

			write (*,'(/, X, "-- Mean envelope SD = ", E10.4)') mean_stdv

			if (j >= 2 .and. mean_stdv >= last_mean_stdv) exit

			last_mean_stdv = mean_stdv

			! Update IMF (vector)
			imf(:,:,i) = inp - mean_env

			j = j + 1

		end do
	end do

	deallocate(inp, proj, Emax, Emin, mean_env, dir)

end subroutine memd

! Build a string containing the entries of a 1D array
pure function vec_str(v)
	real(DP), intent(in)      :: v(:)
	character(len=32*size(v)) :: vec_str
	character(len=16)         :: tmp
	integer                   :: k
	
	write(tmp, '(ES10.4E1)') v(1)
	
	vec_str = '(' // trim(adjustl(tmp))

	do k = 2, size(v)
		write(tmp, '(ES10.4E1)') v(k)
		vec_str = trim(vec_str) // ', ' // trim(adjustl(tmp))
	end do

	vec_str = trim(vec_str) // ')'

end function vec_str

! Calculation of the elements of the interpolation matrix (unstable)
pure function G(cth, bl, lmin, lmax)
	integer, intent(in)  :: lmin, lmax
	real(DP), intent(in) :: cth, bl(0:lmax,1)
	real(DP)             :: P(0:2), G
	integer              :: l
	
	G = 0.0
	P(0) = 1.0
	P(1) = cth
	
	! Initial terms of summation
	if (lmin == 0) G = G + P(0) * bl(0,1) / (4.0*pi)
	if (lmin <= 1) G = G + P(1) * bl(1,1) * 3.0 / (4.0*pi)
	
	! Summation
	do l = 2, lmax
		P(mod(l,3)) = (dble(2*l-1) * cth * P(mod(l-1,3)) - dble(l-1) * P(mod(l-2,3))) / dble(l)
		if (l >= lmin) G = G + P(mod(l,3)) * bl(l,1) * dble(2*l+1) / (4.0*pi)
	end do
	
end function G

! Calculation of the elements of the interpolation matrix (stable)
function Aij(cth, bl, lmin, lmax)
	integer, intent(in)  :: lmin, lmax
	real(DP), intent(in) :: cth, bl(0:lmax,1)
	real(DP)             :: sth, mfac(0:0), recfac(0:1,0:lmax), lam_lm(0:lmax), Aij
	integer              :: l
	
	Aij = 0.0
	sth = sqrt(1.0 - cth**2)
	
	! Recursion factor used in "lambda_00" calculation
	call gen_mfactor(0, mfac)
	
	! Generate recursion factors useful for "lambda_l0"
	call gen_recfactor(lmax, 0, recfac)
	
	! Compute "lambda_l0(theta_ij)=sqrt((2*l+1)/4*pi)*Pl(cos(theta_ij))" for all "l" for "m=0"
	if (sth > epsilon(sth)) then
		call do_lambda_lm(lmax, 0, abs(cth), sth, mfac(0), recfac, lam_lm)
	else
		forall (l=0:lmax) lam_lm(l) = sqrt(dble(2*l+1) / (4.0*pi))
	end if
	
	! Southern-hemisphere parity correction
	if (cth < 0.0) forall (l=0:lmax) lam_lm(l) = (-1.0)**l * lam_lm(l)
	
	! Inner summation
	do l = lmin, lmax
		Aij = Aij + bl(l,1) * sqrt(dble(2*l+1) / (4.0*pi)) * lam_lm(l)
	end do
	
end function Aij

! Calculation of the elements of the interpolation matrix (stable and faster)
! Precomputed quantities (computed ONCE outside any loop):
!   mfac0   : scalar returned by gen_mfactor for m=0
!   recfac  : recursion factors from gen_recfactor for m=0, all l up to lmax
!   wgt     : wgt(l) = sqrt((2l+1)/4pi) * bl(l,1), for l = 0..lmax
!             (folds the beam/operator spectrum and the spherical harmonic
!              normalisation into a single array, removing both the sqrt
!              and the bl lookup from the inner summation loop)
function Aij_precomp(cth, lmin, lmax, mfac0, recfac, wgt)
	integer,  intent(in) :: lmin, lmax
	real(DP), intent(in) :: cth, mfac0, recfac(0:1, 0:lmax), wgt(0:lmax)
	real(DP)             :: Aij_precomp, sth, lam_lm(0:lmax)
	integer              :: l
 
	Aij_precomp = 0.0_DP
	sth = sqrt(1.0_DP - cth**2)
 
	! Compute "lambda_l0(theta)" for all "l" (same branch logic as original Aij)
	if (sth > epsilon(sth)) then
		call do_lambda_lm(lmax, 0, abs(cth), sth, mfac0, recfac, lam_lm)
	else
		forall (l = 0:lmax) lam_lm(l) = sqrt(dble(2*l+1) / (4.0_DP*pi))
	end if
 
	! Southern-hemisphere parity correction
	if (cth < 0.0_DP) forall (l = 0:lmax) lam_lm(l) = (-1.0_DP)**l * lam_lm(l)
 
	! Inner summation: "wgt" already contains sqrt((2l+1)/4pi)*bl(l,1)
	do l = lmin, lmax
		Aij_precomp = Aij_precomp + wgt(l) * lam_lm(l)
	end do
 
end function Aij_precomp

! Generates factor used in "lambda_mm" calculation for all "m" in "0<=m<=m_max"
subroutine gen_mfactor(m_max, m_fact)
	integer, intent(in)   :: m_max
	real(DP), intent(out) :: m_fact(0:m_max)
	integer               :: m

	! fact(m) = fact(m-1) * sqrt((2m+1)/(2m))
	m_fact(0) = 1.0
	do m = 1, m_max
		m_fact(m) = m_fact(m-1) * sqrt(dble(2*m+1) / dble(2*m))
	end do

	! Log_2 ( fact(m) / sqrt(4 Pi) )
	do m = 0, m_max
		m_fact(m) = log(SQ4PI_INV * m_fact(m)) * ALN2_INV
	end do

end subroutine gen_mfactor

! Generates recursion factors used to computes the "lambda_lm" of degree "m" for all "l" in "m<=l<=l_max"
subroutine gen_recfactor( l_max, m, recfac)
	integer, intent(in)   :: l_max, m
	real(DP), intent(out) :: recfac(0:1, 0:l_max)
	real(DP)              :: fm2, fl2
	integer               :: l

	recfac(0:1,0:m) = 0.0
	fm2 = dble(m)**2
	do l = m, l_max
		fl2 = dble(l+1)**2
		recfac(0,l) = sqrt((4.0*fl2-1.0) / (fl2-fm2))
	end do
	
	! Put outside the loop because of problem on some compilers
	recfac(1,m:l_max) = 1.0 / recfac(0,m:l_max)

end subroutine gen_recfactor

! Initialize RESCALE_TAB array
subroutine init_rescale_tab()
	integer           :: s, smax
	real(DP)          :: logOVFLOW
	
	logOVFLOW=log(FL_LARGE)
	smax = int(log(MAX_DP) / logOVFLOW)
	rescale_tab(RSMIN:RSMAX) = 0.0
	do s = -smax, smax
		rescale_tab(s) = FL_LARGE**s
	end do
	rescale_tab(0) = 1.0
	
end subroutine init_rescale_tab

! Computes scalar "lambda_lm(theta)" for all "l" in "m<=l<=lmax" for a given "m" and a given "theta"
subroutine do_lambda_lm(lmax, m, cth, sth, mfac, recfac, lam_lm)
	integer, intent(in)   :: lmax,  m
	real(DP), intent(in)  :: cth, sth, mfac, recfac(0:1,0:lmax)
	real(DP), intent(out) :: lam_lm(0:lmax)
	real(DP)              :: log2val, dlog2lg, ovflow, unflow, corfac, lam_mm, lam_0, lam_1, lam_2
	integer               :: scalel, l, l_min
	
	! Define constants
	ovflow = rescale_tab(1)
	unflow = rescale_tab(-1)
	l_min = l_min_ylm(m, sth)
	dlog2lg = real(LOG2LG, kind=DP)

	! Computes "lambda_mm"
	log2val = mfac + m*log(sth) * ALN2_INV     ! "log_2(lambda_mm)"
	scalel = int(log2val / dlog2lg)
	corfac = rescale_tab(max(scalel,RSMIN))
	lam_mm = 2.0**(log2val - scalel * dlog2lg) ! Rescaled "lambda_mm"
	if (IAND(m,1)>0) lam_mm = -lam_mm          ! Negative for odd "m"

	lam_lm(0:lmax) = 0.0
	
	! "l=m"
	lam_lm(m) = lam_mm * corfac

	! "l>m"
	lam_0 = 0.0
	lam_1 = 1.0
	lam_2 = cth * lam_1 * recfac(0,m)
	do l = m+1, lmax
		! Do recursion
		if (l >= l_min) then
			lam_lm(l) = lam_2 * corfac * lam_mm
		end if
		lam_0 = lam_1 * recfac(1,l-1)
		lam_1 = lam_2
		lam_2 = (cth * lam_1 - lam_0) * recfac(0,l)

		! Do dynamic rescaling
		if (abs(lam_2) > ovflow) then
			lam_1 = lam_1 * unflow
			lam_2 = lam_2 * unflow
			scalel = scalel + 1
			corfac = rescale_tab(max(scalel,RSMIN))
		else if (abs(lam_2) < unflow .and. abs(lam_2) /= 0.0) then
			lam_1 = lam_1 * ovflow
			lam_2 = lam_2 * ovflow
			scalel = scalel - 1
			corfac = rescale_tab(max(scalel,RSMIN))
		end if

	end do
end subroutine do_lambda_lm

! Display usage information
subroutine usage()
	write (*,'(/, X, A, /)') "Usage:"
	write (*,*) "For Hilbert-Huang transform: hht IFN IMF NIT STF TNS GSP"
	write (*,'(X, A, /)') "For finding local extrema: hht -lext IFN STF TNS SMP"
	write (*,*) "IFN = Input file name"
	write (*,*) "IMF = Number of Intrinsic Mode Functions"
	write (*,*) "NIT = Maximum number of iterations in the Empirical Mode Decomposition"
	write (*,*) "TNS = Tension parameter"
	write (*,*) "GSP = FWHM parameter for Gaussian smoothing (in arcminutes)"
	write (*,'(X, A, /)') "STF = Stiffness"
	
end subroutine usage

!! Compute the covariance matrix of "nvar" variables with "nobs" observations stacked in a matrix "A" with weights "w"
!subroutine weighted_cov_mtx(nvar, nobs, A, w, C)
!	integer, intent(in)     :: nvar, nobs
!	real(DP), intent(in)    :: A(nobs,nvar), w(nobs)
!	real(DP), intent(out)   :: C(nvar,nvar)
!	real(DP)                :: wsum, mean_i, mean_j
!	integer                 :: i, j
!
!	wsum = sum(w)
!
!	do i = 1, nvar
!		do j = i, nvar
!			mean_i = sum(w * A(:,i)) / wsum
!			mean_j = sum(w * A(:,j)) / wsum
!			C(i,j) = sum( w * (A(:,i)-mean_i) * (A(:,j)-mean_j) ) / wsum
!			C(j,i) = C(i,j)
!		end do
!	end do
!
!end subroutine

! Compute the covariance matrix of "nvar" variables with "nobs" observations stacked in a matrix "A" with weights "w"
subroutine weighted_cov_mtx(nvar, nobs, A, w, C)
	integer, intent(in)     :: nvar, nobs
	real(DP), intent(in)    :: A(nobs,nvar), w(nobs)
	real(DP), intent(out)   :: C(nvar,nvar)
	real(DP)                :: wsum, mean_i, mean_j
	integer                 :: i, j

	wsum = sum(w)

	do i = 1, nvar
		do j = i, nvar
			C(i,j) = sum( w * A(:,i) * A(:,j) ) / wsum
			C(j,i) = C(i,j)
		end do
	end do

end subroutine

! Compute the covariance matrix of "nvar" variables with "nobs" observations stacked in a matrix "A"
subroutine cov_mtx(nvar, nobs, A, C)
	integer, intent(in)     :: nvar, nobs
	real(DP), intent(inout) :: A(nobs,nvar)
	real(DP), intent(out)   :: C(nvar,nvar)
	integer                 :: i, j
	real(DP)                :: mean_i
	real(DP), allocatable   :: DeltaA(:,:)
	
	! Provide storage for temporary map containing deviations from the mean for each element of A
	allocate(DeltaA(nobs,nvar))
	
	! Substract the mean of each variable from each observation
	do i = 1, nvar
!		mean_i = sum(A(:,i)) / nobs
		mean_i = 0.0
		DeltaA(:,i) = A(:,i) - mean_i
	end do
	
	! Compute the covariance matrix
	do i = 1, nvar
		do j = i, nvar
			C(i,j) = dot_product(DeltaA(:,i), DeltaA(:,j)) / nobs
			C(j,i) = C(j,i)
		end do
	end do
	
	! Free memory from allocatable objects
	deallocate(DeltaA)
end subroutine cov_mtx

! Compute the inverse of a matrix "A" of order "n"
subroutine mtx_inv(n, A, A_inv)
	integer, intent(in)   :: n
	real(DP), intent(in)  :: A(n,n)
	real(DP), intent(out) :: A_inv(n,n)
	real(DP), allocatable :: work(:)
	real(DP)              :: get_lwork(1)
	integer               :: ipiv(n), info, lwork, i, j
	
	A_inv = A ! Initialize inverse matrix
	info  = 0 ! Status indicator ("stat/=0" indicates an error)
	
	! Calculate the optimal size of the WORK array as the first entry of the GET_LWORK array
	call dsytrf('U', n, A_inv, n, ipiv, get_lwork, -1, info)

	! Provide storage for the WORK array
	lwork = int(get_lwork(1))
	allocate(work(lwork), source=0.0)
	
	! LAPACK subroutine to compute the factorization of a real symmetric matrix A using the Bunch-Kaufman diagonal pivoting method
	call dsytrf('U', n, A_inv, n, ipiv, work, lwork, info)

	! Stop the program if necessary
	if (info /= 0) call fatal_error('Singular matrix in subroutine MTX_INV()')
	
	! Free memory from Work
	deallocate(work)

	! Resizing WORK array
	allocate(work(n), source=0.0)

	! LAPACK subroutine to get "A^(-1)" for a symmetric matrix "A"
	call dsytri('U', n, A_inv, n, ipiv, work, info)
	
	! Stop the program if necessary
	if (info /= 0) call fatal_error('Singular matrix in subroutine MTX_INV()')
	
	! Off-diagonal elements of the matrix in its lower triangular part
	do i = 2, n
		do j = 1, i-1
			A_inv(i,j) = A_inv(j,i)
		end do
	end do
	
	deallocate(work)
	
end subroutine mtx_inv

end module
