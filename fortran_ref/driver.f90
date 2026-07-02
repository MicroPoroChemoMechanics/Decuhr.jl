! Driver de comparaison : appelle la subroutine Fortran DECUHR de reference
! (Espelid & Genz, 1994) sur les 6 memes cas que le harnais Julia (cases.jl),
! avec exactement les memes parametres :
!   epsabs=1e-8  epsrel=1e-6  maxpts(budget)=1e6  wrksub=50000  emax=20
!   minpts=0  key=0  restar=0  numfun=1
! Imprime, par cas : resultat, erreur relative, IFAIL, NEVAL.
! But : determiner si le Fortran renvoie les memes IFAIL (notamment IFAIL=1,
! "MaxIters") et les memes valeurs que le portage Julia ① — en particulier
! pour C1 (x*y)^-0.5 et C4 (x*y*z)^-1/3.

program cmp_decuhr
  implicit none
  integer, parameter :: nwmax = 1500000, iwmax = 110000
  double precision :: work(nwmax)
  integer :: iwork(iwmax)
  double precision :: pi, half
  external :: f1, f2, f3, f4, f5, f6

  pi   = 4.d0*atan(1.d0)
  half = pi/2.d0

  write(*,'(A)') '================================================================'
  write(*,'(A)') ' DECUHR Fortran de reference - 6 cas communs'
  write(*,'(A)') ' epsabs=1e-8 epsrel=1e-6 maxpts=1000000 wrksub=50000 emax=20'
  write(*,'(A)') '================================================================'
  write(*,'(A)') 'case          result                exact                relerr      ifail neval'

  call run('C1_xy_vertex', f1, 2, [0.d0,0.d0],   [1.d0,1.d0],   2, -0.5d0,        0, 4.d0,                 work, iwork, nwmax)
  call run('C2_radial   ', f2, 2, [0.d0,0.d0],   [1.d0,1.d0],   2, -1.0d0,        0, 2.d0*log(1.d0+sqrt(2.d0)), work, iwork, nwmax)
  call run('C3_smooth   ', f3, 2, [0.d0,0.d0],   [half,half],   1,  0.0d0,        0, 1.d0,                 work, iwork, nwmax)
  call run('C4_3d_vertex', f4, 3, [0.d0,0.d0,0.d0],[1.d0,1.d0,1.d0], 3, -1.d0/3.d0, 0, 27.d0/8.d0,        work, iwork, nwmax)
  call run('C5_poly     ', f5, 2, [0.d0,0.d0],   [1.d0,1.d0],   1,  0.0d0,        0, 2.d0/3.d0,            work, iwork, nwmax)
  call run('C6_log      ', f6, 2, [0.d0,0.d0],   [1.d0,1.d0],   2,  0.0d0,        1, 2.d0,                 work, iwork, nwmax)

end program cmp_decuhr

subroutine run(name, fun, ndim, alo, ahi, singul, alpha, logf, exact, work, iwork, nwmax)
  implicit none
  character(*), intent(in) :: name
  external :: fun
  integer, intent(in) :: ndim, singul, logf, nwmax
  double precision, intent(in) :: alo(ndim), ahi(ndim), alpha, exact
  double precision :: work(*)
  integer :: iwork(*)
  double precision :: a(15), b(15), result(1), abserr(1), relerr
  integer :: j, neval, ifail
  do j = 1, ndim
     a(j) = alo(j)
     b(j) = ahi(j)
  end do
  result(1) = 0.d0
  abserr(1) = 0.d0
  call decuhr(ndim, 1, a, b, 0, 1000000, fun, singul, alpha, logf, &
       1.d-8, 1.d-6, 0, 50000, nwmax, 0, 20, result, abserr, neval, ifail, work, iwork)
  relerr = abs(result(1) - exact) / abs(exact)
  write(*,'(A,2X,ES22.15,2X,ES20.13,2X,ES10.3,2X,I3,2X,I9)') &
       name, result(1), exact, relerr, ifail, neval
end subroutine run

! ----- integrands : FUNSUB(NDIM, X, NUMFUN, FUNVLS) -----
subroutine f1(ndim, x, numfun, funvls)   ! (x*y)^(-1/2)
  integer ndim, numfun
  double precision x(ndim), funvls(numfun)
  funvls(1) = (x(1)*x(2))**(-0.5d0)
end subroutine
subroutine f2(ndim, x, numfun, funvls)   ! 1/sqrt(x^2+y^2)
  integer ndim, numfun
  double precision x(ndim), funvls(numfun)
  funvls(1) = 1.d0/sqrt(x(1)**2 + x(2)**2)
end subroutine
subroutine f3(ndim, x, numfun, funvls)   ! sin(x)cos(y)
  integer ndim, numfun
  double precision x(ndim), funvls(numfun)
  funvls(1) = sin(x(1))*cos(x(2))
end subroutine
subroutine f4(ndim, x, numfun, funvls)   ! (x*y*z)^(-1/3)
  integer ndim, numfun
  double precision x(ndim), funvls(numfun)
  funvls(1) = (x(1)*x(2)*x(3))**(-1.d0/3.d0)
end subroutine
subroutine f5(ndim, x, numfun, funvls)   ! x^2+y^2
  integer ndim, numfun
  double precision x(ndim), funvls(numfun)
  funvls(1) = x(1)**2 + x(2)**2
end subroutine
subroutine f6(ndim, x, numfun, funvls)   ! -log(x*y)
  integer ndim, numfun
  double precision x(ndim), funvls(numfun)
  funvls(1) = -log(x(1)*x(2))
end subroutine
