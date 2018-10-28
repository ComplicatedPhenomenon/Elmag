!=============================================================================!
!=============================================================================!
module user_result
  implicit none
  save
  integer, parameter ::            &
     n_bin = 81,                   &  ! # of bins in agamma
     n_bin2 = 2,                   &  ! # angular bin
     n_bin3 = 4,                   &  ! # time bins
     n_bin4 = 7                       ! # time bins
  double precision spec(n_bin,0:1),spec_tot(n_bin,0:1)     ! diffuse spectrum
  double precision gam_th(n_bin,n_bin2),gam_th_tot(n_bin,n_bin2) ! spectrum
  double precision spec_t(n_bin,n_bin4),spec_t_tot(n_bin,n_bin4) ! spectrum
  double precision cum_t(n_bin,n_bin4)                      ! cum. spectrum
  integer n_reg
  character*60 ::  filename = '_z15'         ! name in output
end module user_result
!=============================================================================!
!=============================================================================!
module user_variables
  implicit none
  save
  integer,parameter :: model = 4   ! EBL model: 
                                   ! 1: Kneiske Dole best fit, 2: lower limit,
                                   ! 3: Franceschini arXiv:0805.1841
                                   ! 4: Model C Finke et al. arXiv:0905.1115
                                   ! 5: Gilmore et al. arXiv:1104.0671
                                   ! 6: Dominguez et al. arXiv:1007.1459

  integer,parameter :: nmax=1*10**3      ! # injected photons

  double precision,parameter ::        &  
  ethr=1.d5,                           & ! energy threshold, eV
  egmax=2.d15,                         & ! maximal gamma energy, eV
  ir_rescale=1.d0,                     & ! multipl. recaling of IR background
  cohlnth=3.086d21*1.d3,               & ! IGMF coherence length
  th_jet=6.d0,                         & ! jet opening angle/degrees
  a_smp=0.9d0                            ! 1 max. weighted sampling, 0-no sampling
end module user_variables
!=============================================================================!
!=============================================================================!
program Elmagcasc
  use mpi
  use user_variables, only : nmax
  use user_result
  implicit none
  integer myid,ierr,n_proc,n_array

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,ierr)

  call init(myid,n_proc)
  
  call user_main(myid,nmax)

  n_array = 2*n_bin
  call MPI_REDUCE(spec,spec_tot,n_array,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                  MPI_COMM_WORLD,ierr)        ! sum individal arrays spec
  call MPI_REDUCE(gam_th,gam_th_tot,n_array,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                  MPI_COMM_WORLD,ierr)        ! sum individal arrays gam_th
  n_array = n_bin*n_bin4
  call MPI_REDUCE(spec_t,spec_t_tot,n_array,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                  MPI_COMM_WORLD,ierr)        ! sum individal arrays spec


  if (myid==0) then
     call user_output(nmax,n_proc)  ! make plots
     call banner(n_proc,1)
  end if
  close(99)               

  call MPI_Finalize(ierr)

end program Elmagcasc
!=============================================================================!
!=============================================================================!
subroutine user_main(myid,nmax)
!-----------------------------------------------------------------------------
! calls: plots,initial_particle,cascade
!-----------------------------------------------------------------------------
  use internal, only : bal
  use user_result
  implicit none
  integer myid,nl,icq,nmax
  double precision z,e0,weight

  if (myid==0) then
     call plots
     write(*,*)'start with nmax = ',nmax
  end if

  z = 0.15d0                                         ! initial redshift 
  do nl=1,nmax
     call initial_particle(e0,weight)                ! generate initial energy
     icq = 0                                         ! (0 - gamma, +-1 - e+-)
     call cascade(icq,e0,weight,z)                   ! starts e/m cascade
     if (myid==0 .and. mod(nl*100,nmax)==0) write(*,*)'nl=',nl
  enddo

end subroutine user_main
!=============================================================================!
!=============================================================================!
subroutine user_output(n_max,n_proc) 
  use user_variables
  use user_result
  implicit none
  integer n_max,n_proc,i
  
  spec_tot = spec_tot/dble(n_proc*n_max)
  gam_th_tot = gam_th_tot/dble(n_proc*n_max)
  spec_t_tot = spec_t_tot/dble(n_proc*n_max)

  cum_t(:,1) = spec_t_tot(:,1)
  do i=2,n_bin4
     cum_t(:,i) = cum_t(:,i-1)+spec_t_tot(:,i)
  end do

  open(unit=31,file='Output/spec_diff'//filename) ! energy spectrum g,e
  open(unit=32,file='Output/spec_95'//filename)   ! energy spectrum g in PSF
  open(unit=33,file='Output/spec_t'//filename)    ! energy-time spectra g in PSF
  open(unit=34,file='Output/spec_c'//filename)    ! energy-time spectra g in PSF
  do i=1,n_bin-1
     write(31,*) real(ethr*(egmax/ethr)**((i-.5d0)/(n_bin-1))), &
          real(spec_tot(i,0)),real(spec_tot(i,1))
     write(32,*) real(ethr*(egmax/ethr)**((i-.5d0)/(n_bin-1))), &
          real(gam_th_tot(i,1)),real(gam_th_tot(i,2))
     write(33,*) real(ethr*(egmax/ethr)**((i-.5d0)/(n_bin-1))), &
          real(spec_t_tot(i,:))
     write(34,*) real(ethr*(egmax/ethr)**((i-.5d0)/(n_bin-1))), &
          real(cum_t(i,:))
  end do
  close(31)
  close(32)
  close(33)
  close(34)

end subroutine user_output
!=============================================================================!
!=============================================================================!
subroutine register(e0,theta,dt,weight,icq)
  use user_variables
  use user_result
  implicit none
  integer icq,i,j,k
  double precision e0,weight,theta,dt,theta95
  double precision thereg_en

! diffuse energy spectrum:
  i=min(n_bin,int(log(e0/ethr)/log(egmax/ethr)*(n_bin-1))+1)
  i=max(i,1)
  spec(i,abs(icq))=spec(i,abs(icq))+weight*e0/log(egmax/ethr)*(n_bin-1) 

  if (icq.ne.0) return                                 ! forget electrons

! gamma energy spectrum, inside/outside 95% PSF:
  theta95 = thereg_en(e0)
  if (theta<0.d0) write(*,*) 'theta',theta
  if (theta<theta95) then  
     j=1
  else
     j=2
  end if
  gam_th(i,j)=gam_th(i,j)+weight*e0/log(egmax/ethr)*(n_bin-1)

  if (j==1) then                                   ! only gamma's inside PSF
     if (dt==0.d0) then
        k=1
     else
        k = log10(dt)+1
     end if
     if (k<1) k=1
     if (k>n_bin4) k=n_bin4
     spec_t(i,k)=spec_t(i,k)+weight*e0/log(egmax/ethr)*(n_bin-1)
  end if

! angular profile:

  n_reg = n_reg+1

end subroutine register
!=============================================================================!
!=============================================================================!
! sample initial photon/electron from the input spectra                       !
!=============================================================================!
subroutine initial_particle(e0,weight)
!-----------------------------------------------------------------------------
! input:
!        e0 - particle energy;
!        weight - initial particle weight
!-----------------------------------------------------------------------------
  use user_variables
  implicit none
  double precision e0,weight,emin,ebreak,gam1,gam2
  double precision psran
        
  emin=1d9                             !low energy cutoff of injected photons
  ebreak=2.d12
  gam1=-2.d0
  gam2=-2.d0                               ! primary slope (dE/dE ~ E0^gam)
      
  e0=emin*(egmax/emin)**psran()                      !energy: uniform in ln E 
  if (e0<ebreak) then
     weight=(e0/ebreak)**(gam1+1.d0)*log(egmax/emin)
  else
     weight=(e0/ebreak)**(gam2+1.d0)*log(egmax/emin) 
  endif

end subroutine initial_particle
!=============================================================================!
!=============================================================================!
!           strength of EGMF/Gauss                                            !
!=============================================================================!
double precision function bemf(r)
!-----------------------------------------------------------------------------
! input:
!        r - distance to Earth/cm
!-----------------------------------------------------------------------------
  implicit none
  double precision r
  bemf = 1.d-17 ! G
end function bemf
!=============================================================================!
!=============================================================================!
!         theta_PSF/degree (95%) from Fermi                                   !
!=============================================================================!
double precision function thereg_en(en)
  implicit none
  double precision en!,thereg_en

  if (en<1.d6) then                                 ! avoid theta>360 degrees 
     thereg_en = 359.d0
     return
  end if

  if (en .le. 3e11) then                                             ! Fermi
     thereg_en = 10.d0**7.12d0*en**(-0.766d0)             
     if (en>1d9) thereg_en = thereg_en + 0.2d0*exp(-1.d10/en)
  else                                                               ! HESS
     thereg_en = 0.11d0                                      
  endif

end function thereg_en
!=============================================================================!
!=============================================================================!
subroutine plots
  implicit none
  integer icq
  double precision E,x(3),z0,rate(2),psf95
  double precision rate_EBL_tab,thereg_en

  x(1) = 125.*3.086d24
  x(2) = 808.*3.086d24
  x(3) = 1500.*3.086d24

  z0 = 0.d0
  icq = 0 

  E = 1.d8
  open(41,file='Output/lint')
  open(42,file='Output/tau')
  open(43,file='Output/psf95')
  do
     E=E*1.1d0
     rate(1) = rate_EBL_tab(E,z0,icq)
     rate(2) = rate_EBL_tab(E,z0,1)
     write(41,*) real(E),rate(1)*3.086d24,rate(2)*3.086d24
     psf95 = thereg_en(E)
     write(43,*) real(E),psf95
     if (E>1d21) exit
  end do
  close(41)
  close(42)
  close(43)

end subroutine plots
!=============================================================================!
!=============================================================================!
