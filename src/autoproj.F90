!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the Wannier90 code and !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the Wannier90      !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the Wannier90 code is www.wannier.org       !
!                                                            !
! The Wannier90 code is hosted on GitHub:                    !
!                                                            !
! https://github.com/wannier-developers/wannier90            !
!------------------------------------------------------------!

module w90_autoproj
  !! Routines to calculate automatically the starting projections

  use w90_constants, only : cmplx_1
  use w90_io,        only : io_error
  use w90_parameters, only: auto_proj, auto_proj_flavour, u_matrix,&
       AUTOPROJ_FLAVOUR_BLOCH_PHASES, AUTOPROJ_FLAVOUR_OBSTRUCTION_MATRIX 
  
  implicit none

  private

  public :: autoproj_calc_u_matrix
  
contains
  
  !==================================================================!
  subroutine autoproj_calc_u_matrix()
  !==================================================================!
    implicit none

    if (.not. auto_proj) then
       call io_error('Cannot call autoproj_calc_u_matrix if auto_proj is .false.')
    end if
    
    if (auto_proj_flavour .eq. AUTOPROJ_FLAVOUR_BLOCH_PHASES) then
       call autoproj_calc_u_matrix_bloch_phases()
    elseif (auto_proj_flavour .eq. AUTOPROJ_FLAVOUR_OBSTRUCTION_MATRIX) then
       call autoproj_calc_u_matrix_obstruction()
       
       
    else
       call io_error('Unknown auto_proj_flavour value in autoproj_calc_u_matrix')
    end if
    

  end subroutine autoproj_calc_u_matrix


  !==================================================================!
  subroutine autoproj_calc_u_matrix_bloch_phases()
  !==================================================================!

    use w90_parameters, only :num_wann, num_kpts

    implicit none

    integer :: n, m
    
    do n=1,num_kpts
       do m=1,num_wann
          u_matrix(m,m,n)=cmplx_1
       end do
    end do
    
  end subroutine autoproj_calc_u_matrix_bloch_phases

  !==================================================================!
  subroutine autoproj_calc_u_matrix_obstruction()
  !==================================================================!

    use w90_parameters, only :kpt_latt, m_matrix, nntot, num_wann, num_kpts
    
    implicit none

    !call io_error('Automatic projections using the obstruction matrix method not implemented yet')
    
    integer :: n, m, K, N1, N2, N3
    integer, dimension(1:3) :: array_ijk = [0, 0, 0]

    N1 = mp_grid(1)
    N2 = mp_grid(2)
    N3 = mp_grid(3)

    do n=1,num_kpts
       do m=1,num_wann
          u_matrix(m,m,n)=cmplx_1
       end do
    end do
    
    do i = 1,N1    
       ! " u_matrix(:,:,ijk_to_K(i,1,1)) = propagate(Id_num_wann,ijk_to_K(i,1,1))
    end do

    !call K_to_ijk(num_kpts, array_ijk)
    !call ijk_to_K(array_ijk(1), array_ijk(2), array_ijk(3), K)
    !
    !print *, "num_kpts = ", num_kpts
    !print *, "(i,j,k) = ", array_ijk(:)
    !print *, "kpt_latt = ", kpt_latt(:,K)

  end subroutine autoproj_calc_u_matrix_obstruction

  
  !==================================================================!
  subroutine K_to_ijk(n, array_ijk)
  !==================================================================!
    
    use w90_parameters, only :num_kpts, mp_grid

    implicit none

    integer :: n, N1, N2, N3, temp, i, j, k
    integer, dimension(1:3), intent(out) :: array_ijk

    N1 = mp_grid(1)
    N2 = mp_grid(2)
    N3 = mp_grid(3)

    if ((n .le. 0) .or. (n .ge. num_kpts+1)) call io_error('In K_to_ijk(): input n not within bounds')

    if (N1*N2*N3 .ne. num_kpts) then
        call io_error('Size of k-point mesh does not fit with total size')
    endif
   
    ! Assuming that n = (i-1)*(N2*N3) + (j-1)*N3 + k
    ! TODO check that they correspond to the correct k-points
    
    i = ((n-1) / (N2*N3)) + 1 
    temp = (n-1) - (i-1)*N2*N3
    j = (temp / N3) + 1
    k = temp - (j-1)*N3 + 1

    if (n .ne. (i-1)*N2*N3 + (j-1)*N3 + k) call io_error('Wrong computation of i, j, k indices')

    array_ijk = [i, j, k]

    return
        
  end subroutine K_to_ijk
  
  !==================================================================!
  subroutine ijk_to_K(i,j,k,n)
  !==================================================================!
    
    use w90_parameters, only :num_kpts, mp_grid

    implicit none

    integer :: i, j, k, N1, N2, N3
    integer, intent(out) :: n 

    if ((i .le. 0) .or. (i .ge. N1+1)) call io_error('In K_to_ijk(): input i not within bounds')
    if ((j .le. 0) .or. (j .ge. N2+1)) call io_error('In K_to_ijk(): input j not within bounds')
    if ((k .le. 0) .or. (k .ge. N3+1)) call io_error('In K_to_ijk(): input k not within bounds')
    
    N1 = mp_grid(1)
    N2 = mp_grid(2)
    N3 = mp_grid(3) 
    
    ! Assuming that n = (i-1)*(N2*N3) + (j-1)*N3 + k
    ! TODO check that they correspond to the correct k-points
    
    n = (i-1)*N2*N3 + (j-1)*N3 + k

    if (N1*N2*N3 .ne. num_kpts) then
        call io_error('Size of k-point mesh does not fit with total size')
    endif

    return

  end subroutine ijk_to_K




end module w90_autoproj
  
