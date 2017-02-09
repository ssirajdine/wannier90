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
    
    do n=1,num_kpts
       do m=1,num_wann
          u_matrix(m,m,n)=cmplx_1
       end do
    end do
    
  end subroutine autoproj_calc_u_matrix_bloch_phases

  !==================================================================!
  subroutine autoproj_calc_u_matrix_obstruction()
  !==================================================================!

    call io_error('Automatic projections using the obstruction matrix method not implemented yet')
    
  end subroutine autoproj_calc_u_matrix_obstruction

  
end module w90_sitesym
  
