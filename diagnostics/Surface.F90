!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module surface_diagnostics

  use fldebug
  use futils, only: int2str
  use global_parameters, only : OPTION_PATH_LEN
  use spud
  use fields
  use state_module
  use field_options
  use diagnostic_source_fields
  use bedload_diagnostics, only: surface_horizontal_divergence, surface_elevation_smoothing
  use surface_integrals
  
  implicit none
  
  private
  
  public :: calculate_grad_normal, calculate_surface_horizontal_divergence, calculate_surface_elevation_smoothing
  
contains
  
  subroutine calculate_grad_normal(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    
    type(scalar_field), pointer :: source_field
    type(vector_field), pointer :: positions
    
    character(len = OPTION_PATH_LEN) :: base_path
    integer, dimension(2) :: nsurface_ids
    integer, dimension(:), allocatable :: surface_ids
        
    source_field => scalar_source_field(state, s_field)
    positions => extract_vector_field(state, "Coordinate")
    
    base_path = trim(complete_field_path(s_field%option_path)) // "/algorithm"    
    if(have_option(trim(base_path) // "/surface_ids")) then
      nsurface_ids = option_shape(trim(base_path) // "/surface_ids")
      assert(nsurface_ids(1) >= 0)
      allocate(surface_ids(nsurface_ids(1)))
      call get_option(trim(base_path) // "/surface_ids", surface_ids)
      
      call surface_gradient_normal(source_field, positions, s_field, surface_ids = surface_ids)
      
      deallocate(surface_ids)
    else
      call surface_gradient_normal(source_field, positions, s_field)
    end if
    
  end subroutine calculate_grad_normal

subroutine calculate_surface_horizontal_divergence(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    
    type(vector_field), pointer :: source_field
    type(vector_field), pointer :: positions
    
    character(len = OPTION_PATH_LEN) :: base_path
    integer, dimension(2) :: nsurface_ids, bc_nsurface_ids
    integer, dimension(:), allocatable :: surface_ids, bc_surface_ids, faceglobalnodes
    integer, dimension(surface_element_count(s_field)) :: bc_surface_elements
    integer :: i, idx, i_ele, i_bc, i_bc_surf, i_surf, nbcs, snloc, globnod

    real :: beta, smoothing_length
        
    ewrite(1,*) "JN - In calculate_surface_horizontal_divergence"

    source_field => vector_source_field(state, s_field)
    positions => extract_vector_field(state, "Coordinate")
    
    base_path = trim(complete_field_path(s_field%option_path)) // "/algorithm"    

    ! get surface ids
    nsurface_ids = option_shape(trim(base_path) // "/surface_ids")
    assert(nsurface_ids(1) >= 0)
    allocate(surface_ids(nsurface_ids(1)))
    call get_option(trim(base_path) // "/surface_ids", surface_ids)
    call get_option(trim(base_path) // "/beta", beta, default=1.0)
    call get_option(trim(base_path) // "/smoothing_length", smoothing_length, default=0.0)

    call surface_horizontal_divergence(source_field, positions, s_field,&
    beta, smoothing_length, surface_ids = surface_ids,&
    option_path=base_path)

    ! set boundary id's to zero for surface_horizontal_div
    ! get boundary condition path and number of boundary conditions
    !nbcs=option_count(trim(base_path)//'/boundary_conditions')
    !if (nbcs >= 1) then

    if (option_count(trim(base_path)//'/set_to_zero') == 1) then

    ! Loop over boundary conditions for field
    !boundary_conditions: do i_bc=0, nbcs-1

        ! get boundary conditions surface ids
        !bc_nsurface_ids=option_shape(trim(base_path)//'/boundary_conditions['//int2str(i_bc)//']/surface_ids')
        bc_nsurface_ids=option_shape(trim(base_path)//'/set_to_zero[0]/surface_ids')
        allocate(bc_surface_ids(bc_nsurface_ids(1)))
        !call get_option(trim(base_path)//'/boundary_conditions['//int2str(i_bc)//']/surface_ids', bc_surface_ids)
        call get_option(trim(base_path)//'/set_to_zero[0]/surface_ids', bc_surface_ids)

        !bc_surface_id: do i_bc_surf=1, bc_nsurface_ids(1)
            !surface_id: do i_surf=1, nsurface_ids(1)

                ! get the relevant surface elements
                idx=0
                do i = 1, surface_element_count(s_field)
                    if (any(bc_surface_ids == surface_element_id(s_field, i))) then
                    idx=idx+1
                    bc_surface_elements(idx)=i
                    end if
                end do

                do i_ele=1, idx
                    snloc = face_loc(positions, 1)

                    allocate(faceglobalnodes(1:snloc))
                    faceglobalnodes = face_global_nodes(positions, bc_surface_elements(i_ele))

                    do i = 1, snloc
                        globnod = faceglobalnodes(i)
                        call set(s_field,globnod,0.0)
                    end do
                    deallocate(faceglobalnodes)

                end do

                !end do surface_id
                !end do bc_surface_id
                !end do boundary_conditions

                deallocate(bc_surface_ids)

    end if

    deallocate(surface_ids)
    
    ewrite_minmax(s_field)
    !!ewrite(1,*) 'JN - GOT HERE'

  end subroutine calculate_surface_horizontal_divergence
  
  subroutine calculate_surface_elevation_smoothing(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field

    type(scalar_field), pointer :: source_field
    type(vector_field), pointer :: positions

    character(len = OPTION_PATH_LEN) :: base_path
    integer, dimension(2) :: nsurface_ids
    integer, dimension(:), allocatable :: surface_ids
    real :: smoothing_length, theta
    logical :: lump_mass

    ewrite(1,*) "JN - In calculate_surface_elevation_smoothing"

    source_field => scalar_source_field(state, s_field)
    positions => extract_vector_field(state, "Coordinate")

    base_path = trim(complete_field_path(s_field%option_path)) // "/algorithm"

    nsurface_ids = option_shape(trim(base_path) // "/surface_ids")
    assert(nsurface_ids(1) >= 0)
    allocate(surface_ids(nsurface_ids(1)))
    call get_option(trim(base_path) // "/surface_ids", surface_ids)
    call get_option(trim(base_path) // "/smoothing_length", smoothing_length, default=0.0)
    lump_mass = have_option(trim(base_path) // "/lump_mass")
    call get_option(trim(base_path) // "/theta", theta, default=theta)

    call surface_elevation_smoothing(source_field, positions, s_field,&
         smoothing_length, surface_ids = surface_ids,&
         option_path=base_path, lump_mass = lump_mass, theta = theta)

    deallocate(surface_ids)

    ewrite_minmax(s_field)
    !!ewrite(1,*) 'JN - GOT HERE'

  end subroutine calculate_surface_elevation_smoothing

end module surface_diagnostics
