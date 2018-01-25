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

module bedload_diagnostics

  use fldebug
  use global_parameters, only : OPTION_PATH_LEN, FIELD_NAME_LEN, timestep, dt
  use vector_tools
  use quadrature
  use elements
  use spud
  use fetools
  use fields
  use transform_elements
  use state_module
  use field_options
  use sediment, only: get_n_sediment_fields, get_sediment_item
  use sparsity_patterns, only: make_sparsity
  use sparse_tools, only: csr_matrix, csr_sparsity
  use solvers, only: petsc_solve
  use boundary_conditions
  use halos, only: generate_surface_mesh_halos
  use field_derivatives
  use sparse_matrices_fields
  use state_fields_module
  
  implicit none 
  
  private

  public :: calculate_bedload_flux, surface_horizontal_divergence, surface_elevation_smoothing
  
contains 
 
    subroutine calculate_bedload_flux(state, v_field)

          type(state_type), intent(inout)           :: state
          type(mesh_type)                           :: surface_mesh
          type(scalar_field), pointer               :: sediment_field, evisc
          type(tensor_field), pointer               :: tvisc
          type(vector_field), pointer               :: x, bss
          type(vector_field), intent(inout)         :: v_field
          type(vector_field)                        :: bedload_flux_surface

          integer, dimension(surface_element_count(v_field))    :: surface_elements

          integer                                   :: i, j, k, idx, n_sediment_fields, i_field,&
                                                       & stat, i_ele, snloc, sndim, sele, globnod, globnod_prev, globnod_next, globnod_visc

          integer, dimension(:), pointer            :: surface_node_list, s_nodes, ele
          integer, dimension(2)                     :: surface_id_count
          integer, dimension(:), allocatable        :: surface_ids, faceglobalnodes, faceglobalnodes_prev, faceglobalnodes_next, faceglobalnodes_visc
          real, dimension(:), allocatable           :: t_crit, q_star, q_aval, tvisc_val
          real                                      :: density, R, d, g, zevisc, d_star
          real                                      :: b, repose, q_aval_x, q_aval_y, acc_aval
          real, parameter                           :: PI = 4.D0*DATAN(1.D0)
          character(len = OPTION_PATH_LEN)          :: base_path

          ! Declarations for future development
          !type(element_type), pointer               :: s_shape
          !integer, dimension(:), pointer            :: surface_element_list, surface_ids, bedload_flux_surface_nodes
          !type(scalar_field)                        :: masslump
          !real, dimension(:,:)                      :: invmass
          !real, dimension(8)                        :: detwei
          !real, dimension(ele_loc(bedload_flux_surface, i_ele), ele_loc(bedload_flux_surface, i_ele)) :: invmass
          !real, dimension(ele_ngi(bedload_flux_surface, i_ele)) :: detwei

          ewrite(1,*) "JN - In calculate_sediment_bedload_flux"

          ! extract turbulent viscocity field
          !evisc => extract_scalar_field(state, "ScalarEddyViscosity", stat)
          !tvisc => extract_tensor_field(state, "Viscosity", stat)
          tvisc => extract_tensor_field(state, "BackgroundViscosity", stat)

          ! extract turbulent viscocity field in zero dump
          !if (timestep == 0) then
          !  call get_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/tensor_field::BackgroundViscosity/prescribed/value::WholeMesh/isotropic/constant", zevisc)
          !end if
          !ewrite(1,*) 'JN - EDDY VISC IN NODE 10:', zevisc
          !ewrite(1,*) 'JN - EDDY VISC IN NODE 10:', node_val(evisc, 10)

          ! extract coordinate field
          x => extract_vector_field(state, "Coordinate")

          ! extract bed shear stress field
          bss => extract_vector_field(state, "BedShearStress", stat)
          if (stat /= 0) then
            FLExit("A bed shear stress must be specified to calculate bedload flux")
          end if

          n_sediment_fields = get_n_sediment_fields()
          if (n_sediment_fields == 0) return

          ! loop over sediment fields
          sediment_fields: do i_field = 1, n_sediment_fields

            ! get sediment field from state
            call get_sediment_item(state, i_field, sediment_field)

            ! get R
            call get_sediment_item(state, i_field, 'submerged_specific_gravity', R)

            ! get d
            call get_sediment_item(state, i_field, 'diameter', d)

            ! get gravity
            call get_option("/physical_parameters/gravity/magnitude", g)

            ! get constant density from bed shear stress
            call get_option("/material_phase[0]/vector_field::BedShearStress/diagnostic/density", density)

            ! define base path
            base_path = trim(complete_field_path(v_field%option_path))//&
            &"/algorithm"
            !ewrite(1,*) 'JN - PATH:', base_path

            ! obtain sediment bedload surface ids
            surface_id_count=option_shape(trim(base_path)//&
            &"/surface_ids")

            allocate(surface_ids(surface_id_count(1)))

            call get_option(trim(base_path) // "/angle_of_repose", repose, default=0.0)

            call get_option(trim(base_path) // "/avalanche_factor", acc_aval, default=1.0)

            call get_option(trim(base_path)//"/surface_ids", surface_ids)
            !ewrite(1,*) 'JN - SURFACE ID:', surface_ids

            ! get the relevant surface elements
            idx=1
            do i = 1, surface_element_count(v_field)
                if (any(surface_ids == surface_element_id(v_field, i))) then
                surface_elements(idx)=i
                idx=idx+1
                end if
            end do

            ! generate surface_mesh for calculation of bedload flux and create surface field
            call create_surface_mesh(surface_mesh, surface_node_list, &
                & mesh=v_field%mesh, surface_elements=surface_elements(:idx-1), name='SurfaceMesh')

            call allocate(bedload_flux_surface, v_field%dim, surface_mesh, name='BedloadRateSurface')

            call zero(bedload_flux_surface)

            call zero(v_field)

            !ewrite(1,*) 'JN - SURF ELEMENTS:', element_count(bedload_flux_surface)
            !ewrite(1,*) 'JN - SURF NODES:', surface_node_list

            ! calculate bedload flux
            ! loop through elements in surface field
            elements: do i_ele=1, element_count(bedload_flux_surface)
            !elements: do i_ele=1, element_count(v_field)
            !elements: do i_ele=1, size(surface_elements)

                ! check if element is on prescribed surface
                !if (.not.any(surface_element_id(v_field, i_ele) == surface_ids)) then
                     !cycle elements
                !end if

                !s_nodes => ele_nodes(bedload_flux_surface, i_ele)
                !s_nodes => ele_nodes(v_field, i_ele)
                !ele => face_ele(bss, i_ele)

                sndim = mesh_dim(bss)
                !ewrite(1,*) 'JN - BSS DIMENSION:', sndim

                snloc = face_loc(bss, 1)
                !ewrite(1,*) 'JN - BSS NODES IN FACE ONE:', snloc

                !!ewrite(1,*) '******'
                !!ewrite(1,*) 'JN - BSS ELEMENT:', i_ele
                !!ewrite(1,*) 'JN - V_FIELD SURF ELE:', surface_elements((i_ele))
                !ewrite(1,*) 'JN - BSS ELEMENTS:', surface_element_count(bss)

                allocate(tvisc_val(1:sndim))
                allocate(t_crit(1:sndim))
                allocate(q_star(1:sndim))
                allocate(q_aval(1:sndim))
                !allocate(faceglobalnodes_visc(1:snloc))
                allocate(faceglobalnodes(1:snloc))
                !allocate(faceglobalnodes_next(1:snloc))
                !allocate(faceglobalnodes_prev(1:snloc))

                ! Reading total viscosity in nodes
                where (face_val(tvisc, 1, 1, surface_elements(i_ele)) > 0.0)
                    tvisc_val = face_val(tvisc, 1, 1, surface_elements(i_ele))
                elsewhere
                    tvisc_val = 0.0
                end where

                !faceglobalnodes_visc = face_global_nodes(evisc, surface_elements(i_ele))
                !ewrite(1,*) 'JN - EDDY VISC GLOBAL NODES IN FACE:', faceglobalnodes_visc

                !do sele=1, surface_element_count(bss)
                !faceglobalnodes = face_global_nodes(bss, sele)
                faceglobalnodes = face_global_nodes(bss, surface_elements(i_ele))
                !!ewrite(1,*) 'JN - BSS GLOBAL NODES IN FACE:', faceglobalnodes

                !if (i_ele < element_count(bedload_flux_surface) - 1) then
                !faceglobalnodes_next = face_global_nodes(bss, surface_elements(i_ele+1))
                !ewrite(1,*) 'JN - BSS GLOBAL NODES IN NEXT FACE:', faceglobalnodes_next
                !end if

                !if (i_ele > 1) then
                !faceglobalnodes_prev = face_global_nodes(bss, surface_elements(i_ele-1))
                !ewrite(1,*) 'JN - BSS GLOBAL NODES IN PREV FACE:', faceglobalnodes_prev
                !end if

                !!ewrite(1,*) '******'

                do i = 1, snloc

                    !globnod_visc = faceglobalnodes_visc(i)
                    !ewrite(1,*) 'JN - EDDY VISC GLOBAL NODE:', globnod_visc

                    globnod = faceglobalnodes(i)
                    !!ewrite(1,*) 'JN - BSS GLOBAL NODE:', globnod
                    !globnod_next = faceglobalnodes_next(i)
                    !ewrite(1,*) 'JN - BSS GLOBAL NODE NEXT FACE:', globnod_next
                    !globnod_prev = faceglobalnodes_prev(i)
                    !ewrite(1,*) 'JN - BSS GLOBAL NODE PREV FACE:', globnod_prev

                    ! compute dimensionless particle diameter
                    !!ewrite(1,*) 'JN - EDDY VISC:', node_val(evisc, globnod_visc)
                    !d_star = 9.10345 !TEMPORARY CHANGE FOR MATT for 2D case
                    d_star = 6.57471 !TEMPORARY CHANGE FOR MATT for 3D case
                    !d_star = d * ((R*g/(zevisc)**2)**1/3)
                    !d_star = d * ((R*g/(node_val(evisc, globnod_visc))**2)**1/3)
                    !d_star = d * ((R*g/(tvisc_val(i)**2))**1/3)

                    !!ewrite(1,*) 'JN - DIMENSIONLESS PARTICLE DIAMETER:', d_star

                    if (have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon")) then
                        !if (timestep == 0) then
                        !    t_crit = 0.00
                        !else
                            t_crit(1:sndim) = (0.30 / (1 + 1.2*d_star)) + 0.055*(1-exp(-0.020*d_star)) !Soulsby (1997)
                        !end if
                    else
                        t_crit(1:sndim) = 0.05
                        !t_crit = 0.00
                    end if

                    !!ewrite(1,*) 'JN - T_CRIT:', t_crit

                    if ((norm2(node_val(bss, globnod) / (density * R * g * d)) - norm2(t_crit)) > 0) then

                    ! calculate Meyer-Peter and Muller model
                    if (have_option(trim(base_path)//"/transport_model/meyer_peter_muller")) then

                    !q_star = 8.0 * ((node_val(bss, globnod) / (density * R * g * d)) - t_crit) * (sqrt(norm2(node_val(bss, globnod) / (density * R * g * d)) - norm2(t_crit)))

                    q_star = 8.0 * (norm2(node_val(bss, globnod) / (density * R * g * d)) - norm2(t_crit)) ** (3./2.) * (node_val(bss, globnod)/norm2(node_val(bss, globnod)))

                    ! calculate Nielsen model
                    elseif (have_option(trim(base_path)//"/transport_model/nielsen")) then

                    q_star = 12.0 * ((node_val(bss, globnod) / (density * R * g * d)) - t_crit) * (sqrt(norm2(node_val(bss, globnod) / (density * R * g * d))))

                    ! calculate Van Rijn model
                    elseif (have_option(trim(base_path)//"/transport_model/van_rijn")) then

                    q_star = (0.053 / (d_star ** 0.3)) * (((node_val(bss, globnod) / t_crit) - 1) ** 2.1)

                    ! calculate Engelund and Fresoe model
                    elseif (have_option(trim(base_path)//"/transport_model/engelund_fredsoe")) then

                    q_star = 18.74 * ((node_val(bss, globnod) / (density * R * g * d)) - t_crit) * ((sqrt(norm2(node_val(bss, globnod) / (density * R * g * d)))) - (0.7 * sqrt(norm2(t_crit))))

                    end if

                    else

                    q_star = 0.0

                    end if

                    !ewrite(1,*) 'JN - QSTAR:', q_star

                    call set(v_field, globnod, q_star * sqrt(R * g * (d ** 3)))
                    !!ewrite(1,*) 'JN - BEDLOAD RATE:', node_val(v_field, globnod)
                    !!ewrite(1,*) '******'

                    ! AVALANCHE MODEL

                    if (sndim == 2) then

                        call avalanche_2d(state, i_ele, repose, x, v_field, surface_elements, b, q_aval_x, q_aval_y)

                        !if (tan(abs(b)) > tan(repose*PI/180)) then
                        !ewrite(1,*) 'JN - AVALANCHE_X IN BED FLUX:', q_aval_x
                        !ewrite(1,*) 'JN - AVALANCHE_Y IN BED FLUX:', q_aval_y
                        !end if

                    elseif (sndim == 3) then

                        call avalanche_2d(state, i_ele, repose, x, v_field, surface_elements, b, q_aval_x, q_aval_y)

                    end if

                    if ((repose /= 0.0) .and. (dt /= 0.0)) then
                        if (tan(abs(b)) > tan(repose*PI/180)) then
                            q_aval(1:sndim) = (/ q_aval_x, q_aval_y /)
                            !!ewrite(1,*) 'JN - Q_AVAL:', q_aval

                            !q_aval(1:sndim) = (/ abs(q_aval_x), abs(q_aval_y) /)
                            !ewrite(1,*) 'JN - Q_AVAL ABS:', q_aval

                            !!ewrite(1,*) 'JN - BED FLUX PLUS Q_AVAL:', node_val(v_field, globnod) + q_aval/dt
                            !!ewrite(1,*) 'JN - BED FLUX MINUS Q_AVAL:', node_val(v_field, globnod) - q_aval/dt

                            if (b < 0.0) then

                                call addto(v_field, 1, globnod, acc_aval*q_aval_x/dt)
                                call addto(v_field, 2, globnod, acc_aval*q_aval_y/dt)

                                !ewrite(1,*) 'JN - SIGNS:', sign(1.0, node_val(v_field,1,globnod))
                                !ewrite(1,*) 'JN - SIGNS:', sign(1.0, node_val(v_field,2,globnod))
                                !ewrite(1,*) 'JN - SIGNS:', sign(1.0, q_aval(1))
                                !ewrite(1,*) 'JN - SIGNS:', sign(1.0, q_aval(2))

                                !do k = 1, sndim
                                !if (sign(1.0, node_val(v_field, k, globnod)) > 0.0) then
                                    !call addto(v_field, k, globnod, q_aval(k))
                                !else if (sign(1.0, node_val(v_field, k, globnod)) < 0.0) then
                                    !call addto(v_field, k, globnod, -q_aval(k))
                                !end if
                                !end do

                                !if (i_ele < element_count(bedload_flux_surface) - 1) then
                                !call addto(v_field, globnod_next, q_aval)
                                !elseif (i_ele == element_count(bedload_flux_surface)) then
                                !call addto(v_field, globnod, q_aval)
                                !endif

                            elseif (b > 0.0) then

                                call addto(v_field, 1, globnod, -acc_aval*q_aval_x/dt)
                                call addto(v_field, 2, globnod, -acc_aval*q_aval_y/dt)

                                !do k = 1, sndim
                                !if (sign(1.0, node_val(v_field, k, globnod)) > 0.0) then
                                    !call addto(v_field, k, globnod, q_aval(k))
                                !else if (sign(1.0, node_val(v_field, k, globnod)) < 0.0) then
                                    !call addto(v_field, k, globnod, -q_aval(k))
                                !end if
                                !end do

                                !if (i_ele > 1) then
                                !call addto(v_field, globnod_prev, -1 * q_aval)
                                !elseif (i_ele == 1) then
                                !call addto(v_field, globnod, q_aval)
                                !endif

                            end if
                        end if
                    end if

                    !!ewrite(1,*) 'JN - BEDLOAD ADDITIONAL RATE:', node_val(v_field, globnod)
                    !!ewrite(1,*) '******'

                    !if (i_ele < element_count(bedload_flux_surface) - 1) then
                    !ewrite(1,*) 'JN - BEDLOAD ADDITIONAL RATE NEXT:', node_val(v_field, globnod_next)
                    !end if

                    !if (i_ele > 1) then
                    !ewrite(1,*) 'JN - BEDLOAD ADDITIONAL RATE PREV:', node_val(v_field, globnod_prev)
                    !end if

                    ! For future FE development
                    !call set(bedload_flux_surface, s_nodes(i), node_val(v_field, globnod))
                end do
                !end do

                !do j=1, size(surface_node_list)
                    !ewrite(1,*) 'JN - NODE ID:', surface_node_list(j)
                    !ewrite(1,*) 'JN - NODE VAL ID:', node_val(v_field, j)
                !end do

                !do j=1,size(surface_node_list)
                    !call set(bedload_flux_surface,surface_node_list(j),node_val(v_field,j))
                !end do

                deallocate(q_aval)
                deallocate(q_star)
                deallocate(t_crit)
                deallocate(tvisc_val)
                !deallocate(faceglobalnodes_visc)
                deallocate(faceglobalnodes)
                !deallocate(faceglobalnodes_prev)
                !deallocate(faceglobalnodes_next)

            end do elements

          ewrite_minmax(v_field)
          !!ewrite(1,*) 'JN - GOT HERE'

          call deallocate(bedload_flux_surface)
          call deallocate(surface_mesh)
          deallocate(surface_ids)

          end do sediment_fields
    
    end subroutine calculate_bedload_flux

    subroutine avalanche_2d(state, i_ele, repose, X, v_field, surface_elements, b, q_aval_x, q_aval_y)

        type(state_type), intent(in)                       :: state
        type(vector_field), pointer                        :: X
        type(vector_field), intent(in)                     :: v_field
        integer                                            :: i, j, i_ele, sndim, snloc, globnod
        integer, dimension(:), allocatable                 :: faceglobalnodes
        integer, dimension(surface_element_count(v_field)) :: surface_elements
        real                                               :: dx, dy, nl, nf, x1, x2, y1, y2, length, tx, ty,&
                                                              & q_aval, repose
        real, intent(out)                                  :: b, q_aval_x, q_aval_y
        real, parameter                                    :: PI = 4.D0*DATAN(1.D0)

        sndim = mesh_dim(X)
        !ewrite(1,*) 'JN - X DIMENSION:', sndim

        snloc = face_loc(X, 1)
        !ewrite(1,*) 'JN - X NODES IN FACE ONE:', snloc

        !ewrite(1,*) 'JN - BSS ELEMENT:', i_ele
        !ewrite(1,*) 'JN - V_FIELD SURF ELE:', surface_elements((i_ele))

        allocate(faceglobalnodes(1:snloc))
        faceglobalnodes = face_global_nodes(X, surface_elements(i_ele))
        !ewrite(1,*) 'JN - X GLOBAL NODES IN FACE:', faceglobalnodes
        !ewrite(1,*) 'JN - MINVAL:', minval(faceglobalnodes)
        !ewrite(1,*) 'JN - MAXVAL:', maxval(faceglobalnodes)
        !ewrite(1,*) 'JN - MINLOC:', minloc(faceglobalnodes)
        !ewrite(1,*) 'JN - MAXLOC:', maxloc(faceglobalnodes)
        !ewrite(1,*) 'JN - MIN:', faceglobalnodes(minloc(faceglobalnodes))
        !ewrite(1,*) 'JN - MAX:', faceglobalnodes(maxloc(faceglobalnodes))

!        do j = 1, snloc
!            do i = 1, sndim
!                globnod = faceglobalnodes(j)
!                ewrite(1,*) 'JN - COORDINATES:', node_val(X, i, globnod)
!            end do
!        end do

        !nl = faceglobalnodes(maxloc(faceglobalnodes))
        !ewrite(1,*) 'JN - LAST NODE:', nl
        !nf = faceglobalnodes(minloc(faceglobalnodes))
        !ewrite(1,*) 'JN - FIRST NODE:', nf

        !x2 = node_val(X, 1, maxval(faceglobalnodes))
        !x1 = node_val(X, 1, minval(faceglobalnodes))
        x2 = node_val(X, 1, faceglobalnodes(2))
        x1 = node_val(X, 1, faceglobalnodes(1))
        !x2 = node_val(X, 1, nl)
        !x1 = node_val(X, 1, nf)
        !!ewrite(1,*) 'JN - X2:', x2
        !!ewrite(1,*) 'JN - X1:', x1

        if (x2 > x1) then

            y2 = node_val(X, 2, faceglobalnodes(2))
            y1 = node_val(X, 2, faceglobalnodes(1))

        else

            x2 = x2 + x1
            x1 = x2 - x1
            x2 = x2 - x1

            y2 = node_val(X, 2, faceglobalnodes(1))
            y1 = node_val(X, 2, faceglobalnodes(2))

        end if

        !!ewrite(1,*) 'JN - NEW X2:', x2
        !!ewrite(1,*) 'JN - NEW X1:', x1

        !y2 = node_val(X, 2, maxval(faceglobalnodes))
        !y1 = node_val(X, 2, minval(faceglobalnodes))
        !y2 = node_val(X, 2, faceglobalnodes(2))
        !y1 = node_val(X, 2, faceglobalnodes(1))
        !y2 = node_val(X, 2, nl)
        !y1 = node_val(X, 2, nf)
        !!ewrite(1,*) 'JN - Y2:', y2
        !!ewrite(1,*) 'JN - Y1:', y1

        dx = x2 - x1
        !!ewrite(1,*) 'JN - DX:', dx

        dy = y2 - y1
        !!ewrite(1,*) 'JN - DY:', dy

        length = sqrt(dx**2 + dy**2)
        !!ewrite(1,*) 'JN - LENGTH:', length

        tx = dx / length
        !!ewrite(1,*) 'JN - TX:', tx

        ty = dy / length
        !!ewrite(1,*) 'JN - TY:', ty

        !b = atan2(ty,tx)
        b = atan(ty/tx)
        !!ewrite(1,*) 'JN - BOTTOM SLOPE RAD:', b
        !!ewrite(1,*) 'JN - BOTTOM SLOPE DEG:', b * 180/PI

        if (tan(abs(b)) > tan(repose*PI/180)) then

            q_aval = 0.5 * length**2 * (1 - 0.0) * ((tan(abs(b)) - tan(repose*PI/180))/cos(abs(b)))
            !!ewrite(1,*) 'JN - AVALANCHE:', q_aval

            q_aval_x = q_aval * tx
            !return
            !!ewrite(1,*) 'JN - AVALANCHE_X:', q_aval_x

            q_aval_y = q_aval * ty
            !return
            !!ewrite(1,*) 'JN - AVALANCHE_Y:', q_aval_y

        end if

        deallocate(faceglobalnodes)

    end subroutine avalanche_2d

    subroutine surface_horizontal_divergence(source, positions, output, beta, smoothing_length, surface_ids, option_path)
    !!< Return a field containing:
    !!< div_HS source
    !!< where div_hs is a divergence operator restricted to the surface and
    !!< of spatial dimension one degree lower than the full mesh.

        type(vector_field), intent(in) :: source
        type(vector_field), intent(in) :: positions
        type(scalar_field), intent(inout) :: output
        real, intent(in) :: beta, smoothing_length
        integer, dimension(:), intent(in) :: surface_ids
        character(len=*), optional, intent(in) :: option_path

        type(mesh_type) :: surface_mesh, source_surface_mesh, output_surface_mesh

        integer :: i, idx

        integer, dimension(surface_element_count(output)) :: surface_elements
        integer, dimension(:), pointer :: surface_nodes, source_surface_nodes, output_surface_nodes
        real :: face_area, face_integral

        type(vector_field) :: surface_source, surface_positions
        type(scalar_field) :: surface_output
        real, dimension(mesh_dim(positions))  :: val

        !interface
        !    subroutine generate_free_surface_halos(full_mesh, surface_mesh, surface_nodes)
        !        use fields
        !        type(mesh_type), intent(in):: full_mesh
        !        type(mesh_type), intent(inout):: surface_mesh
        !        integer, dimension(:), intent(in):: surface_nodes
        !    end subroutine
        !end interface

        call zero(output)

        !!! get the relevant surface. This wastes a bit of memory.
        idx=1
        do i = 1, surface_element_count(output)
            if (any(surface_ids == surface_element_id(positions, i))) then
                surface_elements(idx)=i
                idx=idx+1
            end if
        end do
        !!ewrite(1,*) 'JN - SURFACE ID:', surface_ids

        !!! make the surface meshes
        call create_surface_mesh(surface_mesh,  surface_nodes, &
         positions%mesh, surface_elements=surface_elements(:idx-1), &
         name='CoordinateSurfaceMesh')
        call create_surface_mesh(source_surface_mesh,  source_surface_nodes, &
         source%mesh, surface_elements=surface_elements(:idx-1), &
         name='SourceSurfaceMesh')
        call create_surface_mesh(output_surface_mesh,  output_surface_nodes, &
         output%mesh, surface_elements=surface_elements(:idx-1), &
         name='OutputSurfaceMesh')

        call generate_surface_mesh_halos(output%mesh, output_surface_mesh, output_surface_nodes)

        call add_faces(surface_mesh)
        call add_faces(source_surface_mesh, surface_mesh)
        call add_faces(output_surface_mesh, surface_mesh)

        call allocate(surface_positions,mesh_dim(surface_mesh),surface_mesh,&
         "Coordinates")
        call allocate(surface_source,mesh_dim(surface_mesh),source_surface_mesh,&
         "Source")
        call allocate(surface_output,output_surface_mesh,"Divergence")

        !!ewrite(1,*) 'JN - POSITIONS (FIELD) DIM:', mesh_dim(positions)
        !!ewrite(1,*) 'JN - POSITIONS (SURFACE) DIM:', mesh_dim(surface_mesh)
        !!ewrite(1,*) 'JN - OUTPUT FIELD DIM:', mesh_dim(output)
        !!ewrite(1,*) 'JN - OUTPUT (SURFACE) DIM:', mesh_dim(surface_output)
        !!ewrite(1,*) 'JN - ELEMENTS OF OUTPUT MESH:', element_count(output_surface_mesh)
        !!ewrite(1,*) 'JN - NODES OF OUTPUT MESH:', output_surface_nodes
        !!ewrite(1,*) '******'

        do i=1,size(surface_nodes)
            val=node_val(positions,surface_nodes(i))
            call set(surface_positions,i,val(:mesh_dim(surface_mesh)))
        end do

        do i=1,size(source_surface_nodes)
            val=node_val(source,source_surface_nodes(i))
            call set(surface_source,i,val(:mesh_dim(surface_mesh)))
        end do

        if (present(option_path)) then
            surface_output%option_path = option_path
        else
            surface_output%option_path = output%option_path
        end if

        !!! now do the low dimensional divergence operation

        call get_div(surface_source,surface_positions,surface_output,beta,smoothing_length)

        do i=1,size(output_surface_nodes)
            call set(output,output_surface_nodes(i),node_val(surface_output,i))
        end do

        !!! I *think* this surfices for parallel, due to the relatively locality of the divergence operator

        call halo_update(output)

        call deallocate(surface_positions)
        call deallocate(surface_source)
        call deallocate(surface_output)

        call deallocate(surface_mesh)
        call deallocate(source_surface_mesh)
        call deallocate(output_surface_mesh)

        deallocate(surface_nodes, source_surface_nodes, output_surface_nodes)

   end subroutine surface_horizontal_divergence

   subroutine get_div(source,positions,output,beta,smoothing_length)

        type(vector_field), intent(in) :: source, positions
        type(scalar_field), intent(inout) :: output
        real, intent(in) ::  beta, smoothing_length

        type(element_type), pointer :: shape
        integer, dimension(:), pointer :: nodes

        type(csr_matrix) :: M
        type(csr_sparsity) :: sparsity
        type(scalar_field) :: rhs
        !type(scalar_field) :: masslump
        integer, dimension(:), allocatable :: faceglobalnodes
        integer :: i, j, ele, face, dim, globnod

        real, dimension(ele_ngi(output,1)) :: detwei
        real, dimension(face_ngi(positions,1)) :: sdetwei, flux
        real, dimension(mesh_dim(positions), face_ngi(positions,1)) :: normal
        real, dimension(ele_loc(positions,1), ele_ngi(positions,1), mesh_dim(positions)) :: do_t

        sparsity = make_sparsity(output%mesh, output%mesh, "SparsityMat")

        call allocate(M, sparsity, name = "MassMat")
        call allocate(rhs, output%mesh, "RHS")
        !call allocate(masslump, output%mesh, 'MassLump')

        call zero(M)
        call zero(rhs)
        call zero(output)
        !call zero(masslump)

        !!ewrite(1,*) 'JN - SOURCE DIM:', mesh_dim(source)
        !!ewrite(1,*) 'JN - POSITIONS DIM:', mesh_dim(positions)
        !!ewrite(1,*) 'JN - OUTPUT DIM:', mesh_dim(output)
        !!ewrite(1,*) 'JN - POSITIONS ELE LOC:', ele_loc(positions,1)
        !!ewrite(1,*) 'JN - POSITIONS ELE GAUSS POINTS:', ele_ngi(positions,1)
        !!ewrite(1,*) 'JN - POSITIONS FACE LOC:', face_loc(positions,1)
        !!ewrite(1,*) 'JN - POSITIONS FACE GAUSS POINTS:', face_ngi(positions,1)

        do ele = 1, element_count(output)
        !!! element loop

        !!ewrite(1,*) 'JN - ELEMENT:', ele
        !!ewrite(1,*) 'JN - VAL IN ELEMENT:', ele_val_at_quad(source, ele)
        !!ewrite(1,*) 'JN - ELEMENT LOC:', ele_loc(output,ele)

            shape => ele_shape(output,ele)
            !ewrite(1,*) 'JN - SHAPES IN ELEMENT:', shape
            nodes => ele_nodes(output,ele)
            !!ewrite(1,*) 'JN - NODES IN ELEMENT:', nodes

            call transform_to_physical(positions, ele, shape=shape, dshape=do_t, detwei=detwei)
            !! build the mass matrix
            call addto(M,nodes,nodes,shape_shape(shape,shape,detwei)+((smoothing_length)**2)*dshape_dot_dshape(do_t,do_t,detwei))
            !ewrite(1,*) 'JN - MASS MATRIX:', node_val(M, ele)

            !! generate element rhs
            call addto(rhs, nodes, get_div_element(ele_loc(output,ele)))
            !!ewrite(1,*) 'JN - RHS VECTOR:', node_val(RHS, ele)

            !! build the lumped mass in the elements
            !call addto(masslump,nodes,shape_shape(shape,shape,detwei))
            !do i = 1,size(nodes)
                !call addto(masslump,nodes(i),sum(detwei)/ele_loc(output,ele))
            !end do
        end do

        do face = 1, face_count(output)
        !!! face loop

        !!ewrite(1,*) 'JN - FACE:', face
        !!ewrite(1,*) 'JN - VAL IN FACE:', face_val_at_quad(source, face)
        !!ewrite(1,*) 'JN - FACE LOC:', face_loc(output,face)

            shape => face_shape(positions,face)
            !ewrite(1,*) 'JN - SHAPES IN FACE:', shape

            !!ewrite(1,*) 'JN - NODES IN FACE:', face_global_nodes(source,face)

            call transform_facet_to_physical(positions, face, detwei_f=sdetwei, normal=normal)
            !! this could be cleverer
            !!ewrite(1,*) 'JN - NORMAL:', normal
            flux = sum(face_val_at_quad(source, face) * normal, dim=1)
            !! generate_element_rhs
            call addto(rhs, face_global_nodes(output,face), get_div_face(face_loc(output,face)))
            !!ewrite(1,*) 'JN - RHS VECTOR PLUS FACE CONTRIBUTION:', node_val(RHS, ele)

            !! build the lumped mass in the faces
            !call addto(masslump,face_global_nodes(source,face),shape_shape(shape,shape,detwei))
            !call addto(masslump,face_global_nodes(source,face),sum(shape_shape(shape,shape,detwei),1))
        end do

        call petsc_solve(output,M,rhs)
        !call petsc_solve(output,lumped(M),rhs)
        !call petsc_solve(output,masslump,rhs)

        !!check continuity in output mesh and divide by the inverse global lumped mass
        !if(continuity(output)>=0) then
        !    where (masslump%val/=0.0)
        !        masslump%val=1./masslump%val
        !    end where
        !    call scale(rhs, masslump)
        !    call deallocate(masslump)
        !end if

        !call set(output,rhs)

        !allocate(faceglobalnodes(1:size(nodes)))
        !do j = 1, element_count(output)
        !    faceglobalnodes = face_global_nodes(output,j)
            !!do i = 1, size(nodes)
            !!ewrite(1,*) 'JN - BEDLOAD DIV:', node_val(output, nodes(i))
        !        globnod = faceglobalnodes(i)
        !        ewrite(1,*) 'JN - BEDLOAD DIV:', node_val(output, globnod)
            !!end do
        !end do
        !deallocate(faceglobalnodes)

        call deallocate(M)
        call deallocate(sparsity)
        call deallocate(rhs)

        contains

        function get_div_element(nloc) result(local_rhs)

            integer ::nloc
            real, dimension(nloc) :: local_rhs
            real, dimension(mesh_dim(output), mesh_dim(output), ele_ngi(output, ele)) :: source_grad

            source_grad = ele_grad_at_quad(source, ele, do_t)
            local_rhs = -beta * dshape_dot_vector_rhs(do_t, ele_val_at_quad(source,ele), detwei)
            do dim=1, mesh_dim(output)
                local_rhs = local_rhs + (1.0-beta) * shape_rhs(shape, source_grad(dim, dim, :) * detwei)
            end do

        end function get_div_element

        function get_div_face(nloc) result(local_rhs)

            integer:: nloc
            real, dimension(nloc) :: local_rhs

            local_rhs = beta * shape_rhs(shape, flux * sdetwei)

        end function get_div_face

    end subroutine get_div

    subroutine surface_elevation_smoothing(source, positions, output, smoothing_length, surface_ids, option_path, lump_mass, theta)
    !!< Return a field containing:
    !!< w = zn+1 - zn / dt
    !!< where w is a grid velocity of the smoothed coordinate restricted to the surface and
    !!< of spatial dimension one degree lower than the full mesh.

        type(scalar_field), intent(in) :: source
        type(vector_field), intent(in) :: positions
        type(scalar_field), intent(inout) :: output
        real, intent(in) :: smoothing_length
        integer, dimension(:), intent(in) :: surface_ids
        character(len=*), optional, intent(in) :: option_path
        logical, intent(in) :: lump_mass
        real, intent(in) :: theta

        type(mesh_type) :: surface_mesh, source_surface_mesh, output_surface_mesh

        integer :: i, idx

        integer, dimension(surface_element_count(output)) :: surface_elements
        integer, dimension(:), pointer :: surface_nodes, source_surface_nodes, output_surface_nodes
        real :: face_area, face_integral

        type(vector_field) :: surface_positions
        type(scalar_field) :: surface_source, surface_elevation, surface_output
        real, dimension(mesh_dim(positions)) :: val
        real :: sval

        call zero(output)

        !!! get the relevant surface. This wastes a bit of memory.
        idx=1
        do i = 1, surface_element_count(output)
            if (any(surface_ids == surface_element_id(positions, i))) then
                surface_elements(idx)=i
                idx=idx+1
            end if
        end do
        !!ewrite(1,*) 'JN - SURFACE ID:', surface_ids

        !!! make the surface meshes
        call create_surface_mesh(surface_mesh,  surface_nodes, &
         positions%mesh, surface_elements=surface_elements(:idx-1), &
         name='CoordinateSurfaceMesh')
        call create_surface_mesh(source_surface_mesh,  source_surface_nodes, &
         source%mesh, surface_elements=surface_elements(:idx-1), &
         name='SourceSurfaceMesh')
        call create_surface_mesh(output_surface_mesh,  output_surface_nodes, &
         output%mesh, surface_elements=surface_elements(:idx-1), &
         name='OutputSurfaceMesh')

        call generate_surface_mesh_halos(output%mesh, output_surface_mesh, output_surface_nodes)

        call add_faces(surface_mesh)
        call add_faces(source_surface_mesh, surface_mesh)
        call add_faces(output_surface_mesh, surface_mesh)

        call allocate(surface_positions,mesh_dim(surface_mesh),surface_mesh,"Coordinates")
        call allocate(surface_elevation,surface_mesh,"Elevation")
        call allocate(surface_source,source_surface_mesh,"Source")
        call allocate(surface_output,output_surface_mesh,"GridVelocity")

        !!ewrite(1,*) 'JN - POSITIONS (FIELD) DIM:', mesh_dim(positions)
        !!ewrite(1,*) 'JN - POSITIONS (SURFACE) DIM:', mesh_dim(surface_mesh)
        !!ewrite(1,*) 'JN - OUTPUT FIELD DIM:', mesh_dim(output)
        !!ewrite(1,*) 'JN - OUTPUT (SURFACE) DIM:', mesh_dim(surface_output)
        !!ewrite(1,*) 'JN - ELEMENTS OF OUTPUT MESH:', element_count(output_surface_mesh)
        !!ewrite(1,*) 'JN - NODES OF OUTPUT MESH:', output_surface_nodes
        !!ewrite(1,*) '******'

        do i=1,size(surface_nodes)
            val=node_val(positions,surface_nodes(i))
            call set(surface_positions,i,val(:mesh_dim(surface_mesh)))
        end do

        do i=1,size(surface_nodes)
            val=node_val(positions,surface_nodes(i))
            call set(surface_elevation,i,val(mesh_dim(source)))
        end do

        do i=1,size(source_surface_nodes)
            sval=node_val(source,source_surface_nodes(i))
            call set(surface_source,i,sval)
        end do

        if (present(option_path)) then
            surface_output%option_path = option_path
        else
            surface_output%option_path = output%option_path
        end if

        call get_coordinate_smooth(surface_source,surface_positions,surface_elevation,surface_output,smoothing_length, lump_mass, theta)

        do i=1,size(output_surface_nodes)
            call set(output,output_surface_nodes(i),node_val(surface_output,i))
        end do

        call halo_update(output)

        call deallocate(surface_positions)
        call deallocate(surface_elevation)
        call deallocate(surface_source)
        call deallocate(surface_output)

        call deallocate(surface_mesh)
        call deallocate(source_surface_mesh)
        call deallocate(output_surface_mesh)

        deallocate(surface_nodes, source_surface_nodes, output_surface_nodes)

  end subroutine surface_elevation_smoothing

  subroutine get_coordinate_smooth(source, positions, elevation, output, smoothing_length, lump_mass, theta)

        type(vector_field), intent(in) :: positions
        type(scalar_field), intent(in) :: source, elevation
        type(scalar_field), intent(inout) :: output
        real, intent(in) :: smoothing_length
        logical, intent(in) :: lump_mass
        real, intent(in) :: theta

        type(element_type), pointer :: shape
        integer, dimension(:), pointer :: nodes

        type(csr_matrix) :: M
        type(csr_sparsity) :: sparsity
        type(scalar_field) :: rhs
        !integer, dimension(:), allocatable :: faceglobalnodes
        !integer :: face, globnod
        integer :: ele, dim, i

        real, dimension(ele_ngi(output,1)) :: detwei
        !real, dimension(face_ngi(positions,1)) :: sdetwei, flux
        !real, dimension(mesh_dim(positions), face_ngi(positions,1)) :: normal
        real, dimension(ele_loc(positions,1), ele_ngi(positions,1), mesh_dim(positions)) :: do_t

        type(scalar_field) :: lumpedmass, inverse_lumpedmass

        sparsity = make_sparsity(output%mesh, output%mesh, "SparsityMat")

        call allocate(M, sparsity, name = "MassMat")
        call allocate(rhs, output%mesh, "RHS")

        call zero(M)
        call zero(rhs)
        call zero(output)

        if (lump_mass) then
            call allocate(lumpedmass, output%mesh, "LumpedMass")
            call allocate(inverse_lumpedmass, output%mesh, "InverseLumpedMass")
            call zero(lumpedmass)
        end if

        !!ewrite(1,*) 'JN - SOURCE DIM:', mesh_dim(source)
        !!ewrite(1,*) 'JN - POSITIONS DIM:', mesh_dim(positions)
        !!ewrite(1,*) 'JN - OUTPUT DIM:', mesh_dim(output)
        !!ewrite(1,*) 'JN - POSITIONS ELE LOC:', ele_loc(positions,1)
        !!ewrite(1,*) 'JN - POSITIONS ELE GAUSS POINTS:', ele_ngi(positions,1)
        !!ewrite(1,*) 'JN - POSITIONS FACE LOC:', face_loc(positions,1)
        !!ewrite(1,*) 'JN - POSITIONS FACE GAUSS POINTS:', face_ngi(positions,1)

        do ele = 1, element_count(output)
        !!! element loop

        !!ewrite(1,*) 'JN - ELEMENT:', ele
        !!ewrite(1,*) 'JN - VAL IN ELEMENT:', ele_val_at_quad(source, ele)
        !!ewrite(1,*) 'JN - ELEMENT LOC:', ele_loc(output,ele)

            shape => ele_shape(output,ele)
            !ewrite(1,*) 'JN - SHAPES IN ELEMENT:', shape
            nodes => ele_nodes(output,ele)
            !!ewrite(1,*) 'JN - NODES IN ELEMENT:', nodes

            call transform_to_physical(positions, ele, shape=shape, dshape=do_t, detwei=detwei)
            !! build the mass matrix
            if (.not. lump_mass) then
                call addto(M,nodes,nodes,shape_shape(shape,shape,detwei)+((smoothing_length)**2)*dshape_dot_dshape(do_t,do_t,detwei))
            end if
            !ewrite(1,*) 'JN - MASS MATRIX:', node_val(M, ele)

            !! generate element rhs
            call addto(rhs, nodes, get_coordinate_smooth_element(ele_loc(output,ele)))
            if (lump_mass) then
                call addto(lumpedmass, nodes, shape_rhs(shape,detwei))
            end if
            !!ewrite(1,*) 'JN - RHS VECTOR:', node_val(RHS, ele)

        end do

        !!do face = 1, face_count(source)
        !!! face loop

        !!ewrite(1,*) 'JN - FACE:', face
        !!ewrite(1,*) 'JN - VAL IN FACE:', face_val_at_quad(source, face)
        !!ewrite(1,*) 'JN - FACE LOC:', face_loc(output,face)

            !!shape => face_shape(positions,face)
            !ewrite(1,*) 'JN - SHAPES IN FACE:', shape

            !!ewrite(1,*) 'JN - NODES IN FACE:', face_global_nodes(source,face)

            !!call transform_facet_to_physical(positions, face, detwei_f=sdetwei, normal=normal)
            !! this could be cleverer
            !!ewrite(1,*) 'JN - NORMAL:', normal
            !!flux = sum(face_val_at_quad(source, face) * normal, dim=1)
            !! generate_element_rhs
            !!call addto(rhs, face_global_nodes(source,face), get_coordinate_smooth_face(face_loc(output,face)))
            !!ewrite(1,*) 'JN - RHS VECTOR PLUS FACE CONTRIBUTION:', node_val(RHS, ele)

        !!end do

        if (lump_mass) then
            call invert(lumpedmass, inverse_lumpedmass)
            call set(output, rhs)
            call scale(output, inverse_lumpedmass)

            call deallocate(lumpedmass)
            call deallocate(inverse_lumpedmass)
        else
            call petsc_solve(output,M,rhs)
        end if

        call addto(output,elevation,-1.0)

        if (dt /= 0.0) then
            call scale(output, 1.0/dt)
        else
            call scale(output, 1.0)
        end if

        call set(output, output, source, theta)

        call deallocate(M)
        call deallocate(sparsity)
        call deallocate(rhs)

        contains

        function get_coordinate_smooth_element(nloc) result(local_rhs)

        integer ::nloc
        real, dimension(nloc) :: local_rhs

        local_rhs = shape_rhs(shape, (ele_val_at_quad(elevation,ele) + ele_val_at_quad(source,ele) * dt) * detwei)

        end function get_coordinate_smooth_element

        !!function get_coordinate_smooth_face(nloc) result(local_rhs)

        !!integer:: nloc
        !!real, dimension(nloc) :: local_rhs

        !!local_rhs = beta * shape_rhs(shape, flux * sdetwei)

        !!end function get_coordinate_smooth_face

   end subroutine get_coordinate_smooth

end module bedload_diagnostics
