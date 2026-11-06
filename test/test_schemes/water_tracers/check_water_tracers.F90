module check_water_tracers
  !    Test scheme which checks that each water tracer
  !    constituent still equals the bulk water species
  !    it tracks multiplied by that water tracer's
  !    "prescribed_ratio" property, and
  !    returns an error if it does not.
  !
  !    For now this scheme should only be included in
  !    test suites if the goal is to ensure that water
  !    tracers are properly conserved.
  !
  !    The check uses the relative difference between each water tracer and
  !    its expected value (bulk water times prescribed ratio):
  !
  !       |tracer - expected| / |expected| <= ratio_rel_tolerance
  !
  !    Points where the expected value is at or below the water tracer's
  !    minimum value are skipped, as clipping at that minimum means the
  !    prescribed ratio cannot be held exactly there.
  !
  !    If no water tracer constituents are present, then a message saying so
  !    is written to the log (every time step), as in that case there is
  !    nothing to check.  This is not treated as an error.
  !
  !    The "bulk_water_index" property of each water tracer is expected to
  !    already have been set by the host model.  The "run" phase is a no-op.

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: check_water_tracers_timestep_final
  public :: check_water_tracers_run

  ! Maximum allowed relative difference between a water tracer and its
  ! bulk water species multiplied by the prescribed ratio:
  real(kind_phys), parameter :: ratio_rel_tolerance = 1.0e-10_kind_phys

contains

!> \section arg_table_check_water_tracers_timestep_final  Argument Table
!! \htmlinclude check_water_tracers_timestep_final.html
  subroutine check_water_tracers_timestep_final(const_props, &
       const_array, amIRoot, iulog, errmsg, errcode)

     ! Use statements
     use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
     use ccpp_constituent_prop_mod, only: stdname_len
     use ccpp_constituent_prop_mod, only: int_unassigned, kphys_unassigned

     ! Input arguments
     type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)
     real(kind_phys),                   intent(in)  :: const_array(:,:,:)
     logical,                           intent(in)  :: amIRoot ! MPI root task?
     integer,                           intent(in)  :: iulog   ! Log output unit

     ! Output arguments
     character(len=*),                  intent(out) :: errmsg
     integer,                           intent(out) :: errcode

     ! Local variables
     integer                     :: tracer_idx   ! Water tracer constituent index
     integer                     :: bulk_idx     ! Bulk water species constituent index
     integer                     :: icol, ilev   ! Column and level indices
     integer                     :: worst_col    ! Column of largest relative difference
     integer                     :: worst_lev    ! Level of largest relative difference
     logical                     :: is_tracer    ! Is this constituent a water tracer?
     logical                     :: has_tracers  ! Are any water tracers present?
     real(kind_phys)             :: ratio_val    ! Water tracer prescribed ratio
     real(kind_phys)             :: min_val      ! Water tracer minimum value
     real(kind_phys)             :: expected     ! Expected water tracer value
     real(kind_phys)             :: rel_diff     ! Relative difference from expected
     real(kind_phys)             :: worst_diff   ! Largest relative difference found
     character(len=stdname_len)  :: tracer_name  ! Water tracer standard name
     character(len=stdname_len)  :: bulk_name    ! Bulk water species standard name
     character(len=16)           :: col_str      ! Column index, as a string
     character(len=16)           :: lev_str      ! Level index, as a string
     character(len=32)           :: tracer_str   ! Water tracer value, as a string
     character(len=32)           :: bulk_str     ! Bulk water value, as a string
     character(len=32)           :: ratio_str    ! Prescribed ratio, as a string
     character(len=16)           :: diff_str     ! Relative difference, as a string
     character(len=*), parameter :: subname = 'check_water_tracers_timestep_final: '

     errcode = 0
     errmsg = ''

     has_tracers = .false.

     do tracer_idx = 1, size(const_props)

        ! Only water tracers are checked here:
        call const_props(tracer_idx)%is_water_tracer(is_tracer, errcode, errmsg)
        if (errcode /= 0) then
           return
        end if
        if (.not. is_tracer) then
           cycle
        end if
        has_tracers = .true.

        call const_props(tracer_idx)%standard_name(tracer_name, errcode, errmsg)
        if (errcode /= 0) then
           return
        end if

        ! Find the bulk water species this water tracer is tracking:
        call const_props(tracer_idx)%bulk_water_index(bulk_idx, errcode, errmsg)
        if (errcode /= 0) then
           return
        end if
        if (bulk_idx == int_unassigned) then
           errcode = 1
           write(errmsg,*) subname, 'Water tracer "', trim(tracer_name), &
                '" has no "bulk_water_index" property set'
           return
        end if
        if ((bulk_idx < 1) .or. (bulk_idx > size(const_array, 3))) then
           errcode = 1
           write(errmsg,*) subname, 'Water tracer "', trim(tracer_name), &
                '" has an out-of-range "bulk_water_index" of ', bulk_idx
           return
        end if

        call const_props(tracer_idx)%prescribed_ratio(ratio_val, errcode, errmsg)
        if (errcode /= 0) then
           return
        end if
        if (ratio_val == kphys_unassigned) then
           errcode = 1
           write(errmsg,*) subname, 'Water tracer "', trim(tracer_name), &
                '" has no "prescribed_ratio" property set'
           return
        end if

        call const_props(tracer_idx)%minimum(min_val, errcode, errmsg)
        if (errcode /= 0) then
           return
        end if

        ! Find the largest relative difference between the water tracer and
        ! its expected value, skipping points where the expected value would
        ! have been clipped at the water tracer's minimum value:
        worst_diff = 0.0_kind_phys
        worst_col  = 0
        worst_lev  = 0
        do ilev = 1, size(const_array, 2)
           do icol = 1, size(const_array, 1)
              expected = const_array(icol,ilev,bulk_idx) * ratio_val
              if (expected <= min_val) then
                 cycle
              end if
              if (abs(expected) <= tiny(expected)) then
                 cycle
              end if
              rel_diff = abs(const_array(icol,ilev,tracer_idx) - expected) / &
                   abs(expected)
              if (rel_diff > worst_diff) then
                 worst_diff = rel_diff
                 worst_col  = icol
                 worst_lev  = ilev
              end if
           end do
        end do

        if (worst_diff > ratio_rel_tolerance) then
           call const_props(bulk_idx)%standard_name(bulk_name, errcode, errmsg)
           if (errcode /= 0) then
              return
           end if
           ! Write each number separately, then build the message by
           ! assignment, which (unlike an internal write) safely truncates
           ! if the message is longer than <errmsg>:
           write(col_str,   '(i0)')     worst_col
           write(lev_str,   '(i0)')     worst_lev
           write(tracer_str,'(es24.16)') const_array(worst_col,worst_lev,tracer_idx)
           write(bulk_str,  '(es24.16)') const_array(worst_col,worst_lev,bulk_idx)
           write(ratio_str, '(es24.16)') ratio_val
           write(diff_str,  '(es10.3)')  worst_diff
           errcode = 1
           errmsg = subname//'Water tracer "'//trim(tracer_name)//         &
                '" does not match bulk water species "'//trim(bulk_name)// &
                '" times its prescribed ratio.  Largest relative '//       &
                'difference ('//trim(adjustl(diff_str))//') is at column '// &
                trim(col_str)//', level '//trim(lev_str)//': tracer = '//  &
                trim(adjustl(tracer_str))//', bulk = '//                   &
                trim(adjustl(bulk_str))//', ratio = '//                    &
                trim(adjustl(ratio_str))
           return
        end if

     end do

     ! Let the user know if there was nothing to check:
     if ((.not. has_tracers) .and. amIRoot) then
        write(iulog, *) 'check_water_tracers: No water tracer constituents ', &
             'are present, so no water tracer ratios were checked.'
     end if

  end subroutine check_water_tracers_timestep_final

!> \section arg_table_check_water_tracers_run  Argument Table
!! \htmlinclude check_water_tracers_run.html
  subroutine check_water_tracers_run(errmsg, errcode)

     ! Dummy arguments
     character(len=*),   intent(out) :: errmsg
     integer,            intent(out) :: errcode

     errcode = 0
     errmsg = ''

  end subroutine check_water_tracers_run

end module check_water_tracers
