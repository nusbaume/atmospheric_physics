module prescribe_water_tracers
  !    Test scheme which, at the start of every physics time step, sets each
  !    water tracer constituent to the values of the bulk water species it
  !    tracks multiplied by that water tracer's "prescribed_ratio" property
  !    (clipped at the water tracer's minimum value).
  !
  !    This scheme is meant for physics snapshot tests (i.e. the null dycore),
  !    where the bulk water species are re-read from the snapshot file at the
  !    start of every time step, so that every water tracer also starts each
  !    time step at exactly its prescribed ratio to its bulk water species.
  !    This allows a physics scheme's conservation of water tracer ratios to
  !    be tested.  Any water tracer values already present (e.g. read from the
  !    snapshot file) are overwritten.
  !
  !    WARNING: This scheme resets the water tracers every time step, and so
  !    should only be added to suites used for physics snapshot testing.
  !
  !    The "bulk_water_index" property of each water tracer is expected to
  !    already have been set by the host model.  The "run" phase is a no-op.

  implicit none
  private

  public :: prescribe_water_tracers_timestep_init
  public :: prescribe_water_tracers_run

contains

!> \section arg_table_prescribe_water_tracers_timestep_init  Argument Table
!! \htmlinclude prescribe_water_tracers_timestep_init.html
  subroutine prescribe_water_tracers_timestep_init(const_props, const_array, &
       errmsg, errcode)

     ! Use statements
     use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
     use ccpp_constituent_prop_mod, only: stdname_len
     use ccpp_constituent_prop_mod, only: int_unassigned, kphys_unassigned
     use ccpp_kinds,                only: kind_phys

     ! Input arguments
     type(ccpp_constituent_prop_ptr_t), intent(in)    :: const_props(:)

     ! Input/output arguments
     real(kind_phys),                   intent(inout) :: const_array(:,:,:)

     ! Output arguments
     character(len=*),                  intent(out)   :: errmsg
     integer,                           intent(out)   :: errcode

     ! Local variables
     integer                     :: tracer_idx  ! Water tracer constituent index
     integer                     :: bulk_idx    ! Bulk water species constituent index
     logical                     :: is_tracer   ! Is this constituent a water tracer?
     real(kind_phys)             :: ratio_val   ! Water tracer prescribed ratio
     real(kind_phys)             :: min_val     ! Water tracer minimum value
     character(len=stdname_len)  :: tracer_name ! Water tracer standard name
     character(len=*), parameter :: subname = 'prescribe_water_tracers_timestep_init: '

     errcode = 0
     errmsg = ''

     ! Note that bulk water species are never water tracers themselves, so
     ! the order in which the water tracers are set does not matter:
     do tracer_idx = 1, size(const_props)

        ! Only water tracers are set here:
        call const_props(tracer_idx)%is_water_tracer(is_tracer, errcode, errmsg)
        if (errcode /= 0) then
           return
        end if
        if (.not. is_tracer) then
           cycle
        end if

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

        ! Set the water tracer from its bulk water species:
        const_array(:,:,tracer_idx) = max(min_val, &
             const_array(:,:,bulk_idx) * ratio_val)

     end do

  end subroutine prescribe_water_tracers_timestep_init

!> \section arg_table_prescribe_water_tracers_run  Argument Table
!! \htmlinclude prescribe_water_tracers_run.html
  subroutine prescribe_water_tracers_run(errmsg, errcode)

     ! Dummy arguments
     character(len=*),   intent(out) :: errmsg
     integer,            intent(out) :: errcode

     errcode = 0
     errmsg = ''

  end subroutine prescribe_water_tracers_run

end module prescribe_water_tracers
