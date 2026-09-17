module initialize_water_tracers
  !    Scheme which initializes the water tracer by doing the following:
  !
  !    A. Setting the "bulk_water_index" property of each water tracer
  !       constituent, which is found by looking up the constituent index of
  !       the tracer's "bulk_water_name" property, and
  !
  !    B. Setting the initial values of each water tracer constituent to the
  !       values of the bulk water species it tracks multiplied by that water
  !       tracer's "prescribed_ratio" property.
  !
  !    This scheme should be called after all other physics schemes which
  !    initialize a water-species constituent have been initalized, as
  !    otherwise the water tracer initial value will be set incorrectly.

  implicit none
  private

  public :: initialize_water_tracers_init

contains

!> \section arg_table_initialize_water_tracers_init  Argument Table
!! \htmlinclude initialize_water_tracers_init.html
  subroutine initialize_water_tracers_init(const_props, const_array, &
       errmsg, errcode)

     ! Use statements
     use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
     use ccpp_constituent_prop_mod, only: stdname_len, kphys_unassigned
     use ccpp_constituent_prop_mod, only: int_unassigned
     use ccpp_scheme_utils,         only: ccpp_constituent_index
     use ccpp_kinds,                only: kind_phys

     ! Input/output arguments
     type(ccpp_constituent_prop_ptr_t), intent(inout) :: const_props(:)
     real(kind_phys),                   intent(inout) :: const_array(:,:,:)

     ! Output arguments
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errcode

     ! Local variables
     integer                    :: tracer_idx    ! Water tracer constituent index
     integer                    :: bulk_idx      ! Bulk water species constituent index
     logical                    :: is_tracer     ! Is this constituent a water tracer?
     real(kind_phys)            :: ratio_val     ! Water tracer prescribed ratio
     character(len=stdname_len) :: tracer_name   ! Water tracer standard name
     character(len=stdname_len) :: bulk_name     ! Bulk water species standard name
     character(len=*), parameter :: subname = 'initialize_water_tracers_init: '

     errcode = 0
     errmsg = ''

     do tracer_idx = 1, size(const_props)

        ! Only water tracers need to be initialized here:
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

        !----------------------------------------------------------------
        ! A. Find and set the index of the bulk water species being tracked
        !----------------------------------------------------------------

        ! Each water tracer is registered with the standard name of the bulk
        ! water species it tracks, so simply look up the constituent index of
        ! that standard name:
        call const_props(tracer_idx)%bulk_water_name(bulk_name, errcode, errmsg)
        if (errcode /= 0) then
           return
        end if

        if (len_trim(bulk_name) == 0) then
           errcode = 1
           write(errmsg,*) subname, 'Water tracer "', trim(tracer_name), &
                '" has no "bulk_water_name" property set'
           return
        end if

        call ccpp_constituent_index(bulk_name, bulk_idx, errcode, errmsg)
        if (errcode /= 0) then
           return
        end if

        ! An unregistered standard name is not flagged as an error, so check
        ! the returned index directly:
        if (bulk_idx == int_unassigned) then
           errcode = 1
           write(errmsg,*) subname, 'No registered constituent has the ', &
                'standard name "', trim(bulk_name), '", which is the bulk ', &
                'water species tracked by water tracer "', &
                trim(tracer_name), '"'
           return
        end if

        ! Record which bulk water species this water tracer is tracking.  Note
        ! that this will fail if the index has already been set:
        call const_props(tracer_idx)%set_bulk_water_index(bulk_idx, errcode, errmsg)
        if (errcode /= 0) then
           return
        end if

        !----------------------------------------------------------------
        ! B. Set the water tracer's initial values from the bulk species
        !----------------------------------------------------------------

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

        const_array(:,:,tracer_idx) = const_array(:,:,bulk_idx) * ratio_val

     end do

  end subroutine initialize_water_tracers_init

end module initialize_water_tracers
