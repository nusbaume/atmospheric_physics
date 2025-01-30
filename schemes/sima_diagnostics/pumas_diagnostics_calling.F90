module pumas_diagnostics_calling

   use ccpp_kinds, only:  kind_phys

   implicit none
   private
   save

   public :: pumas_diagnostics_init ! init routine                                                
   public :: pumas_diagnostics_run  ! main routine
      
CONTAINS



   !> \section arg_table_pumas_diagnostics_init  Argument Table                                   
   !! \htmlinclude pumas_diagnostics_init.html                                                    
   subroutine pumas_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:                                                                                      

      errmsg = ''
      errflg = 0
      
   end subroutine pumas_diagnostics_init

   !> \section arg_table_pumas_diagnostics_run  Argument Table                                    
   !! \htmlinclude pumas_diagnostics_run.html                                                     
   subroutine pumas_diagnostics_run(dudt, dvdt, dsdt, errmsg, errflg)

      use cam_history, only: history_out_field
      !------------------------------------------------                                                       
      !   Input / output parameters                                                                           
      !------------------------------------------------                                                       
      real(kind_phys), intent(in) :: dudt(:,:) !tendency_of_x_wind due to RF                                  
      real(kind_phys), intent(in) :: dvdt(:,:) !tendency_of_y_wind due to RF                                  
      real(kind_phys), intent(in) :: dsdt(:,:) !tendency_of_y_wind due to RF                                  

      ! CCPP error handling variables                                                                         
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      end subroutin pumas_diagnostics_run 

end module pumas_diagnostics_calling

      
      
