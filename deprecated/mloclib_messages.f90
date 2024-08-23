!> Messaging procedures

module mloclib_messages

   implicit none
   save
   
contains

!*****************************************************************************************
   subroutine oops (msg)
   
   ! Error report and stop
   
      character(len=*) :: msg
   
      write (*,'(/a)') 'Oops! Fatal error in '//trim(msg)
   
      stop
      
   end subroutine oops

      
!*****************************************************************************************
   subroutine warnings (msg)
   
   ! Warnings
   
      character(len=*) :: msg
         
      write (*,'(t3,a)') 'Warning from '//trim(msg)
         
      return
      
   end subroutine warnings

      
!*****************************************************************************************
   subroutine fyi (msg)
   
   ! Informational messages
   
      character(len=*) :: msg
         
      write (*,'(t3,a)') 'FYI from '//trim(msg)
         
      return
      
   end subroutine fyi
      
!*****************************************************************************************
   subroutine allocation_error (procedure_name, array_name, error)
   
      integer, intent(in) :: error
      character(len=*), intent(in) :: array_name
      character(len=*), intent(in) :: procedure_name
      character(len=132) :: msg
      
      write (msg,'(a,i3,3a)') trim(procedure_name)//': allocation error (', error,&
       ') for array "', trim(array_name), '"'
      call oops (trim(msg))
      
      return
   
   end subroutine allocation_error


!*****************************************************************************************
end module mloclib_messages

