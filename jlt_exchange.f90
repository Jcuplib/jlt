module jlt_exchange
  use jlt_exchange_class, only : exchange_class
  implicit none
  private

!--------------------------------   public   ---------------------------------!

  public :: init_exchange                     ! subroutine ()
  public :: set_mapping_table                 ! subroutine (my_name, send_comp_name, send_grid_name, recv_comp_name, recv_grid_name,
                                              !             map_tag, intpl_flag, send_grid, recv_grid, coef)
  public :: is_exchange_assigned              ! logical function (send_comp_name, send_grid_name, recv_comp_name, recv_grid_name)
  public :: get_exchange_ptr                  ! type(exchange_type), pointer, function (send_comp_name, send_grid_name, recv_comp_name, recv_grid_name)
  
!--------------------------------   private  ---------------------------------!

  integer, parameter :: MAX_EXCHANGE = 8
  integer :: num_of_exchange

  type(exchange_class), pointer :: exchange(:)

  
contains

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine init_exchange()
  implicit none

  num_of_exchange = 0
  allocate(exchange(MAX_EXCHANGE))
  
end subroutine init_exchange

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_mapping_table(my_name, send_comp_name, send_grid_name, recv_comp_name, recv_grid_name, &
                             map_tag, intpl_flag, send_grid, recv_grid, coef)
  use jlt_grid, only : get_grid_ptr
  implicit none
  character(len=*), intent(IN) :: my_name
  character(len=*), intent(IN) :: send_comp_name, send_grid_name
  character(len=*), intent(IN) :: recv_comp_name, recv_grid_name
  integer, intent(IN)          :: map_tag
  logical, intent(IN)          :: intpl_flag
  integer, intent(IN)          :: send_grid(:), recv_grid(:)
  real(kind=8), intent(IN)     :: coef(:)

  num_of_exchange = num_of_exchange + 1

  exchange(num_of_exchange) = exchange_class(trim(my_name), intpl_flag)

  call exchange(num_of_exchange)%set_mapping_table(trim(send_comp_name), trim(send_grid_name), &
                                                   trim(recv_comp_name), trim(recv_grid_name), &
                                                   map_tag, send_grid, recv_grid, coef)

end subroutine set_mapping_table

!=======+=========+=========+=========+=========+=========+=========+=========+

function is_exchange_assigned(send_comp_name, send_grid_name, &
                              recv_comp_name, recv_grid_name) result (res)
  implicit none
  character(len=*), intent(IN) :: send_comp_name, send_grid_name
  character(len=*), intent(IN) :: recv_comp_name, recv_grid_name
  logical :: res
  integer :: i

  do i = 1, num_of_exchange
     if (exchange(i)%is_my_exchange(send_comp_name, send_grid_name, &
                                    recv_comp_name, recv_grid_name)) then
        res = .true.
        return
     end if
  end do

  res = .false.
  
end function is_exchange_assigned

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_exchange_ptr(send_comp_name, send_grid_name, &
                          recv_comp_name, recv_grid_name) result (res)
  use jlt_utils, only : error
  implicit none
  character(len=*), intent(IN) :: send_comp_name, send_grid_name
  character(len=*), intent(IN) :: recv_comp_name, recv_grid_name
  type(exchange_class), pointer :: res
  integer :: i

  do i = 1, num_of_exchange
     if (exchange(i)%is_my_exchange(send_comp_name, send_grid_name, &
                                    recv_comp_name, recv_grid_name)) then
        res => exchange(i)
        return
     end if
  end do

  call error("[jlt_exchange:get_exchange_ptr]","no such name:"//trim(send_comp_name)//":"// &
             trim(send_grid_name)//":"//trim(recv_comp_name)//":"//trim(recv_grid_name))
  
  res => null()
  
end function get_exchange_ptr

!=======+=========+=========+=========+=========+=========+=========+=========+

end module jlt_exchange
