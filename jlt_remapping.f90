module jlt_remapping
  use jlt_constant, only : STR_SHORT, STR_MID, STR_LONG
  use jlt_grid_class, only : grid_class
  
!--------------------------------   public  ----------------------------------!

public :: init_remapping                   ! subroutine (local_comm)
public :: make_grid_rank_file              ! subroutine (my_comp, my_grid, grid_index)
public :: make_grid_index_file             ! subroutine (my_comp, my_grid, my_rank, grid_index)

!--------------------------------   private  ---------------------------------!

integer :: my_comm
integer :: my_size
integer :: my_rank

contains

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine init_remapping(local_comm)
  use mpi
  implicit none
  integer, intent(IN) :: local_comm
  integer :: ierror
  
  my_comm = local_comm

  call mpi_comm_size(my_comm, my_size, ierror)
  call mpi_comm_rank(my_comm, my_rank, ierror)
  
  
end subroutine init_remapping

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine make_grid_rank_file(my_comp, my_grid, grid_index)
  use mpi
  use jlt_utils, only : set_fid, sort_int_1d, set_fid
  implicit none
  character(len=*), intent(IN) :: my_comp, my_grid    ! comp name and grid name
  integer, intent(IN)          :: grid_index(:)       ! local grid index
  character(len=STR_LONG)      :: file_name           ! grid_rank file name
  integer, allocatable         :: grid_size(:)        ! grid index array size of each process
  integer                      :: global_size         ! number of global grid points
  integer, allocatable         :: global_index(:)     ! global grid index
  integer, allocatable         :: global_rank(:)      ! global rank of grid points
  integer, allocatable         :: displs(:)           ! displacement
  integer                      :: int_array(1)
  integer                      :: local_size          ! local grid index size
  integer                      :: ierror
  integer                      :: i, j, counter
  integer                      :: fid = 222
  
  local_size   = size(grid_index)
  int_array(1) = local_size

  if (my_rank == 0) then
     allocate(grid_size(my_size))
     allocate(displs(my_size))
  else
     allocate(grid_size(1))
     allocate(displs(1))
  end if
  
  call mpi_gather(int_array, 1, MPI_INTEGER, grid_size, 1, MPI_INTEGER, 0, my_comm, ierror)

  global_size = 0
   
  do i = 1, my_size
     global_size = global_size + grid_size(i)
  end do

  if (my_rank == 0) then
     allocate(global_index(global_size))
     allocate(global_rank(global_size))
  else
     allocate(global_index(1))
     allocate(global_rank(1))
  end if

  if (my_rank == 0) then
     displs(1) = 0
     do i = 2, my_size
        displs(i) = displs(i-1) + grid_size(i)
     end do
     
     counter = 0
     do i = 1, my_size
        do j = 1, grid_size(i)
           counter = counter + 1
           global_rank(counter) = i - 1
        end do
     end do
  end if

  call mpi_gatherv(grid_index, local_size, MPI_INTEGER, global_index, grid_size, displs, MPI_INTEGER, 0, my_comm, ierror)

  if (my_rank == 0) then
    call sort_int_1d(global_size, global_index, global_rank)
    do i = 1, global_size
       write(*,*) global_index(i), global_rank(i)
    end do

    write(file_name, '("jlt.",A,".",A,".GRID_INDEX")') trim(my_comp), trim(my_grid)
    call set_fid(fid)
    open(fid, file=trim(file_name), form="unformatted", access="stream", action = 'write', status="replace")
    write(fid) global_index
    close(fid)

    write(file_name, '("jlt.",A,".",A,".GRID_RANK")') trim(my_comp), trim(my_grid)
    call set_fid(fid)
    open(fid, file=trim(file_name), form="unformatted", access="stream", action = 'write', status="replace")
    write(fid) global_rank
    close(fid)

  end if
 
end subroutine make_grid_rank_file

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine make_grid_index_file(my_name, my_grid, my_rank, grid_index)
  use jlt_utils, only : set_fid, sort_int_1d
  implicit none
  character(len=*), intent(IN) :: my_name, my_grid
  integer, intent(IN)          :: my_rank
  integer, intent(IN)          :: grid_index(:)
  character(len=STR_LONG)      :: file_name
  character(len=5)             :: pe_num = "00000"
  integer, allocatable         :: sorted_index(:)
  integer, allocatable         :: sorted_pos(:)
  integer                      :: grid_size
  integer                      :: fid = 201
  integer                      :: i
  
  grid_size = size(grid_index)

  allocate(sorted_index(grid_size))
  allocate(sorted_pos(grid_size))

  do i = 1, grid_size
     sorted_index(i) = grid_index(i)
     sorted_pos(i) = i
  end do

  call sort_int_1d(grid_size, sorted_index, sorted_pos)
  
  write(pe_num, '(I5.5)') my_rank
  
  write(file_name, '("jlt.",A,".",A,".PE",A)') trim(my_name), trim(my_grid), pe_num

  call set_fid(fid)

  open(fid, file=trim(file_name), form="unformatted", action = 'write')


  write(fid) grid_size, sorted_index, sorted_pos

  close(fid)

  deallocate(sorted_index)
  deallocate(sorted_pos)
  
end subroutine make_grid_index_file

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine make_local_mapping_table(global_index, global_target, global_coef, &
                                    grid_index, &
                                    local_size, local_index,  local_target,  local_coef)
  use jlt_utils, only : sort_int_1d, binary_search, error
  implicit none
  integer, intent(IN)      :: global_index(:), global_target(:)
  real(kind=8), intent(IN) :: global_coef(:)
  integer, intent(IN)      :: grid_index(:)
  integer, intent(OUT)     :: local_size
  integer, pointer         :: local_index(:), local_target(:)
  real(kind=8), pointer    :: local_coef(:)
  integer, pointer         :: sorted_index(:), sorted_pos(:)
  integer                  :: table_size
  integer                  :: current_pos, counter
  integer                  :: res
  integer                  :: i
  
  table_size = size(global_index)
  
  allocate(sorted_index(table_size))
  allocate(sorted_pos(table_size))

  do i = 1, table_size
     sorted_index(i) = global_index(i)
     sorted_pos(i)   = i
  end do

  call sort_int_1d(table_size, sorted_index, sorted_pos)

  counter = 0
  current_pos = 1
  do i = 1, size(grid_index)
     do while (sorted_index(current_pos) <= grid_index(i))
        if (sorted_index(current_pos) == grid_index(i)) counter = counter + 1
        current_pos = current_pos + 1
     end do
  end do

  local_size = counter
  allocate(local_index(local_size))
  allocate(local_target(local_size))
  allocate(local_coef(local_size))

  counter = 0
  current_pos = 1
  do i = 1, size(grid_index)
     do while (sorted_index(current_pos) <= grid_index(i))
        if (sorted_index(current_pos) == grid_index(i)) then
           counter = counter + 1
           local_index(counter)  = sorted_index(current_pos)
           local_target(counter) = global_target(sorted_pos(current_pos))
           local_coef(counter)   = global_coef(sorted_pos(current_pos))
        end if
        current_pos = current_pos + 1
     end do
  end do

  do i = 1, size(grid_index)
     res = binary_search(local_index, grid_index(i))
     if (res <= 0) then
        call error("jlt_remapping:make_local_mapping_table", "grid index not listed in mapping table")
     end if
  end do
  
end subroutine make_local_mapping_table

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine get_target_rank(target_name, num_of_target, target_index, target_rank)
  implicit none
  character(len=*), intent(IN) :: target_name
  integer, intent(IN)          :: num_of_target
  integer, intent(IN)          :: target_index(:)
  integer, intent(OUT)         :: target_rank(:)

  target_rank(:) = 0
  
end subroutine get_target_rank

!=======+=========+=========+=========+=========+=========+=========+=========+
end module jlt_remapping

