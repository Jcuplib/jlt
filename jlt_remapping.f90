module jlt_remapping
  use jlt_constant, only : STR_SHORT, STR_MID, STR_LONG
  use jlt_grid_class, only : grid_class
  
!--------------------------------   public  ----------------------------------!

public :: init_remapping                   ! subroutine (local_comm)
public :: make_grid_rank_file              ! subroutine (my_comp, my_grid, grid_index)
public :: delete_grid_rank_file            ! subroutine (my_comp, my_grid)
public :: make_local_mapping_table_no_sort ! subroutine (my_index_global, target_index_global, coef_global,
                                           !             grid_index, my_index_local, target_index_local, coef_local) 
public :: make_local_mapping_table         ! subroutine (my_index_global, target_index_global, coef_global,
                                           !             grid_index, my_index_local, target_index_local, coef_local) 
public :: get_target_grid_rank             ! subroutine (target_comp, target_grid, target_index, target_rank)
public :: reorder_index_by_target_rank     ! subroutine (send_index, recv_index, target_rank)
public :: delete_same_index                ! subroutine (index_in, index_out)
public :: make_exchange_table              ! subroutine (target_rank, local_index, exchange_index, num_of_index)
public :: make_conversion_table            ! subroutine (mapping_table, grid_index, conv_table)

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
  use jlt_utils, only : set_fid, sort_int_1d, set_fid, error, put_log
  implicit none
  character(len=*), intent(IN) :: my_comp, my_grid    ! comp name and grid name
  integer, intent(IN)          :: grid_index(:)       ! local grid index
  character(len=STR_LONG)      :: file_name           ! grid_rank file name
  integer, allocatable         :: sorted_index(:)     ! sorted local grid index
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
  type local_grid_type
     integer :: current_pos = 1
     integer :: local_size
     integer, pointer :: local_index(:)
  end type local_grid_type

  type(local_grid_type), pointer :: local_grid(:)
  integer :: ista(MPI_status_size)
  integer :: current_index
  integer :: current_rank


  if (my_rank == 0) then
     allocate(local_grid(my_size))
  end if
  
  local_size   = size(grid_index)
  int_array(1) = local_size

  if (my_rank == 0) then
     allocate(grid_size(my_size))
  else
     allocate(grid_size(1))
  end if
  
  call mpi_gather(int_array, 1, MPI_INTEGER, grid_size, 1, MPI_INTEGER, 0, my_comm, ierror)

  allocate(sorted_index(local_size))

  call put_log("make_grid_rank_file : local sorting start")

  sorted_index(:) = grid_index(:)
  call sort_int_1d(local_size, sorted_index)
  
  call put_log("make_grid_rank_file : local sorting finish")

  call put_log("make_grid_rank_file : send/recv staart")

  if (my_rank == 0) then
     global_size = 0
     do i = 1, my_size
        local_grid(i)%local_size = grid_size(i)
        allocate(local_grid(i)%local_index(grid_size(i)))
        global_size = global_size + grid_size(i)
     end do
     local_grid(1)%local_index(:) = sorted_index(:)
     do i = 2, my_size
        call mpi_recv(local_grid(i)%local_index, grid_size(i), MPI_INTEGER, i-1, 0, my_comm, ista, ierror)
        if (ierror /= 0) then
           call error("make_grid_rank_file2", "mpi_recv error")
        end if
     end do
     allocate(global_index(global_size))
     allocate(global_rank(global_size))
  else
     call mpi_send(sorted_index, local_size, MPI_INTEGER, 0, 0, my_comm, ierror)
     if (ierror /= 0) then
        call error("make_grid_rank_file2", "mpi_send error")
     end if
  end if

  call put_log("make_grid_rank_file : send/recv finish")

  if (my_rank /= 0) then
     deallocate(grid_size)
     return
  end if
  

  ! root rank only 

  call put_log("make_grid_rank_file : global sorting staart")

  ! sorting
  current_index = 99999999 ! huge number
  current_rank  = 0
  do i = 1, global_size
     do j = 1, my_size
        if (local_grid(j)%current_pos <= local_grid(j)%local_size) then
           if (local_grid(j)%local_index(local_grid(j)%current_pos) < current_index) then
              current_index = local_grid(j)%local_index(local_grid(j)%current_pos)
              current_rank  = j
           end if
        end if
     end do
     global_index(i) = current_index
     global_rank(i)  = current_rank - 1
     local_grid(current_rank)%current_pos = local_grid(current_rank)%current_pos + 1
     current_index = 99999999 ! huge number
  end do

  call put_log("make_grid_rank_file : global sorting finish")

  call put_log("make_grid_rank_file : file output start")

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

  call put_log("make_grid_rank_file : file output finish")

  deallocate(global_index)
  deallocate(global_rank)
  
end subroutine make_grid_rank_file

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine delete_grid_rank_file(my_comp, my_grid)
  use jlt_utils, only : set_fid, sort_int_1d, set_fid, error, put_log
  implicit none
  character(len=*), intent(IN) :: my_comp, my_grid    ! comp name and grid name
  character(len=STR_LONG)      :: file_name           ! grid_rank file name
  integer                      :: stat
  integer                      :: fid = 222

  if (my_rank /= 0) return ! root only

  call put_log("delete_grid_rank_file : file delete start")

  write(file_name, '("jlt.",A,".",A,".GRID_INDEX")') trim(my_comp), trim(my_grid)
  call set_fid(fid)
  open(fid, file=trim(file_name),  iostat = stat)
  close(fid, status="delete", iostat = stat)

  write(file_name, '("jlt.",A,".",A,".GRID_RANK")') trim(my_comp), trim(my_grid)
  call set_fid(fid)
  open(fid, file=trim(file_name), iostat = stat)
  close(fid, status="delete", iostat = stat)

  call put_log("delete_grid_rank_file : file delete finish")

  
end subroutine delete_grid_rank_file

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine make_grid_rank_file_org(my_comp, my_grid, grid_index)
  use mpi
  use jlt_utils, only : set_fid, sort_int_1d, set_fid, put_log
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

  if (my_rank == 0) then
     do i = 1, my_size
       global_size = global_size + grid_size(i)
     end do
     allocate(global_index(global_size))
     allocate(global_rank(global_size))
  else
     allocate(global_index(1))
     allocate(global_rank(1))
  end if

  if (my_rank == 0) then
     displs(1) = 0
     do i = 2, my_size
        displs(i) = displs(i-1) + grid_size(i-1)
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

    call put_log("make_grid_rank_file : sorting start")
     
    call sort_int_1d(global_size, global_index, global_rank)

    call put_log("make_grid_rank_file : sorting finish")

    call put_log("make_grid_rank_file : file output start")

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

    call put_log("make_grid_rank_file : file output finish")

 end if

  deallocate(grid_size)
  deallocate(displs)
  deallocate(global_index)
  deallocate(global_rank)
  
end subroutine make_grid_rank_file_org

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine make_grid_rank_file_local(my_comp, my_grid, grid_index)
  use mpi
  use jlt_utils, only : set_fid, sort_int_1d, set_fid
  implicit none
  character(len=*), intent(IN) :: my_comp, my_grid    ! comp name and grid name
  integer, intent(IN)          :: grid_index(:)       ! local grid index
  character(len=STR_LONG)      :: file_name           ! grid_rank file name
  integer, allocatable         :: grid_info(:,:)      ! grid index array size of each process
  integer                      :: global_size         ! number of global grid points
  integer, allocatable         :: global_index(:)     ! global grid index
  integer, allocatable         :: global_rank(:)      ! global rank of grid points
  integer, allocatable         :: displs(:)           ! displacement
  integer                      :: int_array(3)
  integer                      :: local_size          ! local grid index size
  integer                      :: ierror
  integer                      :: i, j, counter
  integer                      :: fid = 222
  

  local_size   = size(grid_index)
  int_array(1) = local_size
  int_array(2) = minval(grid_index)
  int_array(3) = maxval(grid_index)
  
  if (my_rank == 0) then
     allocate(grid_info(3, my_size))
  else
     allocate(grid_info(1,1))
  end if
  
  call mpi_gather(int_array, 3, MPI_INTEGER, grid_info, 3, MPI_INTEGER, 0, my_comm, ierror)

  if (my_rank == 0) then
     write(file_name, '("jlt.",A,".",A,".GRID_INFO")') trim(my_comp), trim(my_grid)
     call set_fid(fid)
     open(fid, file=trim(file_name), form="unformatted", access="stream", action = 'write', status="replace")
     write(fid) my_size
     write(fid) grid_info
     close(fid)
  end if
  
  global_size = size(grid_index)
  
  allocate(global_index(global_size))
  allocate(global_rank(global_size))
     
  call sort_int_1d(global_size, global_index, global_rank)

  write(file_name, '("jlt.",A,".",A,".GRID_INDEX.PE", I5.5)') trim(my_comp), trim(my_grid), my_rank
  call set_fid(fid)
  open(fid, file=trim(file_name), form="unformatted", access="stream", action = 'write', status="replace")
  write(fid) global_index
  close(fid)

  write(file_name, '("jlt.",A,".",A,".GRID_RANK.PE", I5.5)') trim(my_comp), trim(my_grid), my_rank
  call set_fid(fid)
  open(fid, file=trim(file_name), form="unformatted", access="stream", action = 'write', status="replace")
  write(fid) global_rank
  close(fid)


  deallocate(global_index)
  deallocate(global_rank)
  
end subroutine make_grid_rank_file_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine get_target_grid_rank(target_comp, target_grid, target_index, target_rank)
  use jlt_utils, only : set_fid, sort_int_1d, error
  implicit none
  character(len=*), intent(IN) :: target_comp
  character(len=*), intent(IN) :: target_grid
  integer, intent(IN)          :: target_index(:)
  integer, intent(OUT)         :: target_rank(:)
  integer, allocatable         :: index(:)
  integer, allocatable         :: rank(:)
  integer, allocatable         :: sorted_index(:)
  integer, allocatable         :: sorted_pos(:)
  integer                      :: chunk_size = 20
  integer                      :: num_of_read
  integer                      :: mod_size
  character(len=STR_LONG)      :: file_name
  integer                      :: fid1 = 801
  integer                      :: fid2 = 802
  integer                      :: file_size
  integer                      :: current_pos
  integer                      :: i, j


  if (size(target_index) == 0) return
  
  write(file_name, '("jlt.",A,".",A,".GRID_INDEX")') trim(target_comp), trim(target_grid)
  call set_fid(fid1)
  open(fid1, file=trim(file_name), form="unformatted", access="stream", action="read", status="old", err=900)

  write(file_name, '("jlt.",A,".",A,".GRID_RANK")') trim(target_comp), trim(target_grid)
  call set_fid(fid2)
  open(fid2, file=trim(file_name), form="unformatted", access="stream", action="read", status="old", err=900)

  inquire(fid1, size = file_size)

  file_size = file_size/4 ! byte to data
  num_of_read = file_size/chunk_size
  mod_size    = file_size - chunk_size*num_of_read

  allocate(index(chunk_size))
  allocate(rank(chunk_size))

  allocate(sorted_index(size(target_index)))
  allocate(sorted_pos(size(target_index)))

  do i = 1, size(target_index)
     sorted_index(i) = target_index(i)
     sorted_pos(i) = i
  end do

  call sort_int_1d(size(target_index), sorted_index, sorted_pos)

  !if (trim(target_comp) == "ILSIO") then
  !   write(0, *) "get_target_grid_rank, "//trim(target_comp)//", ", chunk_size, file_size, num_of_read, mod_size, size(target_index), target_index(1:10)
  !end if

  current_pos = 1

  do i = 1, num_of_read
     read(fid1) index
     read(fid2) rank

     if (index(chunk_size) < sorted_index(current_pos)) cycle
     
     !if (trim(target_comp) == "ILSIO") write(0, *) index, "/", rank, "/", current_pos, sorted_index(current_pos),"/"

     do j = 1, chunk_size
        !write(0, '(A,I4, I6,I4)') "current_pos, sorted_index, index = ", current_pos, sorted_index(current_pos), index(j)
        
        if (sorted_index(current_pos) == index(j)) then
           !if (trim(target_comp) == "ILSIO") then
           !   write(0, *) "target_rank = ", i, j, current_pos, index(j), rank(j)
           !end if
           target_rank(sorted_pos(current_pos)) = rank(j)
           current_pos = current_pos + 1
           if (size(sorted_index) < current_pos) goto 100
           do while(sorted_index(current_pos) == sorted_index(current_pos-1))
              !if (trim(target_comp) == "ILSIO") then
              !   write(0, *) "target_rank = ", i, j, current_pos, index(j), rank(j)
              !end if
              target_rank(sorted_pos(current_pos)) = rank(j)
              current_pos = current_pos + 1
              if (size(sorted_index) < current_pos) goto 100
           end do
           if (current_pos > size(target_index)) then
              goto 100
           end if
        end if
     end do
  end do

  if (mod_size > 0) then
     read(fid1) index(1:mod_size)
     read(fid2) rank(1:mod_size)

     do j = 1, mod_size
        if (sorted_index(current_pos) == index(j)) then
           target_rank(sorted_pos(current_pos)) = rank(j)
           current_pos = current_pos + 1
           if (size(sorted_index) < current_pos) goto 100
           do while(sorted_index(current_pos) == sorted_index(current_pos-1)) 
              target_rank(sorted_pos(current_pos)) = rank(j)
              current_pos = current_pos + 1
              if (size(sorted_index) < current_pos) goto 100
           end do
           if (current_pos > size(target_index)) then
              goto 100
           end if
        end if
     end do
  end if

100 continue
  
  close(fid1)
  close(fid2)
  
  deallocate(index)
  deallocate(rank)

  return
  
900 continue
  write(0, *) "file open error!!!!, file name = "//trim(file_name)
  call error("jlt_remapping:get_target_grid_rank", "file open errror, file name = "//trim(file_name))   
  
end subroutine get_target_grid_rank

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

subroutine make_local_mapping_table_no_sort(global_index, global_target, global_coef, &
                                           grid_index, &
                                           local_size, local_index,  local_target,  local_coef)
  use jlt_utils, only : sort_int_1d, binary_search, error
  use jlt_mpi_lib, only : jml_GetMyrankGlobal
  implicit none
  integer, intent(IN)      :: global_index(:), global_target(:)
  real(kind=8), intent(IN) :: global_coef(:)
  integer, intent(IN)      :: grid_index(:)
  integer, intent(OUT)     :: local_size
  integer, allocatable, intent(OUT)     :: local_index(:), local_target(:)
  real(kind=8), allocatable, intent(OUT) :: local_coef(:)
  integer, allocatable     :: sorted_index(:), sorted_pos(:)
  integer, allocatable     :: sorted_target(:)
  real(kind=8), allocatable :: sorted_coef(:)
  integer, allocatable      :: sorted_grid_index(:)
  integer                  :: table_size
  integer                  :: current_pos, counter
  integer                  :: res
  integer                  :: i

  !do i = 1, size(global_index)
  !   write(100+jml_GetMyrankGlobal(), *) "global_mapping ", global_index(i), global_target(i), global_coef(i)
  !end do
  
  !do i = 1, size(grid_index)
  !   write(100+jml_GetMyrankGlobal(), *) "grid_index ", grid_index(i)
  !end do

  table_size = size(grid_index)

  
  allocate(sorted_index(table_size)) 
  
  do i = 1, table_size
     sorted_index(i) = grid_index(i)
  end do

  call sort_int_1d(table_size, sorted_index)

  counter = 0
  do i = 1, size(global_index)
     res = binary_search(sorted_index, global_index(i))
     if (res > 0) counter = counter + 1
  end do

  local_size = counter

  !if (local_size == 0) then
  !   local_index  => null()
  !   local_target => null()
  !   local_coef   => null()
  !   return
  !end if

  allocate(local_index(local_size))
  allocate(local_target(local_size))
  allocate(local_coef(local_size))

  counter = 0

  do i = 1, size(global_index)
     res = binary_search(sorted_index, global_index(i))
     if (res > 0) then
        counter = counter + 1
        local_index(counter) = global_index(i)
        local_target(counter) = global_target(i)
        local_coef(counter)  = global_coef(i)
     end if
  end do
     
  !do i = 1, size(local_index)
  !   write(100+jml_GetMyrankGlobal(), *) "local_mapping ", local_index(i), local_target(i), local_coef(i)
  !end do

  deallocate(sorted_index)
  
end subroutine make_local_mapping_table_no_sort

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine make_local_mapping_table_no_sort_org(global_index, global_target, global_coef, &
                                    grid_index, &
                                    local_size, local_index,  local_target,  local_coef)
  use jlt_utils, only : sort_int_1d, binary_search, error
  use jlt_mpi_lib, only : jml_GetMyrankGlobal
  implicit none
  integer, intent(IN)      :: global_index(:), global_target(:)
  real(kind=8), intent(IN) :: global_coef(:)
  integer, intent(IN)      :: grid_index(:)
  integer, intent(OUT)     :: local_size
  integer, pointer         :: local_index(:), local_target(:)
  real(kind=8), pointer    :: local_coef(:)
  integer, pointer         :: sorted_index(:), sorted_pos(:)
  integer, pointer         :: sorted_target(:)
  real(kind=8), pointer    :: sorted_coef(:)
  integer, pointer         :: sorted_grid_index(:)
  integer                  :: table_size
  integer                  :: current_pos, counter
  integer                  :: res
  integer                  :: i

  do i = 1, size(global_index)
     write(100+jml_GetMyrankGlobal(), *) "global_mapping ", global_index(i), global_target(i), global_coef(i)
  end do
  
  do i = 1, size(grid_index)
     write(100+jml_GetMyrankGlobal(), *) "grid_index ", grid_index(i)
  end do
  
  table_size = size(global_index)

  allocate(sorted_index(table_size))
  allocate(sorted_pos(table_size))
  allocate(sorted_target(table_size))
  allocate(sorted_coef(table_size))
  
  do i = 1, table_size
     sorted_index(i) = global_index(i)
     sorted_pos(i)   = i
  end do

  call sort_int_1d(table_size, sorted_index, sorted_pos)

  do i = 1, table_size
     sorted_target(i) = global_target(sorted_pos(i))
     sorted_coef(i)   = global_coef(sorted_pos(i))
  end do
  
  allocate(sorted_grid_index(size(grid_index)))
  sorted_grid_index(:) = grid_index(:)
  call sort_int_1d(size(grid_index), sorted_grid_index)

  counter = 0
  current_pos = 1
  do i = 1, size(grid_index)
     if (current_pos > size(sorted_index)) exit
     do while (sorted_index(current_pos) <= sorted_grid_index(i))
        if (sorted_index(current_pos) == sorted_grid_index(i)) counter = counter + 1
        current_pos = current_pos + 1
        if (current_pos > size(sorted_index)) exit
     end do
  end do

  local_size = counter

  if (local_size == 0) then ! no local mapping table
     allocate(local_index(0))
     allocate(local_target(0))
     allocate(local_coef(0))
     deallocate(sorted_index)
     deallocate(sorted_pos)
     deallocate(sorted_grid_index)
     return
  end if
  
  allocate(local_index(local_size))
  allocate(local_target(local_size))
  allocate(local_coef(local_size))

  counter = 0
  current_pos = 1
  do i = 1, size(grid_index)
     if (current_pos > size(sorted_index)) exit
     do while (sorted_index(current_pos) <= sorted_grid_index(i))
        if (sorted_index(current_pos) == sorted_grid_index(i)) then
           counter = counter + 1
           local_index(counter)  = sorted_index(sorted_pos(current_pos))
           local_target(counter) = sorted_target(sorted_pos(current_pos))
           local_coef(counter)   = sorted_coef(sorted_pos(current_pos))
        end if
        current_pos = current_pos + 1
        if (current_pos > size(sorted_index)) exit
     end do
  end do

  do i = 1, size(local_index)
     write(100+jml_GetMyrankGlobal(), *) "local_mapping ", local_index(i), local_target(i), local_coef(i)
  end do
  

  do i = 1, size(local_index)
     res = binary_search(sorted_grid_index, local_index(i))
     if (res <= 0) then
        write(0, *) "grid_index = ", i, local_index(i), sorted_grid_index
        call error("jlt_remapping:make_local_mapping_table", "grid index not listed in mapping table")
     end if
  end do

  deallocate(sorted_index)
  deallocate(sorted_pos)
  deallocate(sorted_grid_index)
  
end subroutine make_local_mapping_table_no_sort_org

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
  integer, allocatable, intent(OUT) :: local_index(:), local_target(:)
  real(kind=8), allocatable, intent(OUT) :: local_coef(:)
  integer, allocatable     :: sorted_index(:), sorted_pos(:)
  integer, allocatable     :: sorted_grid_index(:)
  integer                  :: table_size
  integer                  :: current_pos, counter
  integer                  :: res
  integer                  :: i

  !do i = 1, 10
  !   write(0, *) "global_mapping ", global_index(i), global_target(i), global_coef(i)
  !end do
  
  table_size = size(global_index)

  allocate(sorted_index(table_size))
  allocate(sorted_pos(table_size))

  do i = 1, table_size
     sorted_index(i) = global_index(i)
     sorted_pos(i)   = i
  end do

  call sort_int_1d(table_size, sorted_index, sorted_pos)

  allocate(sorted_grid_index(size(grid_index)))
  sorted_grid_index(:) = grid_index(:)
  call sort_int_1d(size(grid_index), sorted_grid_index)

  counter = 0
  current_pos = 1
  do i = 1, size(grid_index)
     if (current_pos > size(sorted_index)) exit
     do while (sorted_index(current_pos) <= sorted_grid_index(i))
        if (sorted_index(current_pos) == sorted_grid_index(i)) counter = counter + 1
        current_pos = current_pos + 1
        if (current_pos > size(sorted_index)) exit
     end do
  end do

  local_size = counter

  if (local_size == 0) then ! no local mapping table
     allocate(local_index(0))
     allocate(local_target(0))
     allocate(local_coef(0))
     deallocate(sorted_index)
     deallocate(sorted_pos)
     deallocate(sorted_grid_index)
     return
  end if
  
  allocate(local_index(local_size))
  allocate(local_target(local_size))
  allocate(local_coef(local_size))

  counter = 0
  current_pos = 1
  do i = 1, size(grid_index)
     if (current_pos > size(sorted_index)) exit
     do while (sorted_index(current_pos) <= sorted_grid_index(i))
        if (sorted_index(current_pos) == sorted_grid_index(i)) then
           counter = counter + 1
           local_index(counter)  = sorted_index(current_pos)
           local_target(counter) = global_target(sorted_pos(current_pos))
           local_coef(counter)   = global_coef(sorted_pos(current_pos))
        end if
        current_pos = current_pos + 1
        if (current_pos > size(sorted_index)) exit
     end do
  end do

  do i = 1, size(local_index)
     res = binary_search(sorted_grid_index, local_index(i))
     if (res <= 0) then
        write(0, *) "grid_index = ", i, local_index(i)
        call error("jlt_remapping:make_local_mapping_table", "grid index not listed in mapping table")
     end if
  end do

  !do i = 1, 10
  !   write(0, *) "local_mapping ", local_index(i), local_target(i), local_coef(i)
  !end do
  
  deallocate(sorted_index)
  deallocate(sorted_pos)
  deallocate(sorted_grid_index)
  
end subroutine make_local_mapping_table

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine reorder_index_by_target_rank(send_index, recv_index, coef, target_rank)
  use jlt_utils, only : sort_int_1d
  implicit none
  integer, intent(INOUT)      :: send_index(:)
  integer, intent(INOUT)      :: recv_index(:)
  real(kind=8), intent(INOUT) :: coef(:)
  integer, intent(INOUT)      :: target_rank(:)
  integer, allocatable        :: sorted_pos(:)
  integer, allocatable        :: int_temp(:)
  real(kind=8), allocatable   :: real_temp(:)
  integer                     :: table_size
  integer                     :: i
  
  table_size = size(target_rank)

  if (table_size <= 0) return
  
  allocate(sorted_pos(table_size))
  allocate(int_temp(table_size))
  allocate(real_temp(table_size))
  
  do i = 1, table_size
     sorted_pos(i) = i
  end do
  
  call sort_int_1d(table_size, target_rank, sorted_pos)

  int_temp(:) = send_index(:)
  do i = 1, table_size
     send_index(i) = int_temp(sorted_pos(i))
  end do
  
  int_temp(:) = recv_index(:)
  do i = 1, table_size
     recv_index(i) = int_temp(sorted_pos(i))
  end do

  real_temp(:) = coef(:)
  do i = 1, table_size
     coef(i) = real_temp(sorted_pos(i))
  end do

  deallocate(sorted_pos)
  deallocate(int_temp)
  deallocate(real_temp)

end subroutine reorder_index_by_target_rank

!=======+=========+=========+=========+=========+=========+=========+=========+
subroutine make_exchange_table(target_rank, local_table, num_of_rank, exchange_rank, num_of_index, offset, exchange_index)
  use jlt_utils, only : sort_int_1d
  implicit none
  integer, intent(IN)  :: target_rank(:)        ! ordered target rank array
  integer, intent(IN)  :: local_table(:)        ! local remapping table of exchange grid
  integer, intent(OUT) :: num_of_rank           ! different rank
  integer, allocatable, intent(OUT) :: exchange_rank(:)      ! exchange rank array
  integer, allocatable, intent(OUT) :: num_of_index(:)       ! num of exchange index by each rank
  integer, allocatable, intent(OUT) :: offset(:)             ! offset address
  integer, allocatable, intent(OUT) :: exchange_index(:)     ! exchange index array
  integer, allocatable :: num_of_all_index(:)   ! num of exchange index by each rank
  integer, allocatable :: temp_index(:)
  integer, allocatable :: same_rank_index(:)
  integer, allocatable :: different_index_array(:)
  integer              :: is, ie, counter 
  integer              :: i, j

  if (size(target_rank) <= 0) then
     num_of_rank = 0
     return
  end if
  
  allocate(temp_index(size(local_table)))
  
  num_of_rank = 1
  do i = 2, size(target_rank)
     if (target_rank(i) /= target_rank(i-1)) num_of_rank = num_of_rank + 1
  end do

  allocate(exchange_rank(num_of_rank))
  allocate(num_of_all_index(num_of_rank))
  allocate(num_of_index(num_of_rank))
  allocate(offset(num_of_rank))
  
  counter = 1
  num_of_rank = 1
  exchange_rank(1) = target_rank(1)
  do i = 2, size(target_rank)
     if (target_rank(i) /= target_rank(i-1)) then
        num_of_all_index(num_of_rank) = counter
        num_of_rank = num_of_rank + 1
        exchange_rank(num_of_rank) = target_rank(i)
        counter = 1
     else
        counter = counter + 1
     end if
  end do
  num_of_all_index(num_of_rank) = counter

  counter = 0
  is = 0
  do i = 1, num_of_rank
     allocate(same_rank_index(num_of_all_index(i)))
     do j = 1, num_of_all_index(i)
        same_rank_index(j) = local_table(counter + j)
     end do
     call delete_same_index(same_rank_index, different_index_array)
     num_of_index(i) = size(different_index_array)
     do j = 1, num_of_index(i)
        temp_index(is + j) = different_index_array(j)
     end do
     is = is + num_of_index(i)
     counter = counter + num_of_all_index(i)
     deallocate(different_index_array)
     deallocate(same_rank_index)
  end do
  
  counter = 0
  do i = 1, num_of_rank
     counter = counter + num_of_index(i)
  end do

  offset(1) = 0
  do i = 2, num_of_rank
     offset(i) = offset(i-1) + num_of_index(i-1)
  end do
  
  allocate(exchange_index(counter))
  exchange_index(1:counter) = temp_index(1:counter)

  deallocate(temp_index)
  deallocate(num_of_all_index)
  
end subroutine make_exchange_table

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine delete_same_index(index_in, index_out)
  use jlt_utils, only : sort_int_1d
  implicit none
  integer, intent(IN)  :: index_in(:)
  integer, allocatable, intent(OUT) :: index_out(:)
  integer, allocatable :: sorted_index(:)
  integer              :: table_size
  integer              :: i, counter
  
  table_size = size(index_in)

  if (table_size <= 0) then
     allocate(index_out(0))
     return
  end if
  
  allocate(sorted_index(table_size))

  sorted_index(:) = index_in(:)
  call sort_int_1d(table_size, sorted_index)

  counter = 1
  do i = 2, table_size
     if (sorted_index(i) /= sorted_index(i-1)) counter = counter + 1 ! count deffernt index number
  end do

  allocate(index_out(counter))

  counter = 1
  index_out(1) = sorted_index(1)
  do i = 2, table_size
     if (sorted_index(i) /= sorted_index(i-1)) then
        counter = counter + 1
        index_out(counter) = sorted_index(i)
     end if
  end do

  deallocate(sorted_index)
  
end subroutine delete_same_index

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine make_conversion_table(remapping_index, grid_index, conv_table)
  use jlt_utils, only : sort_int_1d, binary_search
  use mpi
  implicit none
  integer, intent(IN) :: remapping_index(:)  ! local remapping table
  integer, intent(IN) :: grid_index(:)       ! data grid index
  integer, allocatable, intent(OUT) :: conv_table(:)       ! conversion table from interpolation index i to data grid index
  integer, allocatable :: sorted_index(:)
  integer, allocatable :: sorted_pos(:)
  integer             :: res
  integer             :: i
  integer :: my_rank, ierr

  call mpi_comm_rank(MPI_COMM_WORLD, my_rank, ierr)

  !write(300 + my_rank, *) "make_conversion_table"
  !do i = 1, size(remapping_index)
  !   write(300 + my_rank, *) remapping_index(i)
  !end do
  !do i = 1, size(grid_index)
  !   write(300 + my_rank, *) grid_index(i)
  !end do
  
  allocate(sorted_index(size(grid_index)))
  allocate(sorted_pos(size(grid_index)))

  do i = 1, size(grid_index)
     sorted_index(i) = grid_index(i)
     sorted_pos(i) = i
  end do
  
  call sort_int_1d(size(grid_index), sorted_index, sorted_pos)

  allocate(conv_table(size(remapping_index)))
  
  do i = 1, size(remapping_index)
     res = binary_search(sorted_index, remapping_index(i))
     conv_table(i) = sorted_pos(res)
     !write(300 + my_rank, *) conv_table(i)
  end do

  deallocate(sorted_index)
  deallocate(sorted_pos)
  
end subroutine make_conversion_table

!=======+=========+=========+=========+=========+=========+=========+=========+

end module jlt_remapping

