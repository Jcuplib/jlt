!=-----------------------------------------------------------------------------=
! Time measure module using system_clock routine
!                                                         harayama@rist.or.jp 
!=-----------------------------------------------------------------------------=
module palmtime
  use mpi
  ! *** Custom TIMER Routine with system_clock routine ***
  implicit none
  !
  private

!--------------------------------   public  ----------------------------------!

  public  :: palm_TimeInit        ! subroutine (compk_name, commnicator)
  public  :: palm_TimeFinalize    ! subroutine () 
  public  :: palm_TimeStart       ! subroutine (block_name)
  public  :: palm_TimeEnd         ! subroutine (block_name)
  public  :: palm_TimeStart_w_barrier
  public  :: palm_TimeEnd_w_barrier
  !

!--------------------------------  private  ----------------------------------!

  private :: palm_TimeOutput
  private :: palm_TimeError
  !
  integer, save :: palm_comm_world, palm_myrank, palm_nprocs, ierr
  !
  integer, parameter :: max_records       = 500
  integer, parameter :: max_string_length =  32
  integer, parameter :: ounit = 6
  !
  integer, public, save :: icount_rate, icount_max, record_count = 0
  character(len=max_string_length), public, save :: code_name
  integer :: icount
  !
  real(8), allocatable, save :: ratio(:)
  !
  type time_records
    character(len=max_string_length) :: block_name
    integer :: start_count
    integer :: end_count
    logical :: checked
    integer :: calls
    real(8) :: total_time
    real(8) :: mean_time
    real(8) :: max_time
    real(8) :: min_time
    logical :: ofirst
  end type time_records
  type(time_records), dimension(max_records), public, save :: records
  !
  contains


!=======+=========+=========+=========+=========+=========+=========+=========+
!>

    subroutine palm_TimeInit(block,comm)
      implicit none
      character(len=*), optional :: block  
      integer, optional :: comm
      !
      if ( present(comm) ) then
        palm_comm_world = comm
      else
        call mpi_comm_dup(mpi_comm_world,palm_comm_world,ierr)
      end if
      call mpi_comm_size(palm_comm_world,palm_nprocs,ierr)
      ! initialize
      records(:)%block_name   = ''
      records(:)%start_count  = 0
      records(:)%end_count    = 0
      records(:)%checked      = .false.
      records(:)%calls        = 0
      records(:)%total_time   = 0.d0
      records(:)%mean_time    = 0.d0
      records(:)%max_time     = 0.d0
      records(:)%min_time     = 0.d0
      records(:)%ofirst       = .true.
      !
      code_name = 'Unknown'
      if ( present(block) ) code_name = block
      call system_clock(count_rate=icount_rate, count_max=icount_max)
      call system_clock(count=icount)
      !
      records(1)%block_name  = code_name
      records(1)%start_count = icount
      !
      record_count = record_count +1
      !
    end subroutine palm_TimeInit
    !
!=======+=========+=========+=========+=========+=========+=========+=========+
!>

    subroutine palm_TimeFinalize
      implicit none
      integer :: icount1, icount2, icount, i
      !
      call system_clock(count=icount2)
      records(1)%end_count = icount2
      icount1 = records(1)%start_count
      !
      records(1)%checked = .true.
      records(1)%calls   = 1
      icount = icount2 - icount1
      if ( icount2 < icount1 ) icount = icount + icount_max
      records(1)%total_time = dble(icount)/icount_rate
      records(1)%max_time   = records(1)%total_time
      records(1)%min_time   = records(1)%total_time
      records(1)%ofirst     = .false.
      !
      allocate(ratio(record_count))
      do i = 1, record_count
        ratio(i) = records(i)%total_time/records(1)%total_time * 100.d0
        records(i)%mean_time = records(i)%total_time/records(i)%calls
      enddo
      !
      call palm_TimeOutput
      !
    end subroutine palm_TimeFinalize

!=======+=========+=========+=========+=========+=========+=========+=========+
!>
    !
    subroutine palm_TimeStart(block_name)
      implicit none
      character(len=*), intent(in) :: block_name
      !
      integer :: i, iseek, icount
      !
      ! search block_name
      iseek = 0
      do i = 1, record_count
        if ( records(i)%block_name == block_name ) then 
          iseek = i
          exit
        end if
      enddo 
      !
      if ( iseek == 0 ) then
        record_count = record_count +1
        records(record_count)%block_name  = block_name
        call system_clock(count=icount) 
        records(record_count)%start_count = icount 
      else
        if ( records(iseek)%checked ) then
          records(iseek)%checked   = .false.
          records(iseek)%end_count = 0
          call system_clock(count=icount)
          records(iseek)%start_count = icount
        else
          call palm_TimeError(block_name,1)
        end if
      end if 
      !
    end subroutine palm_TimeStart
    !
!=======+=========+=========+=========+=========+=========+=========+=========+
!>
    subroutine palm_TimeEnd(block_name)
      implicit none
      character(len=*), intent(in) :: block_name
      !
      integer :: i, iseek, icount, icount1, icount2
      real(8) :: elapsed_time, min_time, max_time
      !
      call system_clock(count=icount)
      !
      ! search block_name
      iseek = 0
      do i = 1, record_count
        if ( records(i)%block_name == block_name ) then
          iseek = i
          exit
        end if
      enddo
      !
      if ( iseek == 0 ) then
        call palm_TimeError(block_name,2)
      else
        records(iseek)%end_count = icount
        records(iseek)%checked   = .true.

        ! statistics
        records(iseek)%calls      = records(iseek)%calls +1
        icount1 = records(iseek)%start_count 
        icount2 = records(iseek)%end_count 
        icount  = icount2 - icount1
        if ( icount2 < icount1 ) icount = icount + icount_max
        elapsed_time = dble(icount)/icount_rate
        records(iseek)%total_time = records(iseek)%total_time + elapsed_time 
        ! 
        if ( records(iseek)%ofirst ) then
          records(iseek)%max_time = elapsed_time  
          records(iseek)%min_time = elapsed_time  
          records(iseek)%ofirst = .false.
        else
          if ( records(iseek)%max_time < elapsed_time ) records(iseek)%max_time = elapsed_time
          if ( records(iseek)%min_time > elapsed_time ) records(iseek)%min_time = elapsed_time
        end if
        !
      end if
      !
    end subroutine palm_TimeEnd
    !
!=======+=========+=========+=========+=========+=========+=========+=========+
!>

    subroutine palm_TimeStart_w_barrier(block_name)
      implicit none
      character(len=*), intent(in) :: block_name
      !
      character(len=max_string_length) :: sync_block_name
      !
      sync_block_name = 'sync_'//trim(block_name) 
      call palm_TimeStart(sync_block_name)
        call mpi_barrier(palm_comm_world,ierr)
      call palm_TimeEnd (sync_block_name)
      call palm_TimeStart(block_name)
      !
    end subroutine palm_TimeStart_w_barrier
    !
!=======+=========+=========+=========+=========+=========+=========+=========+
!>

    subroutine palm_TimeEnd_w_barrier(block_name)
      implicit none
      character(len=*), intent(in) :: block_name
      !
      character(len=max_string_length) :: sync_block_name
      !
      sync_block_name = 'sync_'//trim(block_name)
      call palm_TimeStart(sync_block_name)
        call mpi_barrier(palm_comm_world,ierr)
      call palm_TimeEnd (sync_block_name)
      call palm_TimeEnd(block_name)
      !
    end subroutine palm_TimeEnd_w_barrier
    !
!=======+=========+=========+=========+=========+=========+=========+=========+
!>
    subroutine palm_TimeError(block_name,start_or_stop)
      implicit none
      character(len=*), intent(in) :: block_name
      integer, intent(in) :: start_or_stop
      !
      select case(start_or_stop)
        case(1)
          write(ounit,'("[error] invalid palm_TimeStart ---> ",a)') trim(adjustl(block_name))
!          call exit(1)
        case(2)
          write(ounit,'("[error] invalid palm_TimeEnd  ---> ",a)') trim(adjustl(block_name))
!          call exit(2)
      end select
      !
    end subroutine palm_TimeError
    !
!=======+=========+=========+=========+=========+=========+=========+=========+
!>
    subroutine palm_TimeOutput
      implicit none
      character(len=110) :: string   = ''
      character(len= 60) :: fmt_head = ''
      integer i, j
      !
      character(len=255) :: filename = ''
      character(len= 5) :: cid = '00001'
      integer :: ounit
      logical :: isOpen
      !
      real(8), allocatable :: sbuff(:)
      real(8), allocatable :: rbuff(:,:), rbuff_max(:), rbuff_min(:), rbuff_sum(:)
      !
      integer, parameter :: idigits = 13
      character(len=7) :: fmt_head2 = '(i13)'
      character(len=idigits*100) :: string2 = ''
      character(len=idigits) :: string3(3) =(/'Max. ','Min. ','Mean.'/)
      character(len=idigits) :: string4
      integer :: sindx = 32
      !
      call mpi_comm_rank(palm_comm_world,palm_myrank,ierr)
      write(cid,'(i5.5)') palm_myrank
!      filename( 1: 5) = 'etime'
!      filename( 6:10) = cid
!      filename(11:14) = '.txt'
!
!      filename( 1: 4) = 'palm'
!      filename( 5: 9) = cid
!      filename(10:13) = '.txt'
!
      filename = 'TM.'//trim(code_name)//'.pe'//cid
      !
      do i = 10, 100
        inquire(i,opened=isOpen)
        if (.not. isOpen) then
          ounit = i
          exit
        end if
      enddo
      !
      fmt_head = '(i3,1x,a32,1x,i7,1x,f13.6,1x,"(",f6.2,"%)",3(1x,f13.6))'
      string(  1:  3) = 'No.'
      string(  5: 14) = 'Procedure'
      string( 39: 45) = 'Count'
      string( 46: 59) = 'Elapsed(sec.)'
      string( 60: 68) = 'Ratio(%)'
      string( 70: 75) = 'Mean.'
      string( 84: 88) = 'Max.'
      string( 98:102) = 'Min.'
      !
      open(ounit,file=filename,status='replace')
      write(ounit,'(a)') string
      write(ounit,'(a)') repeat('-',110)
      do i = 1, record_count
        if ( i ==2 ) write(ounit,'(a)') repeat('-',110)
        write(ounit,fmt_head) &
         i-1, &
         records(i)%block_name, &
         records(i)%calls,      &
         records(i)%total_time, ratio(i), &
         records(i)%mean_time, records(i)%max_time, records(i)%min_time
      enddo
      close(ounit)
      !
      allocate(sbuff(record_count))

      if (palm_myrank == 0) then
        allocate(rbuff(record_count,palm_nprocs))
      else
        allocate(rbuff(1,1))
      end if

      do i = 1, record_count
        sbuff(i) = records(i)%total_time
      enddo
      !
      ! collect data
      call mpi_gather(sbuff,record_count,mpi_real8, &
                      rbuff,record_count,mpi_real8,0,palm_comm_world,ierr)
      ! 
      if (palm_myrank == 0) then
        allocate(rbuff_max(record_count))
        allocate(rbuff_min(record_count))
        allocate(rbuff_sum(record_count))
      else
        allocate(rbuff_max(1))
        allocate(rbuff_min(1))
        allocate(rbuff_sum(1))
      end if

      call mpi_reduce(sbuff,rbuff_max,record_count,mpi_real8,mpi_max,0,palm_comm_world,ierr)
      call mpi_reduce(sbuff,rbuff_min,record_count,mpi_real8,mpi_min,0,palm_comm_world,ierr)
      call mpi_reduce(sbuff,rbuff_sum,record_count,mpi_real8,mpi_sum,0,palm_comm_world,ierr)
      ! 
      if (palm_myrank == 0) then
        !
!        cid = 'AAAAA'
!        filename(6:10) = cid
        filename = 'TM.'//trim(code_name)
        !
        string2(1:sindx) = 'Procedure'
        sindx = sindx +1
        do i = 1, 3
          string2(sindx:) = string3(i)
          sindx = sindx +idigits
        enddo
!        do i = 0, palm_nprocs-1
!          write(string4,fmt_head2) i
!          string2(sindx:) = string4
!          sindx = sindx +idigits
!        enddo
        !
        open(ounit,file=filename,status='replace')
        write(ounit,'(a)') trim(string2)
        do i = 1, record_count
          write(ounit,'(a32)',advance='no') records(i)%block_name
          write(ounit,'(f13.6)',advance='no') rbuff_max(i)
          write(ounit,'(f13.6)',advance='no') rbuff_min(i)
          write(ounit,'(f13.6)',advance='no') rbuff_sum(i)/palm_nprocs
          !write(ounit,'(f13.6)',advance='no') rbuff_sum(i)/palm_nprocs
          do j = 1, palm_nprocs-1
            write(ounit,'(f13.6)',advance='no') rbuff(i,j)
          enddo
          write(ounit,'(f13.6)') rbuff(i,palm_nprocs)
        enddo
        close(ounit)
      end if
      ! 
      !!call mpi_comm_free(palm_comm_world,ierr)
      !
      deallocate(rbuff)
      deallocate(rbuff_max, rbuff_min, rbuff_sum)

    end subroutine palm_TimeOutput
    !
!=======+=========+=========+=========+=========+=========+=========+=========+
!>
    subroutine del_spaces(string)
      implicit none
      character(len=*),intent(inout) :: string
      !
      character(len=len(string)) :: tmp
      integer i, j
      !
      j = 1
      do i = 1, len(string)
        if ( string(i:i) == ' ' ) cycle
        tmp(j:j) = string(i:i)
        j = j + 1
      enddo
      string = ''
      string = tmp(1:j-1) 
      !
    end subroutine del_spaces
    !
!=======+=========+=========+=========+=========+=========+=========+=========+
!>
end module palmtime
