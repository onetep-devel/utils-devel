! Description: Tool to calculate electron density difference for two
!              .cube files.
!              cubeout = cubein1 - cubein2
!
! Written by Max Phipps, 2016.
! Code derived from visual_output_cube.
!
!---------------------------------------------------------------------

module constants

  integer, parameter :: dp = kind(1.d0)

end module constants


program edd_cube
  !---------------------------------------------------------------------

  use constants, only: dp

  implicit none

  interface

     subroutine read_vol(funit,np_in,vdat)
        use constants, only: dp
        integer, intent(in)          :: funit
        integer, intent(in)          :: np_in(3)
        real(kind=DP), intent(inout) :: vdat(:,:,:)
     end subroutine read_vol

     subroutine write_vol(funit,np_in,vdat)
        use constants, only: dp
        integer, intent(in)          :: funit
        integer, intent(in)          :: np_in(3)
        real(kind=DP), intent(in)    :: vdat(:,:,:)
     end subroutine write_vol

  end interface

  character(len=256) :: fn_cube1in, fn_cube2in, fn_cubeout ! file names
  integer :: funit_cube1in, funit_cube2in, funit_cubeout ! file units
  integer :: ii, ij, ik ! counters
  integer :: ierr
  ! file io:
  character(len=256) :: line_a
  character(len=256) :: line_b
  integer :: nat, nat_check
  integer :: np(3), np_check ! number of grid points check
  ! volumetric data:
  real(kind=DP), allocatable, dimension(:,:,:) :: vdat1, vdat2, vdat_edd

  ! mjsp: check correct number of arguments
  if (iargc() .ne. 3) then
     print*, "ERROR: 3 arguments are required!"
     print*, "  Usage: ./edd_cube cubein1 cubein2 cubeout"
     print*, "  where n(cubeout) = n(cubein1) - n(cubein2)"
     goto 101
  end if

  ! mjsp: initialisation
  funit_cube1in = 7
  funit_cube2in = 8
  funit_cubeout = 9

  ! mjsp: get filenames
  call getarg(1,fn_cube1in)
  call getarg(2,fn_cube2in)
  call getarg(3,fn_cubeout)

  ! mjsp: open .cube files
  open(funit_cube1in,file=fn_cube1in, status='old', &
         access='sequential', form='formatted', action='read', iostat=ierr)
  if (ierr.ne.0) then
     print*,'Error opening file',fn_cube1in,'. iostat=',ierr
     goto 100
  end if

  open(funit_cube2in,file=fn_cube2in, status='old', &
         access='sequential', form='formatted', action='read', iostat=ierr)
  if (ierr.ne.0) then
     print*,'Error opening file',fn_cube2in,'. iostat=',ierr
     goto 99
  end if

  open(funit_cubeout,file=fn_cubeout, & ! status='new', &
         access='sequential', form='formatted', action='write', iostat=ierr)
  if (ierr.ne.0) then
     print*,'Error opening file',fn_cubeout,'. iostat=',ierr
     goto 98
  end if

  ! ======DATA CHECKING==============================

  ! mjsp: check header line 1 and 2
  do ii=1,2

    ! mjsp: read line
    read(funit_cube1in,'(a)') line_a
    read(funit_cube2in,'(a)') line_b

    ! mjsp: check lines are identical
    if (line_a .ne. line_b) then
      ! mjsp: Treat as warning
      print*,"WARNING: cube files' header does not match! (line ",ii,")"
!      goto 98  ! to quit (if raising error)
    end if

    ! mjsp: check passed, so write to output cube file
    write(funit_cubeout,'(a)') trim(line_a)

  end do

  ! mjsp: get number of atoms and perform check
  read(funit_cube1in,'(i4,a)') nat, line_a
  read(funit_cube2in,'(i4,a)') nat_check, line_b
  if (nat .ne. nat_check) then
    print*,"ERROR: cube files' atom counts do not match! (line 3)"
    goto 98
  end if
  if (line_a .ne. line_b) then
    print*,"ERROR: cube files' origin data does not match! (line 3)"
    goto 98
  end if
  write(funit_cubeout,'(i4,a)') nat, line_a(1:40)

  ! mjsp: check voxel (grid) data
  do ii=1,3
    read(funit_cube1in,'(i4,a)') np(ii),line_a
    read(funit_cube2in,'(i4,a)') np_check,line_b
    if ((np(ii) .ne. np_check) .or. (line_a .ne. line_b)) then
      print*,"ERROR: cube files' voxel data does not match! (line ",ii+3,")"
      goto 98
    end if
    write(funit_cubeout,'(i4,a)') np(ii), line_a(1:40)
  end do

  ! mjsp: check atoms and coordinates
  do ii=1,nat
    read(funit_cube1in,'(a)') line_a
    read(funit_cube2in,'(a)') line_b
    if (line_a .ne. line_b) then
      print*,"ERROR: cube files' atomic data does not match! (line ",ii+6,")"
      goto 98
    end if
    write(funit_cubeout,'(a)') line_a(1:57)
  end do

  ! ======END DATA CHECKING==========================

  ! mjsp: Allocations
  allocate(vdat1(np(1),np(2),np(3)),stat=ierr)
  if (ierr.ne.0) then
    print*,"ERROR: allocating vdat1 failed with iostat=",ierr
    goto 98
  end if

  allocate(vdat2(np(1),np(2),np(3)),stat=ierr)
  if (ierr.ne.0) then
    print*,"ERROR: allocating vdat2 failed with iostat=",ierr
    goto 97
  end if

  allocate(vdat_edd(np(1),np(2),np(3)),stat=ierr)
  if (ierr.ne.0) then
    print*,"ERROR: allocating vdat_edd failed with iostat=",ierr
    goto 96
  end if

  ! mjsp: read in volumetric data (cube #1)
  call read_vol(funit_cube1in, np, vdat1)

  ! mjsp: read in volumetric data (cube #2)
  call read_vol(funit_cube2in, np, vdat2)

  ! mjsp: calculate EDD via subtraction
  vdat_edd(:,:,:) = vdat1(:,:,:) - vdat2(:,:,:)

  ! mjsp: write out the EDD volumetric data
  call write_vol(funit_cubeout, np, vdat_edd)

  ! mjsp: deallocations/cleanup
  deallocate(vdat_edd,stat=ierr)
  if (ierr.ne.0) print*,"ERROR: deallocating vdat_edd failed with iostat=",ierr
  96 continue
  deallocate(vdat2,stat=ierr)
  if (ierr.ne.0) print*,"ERROR: deallocating vdat2 failed with iostat=",ierr
  97 continue
  deallocate(vdat1,stat=ierr)
  if (ierr.ne.0) print*,"ERROR: deallocating vdat1 failed with iostat=",ierr
  98 continue
  close(funit_cubeout)
  99 continue
  close(funit_cube2in)
  100 continue
  close(funit_cube1in)

  101 continue

end program edd_cube

subroutine read_vol(funit, np_in, vdat)

  ! reads in .cube format volumetric data,
  ! after the header/metadata, from file

  implicit none

  integer, parameter :: dp = kind(1.d0)

  integer, intent(in)          :: funit        ! file unit
  integer, intent(in)          :: np_in(3)     ! grid points array
  real(kind=DP), intent(inout) :: vdat(:,:,:)  ! volumetric data

  integer :: ii, ij, ik, ic ! counters
  integer :: remainder
  character(len=80) :: remchar   ! string buffer

  remainder = mod(np_in(3),6)
  write(remchar,*) remainder
  do ii=1,np_in(1)
     do ij=1,np_in(2)

        do ik=1,np_in(3)-remainder,6
           read(funit, '(6E13.5)') (vdat(ii, ij, ik+ic),&
                ic =0,5  )
        enddo

        if (remainder.gt.0) then
           read(funit, '('//trim(adjustl(remchar))//'E13.5)') &
                (vdat(ii, ij, np_in(3) -remainder +ic), ic =1, &
                remainder  )
        endif

     end do ! y
  end do ! x

end subroutine read_vol

subroutine write_vol(funit, np_in, vdat)

  ! writes out .cube format volumetric data,
  ! after the header/metadata, to file

  implicit none

  integer, parameter :: dp = kind(1.d0)

  integer, intent(in)          :: funit        ! file unit
  integer, intent(in)          :: np_in(3)     ! grid points array
  real(kind=DP), intent(in)    :: vdat(:,:,:)  ! volumetric data

  integer :: ii, ij, ik, ic ! counters
  integer :: remainder
  character(len=80) :: remchar   ! string buffer

  remainder = mod(np_in(3),6)
  write(remchar,*) remainder
  do ii=1,np_in(1)
     do ij=1,np_in(2)

        do ik=1,np_in(3)-remainder,6
           write(funit, '(6E13.5)') (vdat(ii, ij, ik+ic),&
                ic =0,5  )
        enddo

        if (remainder.gt.0) then
           write(funit, '('//trim(adjustl(remchar))//'E13.5)') &
                (vdat(ii, ij, np_in(3) -remainder +ic), ic =1, &
                remainder  )
        endif

     end do ! y
  end do ! x

end subroutine write_vol
