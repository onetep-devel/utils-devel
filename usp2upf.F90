! PLEASE NOTE:
! THE UPF MODULE USED IN THIS CODE IS GPL'd, SO WOULD NOT BE ABLE TO BE
! SOLD BY ACCELRYS. THIS FILE CAN BE DISTRIBUTED ON ONETEP SITE THOUGH
! AND INCLUDED IN THE ACADEMIC VERSION AS FAR AS I KNOW.
! ndmh
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
program usp2upf
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in the CASTEP .usp form
  !     to unified pseudopotential format
  !
  implicit none
  character(len=256) filein, fileout
  !
  !
  call get_file ( filein )
  open (unit = 1, file = filein, status = 'old', form = 'formatted')
  call scan_usp(1)
  call allocate_usp
  rewind (1)
  call read_usp(1)
  close (1)

  ! convert variables read from user-supplied format into those needed
  ! by the upf format - add missing quantities

  call convert_usp

  fileout=trim(filein)//'.UPF'
  print '(''Output PP file in UPF format :  '',a)', fileout

  open(unit=2,file=fileout,status='unknown',form='formatted')
  call write_upf(2)
  close (unit=2)

stop
20 call errore ('psp2upf', 'Reading pseudo file name ', 1)

end program usp2upf
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module upf
  !
  ! All variables to be written into the UPF file
  ! (UPF = unified pseudopotential format)
  !
  ! pp_info
  integer,parameter :: DP8 = kind(1.0d0)
  integer :: rel
  real(DP8) :: rcloc
  integer :: nwfs
  real(DP8), allocatable :: oc(:), rcut(:), rcutus(:), epseu(:)
  character(len=2), allocatable :: els(:)
  integer, allocatable:: lchi (:), nns (:)
  !
  ! pp_header
  character (len=80):: generated, date_author, comment
  character (len=2) :: psd, pseudotype
  integer :: nv = 0
  integer :: iexch, icorr, igcx, igcc
  integer :: lmax, mesh, nbeta, ntwfc
  logical :: nlcc
  real(DP8) :: zp, ecutrho, ecutwfc, etotps
  real(DP8), allocatable :: ocw(:)
  character(len=2), allocatable :: elsw(:)
  integer, allocatable:: lchiw(:)
  !
  ! pp_mesh
  real(DP8), allocatable :: r(:), rab(:)
  !
  ! pp_nlcc
  real(DP8), allocatable :: rho_atc(:)
  !
  ! pp_local
  real(DP8), allocatable ::  vloc0(:)
  !
  ! pp_nonlocal
  ! pp_beta
  real(DP8), allocatable :: betar(:,:)
  integer, allocatable:: lll(:), ikk2(:)  
  ! pp_dij
  real(DP8), allocatable :: dion(:,:)
  ! pp_qij
  integer ::  nqf, nqlc
  real(DP8), allocatable :: rinner(:), qqq(:,:), qfunc(:,:,:)
  ! pp_qfcoef
  real(DP8), allocatable :: qfcoef(:,:,:,:)
  !
  ! pp_pswfc
  real(DP8), allocatable :: chi(:,:)
  !
  ! pp_rhoatom
  real(DP8), allocatable :: rho_at(:)
end module upf
!
subroutine write_upf(ounps)

  use upf, only: nlcc

  integer :: ounps

  call write_pseudo_comment(ounps)  
  call write_pseudo_header(ounps)  
  call write_pseudo_mesh(ounps)
  if (nlcc)  call write_pseudo_nlcc(ounps)  
  call write_pseudo_local(ounps)  
  call write_pseudo_nl(ounps)  
  call write_pseudo_pswfc(ounps)
  call write_pseudo_rhoatom(ounps)  
  !
  print '("*** PLEASE TEST BEFORE USING!!! ***")'
  print '("review the content of the PP_INFO fields")'
  !
end subroutine write_upf

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_comment (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the comments of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  

    integer :: nb, ios  

    write (ounps, '(a9)', err = 100, iostat = ios) "<PP_INFO>"  

    write (ounps, '(a)', err = 100, iostat = ios) generated
    write (ounps, '(a)', err = 100, iostat = ios) date_author
    write (ounps, '(a)', err = 100, iostat = ios) comment
    if (rel==2) then  
       write (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel,& 
            &"The Pseudo was generated with a Full-Relativistic Calculation"
    else if (rel==1) then  
       write (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel,& 
            &"The Pseudo was generated with a Scalar-Relativistic Calculation"
    else if (rel==0) then 
       write (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel, &
            & "The Pseudo was generated with a Non-Relativistic Calculation"
    endif

    if (rcloc > 0.d0) &
       write (ounps, '(1pe19.11,t24,a)', err = 100, iostat = ios) &
              rcloc, "Local Potential cutoff radius"

    if (nwfs>0) &
       write (ounps, '(a2,2a3,a6,3a19)', err = 100, iostat = ios) "nl", &
              &" pn", "l", "occ", "Rcut", "Rcut US", "E pseu"
    do nb = 1, nwfs  
       write (ounps, '(a2,2i3,f6.2,3f19.11)') els (nb) , nns (nb) , &
            lchi (nb) , oc (nb) , rcut (nb) , rcutus (nb) , epseu(nb)

    enddo

    write (ounps, '(a10)', err = 100, iostat = ios) "</PP_INFO>"  
    return
100 call errore ('write_pseudo_comment', 'Writing pseudo file', abs ( &
         ios))   
  end subroutine write_pseudo_comment

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_header (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the header of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    character (len=4) :: shortname
    character (len=20):: dft  
    integer :: nb, ios  
    !
    !
    write (ounps, '(//a11)', err = 100, iostat = ios) "<PP_HEADER>"  

    write (ounps, '(t3,i2,t24,a)', err = 100, iostat = ios) nv, &
         "Version Number"
    write (ounps, '(t3,a,t24,a)', err = 100, iostat = ios) psd , &
         "Element"
    if (pseudotype == 'NC') then  
       write (ounps, '(a5,t24,a)', err = 100, iostat = ios) "NC", &
            "Norm - Conserving pseudopotential"
    else if (pseudotype == 'US') then
       write (ounps, '(a5,t24,a)', err = 100, iostat = ios) "US", &
            "Ultrasoft pseudopotential"
    else
       call errore ('write_pseudo_header',&
            'Unknown PP type: '//pseudotype, 1)
    endif
    write (ounps, '(l5,t24,a)', err = 100, iostat = ios) nlcc , &
         "Nonlinear Core Correction"
    call dftname (iexch, icorr, igcx, igcc, dft, shortname)
    write (ounps, '(a,t24,a4,a)', err = 100, iostat = ios) &
         dft, shortname," Exchange-Correlation functional"
    write (ounps, '(f17.11,t24,a)') zp , "Z valence"  
    write (ounps, '(f17.11,t24,a)') etotps, "Total energy"  
    write (ounps, '(2f11.7,t24,a)') ecutrho, ecutwfc, &
         "Suggested cutoff for wfc and rho"  

    write (ounps, '(i5,t24,a)') lmax, "Max angular momentum component"  
    write (ounps, '(i5,t24,a)') mesh, "Number of points in mesh"
    write (ounps, '(2i5,t24,a)', err = 100, iostat = ios) ntwfc, &
         nbeta  , "Number of Wavefunctions, Number of Projectors"
    write (ounps, '(a,t24,a2,a3,a6)', err = 100, iostat = ios) &
         " Wavefunctions", "nl", "l", "occ"
    do nb = 1, ntwfc
       write (ounps, '(t24,a2,i3,f6.2)') elsw(nb), lchiw(nb), ocw(nb)
    enddo
    !---> End header writing

    write (ounps, '(a12)', err = 100, iostat = ios) "</PP_HEADER>"
    return   
100 call errore ('write_pseudo_header','Writing pseudo file', abs(ios) )

  end subroutine write_pseudo_header

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_mesh (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the atomic charge density to the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  
    !
    write (ounps, '(//a9)', err = 100, iostat = ios) "<PP_MESH>"  

    write (ounps, '(t3,a6)', err = 100, iostat = ios) "<PP_R>"  
    write (ounps, '(1p4e19.11)', err=100, iostat=ios) (r(ir),  ir=1,mesh )
    write (ounps, '(t3,a7)', err = 100, iostat = ios) "</PP_R>"  
    write (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_RAB>"  
    write (ounps, '(1p4e19.11)', err=100, iostat=ios) (rab(ir), ir=1,mesh )
    write (ounps, '(t3,a9)', err = 100, iostat = ios) "</PP_RAB>"  

    write (ounps, '(a10)', err = 100, iostat = ios) "</PP_MESH>"  

    return

100 call errore ('write_pseudo_rhoatom','Writing pseudo file',abs(ios))

  end subroutine write_pseudo_mesh

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_nlcc (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the core charge for the nonlinear core
    !     correction of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  

    write (ounps, '(//a9)', err = 100, iostat = ios) "<PP_NLCC>"  

    write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                 ( rho_atc(ir), ir = 1, mesh )
    write (ounps, '(a10)', err = 100, iostat = ios) "</PP_NLCC>"  
    return

100 call errore ('write_pseudo_nlcc', 'Writing pseudo file', abs (ios))

  end subroutine write_pseudo_nlcc
  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_local (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the local part of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  

    write (ounps, '(//a10)', err = 100, iostat = ios) "<PP_LOCAL>"  
    write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                ( vloc0(ir), ir = 1, mesh )
    write (ounps, '(a11)', err = 100, iostat = ios) "</PP_LOCAL>"  
    return
100 call errore ('write_pseudo_local', 'Writing pseudo file', abs(ios) )  
  end subroutine write_pseudo_local

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_nl (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the non local part of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: nb, mb, n, ir, nd, i, lp, ios  

    write (ounps, '(//a13)', err = 100, iostat = ios) "<PP_NONLOCAL>"  
    do nb = 1, nbeta  
       write (ounps, '(t3,a9)', err = 100, iostat = ios) "<PP_BETA>"  
       write (ounps, '(2i5,t24,a)', err=100, iostat=ios) &
                                    nb, lll(nb), "Beta    L"
       write (ounps, '(i6)', err=100, iostat=ios) ikk2 (nb)  
       write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                    ( betar(ir,nb), ir=1,ikk2(nb) )
       write (ounps, '(t3,a10)', err = 100, iostat = ios) "</PP_BETA>"  
    enddo

    write (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_DIJ>"  
    nd = 0  
    do nb = 1, nbeta  
       do mb = nb, nbeta  
          if ( abs(dion(nb,mb)) .gt. 1.0d-12 )  nd = nd + 1 
       enddo
    enddo
    write (ounps, '(1p,i5,t24,a)', err=100, iostat=ios) &
                                   nd, "Number of nonzero Dij"
    do nb = 1, nbeta
       do mb = nb, nbeta  
          if ( abs(dion(nb,mb)) .gt. 1.0d-12 ) &
             write(ounps,'(1p,2i5,e19.11)', err=100, iostat=ios) &
                                   nb, mb, dion(nb,mb)
       enddo
    enddo
    write (ounps, '(t3,a9)', err=100, iostat=ios) "</PP_DIJ>"  

    if (pseudotype == 'US') then  
       write (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_QIJ>"  
       write (ounps, '(i5,a)',err=100, iostat=ios) nqf,"     nqf.&
          & If not zero, Qij's inside rinner are computed using qfcoef's"
       if (nqf.gt.0) then
          write (ounps, '(t5,a11)', err=100, iostat=ios) "<PP_RINNER>"  
          write (ounps,'(i5,1pe19.11)', err=100, iostat=ios) &
                                        (i, rinner(i), i = 1, nqlc)
          write (ounps, '(t5,a12)', err=100, iostat=ios) "</PP_RINNER>"  
       end if
       do nb = 1, nbeta 
          do mb = nb, nbeta
             write (ounps, '(3i5,t24,a)', err=100, iostat=ios) &
                                          nb, mb, lll(mb) , "i  j  (l(j))"
             write (ounps, '(1pe19.11,t24,a)', err=100, iostat=ios) &
                                          qqq(nb,mb), "Q_int"
             write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                          ( qfunc (n,nb,mb), n=1,mesh )
             if (nqf.gt.0) then
                write (ounps, '(t5,a11)', err=100, iostat=ios) &
                                          "<PP_QFCOEF>"  
                write(ounps,'(1p4e19.11)', err=100, iostat=ios) &
                                 ((qfcoef(i,lp,nb,mb),i=1,nqf),lp=1,nqlc)
                write (ounps, '(t5,a12)', err=100, iostat=ios) &
                                          "</PP_QFCOEF>"
             end if
          enddo
       enddo
       write (ounps, '(t3,a9)', err = 100, iostat = ios) "</PP_QIJ>"  

    endif
    write (ounps, '(a14)', err = 100, iostat = ios) "</PP_NONLOCAL>"  
    return

100 call errore ('write_pseudo_nl', 'Writing pseudo file', abs (ios) )  

  end subroutine write_pseudo_nl

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_pswfc (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the pseudo atomic functions
    !     of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: nb, ir, ios  

    write (ounps, '(//a10)', err = 100, iostat = ios) "<PP_PSWFC>"  
    do nb = 1, ntwfc
       write (ounps,'(a2,i5,f6.2,t24,a)', err=100, iostat=ios) &
            elsw(nb), lchiw(nb), ocw(nb), "Wavefunction"
       write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
            ( chi(ir,nb), ir=1,mesh )
    enddo
    write (ounps, '(a11)', err = 100, iostat = ios) "</PP_PSWFC>"  
    return

100 call errore ('write_pseudo_pswfc', 'Writing pseudo file', abs(ios) )  
  end subroutine write_pseudo_pswfc
  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_rhoatom (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the atomic charge density to the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  

    write (ounps, '(//a12)', err = 100, iostat = ios) "<PP_RHOATOM>"  
    write (ounps, '(1p4e19.11)', err = 100, iostat = ios) &
                               ( rho_at(ir), ir=1,mesh )
    write (ounps, '(a13)', err = 100, iostat = ios) "</PP_RHOATOM>"  
    return

100 call errore('write_pseudo_rhoatom','Writing pseudo file',abs(ios))
  end subroutine write_pseudo_rhoatom

  !---------------------------------------------------------------------
  subroutine dftname(iexch, icorr, igcx, igcc, longname, shortname)
  !---------------------------------------------------------------------
  implicit none
  integer iexch, icorr, igcx, igcc
  character (len=4) :: shortname
  character (len=20):: longname
  !
  ! The data used to convert iexch, icorr, igcx, igcc
  ! into a user-readable string
  !
  integer, parameter :: nxc = 1, ncc = 9, ngcx = 4, ngcc = 5
  character (len=20) :: exc, corr, gradx, gradc  
  dimension exc (0:nxc), corr (0:ncc), gradx (0:ngcx), gradc (0:ngcc)
  data exc / 'NOX ', 'SLA ' /  
  data corr / 'NOC ', 'PZ  ', 'VWN ', 'LYP ', 'PW  ', 'WIG ', 'HL  ',&
              'OBZ ', 'OBW ', 'GL  ' /
  data gradx / 'NOGX', 'B88 ', 'GGX ', 'PBE ', 'TPSS' /  
  data gradc / 'NOGC', 'P86 ', 'GGC ', 'BLYP', 'PBE ', 'TPSS' /  

  if (iexch==1.and.igcx==0.and.igcc==0) then
     shortname = corr(icorr)
  else if (iexch==1.and.icorr==3.and.igcx==1.and.igcc==3) then
     shortname = 'BLYP'
  else if (iexch==1.and.icorr==1.and.igcx==1.and.igcc==0) then
     shortname = 'B88'
  else if (iexch==1.and.icorr==1.and.igcx==1.and.igcc==1) then
     shortname = 'BP'
  else if (iexch==1.and.icorr==4.and.igcx==2.and.igcc==2) then
     shortname = 'PW91'
  else if (iexch==1.and.icorr==4.and.igcx==3.and.igcc==4) then
     shortname = 'PBE'
  else if (iexch==1.and.icorr==4.and.igcx==4.and.igcc==5) then
     shortname = 'TPSS'
  else
     shortname = ' '
  end if
  write(longname,'(4a5)') exc(iexch),corr(icorr),gradx(igcx),gradc(igcc)

  return
end subroutine dftname

module usp
  !
  ! All variables read from user-supplied file format
  ! Must have distinct names from variables in the "upf" module
  !
  
  INTEGER, PARAMETER :: DP  = KIND(1.0d0)          ! Double precision real type
  INTEGER, PARAMETER :: LONG=selected_int_kind(18) ! Long integer type

  REAL(kind=DP), PARAMETER :: PI=3.141592653589793238462643383279502884197_DP
  REAL(kind=dp), PARAMETER :: TWO_PI=2.0_DP * pi

  ! THE VALUE OF ONE ANGSTROM IN TERMS OF BOHRS
  REAL(kind=DP), PARAMETER :: ANGSTROM=1.889726313_DP 

  ! cks: The value of one Hartree in terms of electron-volts
  REAL(kind=DP), PARAMETER :: HARTREE_IN_EVS=27.2116529_DP 

  ! pdh: square root of minus one
  COMPLEX(kind=DP), PARAMETER :: cmplx_i = (0.0_DP,1.0_DP)

  ! pdh: units for standard output and standard error
  INTEGER, PARAMETER :: stdout = 6
  INTEGER, PARAMETER :: stderr = stdout
  
  ! aam: The symbols of the elements in the periodic table
  CHARACTER(len=2), PARAMETER, dimension(109) :: periodic_table_name= (/ &
       & 'H ',                                                                                'He', &
       & 'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne', &
       & 'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar', &
       & 'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
       & 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
       & 'Cs','Ba', &
       & 'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu', &
       & 'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
       & 'Fr','Ra', &
       & 'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr', &
       & 'Rf','Db','Sg','Bh','Hs','Mt' /)
  
  
  ! cks, 22/1/2004: type for containing pseudopotential information
  type PSEUDO_SPECIES

     ! cks: name of pseudopotential file
     character(len=64) :: pseudo_name

     ! cks: atomic number of atom
     integer :: atomic_number

     ! cks: charge of "pseudised" core of atom
     real(kind=DP) :: ion_charge   ! pa: changed from integer to real

     ! ndmh: whether to subtract the coulomb potential -Z/r
     logical :: subtract_coul
     
     ! cks: number of points in the radial grid
     integer :: n_rad_pts
     
     ! ndmh: spacing of points in the radial grid
     real(kind=DP) :: inv_g_spacing

     ! cks: number of angular momentum shells of projectors
     integer :: n_shells

     ! cks: all projectors of pseudopotential
     integer :: n_proj

     ! cks: angular momentum of each shell, this 
     ! cks: array is allocated to be ang_mom(n_shells)
     integer, pointer, dimension(:) :: ang_mom

     ! cks: maximum radial reciprocal k-vector up to which
     ! cks: the potential and projectors are defined
     real(kind=DP) :: ps_gmax

     ! cks: core radius for each projector shell, 
     ! cks: core_radius(n_shells)
     real(kind=DP), pointer, dimension(:) :: core_radius

     ! cks: Kleinman-Bylander denominator for each projector 
     ! cks: shell, kb_denominator(n_shells)
     real(kind=DP), pointer, dimension(:) :: kb_denominator

     ! cks: the radial recip part of the local potential
     ! cks: rad_locpot_recip(n_rad_pts)
     real(kind=DP), pointer, dimension(:) :: rad_locpot_recip

     ! cks: the values of each projector (shell) on the radial 
     ! cks: grid. The last column of this array contains the 
     ! cks: grid points, so the size of the array is 
     ! cks: rad_proj_recip(n_rad_pts, n_shells +1 )
     real(kind=DP), pointer, dimension(:,:) :: rad_proj_recip

     ! cks: each projector in an fftbox in reciprocal space 
     complex(kind=DP), pointer, dimension(:,:,:,:) :: fftbox_proj_recip

     ! ndmh: whether a core charge is present for this species
     logical :: core_charge
     
     ! ndmh: if core charge is present, array to store core charge in 
     ! ndmh: reciprocal space
     real(kind=DP), pointer, dimension(:) :: core_charge_recip
     
  end type PSEUDO_SPECIES

  type(PSEUDO_SPECIES) :: pp_usp ! ONETEP-format data
  real(kind=DP), allocatable :: usp_rlog(:) 
  real(kind=DP), allocatable :: usp_rab(:)
  real(kind=DP), allocatable :: usp_qfunc(:,:,:)
  real(kind=DP), allocatable :: usp_qfcoef(:,:,:,:)
  real(kind=DP), allocatable :: usp_core_charge_real(:)
  real(kind=DP), allocatable :: usp_qpts(:)
  real(kind=DP), allocatable :: dqtemp(:,:,:)
  real(kind=DP), allocatable :: usp_q(:,:)
  real(kind=DP), allocatable :: usp_D0(:,:)
  real(kind=DP), allocatable :: usp_rinner(:)
  integer, allocatable :: usp_nnn(:), usp_nbl(:)
  integer             :: tot_num_points  
  integer             :: tot_num_projectors  
  integer             :: usp_num_projectors
  character(len=80)   :: current_job
  character(len=80)   :: current_file
  character(len=3)    :: shell_string
  character(len=80)   :: line
  character(len=80)   :: usp_date
  character(len=80)   :: usp_author
  character(len=80)   :: usp_theory
  character(len=80)   :: usp_element
  logical             :: count_points
  real(kind=DP)       :: gmax
  real(kind=DP)       :: temp
  real(kind=DP)       :: fact1
  real(kind=DP)       :: ionic_charge
  real(kind=DP)       :: delta
  real(kind=DP)       :: usp_rcloc
  real(kind=DP)       :: usp_ecutfine
  integer             :: usp_lloc
  integer             :: usp_mesh
  integer             :: usp_kkbeta
  integer             :: usp_nqf
  integer             :: usp_nlc

contains

  subroutine internal_get_core_radii(p_species)

    !==================================================================!
    ! This subroutine calculates the core radii of the pseudopotential !
    ! species by examining the projectors in real space.               !
    !------------------------------------------------------------------!
    ! Written by Peter D. Haynes in early 2004.                        !
    ! Tidied up by Peter D. Haynes on 1/7/2004.                        !
    !==================================================================!

    implicit none

    ! Argument
    type(PSEUDO_SPECIES), intent(inout) :: p_species

    ! Local variables
    integer :: ierr                           ! Error flag
    integer :: i_shell                        ! Shell counter
    integer :: l                              ! Angular momentum of shell
    integer :: g_rad_pt,r_rad_pt              ! Grid point counters
    integer :: n_r_rad_pts                    ! Num points on real radial grid
    real(kind=DP), parameter :: rmax = 5.0_DP ! Maximum core radius
    real(kind=DP), parameter :: halfpi = 0.5_DP * PI
    real(kind=DP) :: gnorm,rnorm              ! Norms in recip/real space
    real(kind=DP) :: proj,proj_max            ! Projector (max value) in rspc
    real(kind=DP) :: val,fac                  ! Variables for integration
    real(kind=DP) :: delr,r,delg,g            ! Grid variables          
    real(kind=DP), allocatable :: integrand(:)! gspc integrand for FT

    ! Set parameters for radial grid
    delg = p_species%ps_gmax / (p_species%n_rad_pts - 1)
    delr = pi / p_species%ps_gmax
    n_r_rad_pts = nint(rmax / delr) + 1
    fac = (7.0_DP/17280.0_DP) * delg / sqrt(halfpi)

    ! Allocate workspace for integration
    allocate(integrand(p_species%n_rad_pts),stat=ierr)

    ! Loop over shells: angular momentum is l
    do i_shell=1,p_species%n_shells
       l = p_species%ang_mom(i_shell)

       ! Calculate norm from reciprocal space representation
       gnorm = 0.0_DP
       do g_rad_pt=1,p_species%n_rad_pts
          g = (g_rad_pt-1) * delg
          val = p_species%rad_proj_recip(g_rad_pt,i_shell) * g
          gnorm = gnorm + val*val
       end do
       gnorm = gnorm * delg

       if (gnorm < epsilon(1.0_DP)) then
          p_species%core_radius(i_shell) = 0.0_DP
          cycle
       end if

       ! Calculate projector on real-space radial grid, and cumulative
       ! norm of projector in real space
       rnorm = 0.0_DP
       proj_max = 0.0_DP

       do r_rad_pt=1,n_r_rad_pts
          r = (r_rad_pt-1) * delr
          p_species%core_radius(i_shell) = 1.8_DP
          integrand(1) = 0.0_DP
          do g_rad_pt=2,p_species%n_rad_pts
             g = (g_rad_pt-1) * delg
             integrand(g_rad_pt) = services_sbessj(l,g*r) * g*g * &
                  p_species%rad_proj_recip(g_rad_pt,i_shell)
          end do

          ! Due to the oscillatory nature of the integrand,
          ! an eight point Newton-Cotes integration method is used.
          ! See  Abramowitz and Stegun Eq. 25.4.17

          proj = 0.0_DP
          do g_rad_pt=1,p_species%n_rad_pts-7,7
             proj = proj + &
                  751.0_DP*(integrand(g_rad_pt)+integrand(g_rad_pt+7)) + &
                  3577.0_DP*(integrand(g_rad_pt+1)+integrand(g_rad_pt+6)) + &
                  1323.0_DP*(integrand(g_rad_pt+2)+integrand(g_rad_pt+5)) + &
                  2989.0_DP*(integrand(g_rad_pt+3)+integrand(g_rad_pt+4))
          end do
          proj = proj * fac
          val = r * proj
          rnorm = rnorm + val*val * delr

          ! Set core radius if magnitude of projector is less then 0.1% of
          ! maximum value and cumulative norm is at least 99.99% of
          ! reciprocal space norm
          ! ndmh: tolerances decreased to avoid chopping end of projector off
          proj_max = max(proj_max,abs(proj))
          if (r_rad_pt > 1) then
             if (abs(proj/proj_max) < 0.5e-3_DP .and. &
                  rnorm/gnorm > 0.99995_DP) then
                p_species%core_radius(i_shell) = r
                exit
             end if
          end if

       end do

    end do

    ! Free workspace
    deallocate(integrand,stat=ierr)

  end subroutine internal_get_core_radii

  real(kind=DP) function services_sbessj(l,x)

    !===============================================!
    ! SPHERICAL BESSEL function OF THE FIRST KIND.  !
    !-----------------------------------------------!
    ! Written by Peter D. Haynes in early 2004.     !
    !===============================================!

    ! Arguments
    integer, intent(in) :: l
    real(kind=DP), intent(in) :: x

    ! Local variables
    integer :: j
    real(kind=DP), parameter :: third = 1.0_DP / 3.0_DP
    real(kind=DP), parameter :: ftnth = 1.0_DP / 14.0_DP
    real(kind=DP) :: x2, sb0, sb1, by, bym, byp, ux
    real(kind=DP) :: sbessj

    if (abs(x) > 0.001_DP) then
       sb0 = sin(x)/x
    else
       x2 = 0.5_DP*x*x
       sb0 = 1.0_DP - third*x2*(1.0_DP - 0.1_DP*x2)
    end if
    if (l == 0) then
       sbessj = sb0
    else
       if (abs(x) > 0.001_DP) then
          sb1 = (sb0 - cos(x)) / x
       else
          sb1 = third*x*(1.0_DP - (0.2_DP*x2)*(1.0_DP - ftnth*x2))
       end if
       if (l == 1) then
          sbessj = sb1
       else if (x == 0.0_DP) then
          sbessj = 0.0_DP
       else
          by = sb1
          bym = sb0
          ux = 1.0_DP / x
          do j=1,l-1
             byp = real(2*J+1,DP)*ux*by - bym
             bym = by
             by = byp
          end do
          sbessj = by
       end if
    end if

    services_sbessj =sbessj

  end function services_sbessj
  
  subroutine internal_radial_transform_b(nqpts,qpts,nrpts,rpts,qfunc,qion,l,rfunc)

    !=========================================================================!
    ! Convert reciprocal space function to the log real grid                  !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 19/02/2009                                  !
    !=========================================================================!

    implicit none

    ! Arguments
    integer,intent(in) :: nqpts
    integer,intent(in) :: nrpts
    real(kind=dp),intent(in) :: qpts(nqpts)
    real(kind=dp),intent(in) :: rpts(nrpts)
    real(kind=dp),intent(in) :: qfunc(nqpts)
    real(kind=dp),intent(in) :: qion
    integer,intent(in) :: l
    real(kind=dp),intent(out) :: rfunc(nrpts)
    
    ! Locals
    real(kind=dp) :: q,dq
    real(kind=dp) :: work(nqpts)
    integer :: ir,iq
    
    dq = qpts(2)-qpts(1)
    rfunc(:) = 0.0_DP

    do ir=1,nrpts

       do iq=1,nqpts
          q = qpts(iq)
          work(iq) = (qfunc(iq)*q**2-qion)*services_sbessj(l,q*rpts(ir))
       end do

       rfunc(ir) = work(1)+work(nqpts)
       do iq = 2,nqpts-1,2
          rfunc(ir) = rfunc(ir) + 4.0_dp*work(iq)+2.0_dp*work(iq+1)
       enddo
       rfunc(ir) = rfunc(ir)*dq/3.0_dp

       !if(l==1) print *,rpts(ir),rfunc(ir)*rpts(ir)

    end do

    ! Normalise
    rfunc = 2.0_DP*rfunc/PI

  end subroutine internal_radial_transform_b
  
  subroutine internal_radial_transform(nrpts,rlog,rab,power,nqpts,qmax, &
       rfunc,qfunc)

    !=========================================================================!
    ! Transforms a radial real space function on a logarithmic grid to a      !
    ! regular reciprocal space grid                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1) nrpts  : input  : number of real space points                        !
    ! 2) rlog   : input  : the logarithmic grid                               !
    ! 3) rab    : input  : the "spacings" of the logarithmic grid             !
    ! 4) power  : input  : power of r needed in integral                      !
    ! 4) nqpts  : input  : number of points in the regular reciprocal grid    !
    ! 5) qmax   : input  : the maximum q-vector for the transform             !
    ! 6) rfunc  : input  : the real space radial function                     !
    ! 7) qfunc  : output : the reciprocal space function                      !
    !=========================================================================!

    implicit none

    integer,intent(in) :: nrpts
    real(kind=dp),intent(in) :: rlog(nrpts)
    real(kind=dp),intent(in) :: rab(nrpts)
    integer,intent(in) :: power
    integer,intent(in) :: nqpts
    real(kind=dp),intent(in) :: qmax   
    real(kind=dp),intent(in) :: rfunc(nrpts)
    real(kind=dp),intent(out) :: qfunc(nqpts)

    ! Local Variables
    integer       :: nq,nr,ll
    integer       :: ierr
    real(kind=dp) :: q,qr
    real(kind=dp) :: rwgt(nrpts)
    real(kind=dp) :: lookup(nqpts,nrpts)

    if(nrpts==2*(nrpts/2)) then
       write(stdout,'(a,a)') 'Error in internal_radial_transform : ', &
            'nrpts must be odd'
       stop
    end if

    ! Set the integration weights
    rwgt(1) = 4.0_dp/12.0_dp
    do nr=2,nrpts-1,2
       rwgt(nr)=rab(nr)*16.0_dp/12.0_dp*rlog(nr)**power
    end do
    do nr=3,nrpts-2,2
       rwgt(nr)=rab(nr)*8.0_dp/12.0_dp*rlog(nr)**power
    end do
    rwgt(nrpts) = 4.0_dp/12.0_dp*rlog(nrpts)**power

    ! Loop over the radial reciprocal space grid
    do nq=1,nqpts
       q = real(nq-1,dp)*qmax/real(nqpts-1,dp)
       lookup(nq,1) = services_sbessj(0,0.0_dp)*rwgt(1)
       do nr=2,nrpts
          qr = q*rlog(nr)
          lookup(nq,nr) = services_sbessj(0,qr)*rwgt(nr)
       end do
    end do

    ! Do the matrix-vector multiplication
    call dgemv('N',nqpts,nrpts,1.0_dp,lookup(1,1),nqpts,rfunc,1,0.0_dp,qfunc,1)

  end subroutine internal_radial_transform
  
end module usp

    !   ----------------------------------------------------------
    subroutine convert_usp
    !   ----------------------------------------------------------

      use usp
      use upf
      implicit none
      
      ! Locals
      integer,parameter :: ndmx = 1999 ! PWSCF max grid size
      real(dp) :: meshmax,minmesh,a_grid,b_grid,r2r3
      real(dp) :: ion_rad,C_ion,d_ion
      real(dp) :: usp_qfunc_recip(pp_usp%n_rad_pts)
      integer :: ib,jb,n

      ! PP_INFO
      rel = 0
      rcloc = usp_rcloc
      generated = 'CASTEP .usp generator'
      date_author = usp_date//usp_author
      comment = 'COMMENT'
      
      ! PP_HEADER      
      psd = usp_element
      if (any(usp_q/=0.0_DP)) then
         pseudotype = 'US'
      else
         pseudotype = 'NC'
      end if
      nlcc = pp_usp%core_charge
      select case (usp_theory)
      case ('LDA')
         icorr=1; iexch=1; igcx=0; igcc=0;
      case ('PBE')
         icorr=4; iexch=1; igcx=3; igcc=4;
      case default
         icorr=1; iexch=1; igcx=0; igcc=0;
      end select
      zp = pp_usp%ion_charge
      etotps = 0.0_DP
      
      ! Guesses for ecuts based on 'FINE' energy in USP
      ecutrho = 4.0_DP*usp_ecutfine/HARTREE_IN_EVS
      ecutwfc = usp_ecutfine/HARTREE_IN_EVS
      
      ! Choose grid size and logarithmic grid coeffs
      mesh = min(ndmx,usp_mesh)
      if (usp_rlog(1)==0.0_DP) then
         ! CASTEP-style grid r(i) = a*(exp[(i-1)b]-1)
         r2r3 = usp_rlog(2)/usp_rlog(3)   ! a(e^b - 1) / a(e^2b - 1)
         ! Solve quadratic in e^b
         b_grid = sqrt(1.0_DP-4.0_DP*r2r3*(1.0_DP-r2r3))
         b_grid = ((1.0_DP+b_grid)/(2.0_DP*r2r3))
         b_grid = log(b_grid)
         ! Find prefactor
         a_grid = (usp_rlog(2)-usp_rlog(3))/(exp(b_grid)-exp(2.0_DP*b_grid))
         
         meshmax=a_grid*(exp(real(mesh-1,DP)*b_grid)-1.0_DP)
         minmesh = 100.0_DP
         if (meshmax<minmesh) then
            ! See if we can increase the number of grid points
            mesh = 1 + int(log(minmesh/a_grid+1.0_DP)/b_grid)
            meshmax=a_grid*(exp(real(mesh-1,DP)*b_grid)-1.0_DP)
            ! If this has failed, scale up a_grid to make meshmax larger
            if (meshmax<minmesh) a_grid = a_grid*(minmesh/meshmax)
         end if
      else
         print '(a)','Error: Unknown grid format, quitting'
         stop
      end if
      
      ! Copy array sizes and allocate arrays
      lmax = maxval(pp_usp%ang_mom(:))
      nbeta = pp_usp%n_shells
      allocate(ikk2(nbeta),lll(nbeta),dion(nbeta,nbeta),qqq(nbeta,nbeta))
      ntwfc = 0
      allocate(elsw(ntwfc), lchiw(ntwfc), ocw(ntwfc))
      nqlc = usp_nlc + 1
      allocate(rinner(1:nqlc))
      nqf = usp_nqf
      allocate(qfcoef(nqf,nqlc,nbeta,nbeta))
      ikk2(:) = mesh
      allocate(r(mesh),rab(mesh),rho_atc(mesh),rho_at(mesh))
      allocate(chi(mesh,nbeta))
      allocate(qfunc(mesh,nbeta,nbeta))
      allocate(betar(maxval(ikk2),nbeta))
      allocate(vloc0(mesh))
      nwfs = 0
      allocate(els(nwfs),nns(nwfs),lchi(nwfs),oc(nwfs),rcut(nwfs), &
           rcutus(nwfs),epseu(nwfs))

      do ib=1,mesh
         r(ib) = a_grid*(exp(b_grid*real(ib-1,DP))-1.0_DP)
         rab(ib) = b_grid*(r(ib)+a_grid)
         !r(ib) = a_grid*(exp(b_grid*real(ib,DP)))
         !rab(ib) = b_grid*(r(ib))
      end do

      ! NLCC Core Charge
      if (nlcc) then
         call internal_radial_transform_b(pp_usp%n_rad_pts,usp_qpts,mesh, &
              r(:),pp_usp%core_charge_recip(:),0.0_DP,0,rho_atc)
              
         ! Remove CASTEP scaling
         rho_atc(:) = rho_atc(:) * r(:)**0 / (4.0_DP * PI)
      end if

      ! Neutral atom density
      ion_rad = 1.1_DP ! Default radius 1.1 Bohr
      d_ion = 1.0_DP/ion_rad**2
      C_ion = zp*(d_ion/PI)**(1.5_DP)
      do ib=1,mesh
         rho_at(ib) = C_ion*exp(-d_ion*r(ib)**2)*4.0_DP*PI*r(ib)**2
         if (r(ib)>8.0_DP) rho_at(ib) = 0.0_DP
      end do
      
      ! Set projector angular momenta
      lll(1:nbeta) = pp_usp%ang_mom(1:nbeta)

      ! Fourier transform local potential to real space
      call internal_radial_transform_b(pp_usp%n_rad_pts,usp_qpts,mesh,&
           r(:),pp_usp%rad_locpot_recip(:),zp,0,vloc0)
      
      ! Convert to Rydbergs
      vloc0 = vloc0 * 2.0_DP
           
      ! Convert beta projectors
      do ib=1,nbeta
      
         ! Find limits of beta functions on logarithmic grid
         do n=1,mesh
            ! Check if r(n) is more than 0.5au past cutoff and n is odd
            if ((r(n)>pp_usp%core_radius(ib)+0.5).and.(modulo(n,2)==1)) then
               ikk2(ib) = n
               exit
            end if
         end do
         
         ! Fourier transform projectors to real space
         call internal_radial_transform_b( &
              pp_usp%n_rad_pts,usp_qpts, &
              ikk2(ib),r(1:ikk2(ib)),&
              pp_usp%rad_proj_recip(1:pp_usp%n_rad_pts,ib), &
              0.0_DP,lll(ib),betar(1:ikk2(ib),ib))
         
         ! CASTEP takes transform of beta/r
         betar(1:ikk2(ib),ib) = betar(1:ikk2(ib),ib) * r(1:ikk2(ib))
         
         ! Zero projectors beyond cutoff
         do n=1,ikk2(ib)
            if (r(n)>pp_usp%core_radius(ib)) then
               betar(n,ib) = 0.0_DP
            end if
         end do
         
         ! Convert to Rydbergs
         betar(:,ib) = betar(:,ib) * 2.0_DP
      end do
      
      ! Set D_0 factors and convert to Rydbergs^-1
      dion(1:nbeta,1:nbeta) = usp_D0(1:nbeta,1:nbeta) * 0.5_DP
      
      ! Set qqq and convert to Rydbergs^-2
      qqq(1:nbeta,1:nbeta) = usp_q(1:nbeta,1:nbeta) * 0.25_DP
      
      ! Set L independent function
      do ib=1,nbeta
         do jb=1,nbeta
            ! Convert to reciprocal space
            call internal_radial_transform(usp_kkbeta,usp_rlog,usp_rab, 2, &
                 pp_usp%n_rad_pts,pp_usp%ps_gmax,usp_qfunc(1:usp_kkbeta,ib,jb),&
                 usp_qfunc_recip)
            ! Interpolate back to real space onto UPF logarithmic grid
            call internal_radial_transform_b(pp_usp%n_rad_pts,usp_qpts,mesh,&
                 r(1:mesh),usp_qfunc_recip,0.0_DP,0,qfunc(1:mesh,ib,jb))
            ! Convert to Rydbergs^-2
            qfunc(1:mesh,ib,jb) = qfunc(1:mesh,ib,jb) * 0.25_DP
            ! Zero functions beyond cutoff
            do n=1,mesh
               if (r(n)>max(pp_usp%core_radius(ib),pp_usp%core_radius(jb))) then
                  qfunc(n,ib,jb) = 0.0_DP
               end if
            end do
         end do
      end do
      
      ! Set L dependent coeffs and convert to Rydbergs^-2
      qfcoef(1:nqf,1:nqlc,1:nbeta,1:nbeta) = &
           usp_qfcoef(1:nqf,0:usp_nlc,1:nbeta,1:nbeta) * 0.25_DP
      
      ! Set rinner
      rinner(1:nqlc) = usp_rinner(0:usp_nlc)
      
      !     ----------------------------------------------------------
      write (6,'(a)') 'Pseudopotential successfully converted'
      !     ----------------------------------------------------------
      return

    end subroutine convert_usp

    subroutine scan_usp(iunps)

      !============================================================! 
      ! This subroutine reads the pseudopotential for a .usp       !
      ! and stores the sizes required to allocate the necessary    !
      ! memory for this element in the p_species array.            !
      !------------------------------------------------------------!
      ! Written by Nicholas Hine on 29/01/09.                      !
      !============================================================!
      
      use usp
      
      implicit none
      
      ! Argument
      integer,intent(in) :: iunps
      
      integer :: ierr, i
      
      ! Local Variables
      integer, parameter :: max_proj=huge(1)
      integer :: version(2), usp_nlcc
      real(dp) :: eeig, zele
      
      character(80) :: line_cut

      current_job = 'seeking marker "END COMMENT"'
      line = repeat(' ',80)
      usp_author = ''
      usp_date = ''
      usp_element = ''
      usp_theory = ''
      usp_ecutfine = -1.0_DP
      do while(index(line,'END COMMENT') == 0)

         if (index(line,'Date of generation')>0) then
            i = index(line,'Date of generation') + 18
            usp_date = line(i:i+10)
         end if
         if (index(line,'Author:')>0) then
            i = index(line,'Author:') + 7
            usp_author = line(i:80)
            i = index(usp_author,'|') - 1
            usp_author = trim(usp_author(1:i))
         end if
         if (index(line,'Element:')>0) then
            i = index(line,'Element:') + 9
            usp_element = line(i:i+2)
         end if
         if (index(line,'z =')>0) then
            i = index(line,'z =') + 3
            line_cut = line(i:80)
            read(line_cut,*) zele
            usp_element = periodic_table_name(int(zele))
         end if
         if (index(line,'Level of theory:')>0) then
            i = index(line,'Level of theory:') + 17
            usp_theory = line(i:80)
            i = index(usp_theory,'|') - 1
            usp_theory = trim(usp_theory(1:i))
         end if
         if (index(line,' loc ')>0) then
            i = index(line,'loc') + 4
            line_cut = line(i:80)
            read(line_cut,*) usp_lloc, eeig, usp_rcloc
         end if
         if (index(line,'FINE')>0.and.(usp_ecutfine<0.0_DP)) then
            read(line,*) usp_ecutfine
         end if
         
         read(iunps,'(a80)',err=101,end=200) line
      end do
      
      ! ** Read the version
      read(iunps,*,err=101,end=200) version(1),version(2)
      
      ! ** Read the Ionic Charge
      read(iunps,*,err=101,end=200) ionic_charge
      
      pp_usp%ion_charge = ionic_charge
      
      ! ** Read the gmax, and nlcc flag
      read(iunps,*,err=101,end=200) gmax, usp_nlcc
      pp_usp%ps_gmax = gmax/ANGSTROM
      
      ! ** Determine whether to include nonlinear core corrections
      pp_usp%core_charge = (usp_nlcc==2)
      
      ! ** Find out the number of projectors, and g-grid points
      count_points = .true.
      tot_num_points = 0
      usp_num_projectors = 0
      line = ''
      do while((len_trim(adjustl(line))/=4).and.(usp_num_projectors<max_proj))
         read(iunps,'(a)',err=101,end=200) line
         if(len_trim(adjustl(line))==1) then 
            usp_num_projectors = usp_num_projectors+1 
            count_points = .false.
         end if
         if(count_points) then
            do i=1,len(line)
               if(line(i:i)=='.') tot_num_points = tot_num_points + 1
            end do
         end if
      end do
      
      pp_usp%n_rad_pts = tot_num_points
      pp_usp%n_shells  = usp_num_projectors
      
      ! ndmh: calculate inverse of spacing of recip space points
      pp_usp%inv_g_spacing = (pp_usp%n_rad_pts-1) / &
           pp_usp%ps_gmax
      
      ! ndmh: usps do not include the -Z/r, so flag to not subtract it
      pp_usp%subtract_coul = .false.

      ! ** Close the file
      close(2,iostat=ierr,err=500)
500   if (ierr /= 0) then
         write(stdout,'(3a,i6)') &
              'Error in scan_usp: closing file "', &
              trim(current_file),'" failed with code ',ierr
         stop
      end if

      return

101   write(stdout,'(4a)') 'Error in scan_usp: reading file "', &
           trim(current_file),'" failed while ',trim(current_job)
      stop

200   write(stdout,'(4a)') 'Error in scan_usp: file "', &
           trim(current_file),'" ended unexpectedly while ',trim(current_job)
      stop

    end subroutine scan_usp
    
    ! 
    !----------------------------------------------------------
    subroutine allocate_usp
    !----------------------------------------------------------
    ! 
      use usp

      implicit none
      
      ! Locals
      integer :: ierr
      
      allocate(pp_usp%ang_mom(pp_usp%n_shells),stat=ierr)

      allocate(pp_usp%core_radius(pp_usp%n_shells),stat=ierr)

      allocate(pp_usp%kb_denominator(pp_usp%n_shells),stat=ierr)

      allocate(pp_usp%rad_locpot_recip(pp_usp%n_rad_pts),stat=ierr)

      allocate(pp_usp%rad_proj_recip(pp_usp%n_rad_pts,pp_usp%n_shells+1), &
           stat=ierr)
      
      ! ndmh: allocate memory for core charge if present
      if (pp_usp%core_charge) then
         allocate(pp_usp%core_charge_recip(pp_usp%n_rad_pts),stat=ierr)
      end if
      
      allocate(usp_qpts(pp_usp%n_rad_pts),stat=ierr)
      
      allocate(usp_nbl(0:3),stat=ierr)
      
      allocate(usp_nnn(pp_usp%n_shells),stat=ierr)
      
      allocate(dqtemp(pp_usp%n_shells,pp_usp%n_shells,0:3),stat=ierr)
      
      allocate(usp_q(pp_usp%n_shells,pp_usp%n_shells),stat=ierr)
      
      allocate(usp_D0(pp_usp%n_shells,pp_usp%n_shells),stat=ierr)
       
    end subroutine allocate_usp
    
    ! 
    !----------------------------------------------------------
    subroutine read_usp(iunps)
    !----------------------------------------------------------
    ! 
      use usp

      implicit none

      ! Arguments
      integer,intent(in) :: iunps

      ! Locals
      integer :: ierr, i, j, k, l, lmin, shell
      real(kind=DP) :: rr

      !==============================================================! 
      ! This subroutine reads the pseudopotential for a .usp         !
      ! and stores the local and nonlocal components of the          !
      ! potential in the arrays allocated in pseudopot_alloc_species !
      !--------------------------------------------------------------!
      ! Written by Nicholas Hine on 29/01/09                         !
      !==============================================================!
      
      ! Local Variables
      integer,parameter :: max_proj=huge(1)
      real(kind=DP) :: factor
      
      ! ndmh: get to the top of the pseudopotential proper
      current_job = 'seeking marker "END COMMENT"'
      line = repeat(' ',80)
      do while (index(line,'END COMMENT') == 0)
         read(iunps,'(a80)',err=102,end=202) line  ! To the end of the comments
      end do
      current_job = 'moving to end of header'
      do i=1,3
         read(iunps,'(a80)',err=102,end=202) line  ! To the end of the header
      end do
      
      ! ndmh: read the local component of the pseudopotential
      read(iunps,*,err=102,end=202) (pp_usp%rad_locpot_recip(i), &
           i=1,pp_usp%n_rad_pts)

      ! ndmh: Scale radial locpot with suitable factors to convert it
      ! ndmh: to atomic units
      factor = ((ANGSTROM**3)/HARTREE_IN_EVS) / (4.0_DP*PI) 
      pp_usp%rad_locpot_recip(:) = &
           pp_usp%rad_locpot_recip(:) * factor
      
      ! ndmh: Read the non-local component of the pseudopotential
      do shell=1,pp_usp%n_shells
      
         write(shell_string,'(i3)') shell

         write(current_job,'(a80)') 'reading angular momentum of shell '// &
              adjustl(shell_string)
         current_job = adjustl(current_job)
         read(iunps,*,err=102,end=202) pp_usp%ang_mom(shell)
         usp_nbl(pp_usp%ang_mom(shell)) = usp_nbl(pp_usp%ang_mom(shell)) + 1
         usp_nnn(shell) = usp_nbl(pp_usp%ang_mom(shell))

         write(current_job,'(a80)') 'reading KB denominator(s) of shell '// &
              adjustl(shell_string)
         current_job = adjustl(current_job)
         !read(iunps,*,err=102,end=202) pp_usp%kb_denominator(shell)
         read(iunps,*,end=102,err=202) &
             (dqtemp(usp_nnn(shell),i,pp_usp%ang_mom(shell)),i=1,usp_nnn(shell))

         write(current_job,'(a80)') 'reading projector of shell '// &
              adjustl(shell_string)
         current_job = adjustl(current_job)
         read(iunps,*,err=102,end=202) (pp_usp%rad_proj_recip(i,shell), &
              i=1,pp_usp%n_rad_pts)

         ! cks: convert projector units from A^{3/2}*eV to a0^{3/2}*Eh 
         ! cks: scaling consistent with CASTEP
         pp_usp%rad_proj_recip(:,shell) = &
              pp_usp%rad_proj_recip(:,shell) * sqrt(ANGSTROM**3) 
         
      end do
      
      ! Fill in other half of dqtemp
      do l=0,maxval(pp_usp%ang_mom(1:usp_num_projectors))
         do i=1,usp_nbl(l)
            do j=1,i-1
               dqtemp(j,i,l) = dqtemp(i,j,l)
               !print '(6(a,i2),a)','dqtemp(',j,',',i,',',l,') = dqtemp(',i,',',j,',',l,')'
            end do
         end do
      end do

      ! Convert dqtemp to usp_D0
      do i=1,usp_num_projectors
         do j=1,usp_num_projectors
            l = pp_usp%ang_mom(i)
            if(pp_usp%ang_mom(j).eq.l) then
               usp_D0(i,j) = dqtemp(usp_nnn(i),usp_nnn(j),l)
            else
               usp_D0(i,j) = 0.0_dp
            end if
         end do
      end do
      
      ! Convert to Hartrees
      usp_D0(:,:) = usp_D0(:,:) / HARTREE_IN_EVS

      ! NB: usps store 1/kb_denom rather than kb_denom
      do i=1,usp_num_projectors
         pp_usp%kb_denominator(i) = 1.0_DP / usp_D0(i,i)
      end do
      
      ! set number of q components
      usp_nlc = 2*maxval(pp_usp%ang_mom(1:pp_usp%n_shells))

      ! cks: initialise the radial grid for the interpolation
      ! pdh: inline this call to avoid awkward copy of arguments
      delta = pp_usp%ps_gmax / real(pp_usp%n_rad_pts-1,DP)
      do i=1,pp_usp%n_rad_pts
         pp_usp%rad_proj_recip(i,pp_usp%n_shells+1) = &
              real(i-1,DP) * delta
      end do
      usp_qpts(:) = pp_usp%rad_proj_recip(:,pp_usp%n_shells+1)
      
      ! ndmh: check we are at the right point of the file
      !    -- There is no '1000' flag if there are 2 projectors per
      !       all channels up to l=3 (the assumed maximum)
      write(current_job,'(a80)') 'reading "1000" '
      if(pp_usp%n_shells < max_proj) then
         read(iunps,'(a)',err=102,end=202) line
         if(trim(adjustl(line))/='1000') goto 202
      end if
      
      ! ndmh: Read in augmentation charges
      write(current_job,'(a80)') 'reading augmentation charges'
      do shell=1,pp_usp%n_shells
         read(iunps,*,end=102,err=202) &
              (dqtemp(usp_nnn(shell),i,pp_usp%ang_mom(shell)), &
              i=1,usp_nnn(shell))
      end do
      
       ! Fill other half
       do l=0,maxval(pp_usp%ang_mom(1:usp_num_projectors))
          do i=1,usp_nbl(l)
             do j=1,i-1
                dqtemp(j,i,l) = dqtemp(i,j,l)
             end do
          end do
       end do

       ! Convert dqtemp to usp_q
       do i=1,usp_num_projectors
          do j=1,usp_num_projectors
             l = pp_usp%ang_mom(i)
             if(l.eq.pp_usp%ang_mom(j)) then
                usp_q(i,j) = dqtemp(usp_nnn(i),usp_nnn(j),l)
             else
                usp_q(i,j) = 0.0_dp
             end if
          end do
       end do
      
      ! ndmh: read in grid sizes for usp_kkbeta and usp_mesh
      write(current_job,'(a80)') 'reading grid sizes'
      if(pp_usp%core_charge) then
         read(iunps,*,err=102,end=202) usp_kkbeta,usp_mesh
      else          
         read(iunps,*,err=102,end=202) usp_kkbeta
         usp_mesh = usp_kkbeta
      end if
      
      ! ndmh: allocate temporary storage
      allocate(usp_rlog(usp_mesh),stat=ierr)
      allocate(usp_rab(usp_mesh),stat=ierr)
      allocate(usp_rinner(0:usp_nlc),stat=ierr)
      allocate(usp_qfunc(usp_kkbeta,usp_num_projectors,usp_num_projectors), &
           stat=ierr)

      ! * Allocate temporary storage for usp_core_charge_real if required
      if (pp_usp%core_charge) then
         allocate(usp_core_charge_real(usp_mesh),stat=ierr)
      end if
     
      ! ndmh: read in the grid points
      write(current_job,'(a80)') 'reading grid points'
      read(iunps,*,err=102,end=202) (usp_rlog(i),i=1,usp_mesh)
      read(iunps,*,err=102,end=202) (usp_rab(i),i=1,usp_mesh)
      
      ! ndmh: read the L independent function
      write(current_job,'(a80)') 'reading L independent function'      
      do i=1,pp_usp%n_shells
         do j=1,i
            read(iunps,*,err=102,end=202) (usp_qfunc(k,i,j),k=1,usp_kkbeta)
         end do
      end do

      ! ndmh: the L dependent coeffs
      ! ndmh: these are zero for ncpp's
      write(current_job,'(a80)') 'reading L dependent coeffs'
      read(iunps,*,err=102,end=202) usp_nqf
      
      allocate(usp_qfcoef(1:usp_nqf,0:usp_nlc,1:usp_num_projectors, &
           1:usp_num_projectors),stat=ierr)

      do l=0,usp_nlc
         read(iunps,*,err=102,end=202) usp_rinner(l)
         do i=1,pp_usp%n_shells
            do j=1,i
               read(iunps,*,err=102,end=202) (usp_qfcoef(k,l,i,j),k=1,usp_nqf)
            end do
         end do
      end do

      ! Assuming MSI corruption of usp_qfunc, regenerate
      do i=1,usp_num_projectors
         do j=1,i
            lmin = abs(pp_usp%ang_mom(i)-pp_usp%ang_mom(j))
            k=0
            do while(usp_rlog(k+1).lt.usp_rinner(lmin))
               k=k+1
               rr = usp_rlog(k)**2.0_dp
               usp_qfunc(k,i,j) = usp_qfcoef(1,lmin,i,j)
               do l=2,usp_nqf
                  usp_qfunc(k,i,j) = usp_qfunc(k,i,j) + &
                       usp_qfcoef(l,lmin,i,j)*rr**real(l-1,dp)
               end do
               usp_qfunc(k,i,j)=usp_qfunc(k,i,j)*usp_rlog(k)**real(lmin+2,dp)
            end do
         end do
      end do

      ! Fill in other halves
      do i=1,usp_num_projectors
         do j=1,i-1
            do k=1,usp_kkbeta
               usp_qfunc(k,j,i) = usp_qfunc(k,i,j)
            end do
            do l=0,usp_nlc
               do k=1,usp_nqf
                  usp_qfcoef(k,l,j,i) = usp_qfcoef(k,l,i,j)
               end do
            end do
         end do
      end do

      ! Finally read in the core charge on the log grid
      if (pp_usp%core_charge) then
         read(iunps,'(3f20.10)',err=102,end=202) &
              (usp_core_charge_real(i),i=1,usp_mesh)

         ! Transform the real-space core charge to reciprocal space
         call internal_radial_transform(usp_mesh,usp_rlog,usp_rab, 0, &
              pp_usp%n_rad_pts, pp_usp%ps_gmax,usp_core_charge_real, &
              pp_usp%core_charge_recip)

      end if
      
      ! pdh: get core radii from projectors
      call internal_get_core_radii(pp_usp)
      
      ! ndmh: deallocate temporary storage
      !if (pp_usp%core_charge) then
      !   deallocate(usp_core_charge_real,stat=ierr)
      !end if
      !deallocate(usp_rab,stat=ierr)
      !deallocate(usp_rlog,stat=ierr)

      ! ndmh: close the file
      close(2,iostat=ierr,err=402)
402   if (ierr /= 0) then
         write(stdout,'(3a,i6)') &
              'Error in read_usp :',&
              ' closing file "', trim(current_file),'" failed with code ',ierr
         stop
      end if
      
      return

102   write(stdout,'(4a)') 'Error in read_usp ','reading file "',&
           trim(current_file),'" failed while ',trim(current_job)
      stop
    
202   write(stdout,'(4a)') 'Error in read_usp ','file "',trim(current_file), &
           '" ended unexpectedly while ',trim(current_job)
      stop
  !
  !     ----------------------------------------------------------
  write (6,'(a)') 'Pseudopotential successfully read'
  !     ----------------------------------------------------------
  !
  return
100 call errore ('read_usp', 'Reading pseudo file', 100 )

end subroutine read_usp

!
! Copyright (C) 2002-2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE input_from_file( )
  !
  ! This subroutine checks program arguments and, if input file is present,
  ! attach input unit ( 5 ) to the specified file
  !
  !
  IMPLICIT NONE
  !
  INTEGER             :: unit = 5, &
                         ilen, iiarg, nargs, ierr
  ! do not define iargc as external: g95 does not like it
  INTEGER             :: iargc
  CHARACTER (LEN=80)  :: input_file
  !
  ! ... Input from file ?
  !
  nargs = iargc()
  !
  DO iiarg = 1, ( nargs - 1 )
     !
     CALL getarg( iiarg, input_file )
     !
     IF ( TRIM( input_file ) == '-input' .OR. &
          TRIM( input_file ) == '-inp'   .OR. &
          TRIM( input_file ) == '-in' ) THEN
        !
        CALL getarg( ( iiarg + 1 ) , input_file )
        !
        OPEN ( UNIT = unit, FILE = input_file, FORM = 'FORMATTED', &
               STATUS = 'OLD', IOSTAT = ierr )
        !
        CALL errore( 'input_from_file', 'input file ' // TRIM( input_file ) &
             & // ' not found' , ierr )
        !
     END IF
     !
  END DO

END SUBROUTINE input_from_file
!
!----------------------------------------------------------------------------
SUBROUTINE get_file( input_file )
  !
  ! This subroutine reads, either from command line or from terminal,
  ! the name of a file to be opened
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=*)  :: input_file
  !
  CHARACTER (LEN=256)  :: prgname
  ! do not define iargc as external: g95 does not like it
  INTEGER             :: nargs, iargc
  LOGICAL             :: exst
  !
  nargs = iargc()
  CALL getarg (0,prgname)
  !
  IF ( nargs == 0 ) THEN
10   PRINT *,'("Input file > ",$)'
     READ (5,'(a)', end = 20, err=20) input_file
     IF ( input_file == ' ') GO TO 10
     INQUIRE ( FILE = input_file, EXIST = exst )
     IF ( .NOT. exst) THEN
        PRINT  '(A,": file not found")', TRIM(input_file)
        GO TO 10
     END IF
  ELSE IF ( nargs == 1 ) then
     CALL getarg (1,input_file)
  ELSE
     CALL errore( TRIM(prgname), 'too many arguments' , nargs )
  END IF
  RETURN
20 CALL errore( TRIM(prgname), 'reading file name' , 1 )
  !
END SUBROUTINE get_file
!
! Copyright (C) 2001-2004 PWSCF-FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE errore( calling_routine, message, ierr )
  !----------------------------------------------------------------------------
  !
  ! ... This is a simple routine which writes an error message to output: 
  ! ... if ierr <= 0 it does nothing, 
  ! ... if ierr  > 0 it stops.
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: calling_routine, message
    ! the name of the calling calling_routinee
    ! the output messagee
  INTEGER,          INTENT(IN) :: ierr
    ! the error flag
  LOGICAL                      :: exists
  !
  !
  IF ( ierr <= 0 ) RETURN
  !
  ! ... the error message is written un the "*" unit
  !
  WRITE( UNIT = *, FMT = '(/,1X,78("%"))' )
  WRITE( UNIT = *, &
         FMT = '(5X,"from ",A," : error #",I10)' ) calling_routine, ierr
  WRITE( UNIT = *, FMT = '(5X,A)' ) message
  WRITE( UNIT = *, FMT = '(1X,78("%"),/)' )
  !
  WRITE( *, '("     stopping ...")' )
  !
  RETURN
  !
END SUBROUTINE errore

! -----------------------------------------------------------------------------

