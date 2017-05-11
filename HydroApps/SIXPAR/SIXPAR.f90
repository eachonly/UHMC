!******************************************************************
! (C) Copyright 2007-2020  ---  Youwei QIN  ---  All rights reserved
!******************************************************************
!                   SIXPAR MODEL
!******************************************************************
module SIXPAR
! Purpose: Contains a Fortran-95 implementation of the SIXPAR model
!          It is based on the Fortran-77 implementation of the SIXPAR model by VIJAI K. GUPTA, 1984
! ---
! Programmer: Youwei QIN
! Kreated: circa Nov 2013 AD, Uni of Adelaide.
! ---
! Comments:
! ---
! 1. please refer to the paper Duan et.al 1992 for detail description,
!    doi: 10.1029/91wr02985
! 2. The F95 version of SIXPAR hasn't been fully compared against F77 version 
! 3. The author suggest using F95 version
! ---
implicit none
! type definitions
! variable definitions
!**********************************
! Procedure availability
private ! to protect type-specific procedure names
!-----
! List of public
public:: mik, mrk, evolveSIXPAR,SIXPARVarType,zero,one
!----------------------------------------------------
integer,parameter::mik=8,mrk=8
real(mrk),parameter::undefRN=-999999999._mrk
real(mrk),parameter::zero=0.0_mrk
real(mrk),parameter::one=1.0_mrk
!----------------------------------------------------
! * SIXPAR parameter-state variable structure
type SIXPARVarType
! - parameters
  real(mrk)::um  =undefRN ! parDefSIXPAR(1)
  real(mrk)::uk  =undefRN ! parDefSIXPAR(2)
  real(mrk)::bm  =undefRN ! parDefSIXPAR(3)
  real(mrk)::bk  =undefRN ! parDefSIXPAR(4)
  real(mrk)::z   =undefRN ! parDefSIXPAR(5)
  real(mrk)::x   =undefRN ! parDefSIXPAR(6)
! - states
  real(mrk)::q            =undefRN
  real(mrk)::r            =undefRN
  real(mrk)::s            =undefRN
  real(mrk)::b            =undefRN
  real(mrk)::us           =undefRN
  real(mrk)::bs           =undefRN
endtype SIXPARVarType
!----------------------------------------------------   
contains
pure subroutine evolveSIXPAR(q,r,s,b,us,bs,&
  rainObs,ETpot,&
  um,uk,bm,bk,z,x,&
  err,message)
! Purpose: Evolves a single step of the full sixpar model,
! ---
! Programmer: Youwei QIN
! Kreated: circa Sep 2014 AD, Uni of Newcastle.
! ---
! Performance
! IN:Q,r,s,b,us,bs,rainObs,rainObs,ETpot,um,uk,bm,bk,z,x
! OUT:Q,r,s,b,us,bs,err,message
! ---
! Comments:
!     Six-Parameter (SIXPAR) model based on Burnash Sacramento 
!     Soil Moisture Accounting (SAC-SMA) model structure with
!     option of using reparameterized percolation parameters
! 
!     VARIABLES DEFINED:
!         P = PRECITATION
!         Q = MODEL OUTPUT
!         R = OVERLAND FLOW
!         S = SUBSURFACE FLOW
!         B = BASEFLOW
!         PERC = PERCOLATION RATE
!         PPERC = PERCOATION DEMAND
!         US = UPPER ZONE CONTENT
!         BS = BOTTOM ZONE CONTENT
!     PARAMETERS DEFINED:
!         UM = UPPER ZONE CAPACITY
!         UK = UPPER ZONE STORAGE COEFFICIENT
!         BM = BOTTOM ZONE CAPACITY
!         BK = BOTTOM ZONE STORAGE COEFFICIENT
!         X, Z = COEFFICIENTS IN PERCOLATION EQUATION
!         A, RKU, RKB = PERCOLATION PARAMETERS
implicit none
integer(mik),parameter::iperc=0
! dummies
real(mrk),intent(inout)::Q,r,s,b,us,bs    ! states
real(mrk),intent(in)::rainObs,ETpot       ! forcings
real(mrk),intent(in)::um,uk,bm,bk,z,x     ! parameters
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="sixpar_engine_youwei"
real(mrk):: a,pet,rain,rlzdr,perc,yy,zz,d
!Initialize the local parameters
a=z;pet=ETpot;rain=rainObs !ETpot not used in sixpar
r=zero; s=zero;b=zero
! Start procedure here
err=0; message=procnam//"/ok"
!Illustration
! Part I Calculate Percolation
us=us+rain
rlzdr=(bm-bs)/bm
if (rlzdr <= a) then
   yy=us*bm*bk/um
   zz=(bm-bs)/(bm*a)
   if (zz <= zero) then
      zz=zero
   else
      zz=zz**x
   end if
   perc=yy+zz*(us-yy)
!  ================================================================
!  If iperc = 1, set the max perc(i) equal to available water us(i)
!  If iperc = 0, perc(i) > us(i), causing negative flow
!  ================================================================
   if (iperc == 1) then
      if (perc > us) then
         perc=us
      end if
   end if
   us=us-perc
else
   perc=us
   us=0
end if

! Part II update the bottom layer water storage,and calculate outflow
bs=bs+perc
if (bs > bm) then
   d=bs-bm
   b=bm*bk
   bs=(bm-b)+d
   if (bs>bm) then
      d=bs-bm
      bs=bm
      us=us+d
   end if
else
   b=bs*bk
   bs=bs-b
end if
!Part III update the upper layer water storage, and calculate outflow
if (us> um) then
   r=us-um
   us=um
else
   r=0      !saturation excess runoff r
end if
s=us*uk
us=us-s
q=(r+s)+b
! End procedure here
endsubroutine evolveSIXPAR
!----------------------------------------------------
endmodule SIXPAR