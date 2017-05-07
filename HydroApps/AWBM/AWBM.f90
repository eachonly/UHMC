!******************************************************************
! (C) Copyright 2017-2027  ---  Youwei Qin ---  All rights reserved
!******************************************************************
!                   Australia Water Balance Model(AWBM)
!******************************************************************
module AWBM
! Purpose: Contains a Fortran-95 implementation of the AWBM model
! ---
! Programmer: Youwei QIN
! Kreated: circa Nov 2013 AD, Uni of Adelaide.
! ---
! Comments:
! ---
! 1. Based on the description in paper Boughton(2003),
!    doi: 10.1016/j.envsoft.2003.10.007
! 2. This version of AWBM has been fully compared against
! ---
implicit none
!**********************************
! Procedure availability
private ! to protect type-specific procedure names
!-----
! List of public
public:: mik, mrk, evolveAWBM,awbmVarType,zero,one
!----------------------------------------------------
integer,parameter::mik=8,mrk=8
real(mrk),parameter::undefRN=-999999999._mrk
real(mrk),parameter::zero=0.0_mrk
real(mrk),parameter::one=1.0_mrk
!----------------------------------------------------
! * AWBM parameter-state variable structure
type awbmVarType
! - parameters
  real(mrk)::A1   =undefRN  ! parDefAWBM(1)
  real(mrk)::A2   =undefRN  ! parDefAWBM(2)
  real(mrk)::BFI  =undefRN  ! parDefAWBM(3)
  real(mrk)::C1   =undefRN  ! parDefAWBM(4)
  real(mrk)::C2   =undefRN  ! parDefAWBM(5)
  real(mrk)::C3   =undefRN  ! parDefAWBM(6)
  real(mrk)::kbase=undefRN  ! parDefAWBM(7)
  real(mrk)::ksurf=undefRN  ! parDefAWBM(8)
! - states
  real(mrk)::Q         =undefRN
  real(mrk)::SurfFlow  =undefRN
  real(mrk)::BaseFlow  =undefRN
  real(mrk)::storage1  =undefRN
  real(mrk)::storage2  =undefRN
  real(mrk)::storage3  =undefRN
  real(mrk)::SurfStor  =undefRN
  real(mrk)::BaseStor  =undefRN
endtype awbmVarType
!----------------------------------------------------   
contains
subroutine evolveAWBM(Storage1,Storage2,Storage3,SurfStor,BaseStor,SurfFlow,BaseFlow,Q,&
  rainObs,ETpot,&
  A1,A2,BFI,C1,C2,C3,kbase,ksurf,&
  err,message)
! Purpose: Evolves a single step of the full AWBM model,
! ---
! Programmer: Youwei QIN
! Kreated: circa Nov 2013 AD, Uni of Adelaide.
! ---
! Performance
! IN:Storage1,Storage2,Storage3,SurfStor,BaseStor,SurfFlow,BaseFlow,Q,rainObs,ETpot,A1,A2,BFI,C1,C2,C3,kbase,ksurf
! OUT:Storage1,Storage2,Storage3,SurfStor,BaseStor,SurfFlow,BaseFlow,Q,err,message
! ---
! Comments:
implicit none
! dummies
real(mrk),intent(inout)::Storage1,Storage2,Storage3,SurfStor,BaseStor,SurfFlow,BaseFlow,Q  ! states
real(mrk),intent(in)::rainObs,ETpot                                                        ! forcings
real(mrk),intent(in)::A1,A2,BFI,C1,C2,C3,kbase,ksurf                                       ! parameters
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="evolveAWBM"
real(mrk)::A11,A22,A33,store1,store2,store3,satexcess
!Initialize the local parameters
A11=zero;A22=zero;A33=zero;store1=zero;store2=zero;store3=zero;satexcess=zero
A11=A1
A22=A2
A33=one-A11-A22
! Start procedure here
err=0; message=procnam//"/ok"
!Illustration
!Daily water balance computation for each bucket (three buckets in this model)
!If moisture int he store becomes negative, it is reset to zero
!If the moisture in the store exceeds the capacity of the store, 
!the moisture in excess becomes runoff and the store is reset to the store capacity
!Daily water balance for bucket1
store1=storage1+rainObs-ETpot
if (store1<=zero) then
    storage1=zero
else if (store1<C1) then
    storage1=store1
else
    storage1=C1
    satexcess=satexcess+A11*(store1-storage1)
end if
!Daily water balance for bucket2
store2=storage2+rainObs-ETpot
if(store2<=zero)then
    storage2=zero
else if(store2<C2)then
    storage2=store2
else
    storage2=C2
    satexcess=satexcess+A22*(store2-storage2)
end if
!Daily water balance for bucket3
store3=storage3+rainObs-ETpot
if(store3<=zero)then
    storage3=zero
else if(store3<C3) then
    storage3=store3
else
    storage3=C3
    satexcess=satexcess+A33*(store3-storage3)
end if
!The process to generate surfaceflow and baseflow
!The excess flow is divided into surface runoff recharge and baseflow recharge
!and the attenuate of the arrival of the generated runoff is represented by
!the surfacestorage store and the basestorage store
!The storage is simulated by linear resorvoir, it is solved by forward differencing method
SurfStor=SurfStor+satexcess*BFI
SurfFlow=ksurf*SurfStor
SurfStor=SurfStor-SurfFlow
BaseStor=BaseStor+satexcess*(one-BFI)
BaseFlow=kbase*BaseStor
BaseStor=BaseStor-BaseFlow
Q=SurfFlow+BaseFlow
! End procedure here
endsubroutine evolveAWBM
!----------------------------------------------------
end module AWBM
