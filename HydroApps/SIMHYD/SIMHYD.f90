!******************************************************************
! (C) Copyright 2007-2020  ---  Dmitri Kavetski  ---  All rights reserved
!******************************************************************
!                   SIMHYD MODEL
!******************************************************************
module SIMHYD
! Purpose: Contains a Fortran-95 implementation of the SIMHYD model
! by Youwei QIN
! ---
! Programmer: Youwei QIN
! Kreated: circa Nov 2013 AD, Uni of Adelaide.
! ---
! Comments:
! ---
! 1. Based on the description in paper Chiew(2005),
!    doi: 10.1.1.375.8824
! 2. This version of SIMHYD has been fully compared against
! ---
implicit none
! type definitions
! variable definitions
!**********************************
! Procedure availability
private ! to protect type-specific procedure names
!-----
! List of public
public:: mik, mrk, evolveSIMHYD,simhydVarType,zero,one
!----------------------------------------------------
integer,parameter::mik=8,mrk=8
real(mrk),parameter::undefRN=-999999999._mrk
real(mrk),parameter::zero=0.0_mrk
real(mrk),parameter::one=1.0_mrk
integer(mik),parameter::et_maxpet=0,et_lins=1
integer(mik),parameter::et_def=et_maxpet
!----------------------------------------------------
! * SIMHYD parameter-state variable structure
type simhydVarType
! - parameters
  real(mrk)::intCap  =undefRN ! parDefSIMHYD(1)
  real(mrk)::smsCap  =undefRN ! parDefSIMHYD(2)
  real(mrk)::sub     =undefRN ! parDefSIMHYD(3)
  real(mrk)::crak    =undefRN ! parDefSIMHYD(4)
  real(mrk)::k       =undefRN ! parDefSIMHYD(5)
  real(mrk)::sq      =undefRN ! parDefSIMHYD(6)
  real(mrk)::coeff   =undefRN ! parDefSIMHYD(7)
! - states
  real(mrk)::Q            =undefRN
  real(mrk)::infilExcess  =undefRN
  real(mrk)::satExcess    =undefRN
  real(mrk)::baseFlow     =undefRN
  real(mrk)::sms          =undefRN
  real(mrk)::gw           =undefRN
endtype simhydVarType
!----------------------------------------------------   
contains
!----------------------------------------------------
pure subroutine evolveSIMHYD(Q,infilExcess,satExcess,baseflow,sms,gw,&
  rainObs,ETpot,&
  intCap,smsCap,sub,crak,k,sq,coeff,&
  err,message)
! Purpose: Evolves a single step of the full SIMHYD model,
! ---
! Programmer: Youwei QIN
! Kreated: circa Nov 2013 AD, Uni of Adelaide.
! ---
! Performance
! IN:Q,infilExcess,satExcess,baseflow,sms,gw,rainObs,ETpot,intCap,smsCap,sub,crak,k,sq,coeff
! OUT:Q,infilExcess,satExcess,baseflow,sms,gw,err,message
! ---
! Comments:
implicit none
! dummies
real(mrk),intent(inout)::Q,infilExcess,satExcess,baseflow,sms,gw  ! states
real(mrk),intent(in)::rainObs,ETpot       ! forcings
real(mrk),intent(in)::intCap,smsCap,sub,crak,k,sq,coeff  ! parameters
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="evolveSIMHYD"
real(mrk)::intMax,int,thrufall,sFrac,pot,et,recharge,smf,rmo,rmoPot,infilAftFlow,pet,rain
integer(mik),parameter::imeth=0 ! 0,1=(good), 2=smooth(bad), -1=noETrain
integer(mik)::etMeth
!Initialize the local parameters
intMax=zero;int=zero;thrufall=zero;sFrac=zero;pot=zero;et=zero;recharge=zero;smf=zero;rmo=zero;rmoPot=zero;infilAftFlow=zero
pet=ETpot;rain=rainObs
! Start procedure here
!initial Default method for calculate ET
etMeth=et_def
err=0; message=procnam//"/ok"
call checkFeasSIMHYD(sms,gw,smsCap,err,message)
if(err/=0)then
  err=10; message="f-"//procnam//"/BC/&"//message; return
endif
!Illustration
! Daily water balance computations
! Interception 
intMax = MIN(intCap,pet)
int = MIN(intMax,rain)
thrufall = rain - int
! Infiltration excess runoff
sFrac = sms/smsCap
rmoPot=coeff*EXP(-sq*sFrac)
rmo = MIN(rmoPot, thrufall)
infilExcess = thrufall - rmo
! Saturation excess runoff
satExcess = rmo*sub*sFrac
infilAftFlow=rmo-satExcess
! Soil moisture store
recharge = infilAftFlow*crak*sFrac
!smf is the soil input
smf = infilAftFlow - recharge
sms=sms+smf
sFrac=sms/smsCap
gw=gw+recharge
!check if the soil moistrue exceeds
if(sFrac>1.0_mrk) then
    gw=gw+sms-smsCap
    sms=smsCap
    sFrac=one
end if
!to calculate the ground flow
baseFlow=k*gw
gw=gw-baseFlow
!to calculate the evaportraspritation
pot = pet - int
selectcase(etMeth)          
case(et_maxpet)             
  et = MIN(sms,10.0_mrk*sFrac,pot)
case(et_lins)               
  et = MIN(sms,10.0_mrk*sFrac,sFrac*pot)
case default
  write(message,'(a,i0,a)')"f-"//procnam//"/unknown[etMeth=",etMeth,"]"
  err=20; return
endselect
sms = sms - et
! Runoff
q = infilExcess + baseFlow + satExcess
! End procedure here
endsubroutine evolveSIMHYD
!----------------------------------------------------
pure subroutine checkFeasSIMHYD(sms,gw,smsCap,err,message)
! Purpose: Checks feasibility of SIMHYD inputs.
! Kreated:  2 Nov 2013 AD, Kathelen Lumely Colledge, SA
implicit none
! dummies
real(mrk),intent(in)::sms,gw,smsCap
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="checkFeasSIMHYD"
character(*),parameter::fmt1='(a,es15.8,a)'
character(*),parameter::fmt2='(a,es15.8,a,es15.8,a)'
! Start procedure here
err=0
if(sms>smsCap)then
  write(message,fmt2)"f-"//procnam//"/badIni[sms>smsCap;(",sms,")>(",smsCap,")]"
  err=10; return
elseif(sms<zero)then
  write(message,fmt1)"f-"//procnam//"/badIni[sms(",sms,")<0.0]"
  err=20; return
elseif(gw<zero)then
  write(message,fmt1)"f-"//procnam//"/badIni[gw(",gw,")<0.0]"
  err=30; return
endif
! End procedure here
endsubroutine checkFeasSIMHYD
!----------------------------------------------------
endmodule SIMHYD