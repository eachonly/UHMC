!******************************************************************
! (C) Copyright 2017-2027  ---  Youwei Qin ---  All rights reserved
!******************************************************************
module SnowMod
! Purpose: Contains a Fortran-95 implementation of the SNOW model
! by Youwei QIN
! ---
! Programmer: Youwei QIN
! Kreated: circa Nov 2013 AD, Uni of Adelaide.
! ---
! Comments:
! 1. Based on the description in paper Kavetski,Kuczera (2007),
!    doi: 10.1029/2006wr005195
! ---
implicit none
!**********************************
! Procedure availability
private ! to protect type-specific procedure names
!-----
! List of public
public:: mik, mrk, SnowPar, SnowFlux, evolveSNOW
!----------------------------------------------------
integer,parameter::mik=8,mrk=8
real(mrk),parameter::zero=0._mrk
!----------------------------------------------------
type SnowPar
  real(8) :: tempT               ! temperature related to snow melt
  real(8) :: factorK             ! snow melt factor
end type SnowPar
type SnowFlux
  real(8) :: SnowStor,melt       ! snow storage, and melt discharge
end type SnowFlux
!----------------------------------------------------   
contains
subroutine evolveSNOW(melt,SnowStor,&
  rainObs,TempObs,&
  tempT,factorK,&
  err,message)
! Purpose: Evolves a single step of the full SNOW model,
! ---
! Programmer: Youwei QIN
! Kreated: circa Nov 2013 AD, Uni of Adelaide.
! ---
! Performance
! IN:melt, SnowStor, rainObs, TempObs, tempT, factorK, err, message
! OUT:melt, SnowStor, err, message
! ---
! Comments:
implicit none
! dummies
real(mrk),intent(inout)::melt,SnowStor      ! states
real(mrk),intent(in)::rainObs,TempObs       ! forcings
real(mrk),intent(in)::tempT,factorK         ! parameters
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="evolveSNOW"
real(mrk)::P
!Initialize the local parameters
P=zero
! Start procedure here
err=0; message=procnam//"/ok"
call checkFeasSNOW(melt,SnowStor,err,message)
if(err/=0)then
  err=10; message="f-"//procnam//"/BC/&"//message; return
endif
!Daily Snow balance
if (TempObs<=tempT) THEN            ! For simplicity, tempT is the threshold for melting/snowing
        P=rainObs
        SnowStor=SnowStor+P
        melt=zero
    ELSE
        P=zero
        melt=min(factorK*(TempObs-tempT),SnowStor)     ! melt bounded by melting or snow storage
        SnowStor=max(SnowStor-melt,zero)               ! use max(XX,zero) as safeguard
    END IF
! End procedure here
endsubroutine evolveSNOW
!----------------------------------------------------
pure subroutine checkFeasSNOW(melt,SnowStor,err,message)
! Purpose: Checks feasibility of SNOW inputs.
! Kreated:  2 Nov 2013 AD, Kathelen Lumely Colledge, SA
implicit none
! dummies
real(mrk),intent(in)::melt,SnowStor
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="checkFeasSNOW"
character(*),parameter::fmt1='(a,es15.8,a)'
! Start procedure here
err=0
if(melt<zero)then
  write(message,fmt1)"f-"//procnam//"/badIni[melt(",melt,")<0.0]"
  err=10; return
elseif(SnowStor<zero)then
  write(message,fmt1)"f-"//procnam//"/badIni[SnowStor(",SnowStor,")<0.0]"
  err=20; return
endif
! End procedure here
endsubroutine checkFeasSNOW
!----------------------------------------------------
end module SnowMod
