program fitModel
use SIMHYD
implicit none
character(256) :: dataFileName
type(simhydVarType) :: simhydType
integer(mik)::err
character(256)::message
character(*),parameter::procnam="fitmodel"
real(mrk), allocatable :: lowPar(:), highPar(:),para(:)
real(mrk), allocatable :: rain(:), EPT(:),obsQ(:),simQ(:)
character(256),allocatable::name(:)
integer(mik)::nWarmup,npar,nData
integer(mik)::status,i
real(mrk)::qObsMean,qObsSS,SSE,NSE
err=0
! Get the basic information for modelling
open (unit=1, file='fitModel.txt', status='old')
    read(1,*) nPar ! Number of parameters
    allocate(lowPar(nPar), highPar(nPar), name(nPar),para(nPar))
    do i = 1, nPar    ! Read data
       read(1,*) name(i), lowPar(i), highPar(i)  ! Parameter name, lower bound, upper bound
    end do
    read(1,*) dataFileName   ! Name of file with rainfall, potential evapotranspiration and observed runoff [mm/day]
    read(1,*) nWarmUp        ! Number of warmup steps
close (unit=1)
!
! Get input and observed response data
open (unit=1, file=dataFileName, status='old')
    read(1,*)   ! Skip header
    nData = 0
    do          ! Read to end of file to get number of data
        read(1,*,iostat=status)
        if (status /= 0) exit
        nData = nData + 1
    end do
    allocate(rain(nData), EPT(nData), obsQ(nData), simQ(nData))
    rewind(unit=1); read(1,*)  ! Now rewind file and read data 
    do i = 1, nData
        read(1,*) rain(i), EPT(i), obsQ(i)
    end do
    qObsMean = sum(obsQ(1+nWarmUp:nData))/real(nData-nWarmUp)  ! Get mean of obsrved runoff
    qObsSS = sum((obsQ(1+nWarmUp:nData)-qObsMean)**2)          ! Get sum of squares of obs runoff devaions from mean
close (unit=1)
! Assign parameters
simhydType%intCap=1.2_mrk              ! assign your number
simhydType%smsCap=100.0_mrk            ! assign your number
simhydType%sub=0.5_mrk                 ! assign your 
simhydType%crak=0.2_mrk                ! assign your 
simhydType%k=0.25_mrk                  ! assign your 
simhydType%sq=1.0_mrk                  ! assign your 
simhydType%coeff=500_mrk               ! assign your 
!
! Initialize initial snowstorage
simhydType%Q=zero                  ! initialize storage
simhydType%infilExcess=zero        ! initialize storage
simhydType%satExcess=zero          ! initialize storage
simhydType%baseFlow=zero           ! initialize storage
simhydType%sms=zero                ! initialize storage
simhydType%gw=zero                 ! initialize storage
! Warmup model so that memory of arbitary initial conditions is forgotten
do i = 1, nWarmUp
   call evolveSIMHYD(simhydType%Q,simhydType%infilExcess,simhydType%satExcess,simhydType%baseFlow,simhydType%sms,simhydType%gw,&
                   rain(i),EPT(i),&
                   simhydType%intCap,simhydType%smsCap,simhydType%sub,simhydType%crak,simhydType%k,simhydType%sq,simhydType%coeff,&
                   err,message)
   if(err/=0)then
      message="f-"//procnam//message
      write(*,*) message
      exit
   endif
   simQ(i) = simhydType%Q
  write(*,*) obsQ(i), simQ(i)
end do
! Compute R^2  for remiander of observed record
SSE = 0.0_mrk
do i = nWarmUp+1, nData
   call evolveSIMHYD(simhydType%Q,simhydType%infilExcess,simhydType%satExcess,simhydType%baseFlow,simhydType%sms,simhydType%gw,&
                   rain(i),EPT(i),&
                   simhydType%intCap,simhydType%smsCap,simhydType%sub,simhydType%crak,simhydType%k,simhydType%sq,simhydType%coeff,&
                   err,message)
   if(err/=0)then
      message="f-"//procnam//message
      write(*,*) message
      exit
   endif
   simQ(i) = simhydType%Q
   write(*,*) obsQ(i), simQ(i)
   SSE = SSE + (obsQ(i)-simQ(i))**2
end do
NSE =1.0_mrk - SSE/qObsSS
write(*,*) "NSE:" , NSE
PAUSE
end program fitModel