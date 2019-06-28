module start

implicit none

integer :: i,j,k,ii,jj,kk,iii,iiii,jjj,n,ng,o,p! general counters
integer :: ios,logp,t,et,tmax,etmax            ! more general counters
integer :: mzadhoc,capped,reco,linear          ! Some switchers
integer :: training,replicas                   ! If training=1 -> Training set, starting from W=0. Otherwise W from file
integer :: hillclimber                         ! If set to 1-> Strict hill-climber, deterministic selection.
integer :: positivecues                        ! If set to 1 cues are in range 0:maxepigen, otherwise range -maxepigen,maxepigen
integer :: ret,pp,im                           ! integers for calling external functions or modifyind datafiles.
integer :: lapso,intervals                     ! Intervals=number of intervals for data recording. Lapso: Generations in an interval.
integer :: ncels,replica,PD,EF                 ! PD=phenotypic dimensionality (number of traits),EF=Environmental factors
integer,allocatable :: thresholdsN(:)          ! thresholds in the N environments (for target switchings)
integer,allocatable :: premutWW(:,:)           ! Stores WW matrix before mutation (reversible if unstable GRN)
real*4 :: a,aa,b,c,q,u,v,x,y,z,m,deg           ! Real auxiliar numbers
real*4 :: fmax,sdev,ss,fmaxval,fmaxabs         ! Real numbers variables, maximum fitness in a eneration and in simulation
real*4, allocatable ::  prepattern(:,:),blocke(:,:)
real*4 :: conWW,conMZZ,maxepigen               ! Connectivity parameters for binary matrices. Maximum absolute value for env. cue
real*4 :: ginitial_det, ginitial_rand, initialW! Initial gene concentration (deterministic and stochastic) and W matrix weights
real*4 ,allocatable ::  thresholds(:)          ! thresholds in concentration of EFs (for target switchings)
real*4,allocatable  ::  premutW(:,:)           ! Stores W matrix before mutation (reversible if unstable GRN)
real*4,allocatable  ::  stab(:,:)              ! Stores the gene expresion during developmental time. To check stability.
real*4, parameter   ::  pi=3.1415926535, delta=0.00001, e=2.7182818284 ! some constants
character(len=48)   ::  arxiv,arxifin          ! for writting and reading datafiles
character(len=24)   ::  arxaux                 ! for writting and reading datafiles

type,public :: inds
  integer              :: ncels, ngs, sat      ! number of cells, number of genes and per cel and active genes, saturation time
  real*4               :: fitness(2)           ! individual fitness in t=0, absolute and normalized
  real*4 ,allocatable  :: w(:,:)               ! interaction strenghts
  real*4, allocatable  :: MZ(:,:)              ! determines importance of node for determining final trait(s), continuous (can be negative)
  integer,allocatable  :: ww(:,:)              ! interaction strenghts discreteee (connected or not?)
  integer,allocatable  :: MZZ(:,:)             ! determines if node is connected to output layer, binary, always 0 for nodes connected to input layer
  real*4 ,allocatable  :: g(:,:)               ! concentraction of each gene in each cell (cell,gen)
  real*4 ,allocatable  :: phen(:,:)            ! phenotype
  real*4 ,allocatable  :: epigen(:,:)          ! environmental effects
end type

type,public :: gens
  real*4    :: deg
end type

type(inds), public, allocatable :: ind(:),indt(:)         ! individual based
type(gens), public, allocatable :: gen(:)

contains

subroutine arxivpublic ; end subroutine arxivpublic       ! just to acess this module

subroutine inicial                                        ! allocate the matrices for cells and individuals

p=2                                                       ! number of individuals (must be an EVEN NUMBER !!!)
if(mod(p,2).ne.0)then ; write(*,*)'p must be an EVEN NUMBER' ; end if
logp=1+int(log(real(p))/log(2d0))
tmax=20                                                   ! developmental time
etmax=1.0E5                                               ! evolutionary time
EF=1                                                      ! EF=Number of environmental factors (inputs)
n=6                                                       ! number of different environments
ng=4                                                      ! initial number of genes
PD=1                                                      ! phenotypic dimensionality (number of traits)
sdev=0.005                                                 ! standard deviation for the mutator algorithm
ss=0.2                                                    ! selection strenght
reco=0                                                    ! recombination; 1=yes, 0=no
capped=0                                                  ! If 1, GRN (W-matrix) values are (-1,1); if 0, unconstrained values.
training=0                                                ! If 1 -> Training set, starting from W=0. Otherwise Test set (W from file).
replicas=10                                                ! Number of replicates
conWW=1.0                                                 ! Probability of having non-zero entries in WW  matrix (0,1)
conMZZ=0.5                                                ! Probability of having non-zero entries in MZZ matrix (0,1)
intervals=10                                               ! Number of intervals for data recording.
lapso=int(etmax/intervals)  	            					      ! Lapso: Generations in an interval.
hillclimber=1                                             ! If set to 1-> Strict hill-climber, deterministic selection.If 0->Probabilistic NS.
ginitial_det=0.5                                          ! Gene concentrations at the start of development, deterministic value
ginitial_rand=0.0                                           ! Gene concentrations at the start of development, randomized value
maxepigen=0.5                                             ! Maximum absolute value for env. cue
positivecues=0                                            ! If 1 then cues in range 0:maxepigen, otherwise range -maxepigen,maxepigen
initialW=5.0E-4                                           ! Connection weights at the start of the simulation
mzadhoc=1                                                 ! If 0: Mz and Mzz matrices read/generated normally. If 1: Mz and Mzz matrices uploaded from external file. For all P and Training.
linear=1                                                  ! If 1, use linear activation function. Otherwise, use logit

if(intervals.gt.etmax)then ; write(*,*)'Etmax MUST BE greater than Intervals' ; end if
if ((hillclimber.eq.1).and.(p.gt.2))then
write(*,*)'WARNING! You are using large (p>2) populations with a hill-climber NS !'; end if

if(allocated(ind))then
  deallocate(ind) ; deallocate(indt) ; deallocate(gen)    ! Allocating variables
  deallocate(prepattern) ; deallocate(blocke)             ! Allocating variables
  deallocate(thresholds) ; deallocate(thresholdsN)        ! Allocating variables
  deallocate(premutW)    ; deallocate(premutWW)           ! Allocating variables
  deallocate(stab)                                        ! Allocating variables
end if                                                    ! Allocating variables
allocate(ind(p),indt(p)) ; allocate(gen(ng))              ! Allocating variables
allocate(thresholds(EF)) ; allocate(thresholdsN(EF))      ! Allocating variables
allocate(premutW(ng,ng)) ; allocate(premutWW(ng,ng))      ! Allocating variables
allocate(stab(n,ng))                                      ! Allocating variables

allocate(prepattern(n,ng))
allocate(blocke(PD,n))                                     ! target dimensionality
blocke=0.0
stab=0.0

if((training.eq.0).and.(replica.lt.1))then ; return ; end if

!!!!!!!!!!!!!!!!!!!!!!                                    ! open file for seting a population if we are in TEST SET.
  if(training.ne.1)then                                   ! Introducing manually the filename from where the system uploads the population
           !GRN_1234_R12_T1234                            ! Follow this template
    arxaux='GRN_312452_C_4_R_1_T00'
           !123456789012345678
    arxiv(1:3)='GRN' ; arxiv(4:22)='_' ; arxiv(23:44)=arxaux(1:22)
    arxiv(45:48)='.dat'                                             ! composing filename
    do im=1,44 ; if (arxiv(im:im)==" ") arxiv(im:im)="_" ; end do   ! composing filename
    open(9000,file='files/'//arxiv,action='read',iostat=ios)
    do i=1,13 ;  read(9000,*)  ;  end do                    ! skip first human readable lines about parameters.
  else
    arxiv(23:25)='GRN' ; arxiv(26:44)='_'
  end if
!!!!!!!!!!!!!!!!!!!!!!!

do i=1,p                                                  ! for all individuals in the population

  ind(i)%ncels=n          ; ind(i)%ngs=ng                 ! initial conditions
  ind(i)%fitness=0.0                                      ! initial conditions

  allocate(ind(i)%w(ind(i)%ngs,ind(i)%ngs))               ! it contains all potential genes
  allocate(ind(i)%ww(ind(i)%ngs,ind(i)%ngs))              ! it contains active interactions
  allocate(ind(i)%g(ind(i)%ncels,ind(i)%ngs))             ! it contaisn all potential genes
  allocate(ind(i)%phen(PD,ind(i)%ncels))                  ! phenotype
  allocate(ind(i)%epigen(ng,n))                           ! environmental effects

  allocate(ind (i)%MZ(ind(i)%ngs,PD))  ;  allocate(ind (i)%MZZ(ind(i)%ngs,PD))
  allocate(indt(i)%MZ(ind(i)%ngs,PD))  ;  allocate(indt(i)%MZZ(ind(i)%ngs,PD))

  indt(i)%ncels=n          ; indt(i)%ngs=ng               ! initial conditions
  indt(i)%fitness=0.0      ; ind(i)%epigen=0.0            ! initial conditions

  allocate(indt(i)%w(ind(i)%ngs,ind(i)%ngs))              ! it contains all potential genes
  allocate(indt(i)%ww(ind(i)%ngs,ind(i)%ngs))             ! it contains all potential genes
  allocate(indt(i)%g(ind(i)%ncels,ind(i)%ngs))            ! it contaisn all potential genes
  allocate(indT(i)%phen(PD,ind(i)%ncels))                 ! phenotype

!!!!!!!!!!!!!!!!!!!! WARNING. MANUAL IMPLEMENTATION !!!!!!!!!!!!!!!!!!!!!!! WARNING. MANUAL IMPLEMENTATION !!!!!!!!

  !ind(i)%epigen(1:2,1)=(/-1.0, 1.0/)                     ! ENVIRONMENTAL FACTORS IN ENVIRONMENT 1
  !ind(i)%epigen(1:2,2)=(/ 1.0,-1.0/)                     ! ENVIRONMENTAL FACTORS IN ENVIRONMENT 2
  thresholds(1)=0.0 ; thresholdsN(1)=0                    ! Re-do for EF>1 !!! WARNING !!!!
  do ii=1,n                                               ! Linear dacaying function [1,-1](equivalent to any non-linear function with a threshold)
    ind(i)%epigen(1,ii)=1.0-(real(ii-1)/real(n-1))     ! Create linear input in range 1:0
    if(positivecues.eq.0)then
      ind(i)%epigen(1,ii)=(2.0*ind(i)%epigen(1,ii))-1.0  ! Rescale input to range 1:-1
    end if
    ind(i)%epigen(1,ii)=maxepigen*ind(i)%epigen(1,ii)  ! Rescale input to maxepigen value (maxepigen:-maxepigen or maxepigen:0)
    if(ind(i)%epigen(1,ii).le.thresholds(1))then          ! Finding the "cell index" where the threshold is applied. Re-do for EF>1 !!! WARNING !!!
      if(thresholdsN(1).eq.0)then                         ! Finding the "cell index" where the threshold is applied. Re-do for EF>1 !!! WARNING !!!
      thresholdsN(1)=ii ; end if                          ! Finding the "cell index" where the threshold is applied. Re-do for EF>1 !!! WARNING !!!
    end if                                                ! Finding the "cell index" where the threshold is applied. Re-do for EF>1 !!! WARNING !!!
  end do                                                  !
! Input targets for each trait across all environments
!blocke(1,1:n)=(/ 0.6, 0.1, 0.3, 0.7, 0.9, 0.4/)    ! Problem A
!blocke(1,1:n)=(/ 0.4, 0.9, 0.7, 0.3, 0.1, 0.6/)    ! Problem B (reverse order of A)
!blocke(1,1:n)=(/ 0.3, 0.1, 0.2, 0.4, 0.5, 0.2/)    ! Problem C (A with phenotype/2)
!blocke(1,1:n)=(/ 0.3, 0.4, 0.5, 0.6, 0.7, 0.8/)    ! Problem D (linear1)
!blocke(1,1:n)=(/ 0.8, 0.7, 0.6, 0.5, 0.4, 0.3/)    ! Problem E (D reversed)
blocke(1,1:n)=(/ 0.3, 0.3, 0.8, 0.8, 0.2, 0.2/)    ! Problem F (step funciton w 3 steps)
!blocke(1,1:n)=(/ 0.5, 0.5, 0.8, 0.8, 0.5, 0.5/)    ! Problem G (step funciton w 2 steps)
!blocke(1,1:n)=(/ 0.8, 0.8, 0.5, 0.5, 0.8, 0.8/)    ! Problem H (step funciton w 2 steps, not working due to neg concentration)
!blocke(1,1:n)=(/ 0.5, 0.5, 0.2, 0.2, 0.5, 0.5/)    ! Problem I (step funciton w 2 steps, lower phenotype)
!blocke(1,1:n)=(/ 0.2, 0.2, 0.5, 0.5, 0.2, 0.2/)    ! Problem J (reverse of I)
!blocke(1,1:n)=(/ 0.2, 0.2, 0.5, 0.5, 0.7, 0.7/)    ! Problem K (monotonic)
!blocke(1,1:n)=(/ 0.7, 0.7, 0.5, 0.5, 0.2, 0.2/)    ! Problem L (reverse K)
!blocke(1,1:n)=(/ 0.8, 0.4, 0.5, 0.6, 0.7, 0.3/)    ! Problem M (a with shifted intercept)
!blocke(1,1:n)=(/ 0.95, 0.94, 0.93, 0.92, 0.91, 0.96/)    ! Problem N (nonplastic, high intercept)
!blocke(1,1:n)=(/ 0.6, 0.6, 0.9, 0.9, 0.6, 0.6/)    ! Problem O (same shape as GI, different values)

!!!!!!!!!!!!!!!!!!!! WARNING. MANUAL IMPLEMENTATION !!!!!!!!!!!!!!!!!!!!!!! WARNING. MANUAL IMPLEMENTATION !!!!!!!!

  do ii=1,ind(i)%ncels                                    ! setting cell topology
    ind(i)%g(ii,:)=0.0
    do iii=1,ind(i)%ngs                                   ! setting initial gene concentrations in each cell
       call random_number(x)
       ind(i)%g(ii,iii)=ginitial_det + ginitial_rand*(x-0.5)    ! small or small noisy initial gene expression
       prepattern(ii,iii)=ind(i)%g(ii,iii)
    end do
  end do

  ind(i)%w=0.0  ; ind(i)%ww=0 ; gen(:)%deg=0.0
  ind(i)%MZ=0.0 ; ind(i)%MZZ=0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CREATING / READING POPULATIONS
  if(training.eq.1)then                                                 ! creating random population in t=0 in training phase
    do iii=1,ind(i)%ngs
      do iiii=1,ind(i)%ngs
        call random_number(x)                                           ! Random W matrix in t=0
        ind(i)%w(iii,iiii)=(x-0.5)*initialW                             ! Use small random weights centered around zero
        !call random_number(x)
        !if(x.le.conWW)then ; jjj=1 ; else ; jjj=0 ; end if
        jjj=1
        ind(i)%ww(iii,iiii)=jjj                                          ! Random WW matrix in t=0
      end do
      if(i.eq.1)then
        do iiii=1,pd
          call random_number(x)
          ind(1)%MZ(iii,iiii)=1.0-(2.0*x)                               ! Random MZ matrix in t=0
          call random_number(x)
          if(x.le.conMZZ)then ; jjj=1 ; else ; jjj=0 ; end if
          ind(1)%MZZ(iii,iiii)=jjj                                      ! Random MZZ matrix in t=0
        end do
      end if
    end do
    ind(i)%MZ=ind(1)%MZ		     										! homogeneous MZ matrix for all individuals
    ind(i)%MZZ=ind(1)%MZZ                                               ! homogeneous MZ matrix for all individuals
  else
    do iii=1,ind(i)%ngs
      read(9000,*) ind(i)%w(iii,1:ind(i)%ngs)                           ! uploading population W from file
    end do
    do iii=1,ind(i)%ngs
      read(9000,*) ind(i)%ww(iii,1:ind(i)%ngs)                           ! uploading population Ww from file
    end do
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end do

if(training.ne.1)then                                                   ! Test set.
  do iii=1,ind(1)%ngs
    read(9000,*) ind(1)%MZ(iii,1:pd)                                    ! uploading MZ  for the whole population from file
  end do
  do iii=1,ind(1)%ngs
    read(9000,*) ind(1)%MZZ(iii,1:pd)                                   ! uploading MZZ for the whole population from file
  end do
  do iii=1,p
    ind(iii)%MZ=ind(1)%MZ	    										! homogeneous MZ matrix for all individuals
    ind(iii)%MZZ=ind(1)%MZZ                                             ! homogeneous MZ matrix for all individuals
  end do
  close(9000)
end if

 indt=ind                                                               ! just initializing all the matrices for time steps

 if(mzadhoc.eq.1)then                                                   ! Uploading Mz and Mzz matrices from external file IFF Mzadhoc==1.
   open(1133,file='files/mzadhoc.dat',action='read',iostat=ios)         ! Uploading Mz and Mzz matrices from external file IFF Mzadhoc==1.
   do iii=1,ind(1)%ngs
     read(1133,*) ind(1)%MZ(iii,1:pd)                                    ! uploading MZ  for the whole population from file
   end do
   do iii=1,ind(1)%ngs
     read(1133,*) ind(1)%MZZ(iii,1:pd)                                   ! uploading MZZ for the whole population from file
   end do
   do iii=1,p
     ind(iii)%MZ=ind(1)%MZ	    								 		! homogeneous MZ matrix for all individuals
     ind(iii)%MZZ=ind(1)%MZZ                                            ! homogeneous MZ matrix for all individuals
   end do
   close(1133)
 end if

 do i=1,p
   do ii=1,n                                              ! (ncels)! all individuals must have identical prepatterns
     ind(i)%g(ii,1:ng)=ind(1)%g(ii,1:ng)
     prepattern(ii,1:ng)=ind(1)%g(ii,1:ng)
   end do
 end do

end subroutine inicial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module start
