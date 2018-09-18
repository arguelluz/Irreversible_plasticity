module start

implicit none

integer :: i,j,k,ii,jj,kk,iii,iiii,jjj,n,ng,o,p! general counters
integer :: ios,logp,t,et,tmax,etmax,capped,reco! general counters
integer :: training,replicas                   ! If training=1 -> Training set, starting from W=0. Otherwise W from file
integer :: ret,pp,im                           ! integers for calling external functions or modifyind datafiles. 
integer :: lapso,intervals                     ! Intervals=number of intervals for data recording. Lapso: Generations in an interval.
integer :: ncels,replica,PD,EF                 ! PD=phenotypic dimensionality (number of traits),EF=Environmental factors
real*4 :: a,aa,b,c,q,u,v,x,y,z,m,deg           ! Real auxiliar numbers 
real*4 :: fmax,fmaxx,fmaxabs,sdev,ss           ! Real numbers variables
real*4 :: conWW,conMZZ                        ! Connectivity arameters for binary matrices
real*4, allocatable ::  prepattern(:,:),block(:,:)
real*4, parameter   ::  pi=3.1415926535, delta=0.00001, e=2.7182818284 ! some constants 
character(len=40)   ::  arxiv,arxifin          ! for writting and reading datafiles
character(len=18)   ::  arxaux                 ! for writting and reading datafiles

type,public :: inds   
  integer              :: ncels, ngs, sat      ! number of cells, number of genes and per cel and active genes, saturation time
  real*4               :: fitness              ! individual fitness in t=0
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

                                                          ! in the future it will take the IC from external files                                                         
p=8                                                       ! number of individuals (must be an EVEN NUMBER !!!)
if(mod(p,2).ne.0)then ; write(*,*)'p must be an EVEN NUMBER' ; end if
logp=1+int(log(real(p))/log(2d0))
tmax=200                                                  ! developmental time
etmax=100                                                 ! evolutionary time      
EF=2                                                      ! EF=Number of environmental factors (inputs)    
n=2                                                       ! number of different environmental input nodes
ng=4                                                      ! initial number of genes
PD=2                                                      ! phenotypic dimensionality (number of traits)
sdev=0.1                                                  ! standard deviation for the mutator algorithm
ss=1.0                                                    ! selection strenght
reco=0                                                    ! recombination; 1=yes, 0=no
capped=0                                                  ! If 1, GRN (W-matrix) values are (-1,1); if 0, unconstrained values.
training=1                                                ! If 1 -> Training set, starting from W=0. Otherwise Test set (W from file).
replicas=2                                                ! Number of replicates 
conWW=1.0                                                 ! Probability of having non-zero entries in WW  matrix (0,1) 
conMZZ=1.0                                                ! Probability of having non-zero entries in MZZ matrix (0,1)  
intervals=4                                               ! Number of intervals for data recording. 
lapso=int(etmax/intervals)  						      ! Lapso: Generations in an interval.
if(intervals.gt.etmax)then ; write(*,*)'Etmax MUST BE greater than Intervals' ; end if

if(allocated(ind))then    ; deallocate(ind)   ; deallocate(indt)  ; deallocate(gen)
deallocate(prepattern)    ; deallocate(block) ; end if
allocate (ind(p),indt(p)) ; allocate(gen(ng)) 

allocate(prepattern(n,ng))
allocate(block(PD,2))                                     ! target dimensionality
block=0.0                                   

block(1:2,1)= (/-1.0,-1.0/)                               ! target in Environment 1 (/trait1, trait2/)
block(1:2,2)= (/ 1.0,-1.0/)                               ! target in Environment 2 (/trait1, trait2/)

if((training.eq.0).and.(replica.lt.1))then ; return ; end if

!!!!!!!!!!!!!!!!!!!!!!                                    ! open file for seting a population if we are in TEST SET.
  if(training.ne.1)then                                   ! Introducing manually the filename from where the system uploads the population
           !GRN_1234_R12_T1234                            ! Follow this template
    arxaux='GRN_0010_R01_T0001' 
           !123456789012345678
    arxiv(1:3)='GRN' ; arxiv(4:18)='_' ; arxiv(19:36)=arxaux(1:18) 
    arxiv(37:40)='.dat'                                             ! composing filename
    do im=1,40 ; if (arxiv(im:im)==" ") arxiv(im:im)="_" ; end do   ! composing filename   
    WRITE(*,*)'arxiPRE es:       ',arxiv 
    open(9000,file='files/'//arxiv,action='read',iostat=ios) 
    do i=1,11 ;  read(9000,*)  ;  end do                    ! skip first human readable lines about parameters.
  else
    arxiv(19:21)='GRN' ; arxiv(22:36)='_'  
  end if      
!!!!!!!!!!!!!!!!!!!!!!!

do i=1,p                                                  ! for all individuals in the population  
                   
  ind(i)%ncels=n          ; ind(i)%ngs=ng                 ! initial conditions
  ind(i)%fitness=0.0                                      ! initial conditions

  allocate(ind(i)%w(ind(i)%ngs,ind(i)%ngs))               ! it contains all potential genes 
  allocate(ind(i)%ww(ind(i)%ngs,ind(i)%ngs))              ! it contains active interactions
  allocate(ind(i)%g(ind(i)%ncels,ind(i)%ngs))             ! it contaisn all potential genes
  allocate(ind(i)%phen(PD,ind(i)%ncels))                  ! phenotype
  allocate(ind(i)%epigen(ng,n)) 

  allocate(ind (i)%MZ(ind(i)%ngs,PD))  ;  allocate(ind (i)%MZZ(ind(i)%ngs,PD)) 
  allocate(indt(i)%MZ(ind(i)%ngs,PD))  ;  allocate(indt(i)%MZZ(ind(i)%ngs,PD))  

  indt(i)%ncels=n          ; indt(i)%ngs=ng               ! initial conditions
  indt(i)%fitness=0.0      ; ind(i)%epigen=0.0            ! initial conditions
  
  allocate(indt(i)%w(ind(i)%ngs,ind(i)%ngs))              ! it contains all potential genes 
  allocate(indt(i)%ww(ind(i)%ngs,ind(i)%ngs))               ! it contains all potential genes 
  allocate(indt(i)%g(ind(i)%ncels,ind(i)%ngs))            ! it contaisn all potential genes
  allocate(indT(i)%phen(PD,ind(i)%ncels))                 ! phenotype
  
    ind(i)%epigen(1:2,1)=(/-1.0, 1.0/)                    ! ENVIRONMENTAL FACTORS IN ENVIRONMENT 1
    ind(i)%epigen(1:2,2)=(/ 1.0,-1.0/)                    ! ENVIRONMENTAL FACTORS IN ENVIRONMENT 2
    
  do ii=1,ind(i)%ncels                                    ! setting cell topology    
    ind(i)%g(ii,:)=0.0   
    do iii=1,ind(i)%ngs                                   ! setting initial gene concentrations in each cell
       call random_number(x) 
       ind(i)%g(ii,iii)=0.01!*x                           ! small or small noisy initial gene expression     
       prepattern(ii,iii)=ind(i)%g(ii,iii)
    end do
  end do

  ind(i)%w=0.0  ; ind(i)%ww=0 ; gen(:)%deg=0.2 
  ind(i)%MZ=0.0 ; ind(i)%MZZ=0
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CREATING / READING POPULATIONS 
  if(training.eq.1)then                                                 ! creating random population in t=0 in training phase
    do iii=1,ind(i)%ngs
      do iiii=1,ind(i)%ngs 
        call random_number(x)                                           ! Random W matrix in t=0
        ind(i)%w(iii,iiii)=(1.0-(2.0*x))    
        call random_number(x)        
        if(x.le.conWW)then ; jjj=1 ; else ; jjj=0 ; end if
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

 indt=ind                                                 ! just initializing all the matrices for time steps

 do i=1,p
   do ii=1,n                                              ! (ncels)! all individuals must have identical prepatterns
     ind(i)%g(ii,1:ng)=ind(1)%g(ii,1:ng)
     prepattern(ii,1:ng)=ind(1)%g(ii,1:ng)
   end do 
 end do 

end subroutine inicial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module start





