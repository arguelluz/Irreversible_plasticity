module start

implicit none

integer :: i,j,k,ii,jj,kk,iii,iiii,jjj,n,ng,o,p,logp,t,et,tmax,etmax,capped
integer :: training,replicas                   ! If training=1 -> Training set, starting from W=0. Otherwise W from file
integer :: ncels,replica,PD                    ! PD=phenotypic dimensionality (number of traits)
real*4 :: a,aa,b,c,q,u,v,x,y,z,m,deg,fmax,fmaxx,fmaxabs,sdev
integer, parameter  ::  EF=2                   ! EF=Number of environmental factors (inputs)
real*4, allocatable ::  prepattern(:,:),block(:,:)
real*4, parameter   ::  pi=3.1415926535, delta=0.00001
character(len=20)   ::  arxiv,arxifin          ! for writting and reading datafiles

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
p=4                                                       ! number of individuals (must be an EVEN NUMBER !!!)
if(mod(p,2).ne.0)then ; write(*,*)'p must be an EVEN NUMBER' ; end if
logp=1+int(log(real(p))/log(2d0))
tmax=200                                                  ! developmental time
etmax=1000                                                ! evolutionary time          
n=2                                                       ! number of different environments
ng=4                                                      ! initial number of genes
PD=2                                                      ! phenotypic dimensionality (number of traits)
sdev=0.1                                                  ! standard deviation for the mutator algorithm
capped=1                                                  ! If 1, GRN (W-matrix) values are (-1,1); if 0, unconstrained values.
training=1                                                ! If 1 -> Training set, starting from W=0. Otherwise Test set (W from file).
if(training.eq.1)then ; replicas=3 ; end if               ! Replicates for the training set.

if(allocated(ind))then    ; deallocate(ind)   ; deallocate(indt)  ; deallocate(gen)
deallocate(prepattern)    ; deallocate(block) ; end if
allocate (ind(p),indt(p)) ; allocate(gen(ng)) 

allocate(prepattern(n,ng))
allocate(block(PD,2))                                     ! target dimensionality
block=0.0                                   

block(1:2,1)= (/-1.0,-1.0/)                               ! target in Environment 1
block(1:2,2)= (/ 1.0,-1.0/)                               ! target in Environment 2

if((training.eq.0).and.(replicas.lt.1))then ; return ; end if

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
  
  ind(i)%MZ=0.0 ; ind(i)%MZZ=0.0                            ! MATRICES FOR PHENOTYPING

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

  ind(i)%w=0.0 ; gen(:)%deg=0.0 
  if(i.eq.1)then                                           ! homogeneous population in t=0
    if(training.ne.1)then ; read(369,*)arxiv(1:20) ; write(*,*)arxiv 
    open(9000+replica,file='grnfiles/'//arxiv,action='read') ; end if  
    do iii=1,ind(i)%ngs
      if(training.ne.1)then ; read(9000+replica,*) ind(i)%w(iii,1:ind(i)%ngs)
      else
        do iiii=1,ind(i)%ngs
          call random_number(x)
          ind(i)%w(iii,iiii) =x!*0.01        
        end do
      end if       
      gen(iii)%deg=-0.2
      ind(i)%MZ(iii,1:PD) =1.0  ! Provisional, it can be mutated if de-comment in grns.mod.f90
      ind(i)%MZZ(iii,1:PD)=1    ! Provisional, it can be mutated if de-comment in grns.mod.f90
    end do     
    if(training.ne.1)then ; close(9000+replica) ; end if
    ind(i)%ww=0 
  else
    ind(i)%w=ind(1)%w   ; ind(i)%ww=ind(1)%ww  
    ind(i)%mz=ind(1)%MZ ; ind(i)%mzz=ind(1)%MZZ
  end if  
end do

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





