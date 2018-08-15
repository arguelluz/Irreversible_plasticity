module start

implicit none

integer :: i,j,k,ii,jj,kk,iii,iiii,jjj,n,ng,o,p,t,et,tmax,etmax
integer :: ncels,replica,PD                ! PD=phenotypic dimensionality (number of traits)
real*4 :: a,aa,b,c,q,u,v,x,y,z,m,deg,fmax
integer, parameter  ::  replicas=1 , EF=2      ! EF=Number of environmental factors (inputs)
real*4, allocatable ::  prepattern(:,:),block(:,:)
real*4, parameter   ::  pi=3.1415926535 
character(len=31)   ::  arxiv                  ! for writting and reading datafiles

type,public :: inds   
  integer              :: ncels, ngs, sat      ! number of cells, number of genes and per cel and active genes, saturation time
  real*4               :: fitness              ! individual fitness in t=0
  real*4               :: fp(replicas)         ! partial fitness
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
p=2                                                       ! number of individuals
tmax=20                                                   ! developmental time
etmax=10000                                               ! evolutionary time          
n=2                                                       ! number of different environments
ng=9                                                      ! initial number of genes
PD=2                                                      ! phenotypic dimensionality (number of traits)

if(allocated(ind))then    ; deallocate(ind)   ; deallocate(indt)  ; deallocate(gen)
deallocate(prepattern)    ; deallocate(block) ; end if
allocate (ind(p),indt(p)) ; allocate(gen(ng)) 

allocate(prepattern(n,ng))
allocate(block(PD,2))                                     ! target dimensionality
block=0.0                                   

block(1:2,1)= (/-1.0,-1.0/)                               ! target in Environment 1
block(1:2,2)= (/ 1.0,-1.0/)                               ! target in Environment 2

do i=1,p                                                  ! for all individuals in the population  

!  open(9000+i,file='rd0c.dat',action='read')              ! open file for taking GRNs
!  open(9000+i,file='rd1.dat',action='read')              ! open file for taking GRNs
                   
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
  
  ind(i)%mz=0.0 ; ind(i)%mzz=0                            ! MATRICES FOR PHENOTYPING

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
  
  do iii=1,ind(i)%ngs
  !  read(9000+i,*) THIS WAS VAUX.. BUT UNNECESSARY...
    do iiii=1,ind(i)%ngs
      call random_number(x)
      ind(i)%w(iii,iiii) =0.01*x  
      ind(i)%ww(iii,iiii)=0 
    end do     
    gen(iii)%deg=-0.2
  end do  
 
!  close(9000+i) 

end do

 indt=ind                                                 ! just initializing all the atrices for time steps

 do i=1,p
   do ii=1,n !(ncels)! all individuals must have identical prepatterns
     ind(i)%g(ii,1:ng)=ind(1)%g(ii,1:ng)
     prepattern(ii,1:ng)=ind(1)%g(ii,1:ng)
   end do 
 end do 

end subroutine inicial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module start




