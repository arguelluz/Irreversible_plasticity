!!! program BOMB. It introduces point mutations in the GRNs read from file.
!!! input  :: GRN file
!!! Output :: different GRN files, each one having one mutation from the original one. 
!!! Compilation: write in terminal:                        f95 bomb.f90 -o bomb
!!! To run the program afer compiling, write in terminal : ./bomb

program bomb

implicit none

integer :: p,training,replica,replicas,et,ret,nfiles,mcc
integer :: etmax,lapso,EF,n,ng,PD,tmax,reco,capped,mutations,mutmz
integer :: i,j,k,pp,Mi,Mj,Mk,im,ios,xi,xj,xk
real*4  :: thresholds,x,y,z,Mx,My,M0,sdev,ss,conWW,conMZZ,pi,blocke(6)
real*4, allocatable :: epigen(:),w(:,:),wm(:,:),MZ(:,:),mMZ(:,:)
integer,allocatable :: ww(:,:),wwm(:,:),MZZ(:,:),mMZZ(:,:)
character(len=48)arxfin
character(len=53)arxfin2
character(len=29)arxkk,arxikk,arxikk1,arxikk2,arxikk3,arxikk4,arxikk5 ! auxiliar 
character(len=44)arxco                                                ! auxiliar
character(len=1)arxuno                                                ! auxiliar
pi=3.1415926535

mutations=5!105 ! How many mutant GRNs will be created ??
mutmz=0       ! mutmz=0 -> MZ Matrix does not mutate ; mutmz=1 -> MZ Matrix mutates

  open(676,file='GRNfiles.txt',action='read',iostat=ios)                ! automatic reading from file
  nfiles=0
  do while(ios.eq.0)
    read(676,*,iostat=ios)arxfin ; nfiles=nfiles+1
  end do
  nfiles=nfiles-1
  rewind(676)

do mcc=1,nfiles
    read(676,*,iostat=ios)arxfin                                 ! AUTOMATIC          
    !arxfin='GRN_312452_C_4_R_3_T10GRN_224411_C_4_R10_T10.dat'   ! MANUAL      ! Needs to be introduced manually .... (I can make it to be read from the command line if necessary)
    !       !123456789012345678901234567890123456789012345678    ! MANUAL      ! it must occupy 48 characters including extension.        
    !       !         1         2         3         4    .dat    ! MANUAL      ! it must occupy 48 characters including extension.
        
     arxco(1:44)=arxfin(1:44)                                                  ! this keeps all filename except the file extension 
     open(7000,file='files/'//arxfin,status='unknown',action='read',iostat=ios)! reading input file from files/* directory 
     !open(7000,file=arxfin,status='unknown',action='read',iostat=ios)         ! reading input file from the folder where bom.f90 is located 

     read(7000,*)arxkk,arxkk,arxkk,arxkk,arxkk,blocke(1:6)              ! 1
     read(7000,*)arxkk,thresholds                                       ! 2
     read(7000,*)arxkk,arxkk,p                                          ! 3
     read(7000,*)arxkk,arxkk,arxkk,ss                                   ! 4
     read(7000,*)arxkk,arxkk,arxkk,reco                                 ! 5
     read(7000,*)arxkk,arxkk,arxkk,arxkk,arxkk,arxkk,training           ! 6
     read(7000,12)arxkk,arxkk,arxkk,arxkk,replica,replicas              ! 7
     read(7000,*)arxkk,arxkk,arxkk,arxkk,et,etmax,lapso                 ! 8
     read(7000,*)arxkk,arxkk,EF,n                                       ! 9
     EF=1 ; n=6 ! this line can not be read because of nonstandard ASCI symbols
     read(7000,*)arxkk,arxkk,arxkk,arxkk,ng,PD                          ! 10
     read(7000,*)arxkk,arxkk,arxkk,arxkk,arxkk,tmax,sdev,ss,reco,capped ! 11     
     read(7000,13)arxikk,ARXUNO,arxuno,arxuno,arxuno,conWW,conMZZ       ! 12 
     if((conww.gt.1.0) .or.(conww.lt.0.0))then  ; conww=1.0  ; write(*,*)'reading error, conWW  set to 1.0' ; end if
     if((conMZZ.gt.1.0).or.(conMZZ.lt.0.0))then ; conMZZ=1.0 ; write(*,*)'reading error, conMZZ set to 0.5' ; end if
     if(allocated(epigen))then ; deallocate(epigen) ; end if
     allocate(epigen(n))
     read(7000,*)arxkk,arxkk,arxkk,epigen(1:n)                          ! 13
          
     12 FORMAT(A6,A1,A5,A18,I20,I20)                                    ! Some required formatting for reading issues     
     13 FORMAT(A29,4A1,F11.1,F19.1)                                     ! Some required formatting for reading issues     
     !!!!!!!!!!!!!!!!!!!!
     
     if(allocated(w))then 
        deallocate(w) ;deallocate(ww) ;deallocate(wm) ;deallocate(wwm)
        deallocate(MZ);deallocate(MZZ);deallocate(mMZ);deallocate(mMZZ)
     end if
     allocate(w(ng,ng),ww(ng,ng),wm(ng,ng),wwm(ng,ng))     ! ALLOCATING WILDTYPE AND MUTANT (M) MATRICES 
     allocate(MZ(ng,PD),MZZ(ng,PD),mMZ(ng,PD),mMZZ(ng,PD)) ! ALLOCATING WILDTYPE AND MUTANT (M) MATRICES 
     
     do pp=1,p
       do i=1,ng
         read(7000,*)w(i,:)
       end do
       do i=1,ng
         read(7000,*)ww(i,:)
       end do      
     end do

     do i=1,ng
       read(7000,*)MZ(i,:)
     end do
     do i=1,ng
       read(7000,*)MZZ(i,:)
     end do
     close(7000)     

if(mutations.gt.999)then ; write(*,*)'WARNING! can not create more then 999 filenames' ; read(*,*) ; end if
do j=1,mutations                                                   ! Here the new mutant GRNs are created.

     wm=w ; wwm=ww ; mMZ=MZ ; mMZZ=MZZ                             ! before mutation, matrices are identical
     
     ! MUTATING THE W MATRIX
     call random_number(Mx)  ; call random_number(My)
     Mi=int(Mx*real(ng)+1) ;  Mj=int(My*real(ng)+1)                ! which element of W will mutate
     call random_number(Mx)  ; call random_number(My)              ! random new value for the mutation
     M0=sdev*sqrt(-2*log(Mx))*cos(2*pi*My)                         ! Box-Muller algotithm. Normal distribtion N(0,sdev)
     wm(Mi,Mj)= w(Mi,Mj)+M0                                        ! adding the new random value to the previous one
     if(capped.eq.1) then                                          ! If we are using capped weights
       if(wm(Mi,Mj).gt.1.0) then                                   ! and if the mutation makes weights greater than 1
         wm(Mi,Mj)=1.0                                             ! Then set them to 1
       else if (wm(Mi,Mj).lt.-1.0) then                            ! if the mutation makes weights lesser than -1
         wm(Mi,Mj)=-1.0                                            ! Then set them to -1
       end if
     end if
     !write(*,*)'In the ',j,'th GRN, the element W(',Mi,Mj,') is changed from ',w(Mi,Mj),' to ',wm(Mi,Mj)

     ! MUTATING THE MZ MATRICES
     if(mutmz.eq.1)then                                            ! If mutmz=0, MZ and MZZ Matrices do not mutate
       call random_number(Mx)  ; call random_number(My)
       Mi=int(Mx*real(ng)+1) ;  Mj=int(My*real(PD)+1)              
       call random_number(M0) ; M0=1.0-2*M0 ; M0=M0*0.2
       mMZ(Mi,Mj)=MZ(Mi,Mj)+M0
       !write(*,*)'In the ',j,'th GRN, the element MZ(',Mi,Mj,') is changed from ',MZ(Mi,Mj),' to ',MZ(Mi,Mj)

       call random_number(Mx)  ; call random_number(My)
       Mi=int(Mx*real(ng)+1) ;  Mj=int(My*real(PD)+1)              
       if(MZZ(Mi,Mj).eq.0)then ; mMZZ(Mi,Mj)=1 ; else 
       mMZZ(Mi,Mj)=1 ; end if                                           ! topological change
       !write(*,*)'In the ',j,'th GRN, the element MZZ(',Mi,Mj,') is changed from ',MZZ(Mi,Mj),' to ',MZZ(Mi,Mj)
     end if

     !! WRITTING THE OUTPUT FILES  
     write(arxfin2,"(A44,A1,I4,A4)")arxco,'_',j,'.dat'                  ! composing filename 
     do im=1,53 ; if (arxfin2(im:im)==" ") arxfin2(im:im)="_" ; end do  ! composing filename
     open(7001,file=arxfin2,status='unknown',action='write',iostat=ios)            ! creating datafile in the folder where bomb.f90 is located
     
     !!!!!!!!!!!!!!   writting the new mutant networks     
     write(7001,*)'TARGETS (E1T1,E1T2,ENT1,ENT2)',blocke(1:6)                 ! 1
     write(7001,*)'THRESHOLDS(CELL).............',thresholds                  ! 2
     write(7001,*)'POPULATON SIZE...............',p                           ! 3
     write(7001,*)'STRENGHT OF SELECTION........',ss                          ! 4
     write(7001,*)'RECOMBINATION; 1=YES; 0=NO   ',reco                        ! 5
     write(7001,*)'TRAINING (1) vs TEST (0) SET ',training                    ! 6
     write(7001,*)'NUMBER /  TOTAL REPLICATES...',replica,replicas            ! 7
     write(7001,*)'CURRENT VS MAXIMUM GENERATION',et,etmax,lapso              ! 8
     write(7001,*)'ENV. FACTORS/ENVIRONMENTS....',EF,n                        ! 9
     write(7001,*)'NUMBER GENES,PHEN. DIMENSIONS',ng,PD                       ! 10
     write(7001,*)'TMAX,SDEV,SS,RECO,CAPPED.....',tmax,sdev,ss,reco,capped    ! 11
     write(7001,*)'CONNECTIVITIES WW / MZZ .....',conWW,conMZZ                ! 12
     write(7001,*)'INPUT VALUES ................',epigen(1:n)                 ! 13     
          
     do pp=1,p
       do i=1,ng
         write(7001,*)wm(i,:)
       end do
       do i=1,ng
         write(7001,*)wwm(i,:)
       end do      
     end do

     do i=1,ng
       write(7001,*)mMZ(i,:)
     end do
     do i=1,ng
       write(7001,*)mMZZ(i,:)
     end do
     close(7001)     
     !!!!!!!!!!!!!!

end do ! from mutations    
end do ! from mcc (number of input files)

      close(676)
      
      !ret=SYSTEM('mv *T??_*.dat files')      ! put the new files in the /files folder

end program
