program startodo  ! just the main program boosting the system ...

use start
use development

 integer             :: whois              ! whois=fittest individual
 character(len=26)   :: phenfile           ! name of the file storing fitness and phenotypes over time

 integer, allocatable :: seed(:)           ! setting the seed for random number generator (not accessible)
 integer size                              ! this way all replicates give the same result
 call random_seed(size=size)               !
 allocate(seed(size))                      !
 call random_seed(put=seed)                !

 lapso=10                                  ! every "lapso"
 call inicial                              ! just to get etmax
  
 call arxivpublic                          ! just to enter the inicial module and set the arxiv variable to public
 ret=SYSTEM('pkill gnuplot')               ! ret=SYSTEM('rm dynamic.dat')
 
 open(267,file='GRNstatus.txt',status='unknown',action='write')           ! just to keep the track of the simulations running
 
do supereplica=1,nfiles                                                   ! runs the program once per filename

do replica=1,replicas!4
  
  call inicial                                                            ! it allocates and inicializes everything ...
  
    if(training.eq.1)then
    write(phenfile,"(A8,I1,I1,I1,I1,I1,I1,A2,I2,A2,I2)")'PHEN_TR_',&            ! PHEN_TR for training
    (int(10.0*blocke(1,1))+1)/2,(int(10.0*blocke(1,2))+1)/2,&                       ! creating datafile for phenotypes and fitnesses over time
    (int(10.0*blocke(1,3))+1)/2,(int(10.0*blocke(1,4))+1)/2,&                       ! creating datafile for phenotypes and fitnesses over time
    (int(10.0*blocke(1,5))+1)/2,(int(10.0*blocke(1,6))+1)/2,&
    '_C',thresholdsN(1),'_R',replica                                      ! creating datafile for phenotypes and fitnesses over time
     phenfile(23:26)='.dat'                                               ! creating datafile for phenotypes and fitnesses over time
   else
    write(phenfile,"(A8,I1,I1,I1,I1,I1,I1,A2,I2,A2,I2)")'PHEN_TE_',&            ! PHEN_TE for test
    (int(10.0*blocke(1,1))+1)/2,(int(10.0*blocke(2,1))+1)/2,&                       ! creating datafile for phenotypes and fitnesses over time
    (int(10.0*blocke(1,3))+1)/2,(int(10.0*blocke(2,4))+1)/2,&                       ! creating datafile for phenotypes and fitnesses over time
    (int(10.0*blocke(1,5))+1)/2,(int(10.0*blocke(1,6))+1)/2,&
    '_C',thresholdsN(1),'_R',replica                                      ! creating datafile for phenotypes and fitnesses over time
     phenfile(23:26)='.dat'                                               ! creating datafile for phenotypes and fitnesses over time
   end if

  do im=1,24 ; if (phenfile(im:im)==" ") phenfile(im:im)="0" ; end do   ! composing filename
  open(20067,file=phenfile,status='unknown',action='write')             ! composing filename

fmax=0                                                                  ! records maximum fitness over evol time
fmaxabs=0.0                                                             ! absolute maximum fitness attained over simulation time (initializing)
!open(666,file='debug.dat',status='unknown',action='write')             ! records maximum fitness over evol time

do et=1,etmax                                                           ! evolutionary time (Main Loop)

   do pp=1,p
        ind(pp)%g(:,:)=0.0 ; do i=1,ind(1)%ngs ; ind(pp)%g(1:n,i)=prepattern(1:n,i) ; end do
        ind(pp)%sat=0 ; call dev(pp)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FITNESS CALCULATION
        ind(pp)%fitness=0.0
        do i=1,pd                                                                         ! individual absolute fitness calculation
          do jjj=1,ind(pp)%ncels
	        !ind(pp)%fitness=ind(pp)%fitness+ind(pp)%phen(i,jjj)*blocke(i,jjj)            ! scalar product (P*S)-based fitness (from the paper)
	        ind(pp)%fitness(1)=ind(pp)%fitness(1)+(ind(pp)%phen(i,jjj)-blocke(i,jjj))**2  ! Euclidean-distance based (Used for exponential)
	      end do
	    end do
	    ind(pp)%fitness(1)=sqrt(ind(pp)%fitness(1))                                        ! Euclidean (**2 in Dragui)
	    ind(pp)%fitness(1)=-1.0*ind(pp)%fitness(1)/(2.0*ss)                                ! Exponential fitness function
        ind(pp)%fitness(1)=e**(ind(pp)%fitness(1))                                         ! Exponential fitness function

        !write(*,*)replica,et,'FIT',pp,ind(pp)%fitness(1),'PHEN',ind(pp)%phen(:,1),ind(pp)%phen(:,2)

   end do                                                                                  ! for each individual
   !!!!!!                                                                natural selection for population-based (super-optimized)

    fmaxval=maxval(ind(:)%fitness(1))
    if(fmaxval.gt.fmaxabs)then ; fmaxabs=fmaxval ; end if               ! records absolute maximum fitness attained over simulation time

    indt(:)%fitness(1)=0.0 ; indt(:)%fitness(2)=0.0
    ind(:)%fitness(2)=ind(:)%fitness(1)-abs(minval(ind(:)%fitness(1))) ! negative fitness may

    do i=2,p
      ind(i)%fitness(2)=ind(i)%fitness(2)+ind(i-1)%fitness(2)           ! ordered sumation
    end do

    fmax=maxval(ind(:)%fitness(2))                                      ! scaling fitness (0-1)
    do i=1,p
      ind(i)%fitness(2)=ind(i)%fitness(2)/fmax
    end do
    if(hillclimber.ne.1)then                                            ! Probabilistic selection
      do i=1,p														    ! super-optimized screening
        if(fmax.le.real(p)*delta)then ; whois=i ; goto 642 ; end if     ! all individuals are equal
        call random_number(x)											! non-deterministic selection
        kk=p/2 ; k=p/2
        do j=1,logp
          k=k/2
          if(j.eq.logp-1)then;k=1;end if                                ! allows for p non congruent log2(p)
          if(x.le.ind(kk)%fitness(2))then
            if(j.eq.logp)then;whois=kk;end if
            kk=kk-k
            if(kk.le.0)then ; whois=1 ; goto 642 ; end if
          else if(x.gt.ind(kk)%fitness(2))then
            if(j.eq.logp)then;whois=kk+1;end if
            kk=kk+k
            if(kk.gt.p)then ; whois=p ; goto 642 ; end if
          end if
        end do
        642 indt(i)=ind(whois)                                          ! it fills up next generation
      end do
      ind=indt
    else
      do i=1,p                                                          ! Strict hill-climber
        if(ind(i)%fitness(1).eq.maxval(ind(:)%fitness(1)))then          ! Strict hill-climber
          whois=i ; exit                                                ! Strict hill-climber
        end if                                                          ! Strict hill-climber
      end do
      do i=1,p
        ind(i)=ind(whois)
      end do
    end if
   !!!!!!!!!!!!!!!!!!!!!!!!
   if((mod(et,lapso).eq.0).or.(et.eq.1))then                            ! writting datafile with final matrix before mutation
       
     write(arxaux,"(A4,I1,I1,I1,I1,I1,I1,A2,I2,A2,I2,A2,I2)")'GRN_',(int(10.0*blocke(1,1))+1)/2,(int(10.0*blocke(1,2))+1)/2,&
     (int(10.0*blocke(1,3))+1)/2,(int(10.0*blocke(1,4))+1)/2,&
     (int(10.0*blocke(1,5))+1)/2,(int(10.0*blocke(1,6))+1)/2,'_C',thresholdsN(1),'_R',replica,'_T',int(et/lapso)
     if(arxaux(11:11)==" ") arxaux(11:11)="0"                           ! composing filename threshold
     if(arxaux(15:15)==" ") arxaux(15:15)="0"                           ! composing filename replicate
     do im=19,22 ; if (arxaux(im:im)==" ") arxaux(im:im)="0" ; end do   ! composing filename
     arxifin(1:22)=arxiv(23:44)                                         ! composing filename
     arxifin(23:44)=arxaux(1:22) ; arxifin(45:48)='.dat'                ! composing filename

     do im=1,44 ; if (arxifin(im:im)==" ") arxifin(im:im)="_" ; end do  ! composing filename
     
     open(7000,file=arxifin,status='unknown',action='write',iostat=ios)                       ! creating datafile

     write(7000,*)'TARGETS (E1T1,E1T2,ENT1,ENT2)',blocke(1,1:6)!, blocke(1:2,n)  ! 1
     write(7000,*)'THRESHOLDS(CELL).............',thresholds(1)               ! 2
     write(7000,*)'POPULATON SIZE...............',p                           ! 3
     write(7000,*)'STRENGHT OF SELECTION........',ss                          ! 4
     write(7000,*)'RECOMBINATION; 1=YES; 0=NO   ',reco                        ! 5
     write(7000,*)'TRAINING (1) vs TEST (0) SET ',training                    ! 6
     write(7000,*)'NUMBER /  TOTAL REPLICATES...',replica,replicas            ! 7
     write(7000,*)'CURRENT VS MAXIMUM GENERATION',et,etmax,lapso              ! 8
     write(7000,*)'ENV. FACTORS/ENVIRONMENTS....',EF,n                        ! 9
     write(7000,*)'NUMBER GENES,PHEN. DIMENSIONS',ng,PD                       ! 10
     write(7000,*)'TMAX,SDEV,SS,RECO,CAPPED.....',tmax,sdev,ss,reco,capped    ! 11
     write(7000,*)'CONNECTIVITIES WW / MZZ .....',conWW,conMZZ                ! 12
     write(7000,*)'INPUT VALUES ................',ind(1)%epigen(1,1:n)        ! 13
     !!!!!!!!!!!!!!!!!!!!
     do pp=1,p
       do i=1,ind(1)%ngs
         write(7000,*)ind(pp)%w(i,:)
       end do
       do i=1,ind(1)%ngs
         write(7000,*)ind(pp)%ww(i,:)
       end do
       do j=1,n
            do i=1,PD                                                                 ! PD, ind(i)%ncels)
               write(20067,*)replica,et,pp,j,i,ind(pp)%phen(i,j),ind(pp)%fitness(1)   ! Fitnesses and phenotypes in the general file ...
               call flush(20067)                                                      ! Fitnesses and phenotypes in the general file ...
            end do
       end do
     end do

     do i=1,ind(1)%ngs
       write(7000,*)ind(1)%MZ(i,:)
     end do
     do i=1,ind(1)%ngs
       write(7000,*)ind(1)%MZZ(i,:)
     end do
     close(7000)
     write(267,*)'supereplica',supereplica,'last created file',arxaux,'iteration',et  ! just to keep the track of the simulations running
   end if
   
   !do pp=1,p                    ! Mutation in the generative matrices for each individual
   !  call mutation(pp)          ! independent subroutine (for stability criteria)
   !end do

end do     ! evolutionary time
 !close(666)                      ! temporary file for debugging
 close(20067)                    ! closes file
end do     ! replicates
end do     ! supereplicates

!ret=SYSTEM('rm fort.*')         ! removes spurious stuffs
ret=SYSTEM('mv GRN_* files/')    ! replaces files into a folder
ret=SYSTEM('mv PHEN_* files/')   ! replaces files into a folder

close(676) ; close(267)

end program startodo
