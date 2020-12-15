program startodo  ! just the main program boosting the system ...

use start
use development

 integer             :: whois              ! whois=fittest individual
 integer             :: npoint             ! which of the three point is to be recorded (just for the datafile)
 integer             :: record,recorda     ! it determines if data will be recorded (different criteria for training and test)
 integer             :: wlen,wpnt,wlap     ! (moving) window leght, number of points and lapse.
 real*4,allocatable  :: mwin(:)            ! moving window. Filled up with fitness values. Criterium for sampling points
 real*4              :: avg0,mmm0,var0     ! statistics for the moving window. Intial (in t<wlen) Average, maxmum and variance
 real*4              :: avg1,mmm1,var1     ! statistics for the moving window. Average, maximum and variance.
 real*4              :: thresvar,avgc      ! threshold in the variance        (to record the moving window)
 real*4              :: thresavg           ! threshold in the average fitness (to record the moving window)
 character(len=26)   :: phenfile           ! name of the file storing fitness and phenotypes over time ! Short, auxiliar file
 character(len=48)   :: phenfileL          ! name of the file storing fitness and phenotypes over time ! Long to know intial conditions in tests.

 integer, allocatable :: seed(:)           ! setting the seed for random number generator (not accessible)
 integer size                              ! this way all replicates give the same result
 call random_seed(size=size)               !
 allocate(seed(size))                      !
 call random_seed(put=seed)                !
 call inicial                              ! just to get etmax

 call arxivpublic                          ! just to enter the inicial module and set the arxiv variable to public
 ret=SYSTEM('pkill gnuplot')               ! ret=SYSTEM('rm dynamic.dat')

 wlen=etmax/100                            ! (moving) window leght=tmax/100. This should be enough.
 wpnt=50                                   ! (moving) window points. They must be few (calculating variance is time-consuming)
 if(wlen.lt.wpnt)then ; wlen=wpnt ; end if ! Number of points must be smaller than window lenght
 wlap=wlen/wpnt                            ! (moving) window lapse. That is the lenght of the intervals to recods fitness.
 thresvar=1d-4                             ! threshold in the variance        (to record the moving window)
 thresavg=1d-5                             ! threshold in the average fitness (to record the moving window)
 if(allocated(mwin))then                   ! allocating the moving window vector.
 deallocate(mwin) ; end if                 ! allocating the moving window vector.
 allocate(mwin(wpnt))                      ! allocating the moving window vector.

do supereplica=1,nfiles                                                   ! runs the program once per filename

do replica=1,replicas!4

  call inicial                                                            ! it allocates and inicializes everything ...
   if(training.eq.1)then
    write(phenfile,"(A8,I1,I1,I1,I1,I1,I1,A2,I2,A2,I2)")'PHEN_TR_',&      ! PHEN_TR for training
    (int(10.0*blocke(1,1))+1)/2,(int(10.0*blocke(1,2))+1)/2,&             ! creating datafile for phenotypes and fitnesses over time
    (int(10.0*blocke(1,3))+1)/2,(int(10.0*blocke(1,4))+1)/2,&             ! creating datafile for phenotypes and fitnesses over time
    (int(10.0*blocke(1,5))+1)/2,(int(10.0*blocke(1,6))+1)/2,&
    '_C',thresholdsN(1),'_R',replica                                      ! creating datafile for phenotypes and fitnesses over time
     phenfile(23:26)='.dat'                                               ! creating datafile for phenotypes and fitnesses over time
     do im=1,24 ; if (phenfile(im:im)==" ") phenfile(im:im)="0" ; end do  ! composing filename
     phenfileL(23:48)=phenfile(1:26) ; do im=1,22 ; phenfileL(im:im)='_' ; end do
     phenfileL(1:3)='PHE'
   else
    write(phenfile,"(A8,I1,I1,I1,I1,I1,I1,A2,I2,A2,I2)")'PHEN_TE_',&      ! PHEN_TE for test
    (int(10.0*blocke(1,1))+1)/2,(int(10.0*blocke(1,2))+1)/2,&             ! creating datafile for phenotypes and fitnesses over time
    (int(10.0*blocke(1,3))+1)/2,(int(10.0*blocke(1,4))+1)/2,&             ! creating datafile for phenotypes and fitnesses over time
    (int(10.0*blocke(1,5))+1)/2,(int(10.0*blocke(1,6))+1)/2,&             ! creating datafile for phenotypes and fitnesses over time.
    '_C',thresholdsN(1),'_R',replica                                      ! creating datafile for phenotypes and fitnesses over time
     phenfile(23:26)='.dat'                                               ! creating datafile for phenotypes and fitnesses over time
     do im=1,24 ; if (phenfile(im:im)==" ") phenfile(im:im)="0" ; end do  ! composing filename
     phenfileL(23:48)=phenfile(1:26) ; phenfileL(1:22)=arxiv(23:44)
     phenfileL(1:3)='PHE'                                                 ! EASY TO GREP
   end if
   open(20067,file=phenfileL,status='unknown',action='write')             ! composing filename

fmax=0                                                                  ! records maximum fitness over evol time
fmaxabs=0.0                                                             ! absolute maximum fitness attained over simulation time (initializing)
recorda=0                                                               ! A counter to ensure Point2 has been taken.
mwin=0.0                                                                ! initializing the moving window vector.
 avg0=0.0 ; mmm0=0.0 ; var0=0.0                                         ! initializing ths moving window statistics.
 avg1=0.0 ; mmm1=0.0 ; var1=0.0                                         ! initializing ths moving window statistics.
do et=1,etmax                                                           ! evolutionary time (Main Loop)
  !if(mod(et,1000).eq.0)then ; write(*,*)'rep,et',replica,et ; end if   ! just a tracker
   do pp=1,p
        ind(pp)%g(:,:)=0.0 ; do i=1,ind(1)%ngs ; ind(pp)%g(1:n,i)=prepattern(1:n,i) ; end do
        call dev(pp)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FITNESS CALCULATION
        ind(pp)%fitness=0.0
        do i=1,pd                                                                         ! individual absolute fitness calculation
          do jjj=1,ind(pp)%ncels
	        ind(pp)%fitness(1)=ind(pp)%fitness(1)+(ind(pp)%phen(i,jjj)-blocke(i,jjj))**2  ! Euclidean-distance based (Used for exponential)
	      end do
	    end do
	    ind(pp)%fitness(1)=sqrt(ind(pp)%fitness(1))                                        ! Euclidean (**2 in Dragui)
	    ind(pp)%fitness(1)=-1.0*ind(pp)%fitness(1)/(2.0*ss)                                ! Exponential fitness function
        ind(pp)%fitness(1)=e**(ind(pp)%fitness(1))                                         ! Exponential fitness function

   end do                                                                                  ! for each individual
   !!!!!!                                                                natural selection for population-based (super-optimized)
    fmaxval=maxval(ind(:)%fitness(1))
    if(fmaxval.gt.fmaxabs)then ; fmaxabs=fmaxval ; end if               ! records absolute maximum fitness attained over simulation time

    indt(:)%fitness(1)=0.0 ; indt(:)%fitness(2)=0.0
    ind(:)%fitness(2)=ind(:)%fitness(1)-abs(minval(ind(:)%fitness(1))) ! negative fitness may

    !!!!!!!!!                                !  (all the following piece of code)
    if(training.eq.0)then                    ! thois piece of code only runs in test simulations
     if(et.lt.wlen)then                      ! it fills up the moving window vector at the beginning
       if(et.eq.1)then                       ! first generation is taken into account
         mwin(1)=fmaxabs!val                 ! maximun historical fitnes is used
       else if(mod(et,wlap).eq.0)then        ! initial moving window
         mwin(1+(et/wlap))=fmaxval           ! initial moving window
       end if
       if(et.eq.wlen-wlap)then               ! it calculates the values of the "fist window" that will be used for later comparisons
         call stati(mwin,wpnt,avg0,mmm0,var0)! it calculates the values of the "fist window" that will be used for later comparisons
       end if
       avg1=avg0 ; mmm1=mmm0 ; var1=var0     ! otherwise stats in the second window would be zero
     else                                    ! it displaces the moving window
       avg0=avg1 ; mmm0=mmm1 ; var0=var1     ! each window would be compared with the previous one (ow with the 1t one this line is commented)
       if(mod(et,wlap).eq.0)then
         do k=1,wpnt-1                       ! this displaces one cell the whole array of values
           mwin(k)=mwin(k+1)
         end do
         mwin(wpnt)=fmaxabs!val              ! this introduces the last value in the moving window
         call stati(mwin,wpnt,avg1,mmm1,var1)! it calculates the values of to be compared with those of the 1st window
       end if
     end if
    end if                                   ! of training=1
    !!!!!!!!!                                !  end

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
   record=0
   if(training.eq.1)then                                                    !
     if((et.eq.lapso(klog)).or.(et.eq.1))then                               ! writting datafile with final matrix before mutation
       record=1                                                             ! In training simulations, it writes the files as always.
     end if                                                                 !
   else                                                                     !
     avgc=(avg1/avg0)-1.0                                                   !  Fitness increases must be significative
     if(et.eq.wlen+1) then                                                  ! criterium for 1_st point.
       npoint=1                                                             ! criterium for 1_st point.
       record=1                      !!!!!!!! COMMENT THIS LINE IF YOU DON'T WANT THE 1ST POINT TO BE RECORDED   !!!!!!!
     goto 997   ;   end if                                                  ! criterium for 1_st point.
     if((mod(et,wlap).eq.0).and.(avgc.lt.thresavg).and.&
     (var1.lt.thresvar).and.(npoint.eq.1))then                              ! criterium for 2_nd point. SIMPLER. THIS WORKS !
     ; recorda=1 ; record=1; ; npoint=2 ; goto 997; end if                  ! criterium for 2_nd point.
     if((et.eq.etmax-1).and.(recorda.eq.0).and.(npoint.eq.1))then           ! If point 2 has not been taken and the simulation is finishing.
       record=1 ; npoint=2 ; goto 997  ; end if                             ! If point 2 has not been taken and the simulation is finishing.
     if((et.eq.etmax).and.(npoint.eq.2)) then     ; npoint=3      ; record=1; goto 997; end if ! criterium for 3_rd point.
   end if

 997 if(record.eq.1)then
     if(training.eq.1)then                                                                                                      !
       write(arxaux,"(A4,I1,I1,I1,I1,I1,I1,A2,I2,A2,I2,A2,I2)")'GRN_',(int(10.0*blocke(1,1))+1)/2,(int(10.0*blocke(1,2))+1)/2,&
       (int(10.0*blocke(1,3))+1)/2,(int(10.0*blocke(1,4))+1)/2,&
       (int(10.0*blocke(1,5))+1)/2,(int(10.0*blocke(1,6))+1)/2,'_C',thresholdsN(1),'_R',replica,'_T',klog!int(et/lapso)
     else                                                                                                                       !
      if(npoint.eq.2)then ; write(*,*)et,'statsP2=',avg1,avg0,var1,var0; end if
      write(arxaux,"(A4,I1,I1,I1,I1,I1,I1,A2,I2,A2,I2,A2,I2)")'GRN_',(int(10.0*blocke(1,1))+1)/2,(int(10.0*blocke(1,2))+1)/2,&  !
       (int(10.0*blocke(1,3))+1)/2,(int(10.0*blocke(1,4))+1)/2,&                                                                !
       (int(10.0*blocke(1,5))+1)/2,(int(10.0*blocke(1,6))+1)/2,'_C',thresholdsN(1),'_R',replica,'_T',npoint                     !
     end if                                                                                                                     !
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
     write(7000,*)'CURRENT VS MAXIMUM GENERATION',et,etmax,klog               ! 8
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
       if(training.eq.1)then                                                        ! This loop must be below. Otherwise only 3 phenotypic points.
         do j=1,n                                                                   ! This loop must be below. Otherwise only 3 phenotypic points.
           do i=1,PD                                                                ! PD, ind(i)%ncels)
             write(20067,*)replica,et,pp,j,i,ind(pp)%phen(i,j),ind(pp)%fitness(1)   ! Fitnesses and phenotypes in the general file ...
             call flush(20067)                                                      ! Fitnesses and phenotypes in the general file ...
           end do
         end do
       end if                                                                        ! This loop must be below. Otherwise only 3 phenotypic points.
     end do

     do i=1,ind(1)%ngs
       write(7000,*)ind(1)%MZ(i,:)
     end do
     do i=1,ind(1)%ngs
       write(7000,*)ind(1)%MZZ(i,:)
     end do
     close(7000)

     if(training.eq.1)then                                                   ! Otherwise only 3 phenotypic points.
       klog=klog+1                 ! specific counter for logaritmic timepoints
     end if                                                                  ! This loop must be below. Otherwise only 3 phenotypic points.
   end if                          ! this "if" closing is to take timepoints ...

   if(training.eq.0)then                                                       ! This loop must be below. Otherwise only 3 phenotypic points.
    if((et.eq.lapso(klog)).or.(et.eq.1))then                                   ! This must be here now. Otherwise only 3 phenotypic points.
   ! if((mod(et,5000).eq.0).or.(et.eq.1))then                                  ! Even (linear) sampling rather than Logtime sampling.
     do pp=1,p                                                                 !
      do j=1,n                                                                 !
        do i=1,PD                                                              ! PD, ind(i)%ncels)
          write(20067,*)replica,et,pp,j,i,ind(pp)%phen(i,j),ind(pp)%fitness(1) ! Fitnesses and phenotypes in the general file ...
          call flush(20067)                                                    ! Fitnesses and phenotypes in the general file ...
        end do                                                                 !
      end do                                                                   !
     end do                                                                    !
     klog=klog+1
     end if
   end if                                                          !

end do     ! evolutionary time

 close(20067)                    ! closes file
end do     ! replicates
end do     ! supereplicates


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !
   contains
   subroutine stati(mwina,snn,sxx,syy,szz)
   real*4 :: sxx,syy,szz
   integer:: snn,sni
   real*4 :: mwina(snn)
     sxx=0.0 ; syy=0.0 ; szz=0.0
     sxx=sum(mwina)/real(snn) ! average
     syy=maxval(mwina)        ! maximum
     do sni=1,snn             ! variance
       szz=szz+(sxx-mwina(sni))**2
     end do
     szz=sqrt(szz/real(snn))
   end subroutine
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !

end program startodo
