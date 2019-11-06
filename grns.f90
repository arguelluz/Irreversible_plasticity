program startodo  ! just the main program boosting the system ...

use start
use development

 integer             :: ret,pp,ios
 integer             :: lapso,whois        ! lapso=intervals for data recording, whois=fittest individual
 real*4,allocatable  :: averagew(:,:)

 integer, allocatable :: seed(:)           ! setting the seed for random number generator
 integer size                              ! this way all replicates give the same result
 call random_seed(size=size)               !
 allocate(seed(size))                      !
 call random_seed(put=seed)                !

 lapso=10
 call inicial                              ! just to get etmax

 if(training.eq.0)then                     ! Set replicate as number of datafiles in Test set.
   pp=0
   open(369,file='GRNames.dat',action='read',iostat=ios)
   do while (ios.eq.0)
     read(369,*,iostat=ios)arxiv ; pp=pp+1
   end do
   replicas=pp-1 ; rewind(369)
 end if

 allocate(averagew(etmax/lapso+1,replicas))
 averagew=0

call arxivpublic                           ! just to enter the inicial module and set the arxiv variable to public
ret=SYSTEM('pkill gnuplot')                ! ret=SYSTEM('rm dynamic.dat')

do replica=1,replicas!4

call inicial                               ! it allocates and inicializes everything ...

open(20067,file='Fitnessfile.dat',status='unknown',action='write')     ! okflush
open(20068,file='Phenotyfile.dat',status='unknown',action='write')     ! okflush

fmax=0  ! records maximum fitness over evol time

do et=1,etmax    ! evolutionary time

   do pp=1,p

   !do replica=1,replicas!4

         ind(pp)%g(:,:)=0.0 ; do i=1,ind(1)%ngs ; ind(pp)%g(1:n,i)=prepattern(1:n,i) ; end do
         ind(pp)%sat=0 ; call dev(pp)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FITNESS CALCULATION
         ind(pp)%fitness=0.0

        do i=1,pd                                                               ! individual absolute fitness calculation
          write(*,*)'pheno',ind(i)%phen
          do jjj=1,ind(pp)%ncels
	        ind(pp)%fitness=ind(pp)%fitness+ind(pp)%phen(i,jjj)*block(i,jjj)    ! scalar product (P*S)-based fitness (from the paper)
	      end do
	    end do

   end do  ! for each individual

   !!!!!! natural selection for Hill-climber
   !if(ind(1)%fitness.ge.ind(2)%fitness)then ! reproduction (the fittest individual is copied)
   !  ind(2)=ind(1) ; whois=1
   !  do jjj=1,PD ; write(20068,*)et,ind(1)%phen(jjj,:) ; end do
   !else
   !  ind(1)=ind(2) ; whois=2
   !  do jjj=1,PD ; write(20068,*)et,ind(1)%phen(jjj,:) ; end do
   !end if
   !!!!!! natural selection for population-based (non-optimized)
   write(*,*)et,'fitnesses',ind(:)%fitness

   do i=2,p                                  ! scaled sumation of fitnesses
     ind(i)%fitness=ind(i)%fitness+ind(i-1)%fitness
   end do
   fmax=maxval(ind(:)%fitness)                  ! first relative fitness (0-1)
   do i=1,p
     ind(i)%fitness=ind(i)%fitness/fmax
   end do


   write(*,*)et,'fitnESFIN',ind(:)%fitness

   !!!!!!!!!!!!!!!!!!!!!!!!

   if(fmax.gt.fmaxx)then   ! data racording (re-check)
     fmaxx=fmax
   end if
   if(mod(et,lapso).eq.0)then
     averagew(et/lapso,replica)=fmaxx/real(n)
   end if
   write(20067,*)et,fmax ; call flush(20067)
   !!!!!!!!!!!!!!!!!!!!!!!

   if(et.eq.etmax)then                       ! writting datafile with final matrix before mutation
     write(arxifin,"(A3,I1,A16)")'GRN',replica,'____________.dat'
     open(7000+replica,file=arxifin,action='write')
     do i=1,ind(1)%ngs
       write(7000+replica,*)ind(whois)%w(i,:)
     end do
     close(7000+replica)
   end if

   do pp=1,1!p      pq solo hay que actualizar el primeroooo ! conserved prepattern, genotype resetting
    ind(pp)%g(:,:)=0.0
    do i=1,ind(1)%ngs
      ind(pp)%g(1:n,i)=prepattern(1:n,i)
    end do
   end do

     !!!! MUTATION IN THE FOUR MATRICES !!!!!!!

     call random_number(x)  ; call random_number(y)
     i=int(x*real(ng)+1) ;  j=int(y*real(ng)+1)               ! which element of W will mutate
     747 call random_number(x)  ; call random_number(y)       ! random new value for the mutation
     z=sdev*sqrt(-2*log(x))*cos(2*pi*y)                       ! Box-Muller algotithm. Normal distribtion N(0,sdev)
     ind(1)%w(i,j)= ind(1)%w(i,j)+z                           ! adding the new random value to the previous one
     if((capped.eq.1).and.((ind(1)%w(i,j).gt.1.0).or.(ind(1)%w(i,j).lt.-1.0)))then ; goto 747 ; end if   ! capped ??

     call random_number(x)  ; call random_number(y)
     i=int(x*real(ng)+1) ;  j=int(y*real(ng)+1)               ! permito mutar un gen mas ...
     if(ind(1)%ww(i,j).eq.0)then ; ind(1)%ww(i,j)=1 ; else ;  ind(1)%ww(i,j)=1 ; end if ! topological change

   !  call random_number(x)  ; call random_number(y)
   !  i=int(x*real(ng)+1) ;  j=int(y*real(PD)+1)              ! permito mutar un gen mas ...
   !  call random_number(z) ; z=1.0-2*z ; z=z*0.2
   !  ind(1)%MZ(i,j)= ind(1)%MZ(i,j)+z

   !  call random_number(x)  ; call random_number(y)
   !  i=int(x*real(ng)+1) ;  j=int(y*real(PD)+1)              ! permito mutar un gen mas ...
   !  if(ind(1)%MZZ(i,j).eq.0)then ; ind(1)%MZZ(i,j)=1 ; else ;  ind(1)%MZZ(i,j)=1 ; end if ! topological change

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end do     ! evolutionary time
end do     ! replicates

close(20067) ; close(20068)


end program startodo
