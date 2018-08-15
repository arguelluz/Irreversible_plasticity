program startodo  ! just the main program boosting the system ...

use start
use development

 integer             :: ret,pp,lapso                 
 real*4              :: maxx,minn,time,xxif,dxif
 real*4,allocatable  :: averagew(:,:)
 
 integer, allocatable :: seed(:)           ! setting the seed for random number generator
 integer size                              ! this way all replicates give the same result
 call random_seed(size=size)               ! 
 allocate(seed(size))                      !
 call random_seed(put=seed)                ! 

 lapso=10
 call inicial                              ! just to get etmax
 allocate(averagew(etmax/lapso+1,replicas))
 averagew=0 

call arxivpublic                           ! just to enter the inicial module and set the arxiv variable to public
ret=SYSTEM('pkill gnuplot') !; ret=SYSTEM('rm dynamic.dat')

do replica=1,replicas!4

call inicial                                              ! it allocates and inicializes everything ...

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
          do jjj=1,ind(pp)%ncels    
	        ind(pp)%fitness=ind(pp)%fitness+ind(pp)%phen(i,jjj)*block(i,jjj) ! scalar product (P*S)-based fitness (from the paper)       
	      end do
	    end do 
         
   end do  ! for each individual    

   if(maxval(ind(:)%fitness).gt.fmax)then       !! data racording
     fmax=maxval(ind(:)%fitness) 
   end if
   if(mod(et,lapso).eq.0)then
     averagew(et/lapso,replica)=fmax/real(n)
   end if
   write(20067,*)et,maxval(ind(:)%fitness) ; call flush(20067)  
  
   !!!!!! natural selection
   if(ind(1)%fitness.ge.ind(2)%fitness)then ! reproduction (the fittest individual is copied) 
     ind(2)=ind(1) 
     do jjj=1,PD ; write(20068,*)et,ind(1)%phen(jjj,:) ; end do
   else 
     ind(1)=ind(2) 
     do jjj=1,PD ; write(20068,*)et,ind(1)%phen(jjj,:) ; end do
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
   !  i=int(x*real(ng)+1) ;  j=int(y*real(PD)+1)           ! permito mutar un gen mas ...
   !  if(ind(1)%MZZ(i,j).eq.0)then ; ind(1)%MZZ(i,j)=1 ; else ;  ind(1)%MZZ(i,j)=1 ; end if ! topological change      
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
end do     ! evolutionary time 
end do     ! replicates

close(20067) ; close(20068)


end program startodo










