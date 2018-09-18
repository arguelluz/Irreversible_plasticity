program startodo  ! just the main program boosting the system ...

use start
use development

 integer             :: whois              ! whois=fittest individual            
 real*4,allocatable  :: averagew(:,:)      ! a matrix for average data
 
 integer, allocatable :: seed(:)           ! setting the seed for random number generator (not accessible)
 integer size                              ! this way all replicates give the same result
 call random_seed(size=size)               ! 
 allocate(seed(size))                      !
 call random_seed(put=seed)                ! 

 lapso=10                                  ! every "lapso"
 call inicial                              ! just to get etmax
 
 allocate(averagew(etmax/lapso+1,replicas))
 averagew=0 

call arxivpublic                           ! just to enter the inicial module and set the arxiv variable to public
ret=SYSTEM('pkill gnuplot')                ! ret=SYSTEM('rm dynamic.dat')

do replica=1,replicas!4

call inicial                               ! it allocates and inicializes everything ...

open(20067,file='Fitnessfile.dat',status='unknown',action='write')     ! okflush  !! WITHIN REPLICATEE !!*****
open(20068,file='Phenotyfile.dat',status='unknown',action='write')     ! okflush !*********************

fmax=0  ! records maximum fitness over evol time

do et=1,etmax    ! evolutionary time
   
   do pp=1,p           
          
         ind(pp)%g(:,:)=0.0 ; do i=1,ind(1)%ngs ; ind(pp)%g(1:n,i)=prepattern(1:n,i) ; end do   
         ind(pp)%sat=0 ; call dev(pp)     
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FITNESS CALCULATION        
        ind(pp)%fitness=0.0     
        do i=1,pd                                                                  ! individual absolute fitness calculation            
          do jjj=1,ind(pp)%ncels                
	        !ind(pp)%fitness=ind(pp)%fitness+ind(pp)%phen(i,jjj)*block(i,jjj)      ! scalar product (P*S)-based fitness (from the paper)       
	        ind(pp)%fitness=ind(pp)%fitness+(ind(pp)%phen(i,jjj)-block(i,jjj))**2  ! Euclidean-distance based (Used for exponential)
	      end do
	    end do 
	    ind(pp)%fitness=sqrt(ind(pp)%fitness)/(2.0*ss)                             ! Exponential fitness function
        ind(pp)%fitness=e**(ind(pp)%fitness)                                       ! Exponential fitness function
   end do  ! for each individual    
    
   !!!!!!                                                                natural selection for population-based (super-optimized)
    fmaxabs=maxval(ind(:)%fitness)
    indt(:)%fitness=0.0
    ind(:)%fitness=ind(:)%fitness-minval(ind(:)%fitness)+delta          ! negative fitness may 
    do i=2,p
      ind(i)%fitness=ind(i)%fitness+ind(i-1)%fitness                    ! ordered sumation
    end do
    fmax=maxval(ind(:)%fitness)                                         ! scaling fitness (0-1)
    do i=1,p 
      ind(i)%fitness=ind(i)%fitness/fmax 
    end do
    do i=1,p														    ! super-optimized screening
      if(fmax.le.real(p)*delta)then ; whois=i ; goto 642 ; end if       ! all individuals are equal
      call random_number(x)												! non-deterministic selection
      kk=p/2 ; k=p/2
      do j=1,logp
        k=k/2   
        if(j.eq.logp-1)then;k=1;end if                                  ! allows for p non congruent log2(p)
        if(x.le.ind(kk)%fitness)then                
          if(j.eq.logp)then;whois=kk;end if
          kk=kk-k    
          if(kk.le.0)then ; whois=1 ; goto 642 ; end if  
        else if(x.gt.ind(kk)%fitness)then
          if(j.eq.logp)then;whois=kk+1;end if
          kk=kk+k     
          if(kk.gt.p)then ; whois=p ; goto 642 ; end if
        end if         
      end do        
      642 indt(i)=ind(whois)                                            ! it fills up next generation 
    end do      
    ind=indt
   !!!!!!!!!!!!!!!!!!!!!!!!

   if(fmaxabs.gt.fmaxx)then                                             ! data racording (re-check)
     fmaxx=fmaxabs 
   end if
   if(mod(et,lapso).eq.0)then
     averagew(et/lapso,replica)=fmaxx/real(n)
   end if
   write(20067,*)et,fmaxabs ; call flush(20067) 
   
   !!!!!!!!!!!!!!!!!!!!!!!   
   if(mod(et,lapso).eq.0)then                                           ! writting datafile with final matrix before mutation      
     write(arxaux,"(A4,I1,I1,I1,I1,A2,I2,A2,I4)")'GRN_',(int(block(1,1))+1)/2,(int(block(2,1))+1)/2,&
     (int(block(1,2))+1)/2,(int(block(2,2))+1)/2,'_R',replica,'_T',int(et/lapso)            
     if(arxaux(11:11)==" ") arxaux(11:11)="0"                           ! composing filename
     do im=15,18 ; if (arxaux(im:im)==" ") arxaux(im:im)="0" ; end do   ! composing filename
     arxifin(1:18)=arxiv(19:36)                                         ! composing filename
     arxifin(19:36)=arxaux(1:18) ; arxifin(37:40)='.dat'                ! composing filename
     do im=1,40 ; if (arxifin(im:im)==" ") arxifin(im:im)="_" ; end do  ! composing filename   
     WRITE(*,*)'arxifin es:       ',7000,arxifin       
     open(7000,file=arxifin,status='unknown',action='write',iostat=ios)                       ! creating datafile 
     write(7000,*)'TARGETS (E1T1,E1T2,E2T1,E2T2)',block(1:2,1), block(1:2,2)  ! 1
     write(7000,*)'POPULATON SIZE...............',p                           ! 2
     write(7000,*)'STRENGHT OF SELECTION........',ss                          ! 3
     write(7000,*)'RECOMBINATION; 1=YES; 0=NO   ',reco                        ! 4
     write(7000,*)'TRAINING (1) vs TEST (0) SET ',training                    ! 5
     write(7000,*)'NUMBER /  TOTAL REPLICATES...',replica,replicas            ! 6
     write(7000,*)'CURRENT VS MAXIMUM GENERATION',et,etmax,lapso              ! 7
     write(7000,*)'ENV. FACTORS/ENVIRONMENTS....',EF,n                        ! 8
     write(7000,*)'NUMBER GENES,PHEN. DIMENSIONS',ng,PD                       ! 9
     write(7000,*)'TMAX,SDEV,SS,RECO,CAPPED.....',tmax,sdev,ss,reco,capped    ! 10   
     write(7000,*)'CONNECTIVITIES WW / MZZ .....',conWW,conMZZ                ! 11
     !!!!!!!!!!!!!!!!!!!!
     do pp=1,p
       do i=1,ind(1)%ngs
         write(7000,*)ind(pp)%w(i,:)
       end do  
       do i=1,ind(1)%ngs
         write(7000,*)ind(pp)%ww(i,:)
       end do       
     end do
     do i=1,ind(1)%ngs
       write(7000,*)ind(1)%MZ(i,:)
     end do
     do i=1,ind(1)%ngs
       write(7000,*)ind(1)%MZZ(i,:)
     end do
     close(7000)    
   end if  

   !!!! MUTATION IN THE FOUR MATRICES !!!!!!!     
   do pp=1,p
     call random_number(x)  ; call random_number(y)     
     i=int(x*real(ng)+1) ;  j=int(y*real(ng)+1)               ! which element of W will mutate
     747 call random_number(x)  ; call random_number(y)       ! random new value for the mutation
     z=sdev*sqrt(-2*log(x))*cos(2*pi*y)                       ! Box-Muller algotithm. Normal distribtion N(0,sdev) 
     ind(pp)%w(i,j)= ind(pp)%w(i,j)+z                         ! adding the new random value to the previous one
     if((capped.eq.1).and.((ind(pp)%w(i,j).gt.1.0).or.(ind(pp)%w(i,j).lt.-1.0)))then ; goto 747 ; end if   ! capped ??
     
     call random_number(x)  ; call random_number(y)     
     i=int(x*real(ng)+1) ;  j=int(y*real(ng)+1)               ! permito mutar un gen mas ...
     if(ind(pp)%ww(i,j).eq.0)then ; ind(pp)%ww(i,j)=1 ; else ;  ind(pp)%ww(i,j)=1 ; end if ! topological change
     
     !call random_number(x)  ; call random_number(y)     
     !i=int(x*real(ng)+1) ;  j=int(y*real(PD)+1)              ! permito mutar un gen mas ...
     !call random_number(z) ; z=1.0-2*z ; z=z*0.2
     !ind(1)%MZ(i,j)= ind(1)%MZ(i,j)+z 
     
     !call random_number(x)  ; call random_number(y)     
     !i=int(x*real(ng)+1) ;  j=int(y*real(PD)+1)              ! permito mutar un gen mas ...
     !if(ind(1)%MZZ(i,j).eq.0)then ; ind(1)%MZZ(i,j)=1 ; else ;  ind(1)%MZZ(i,j)=1 ; end if ! topological change      
   end do
 
end do     ! evolutionary time 
end do     ! replicates

close(20067) ; close(20068)     ! closes files
ret=SYSTEM('rm fort.*')         ! removes spurious stuffs
ret=SYSTEM('mv GRN_* files/')   ! replaces files into a folder

end program startodo










