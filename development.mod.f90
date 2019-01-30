module development

use start

implicit none

contains

!!!!!!!!!!!!!!!!!!!!
subroutine dev(i)                                          ! it runs development for the individual i
integer :: i,ii,j,jj,pp,k,jjj,ret
real*4  :: r1,r2,u,q,y,z,fi,xx,stable                      ! fi=final increment

  do jjj=1,ng
    premutW (jjj,1:ng)=ind(i)%w (jjj,1:ng)                 ! Stores W matrix before mutation (reversible if unstable GRN)
    premutWW(jjj,1:ng)=ind(i)%ww(jjj,1:ng)                 ! Stores W matrix before mutation (reversible if unstable GRN)
  end do

  do jjj=1,ng                                         ! REVERSE
    ind(i)%w (jjj,1:ng)=premutW (jjj,1:ng)                 ! Reverse if unstable GRN)
    ind(i)%ww(jjj,1:ng)=premutWW(jjj,1:ng)                 ! Reverse if unstable GRN)
  end do

  pp=i                                                     ! Mutation in the generative matrices for each individual
  call mutation(pp)                                        ! independent subroutine (for stability criteria)

  ind(i)%sat=0
  indt(i)%g=0.0
  do t=1,tmax                                              ! developmental time
      do j=1,ind(i)%ncels                                  ! for each cell of the individual
        do k=1,ind(i)%ngs                                  ! for each gene ef this cell
          x=ind(i)%g(j,k)                                  ! concentration of gene k in cell j of ind
          q=0.0	                                           ! REACTION
          do jjj=1,ind(i)%ngs
              if(ind(i)%ww(k,jjj).ne.0)then                ! for all active gene interaction
              q=q+ind(i)%w(k,jjj)*ind(i)%g(j,jjj)
            end if
          end do
          q=0.5*q+ind(i)%epigen(k,j)                       ! EPIGENESIS
          y=gen(k)%deg*x                                   ! DEGRADATION

          indt(i)%g(j,k)=x+tanh(q)-y                       ! t+1      ! iterative developmental function Kostas
          if(indt(i)%g(j,k).le.0.0)then ; indt(i)%g(j,k)=0.0 ; end if ! uncommented if POSITIVE STATE VARIABLE. CHOOSE YOURSELF :)
          if(t.eq.tmax-1)then                              !stability criterium
            stab(j,k)=indt(i)%g(j,k)                       !stability criterium
          end if
        end do                                             ! end loop for each gene
      end do                                               ! end loop for each cell

!      if(t.eq.tmax)then                                    ! stability criterium (Same as Dragui) (compares expression in tmax-1,tmax)
!        stable=0.0                                         ! stability criterium (Same as Dragui)
!        do j=1,ind(i)%ncels                                ! for each cell of the individual
!          do k=1,ind(i)%ngs                                ! for each gene ef this cell
!            stable=stable+(stab(j,k)-indt(i)%g(j,k))**2    ! Squared differences gene expression, all cells and environments
!          end do
!        end do
!        if(sqrt(stable).gt.0.1)then                        ! stability criterium (Same as Dragui) ! THRESHOLD CAN BE CHANGED MANUALLY !
!          goto 8881
!        end if                                             ! stability criterium (Same as Dragui)
!      end if

      ind(i)%g=indt(i)%g ; indt(i)%g=0.0                   ! "valid" individuals are updated and the loop closed

  end do                                                   ! end loop developmental time

  ind(i)%phen=0.0                                          ! phenotyping
  do jjj=1,ind(i)%ncels
    do t=1,pd
      do j=EF,ind(i)%ngs                                          ! Excluding environmentally-sensitive genes
        if(ind(i)%MZZ(j,t).ne.0)then
          ind(i)%phen(t,jjj)=ind(i)%phen(t,jjj)+ind(i)%g(jjj,j)*ind(i)%MZ(j,t)
        end if
      end do
    end do
  end do
end subroutine

!!!!!!!!!!!!!!!!!! !!!! MUTATION IN THE FOUR MATRICES !!!!!!!

subroutine mutation(pp)

integer :: pp,Mi,Mj,Mk
real*4  :: Mx,My,Mz

   call random_number(Mx)  ; call random_number(My)
     Mi=int(Mx*real(ng)+1) ;  Mj=int(My*real(ng)+1)                ! which element of W will mutate
     747 call random_number(Mx)  ; call random_number(My)          ! random new value for the mutation
     Mz=sdev*sqrt(-2*log(Mx))*cos(2*pi*My)                         ! Box-Muller algotithm. Normal distribtion N(0,sdev)
     ind(pp)%w(Mi,Mj)= ind(pp)%w(Mi,Mj)+Mz                         ! adding the new random value to the previous one
     if((capped.eq.1).and.((ind(pp)%w(Mi,Mj).gt.1.0).or.(ind(pp)%w(Mi,Mj).lt.-1.0)))then
     goto 747 ; end if   ! capped ??

     call random_number(Mx)  ; call random_number(My)
     i=int(Mx*real(ng)+1) ;  j=int(My*real(ng)+1)                  ! Allow multiple mutations
     if(ind(pp)%ww(Mi,Mj).eq.0)then ; ind(pp)%ww(Mi,Mj)=1 ; else ;  ind(pp)%ww(Mi,Mj)=0 ; end if ! topological change
     
     if(mzadhoc.ne.1)then                                           ! If Mz and Mzz are pre-specified they do not mutate
       call random_number(Mx)  ; call random_number(My)     
       Mi=int(Mx*real(ng)+1) ;  Mj=int(My*real(PD)+1)               ! Allow multiple mutations
       call random_number(Mz) ; Mz=1.0-2*Mz ; Mz=Mz*0.2
       ind(pp)%MZ(Mi,Mj)=ind(pp)%MZ(Mi,Mj)+Mz 
     
       call random_number(Mx)  ; call random_number(My)     
       Mi=int(Mx*real(ng)+1) ;  Mj=int(My*real(PD)+1)               ! Allow multiple mutations
       if(ind(pp)%MZZ(Mi,Mj).eq.0)then ; ind(pp)%MZZ(Mi,Mj)=1 ; else ;  ind(1)%MZZ(Mi,Mj)=1 ; end if ! topological change    
     end if
end subroutine



end module development
