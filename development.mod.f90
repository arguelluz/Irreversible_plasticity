module development

use start

implicit none

contains

!!!!!!!!!!!!!!!!!!!!
subroutine dev(i)                                          ! it runs development for the individual i
integer :: i,ii,j,jj,pp,k,jjj,ret
real*4  :: r1,r2,u,q,y,z,fi,xx,stable,eps                  ! fi=final increment

  !write(*,*)'WDEV1',ind(I)%w(1,:),ind(I)%w(2,:)
  do jjj=1,ng
    premutW (jjj,1:ng)=ind(i)%w(jjj,1:ng)                  ! Stores W matrix before mutation (reversible if unstable GRN)
    !WRITE(*,*)i,jjj,'indi',ind(i)%w(jjj,1:ng)
    premutWW(jjj,1:ng)=ind(i)%ww(jjj,1:ng)                 ! Stores W matrix before mutation (reversible if unstable GRN)
  end do

  8881 do jjj=1,ng                                         ! REVERSE
    ind(i)%w (jjj,1:ng)=premutW (jjj,1:ng)                 ! Reverse if unstable GRN)
    ind(i)%ww(jjj,1:ng)=premutWW(jjj,1:ng)                 ! Reverse if unstable GRN)
    ind(i)%g(1:n,jjj)=prepattern(1:n,jjj)                  ! new Jan-2019
  end do
  !write(*,*)'WDEV2',ind(I)%w(1,:),ind(I)%w(2,:)

  pp=i                                                     ! Mutation in the generative matrices for each individual
  if(hillclimber.ne.1)then
    call mutation(pp)                                      ! independent subroutine (for stability criteria)
  else
    if(pp.gt.1)then
       call mutation(pp)                                   ! only one individual mutates in hill climber
    end if
  end if

  ind(i)%sat=0
  indt(i)%g=0.0
  do t=1,tmax                                              ! developmental time
      !write(*,*)'i,t',et,i,t
      do j=1,ind(i)%ncels                                  ! for each cell of the individual
        do k=1,ind(i)%ngs                                  ! for each gene ef this cell
          x=ind(i)%g(j,k)                                  ! concentration of gene k in cell j of ind
!          q=ind(i)%epigen(k,j)	                           ! baseline gene input values, set as environment
!          call random_number(eps)                          ! sample random noise
!          q=q+((eps-0.5)*0.02)                             ! center, scale and add random noise (-.01,.01) to environment (same noise for all interactions)
          do jjj=1,ind(i)%ngs
            if(ind(i)%ww(k,jjj).ne.0)then                  ! for all active gene interaction
              q=ind(i)%g(j,jjj)                            ! Set gene concentration to same value as t-1
              q=(q+ind(i)%epigen(k,jjj))*0.5                 ! add and average with environmental component
              !q=(q+maxepigen)*0.5                          ! add and average maximum env input (for invariable environments)
              q=q*ind(i)%w(k,jjj)                          ! gene values by activation matrix (w)
              x=x+q                                        ! add activation to gene concentration
            end if
          end do
          !y=gen(k)%deg*x                                   ! DEGRADATION multiplied by current gene concentration (x)
          ! in this loop x=concentration for the gene of interest at t-1, q=change in concentration from reactions and env, y=loss from degradation
          x=tanh(x)                             ! t+1      ! activation function
          indt(i)%g(j,k)=(x+1)*0.5                         ! rescale tanh output to logistic and store as gene value
          if(indt(i)%g(j,k).le.0.0)then ; indt(i)%g(j,k)=0.0 ; end if ! uncommented if POSITIVE STATE VARIABLE. CHOOSE YOURSELF :)
          if(t.eq.tmax-1)then                              !stability criterium
            stab(j,k)=indt(i)%g(j,k)                       !stability criterium
          end if
        end do                                             ! end loop for each gene
      end do                                               ! end loop for each cell

      if(t.eq.tmax)then                                    ! stability criterium (Same as Dragui) (compares expression in tmax-1,tmax)
        !write(*,*)'maxvalG=',maxval(ind(i)%g)             ! maximum values
        stable=0.0                                         ! stability criterium (Same as Dragui)
        do j=1,ind(i)%ncels                                ! for each cell of the individual
          do k=1,ind(i)%ngs                                ! for each gene ef this cell
            stable=stable+(stab(j,k)-indt(i)%g(j,k))**2    ! Squared differences gene expression, all cells and environments
          end do
        end do
        if(sqrt(stable).gt.0.1)then                        ! stability criterium (Same as Dragui) ! THRESHOLD CAN BE CHANGED MANUALLY !
          goto 8881 ; write(*,*)'unstable'
        end if                                             ! stability criterium (Same as Dragui)
      end if

      ind(i)%g=indt(i)%g ; indt(i)%g=0.0                   ! "valid" individuals are updated and the loop closed

  end do                                                   ! end loop developmental time

  ind(i)%phen=0.0                                          ! phenotyping
  do jjj=1,ind(i)%ncels                                    ! For each environment
    do t=1,pd                                              ! For each trait
      do j=1,ind(i)%ngs                                    ! For each gene ! Set 1 to EF to exclude environmentally-sensitive genes.
        if(ind(i)%MZZ(j,t).ne.0)then
          ind(i)%phen(t,jjj)=ind(i)%phen(t,jjj)+ind(i)%g(jjj,j)*ind(i)%MZ(j,t) ! Calculate new phenotype by adding all gene value by MZ matrix
        end if
      end do
    end do
  end do

  !do jjj=1,ng
  ! WRITE(*,*)i,jjj,' indiP',ind(i)%w(jjj,1:ng)
  !end do

end subroutine

!!!!!!!!!!!!!!!!!! !!!! MUTATION IN THE FOUR MATRICES !!!!!!!

subroutine mutation(pp)

integer :: pp,Mi,Mj,Mk
real*4  :: Mx,My,Mz

   call random_number(Mx)  ; call random_number(My)
     Mi=int(Mx*real(ng)+1) ;  Mj=int(My*real(ng)+1)                ! which element of W will mutate
     call random_number(Mx)  ; call random_number(My)          ! random new value for the mutation
     Mz=sdev*sqrt(-2*log(Mx))*cos(2*pi*My)                         ! Box-Muller algotithm. Normal distribtion N(0,sdev)
     !WRITE(*,*)pp,'    PREMUTw',Mi,Mj,ind(pp)%w(Mi,Mj)
     ind(pp)%w(Mi,Mj)= ind(pp)%w(Mi,Mj)+Mz                         ! adding the new random value to the previous one
     if(capped.eq.1) then                                          ! If we are using capped weights
       if(ind(pp)%w(Mi,Mj).gt.1.0) then                            ! and if the mutation makes weights greater than 1
         ind(pp)%w(Mi,Mj)=1.0                                      ! Then set them to 1
       else if (ind(pp)%w(Mi,Mj).lt.-1.0) then                     ! if the mutation makes weights lesser than -1
         ind(pp)%w(Mi,Mj)=-1.0                                     ! Then set them to -1
       end if
     end if

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
