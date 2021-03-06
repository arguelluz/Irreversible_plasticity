module development

use start

implicit none

contains

!!!!!!!!!!!!!!!!!!!!
subroutine dev(i)                                          ! it runs development for the individual i
integer :: i,ii,j,k,jj,kkk,pp                              ! local variables
real*4  :: u,q,y,z,stable                                  ! local variables

  do jjj=1,ng
    premutW (jjj,1:ng)=ind(i)%w(jjj,1:ng)                  ! Stores W matrix before mutation (reversible if unstable GRN)
    premutWW(jjj,1:ng)=ind(i)%ww(jjj,1:ng)                 ! Stores W matrix before mutation (reversible if unstable GRN)
  end do

  8881 do jjj=1,ng                                         ! REVERSE
    ind(i)%w (jjj,1:ng)=premutW (jjj,1:ng)                 ! Reverse if unstable GRN)
    ind(i)%ww(jjj,1:ng)=premutWW(jjj,1:ng)                 ! Reverse if unstable GRN)
    ind(i)%g(1:n,jjj)  =prepattern(1:n,jjj)                ! new Jan-2019
  end do

  pp=i                                                     ! Mutation in the generative matrices for each individual
  if(hillclimber.ne.1)then
    call mutation(pp)                                      ! independent subroutine (for stability criteria)
  else
    if(pp.gt.1)then
       call mutation(pp)                                   ! only one individual mutates in hill climber
    end if
  end if

  indt(i)%g=0.0
  do t=1,tmax                                              ! developmental time
      do j=1,ind(i)%ncels                                  ! for each cell of the individual
        do k=1,ind(i)%ngs                                  ! for each gene ef this cell
          x=ind(i)%g(j,k)                                  ! concentration of gene k in cell j of ind, minus degradation
          do jjj=1,ind(i)%ngs                              ! For all interacting genes jjj
            if(ind(i)%ww(k,jjj).ne.0)then                  ! if the interaction is active
              q=ind(i)%g(j,jjj)                            ! Set concentration of interacting gene (jjj) in same cell (j)
              q=(q+ind(i)%epigen(jjj,j))                   ! add and average with environmental component for interacting gene (jjj) in same cell (j)
              q=q*ind(i)%w(k,jjj)                          ! gene values by activation matrix (w)
              x=x+q                                        ! add activation to gene concentration
            end if
          end do
          ! in this loop x=concentration for the gene of interest at t-1, q=change in concentration from reactions and env, y=loss from degradation
          if(linear.ne.1)then
            x=tanh(x)                                      ! activation function
            x=(x+1.0)*0.5                                  ! rescale tanh output to logistic and store as gene value
          end if
          indt(i)%g(j,k)=x                                 ! concentration of target gene k in cell j at t+1
          if(linear.ne.1)then                              ! gene state variables can be negatives if we use linear activation function
            if(indt(i)%g(j,k).le.0.0)then
              indt(i)%g(j,k)=0.0
            end if ! uncommented if POSITIVE STATE VARIABLE. CHOOSE YOURSELF :)
          end if
          if(t.eq.tmax-1)then                              !stability criterium
            stab(j,k)=indt(i)%g(j,k)                       !stability criterium
          end if
        end do                                             ! end loop for each gene
      end do                                               ! end loop for each cell

      if(t.eq.tmax)then                                    ! stability criterium (Same as Dragui) (compares expression in tmax-1,tmax)
        stable=0.0                                         ! stability criterium (Same as Dragui)
        do j=1,ind(i)%ncels                                ! for each cell of the individual
          do k=1,ind(i)%ngs                                ! for each gene ef this cell
            stable=stable+(stab(j,k)-indt(i)%g(j,k))**2    ! Squared differences gene expression, all cells and environments
          end do
        end do
        if(sqrt(stable).gt.0.1)then                        ! stability criterium (Same as Dragui) ! THRESHOLD CAN BE CHANGED MANUALLY !
          goto 8881 !; write(*,*)'unstable'
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

end subroutine

!!!!!!!!!!!!!!!!!! !!!! MUTATION IN THE FOUR MATRICES !!!!!!!

subroutine mutation(pp)

integer :: pp,Mi,Mj,Mk
real*4  :: Mx,My,Mz

     59 call random_number(Mx) ; call random_number(My)
     Mi=int(Mx*real(ng)+1) ; Mj=int(My*real(ng)+1)                 ! which element of W will mutate
     if((Mi.le.ng).and.(Mi.ge.1).and.(Mj.le.ng).and.(Mj.ge.1))then ; goto 60 ; else ; goto 59 ; end if

     60 call random_number(Mx) ;  call random_number(My)           ! random new value for the mutation
     if((Mx.lt.1.0).and.(Mx.gt.0.0).and.(My.lt.1.0).and.(My.gt.0.0))then ; goto 92 ; else ; goto 60 ; end if  ! random new value for the mutation

     92 Mz=sdev*sqrt(-2*log(Mx))*cos(2*pi*My)                      ! Box-Muller algotithm. Normal distribtion N(0,sdev)
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
       62 call random_number(Mx) ;  if((Mx.lt.1.0).and.(Mx.gt.0.0))then ; goto 63 ; else ; goto 62 ; end if
       63 call random_number(My) ;  if((My.lt.1.0).and.(My.gt.0.0))then ; goto 64 ; else ; goto 63 ; end if
       64 Mi=int(Mx*real(ng)+1) ;  Mj=int(My*real(PD)+1)               ! Allow multiple mutations
       call random_number(Mz) ; Mz=1.0-2*Mz ; Mz=Mz*0.2
       ind(pp)%MZ(Mi,Mj)=ind(pp)%MZ(Mi,Mj)+Mz

       65 call random_number(Mx) ;  if((Mx.lt.1.0).and.(Mx.gt.0.0))then ; goto 66 ; else ; goto 65 ; end if
       66 call random_number(My) ;  if((My.lt.1.0).and.(My.gt.0.0))then ; goto 67 ; else ; goto 66 ; end if
       67 Mi=int(Mx*real(ng)+1) ;  Mj=int(My*real(PD)+1)               ! Allow multiple mutations
       if(ind(pp)%MZZ(Mi,Mj).eq.0)then ; ind(pp)%MZZ(Mi,Mj)=1 ; else ;  ind(1)%MZZ(Mi,Mj)=1 ; end if ! topological change
     end if
end subroutine

end module development
