module development

use start

implicit none

contains

!!!!!!!!!!!!!!!!!!!!
subroutine dev(i)    ! it runs development for the individual i 
integer :: i,ii,j,jj,k,jjj,ret,satur(ind(i)%ncels),saturado
real*4  :: r1,r2,u,q,y,z,fi,maxx,minn,xx ! fi=final increment
  
  saturado=0 ; ind(i)%sat=0	   
  do t=1,tmax        ! developmental time
      do j=1,ind(i)%ncels       ! for each cell of the individual
        do k=1,ind(i)%ngs!act     ! for each gene ef this cell       
          x=ind(i)%g(j,k)       ! concentration of gene k in cell j of ind
          q=0.0	                                           ! REACTION
          do jjj=1,ind(i)%ngs                                   
              if(ind(i)%ww(k,jjj).ne.0)then                 ! for all active gene interaction   
              q=q+ind(i)%w(k,jjj)*ind(i)%g(j,jjj)
            end if
          end do 
          q=0.5*q+ind(i)%epigen(k,j)                       ! EPIGENESIS         
          y=gen(k)%deg*x                                   ! DEGRADATION 
                     
          indt(i)%g(j,k)=x+tanh(q)-y                       ! t+1      ! iterative developmental function Kostas
          if(indt(i)%g(j,k).le.0.0)then ; indt(i)%g(j,k)=0.0 ; end if ! uncommented if POSITIVE STATE VARIABLE ?? CHOOSE YOURSELF :)
                   
        end do ! for each gene
      end do   ! for each cell      
      ind(i)%g=indt(i)%g ; indt(i)%g=0.0 

  end do       ! developmental time  
  
  ind(i)%phen=0.0   ! phenotyping
  do jjj=1,ind(i)%ncels 
    do t=1,pd 
      do j=1,ind(i)%ngs
        if(ind(i)%MZZ(j,t).ne.0)then
          ind(i)%phen(t,jjj)=ind(i)%phen(t,jjj)+ind(i)%g(jjj,j)*ind(i)%MZ(j,t)
        end if
      end do
    end do
  end do

end subroutine
!!!!!!!!!!!!!!!!!!

end module development
