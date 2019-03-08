program integral ! auxiliary program, calculates the maximum gene value of an individual in the Irreversible paslticity (Dragui) model
! Solves the maximum final gene concentration for the hidden layer genes
implicit none

integer :: i,j,k,ii,jj,kk,n,ng,t,tmax,ngs
real*4  :: q,u,v,x,y,z,g(2),gt(2),w,env,decay

tmax=20   ! developmental time
ngs=2     ! numnber of genes
w=1.0     ! maxval W matrix
g=1.0     ! Initial conditions
env=0.5   ! max environmental input
decay=0.0 ! gene degradation per iteration (proportion)

open(1, file="MaxVal.tsv", status="new", action="write") ! Create output file
write(1,*)'tMax',tmax,'nGenes',ngs,'maxMatrix',w,'initialConc',g ! Write as tab separated
close(1)

 do t=1,tmax                                               ! developmental time
        do k=1,ngs                                         ! for each gene ef this cell
          x=g(k)                                           ! concentration of gene k in cell j of ind
          q=(x+env)*0.5                                    ! add and average concentration and environment
          do jj=1,ngs
              q=q+w*g(jj)
          end do
          y=decay*x                                        ! DEGRADATION
          gt(k)=tanh(x-y)                                  ! t+1      ! iterative developmental function Kostas
          gt(k)=(gt(k)+1)*0.5                              ! rescale tanh to logistic
          if(gt(k).le.0.0)then ; gt(k)=0.0 ; end if        ! POSITIVE STATE VARIABLE. CHOOSE YOURSELF :)
      end do

      open(1, file="MaxVal.tsv", status="old", position="append", action="write") ! Open output file in append mode
      write(1,*)t,'initial',x,'react',q,'degr',y,'FINAL',Gt(K)                    ! Write as tab separated
      close(1)                                                                    ! Close output file

      g=gt
 end do



end program
