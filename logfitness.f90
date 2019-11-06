program logfitness
implicit none

integer :: i,j,k,kk,n,p,logp,whois
real*4  :: x,y,z,w,fmax,delta
real*4, allocatable :: fitness(:), fitnesst(:)

delta=0.00001
p=10                              ! population
logp=1+int(log(real(p))/log(2d0))
write(*,*)'logp',logp

allocate(fitness(p),fitnesst(p))
fitnesst=0.0

do i=1,p                              ! random fitness assignment
  call random_number(x)
  fitness(i)=real(int(10.0*x))-5.0
end do
fitness(1:p)=(/1000001,1000000,1000000,1000000,1000001,1000000,1000000,1000000,1000000,1000000/)
write(*,*)'fit1',fitness(:)

fitness=fitness-minval(fitness)+delta
write(*,*)fitness(:)/maxval(fitness)

do i=2,p                               ! suma
  fitness(i)=fitness(i)+fitness(i-1)
end do
write(*,*)'fit2',fitness(:)

fmax=maxval(fitness(:))               ! first relative fitness (0-1)
do i=1,p
  fitness(i)=fitness(i)/fmax
end do
write(*,*)'fit3',fitness(:)

do i=1,1000!p
  !read(*,*)
  call random_number(x)
  write(*,*)fitness(:)
  write(*,*)'x',x
  kk=p/2 ; k=p/2
  do j=1,logp
    k=k/2
    if(j.eq.logp-1)then;k=1;end if     ! allows for p non congruent log2(p)
    write(*,*)'engancha',kk,fitness(kk)
    if(x.le.fitness(kk))then
      write(*,*)'kk',k,kk,fitness(kk),'menos'
      if(j.eq.logp)then;whois=kk;end if
      kk=kk-k
      if(kk.eq.0)then ; whois=1 ; goto 754 ; end if
    else if(x.gt.fitness(kk))then
      write(*,*)'kk',k,kk,fitness(kk),'mas'
      if(j.eq.logp)then;whois=kk+1;end if
      kk=kk+k
      if(kk.gt.p)then ; whois=p ; goto 754 ; end if
    end if
    write(*,*)'kk va por',kk
  end do
  754 write(*,*)'final',whois

  kk=0
  do j=1,p
    kk=kk+1
    if(x.lt.fitness(j))then ; write(*,*)x,j,fitness(j); exit ; endif
  end do
  write(*,*)'check',kk

  if(kk.ne.whois)then ; write(*,*)fitness(:); read(*,*) ; end if
  !fitnesst(i)=whois! fill up next generation

  fitnesst(whois)=fitnesst(whois)+1
end do
  fitnesst=fitnesst/maxval(fitnesst)
  write(*,*)fitnesst(:)

end program
