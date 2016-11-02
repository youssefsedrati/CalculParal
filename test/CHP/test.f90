program test
  use grad
  implicit none
  integer,parameter::Nx=3,Ny=3
  real*8,dimension(Nx*Ny)::X,AX
  real*8::B,Cx,Cy,D=1.,dt,dx,dy,Lx=1.,Ly=1.
  integer::i
  
  dx=Lx/(1+Nx)
  dy=Ly/(1+Ny)
  B=3
  Cx=1
  Cy=1

  do i=1,size(X)
     X(i)=i*1.
  end do

  call prodAx(B,Cx,Cy,X,Nx,Ny,AX)

  print*,AX



end program test
