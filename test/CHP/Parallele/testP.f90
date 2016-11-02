program parallele
  !use grad
  use MPI
  implicit none
  !include "mpif.h"
  integer::Nx,Ny
  integer::Me,statinfo,tag,i,j,i1,iN,Np,k
  integer,dimension(MPI_STATUS_SIZE)::status
  real*8::B,Cx,Cy,D,Lx,Ly,dt,dx,dy
  real*8,dimension(:),allocatable::U,UP,U0,AX,AXP,ResultatP,Xn,Xp
  character(len=20)::parametres='param.txt'

  call MPI_INIT(statinfo)
  call MPI_COMM_RANK(MPI_COMM_WORLD,Me,statinfo)
  tag=100
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)

  open (unit = 10, file =parametres , form = 'formatted', status = 'old', action = 'read')
  read(10,*)Nx,Ny,Lx,Ly,D
  close(10)

  call charge(Nx*Ny,Np,me,i1,iN)
  allocate(Xn(Nx),Xp(Nx))
  allocate(U(Nx*Ny))
  allocate(UP(iN-i1+1))
  allocate(AX(Nx*Ny))
  allocate(AXP(iN-i1+1))
  allocate(ResultatP(Nx*Ny))



  !stationnaire
  dx=Lx/(1+Nx)
  dy=Ly/(1+Ny)
  Cx=-D/dx**2
  Cy=-D/dy**2
  B=2*D/dx**2+2*D/dy**2
  U=1.
  UP=1.
  B=3.
  Cx=1.
  Cy=1.
  Xn=1.
  Xp=1.
  call prodAx(B,Cx,Cy,U,Nx,Ny,AX)
  call MatMulP(AXP,UP,Xp,Xn,B,Cx,Cy,Nx,Ny,i1,iN)
  if(me==0)then
     print*, me,i1,iN
     print*, AXP,me
     print*, AX(i1:iN),me
  end if



  deallocate(U)
  deallocate(UP)
  deallocate(AX)
  deallocate(AXP)
  deallocate(ResultatP)



  call MPI_FINALIZE(statinfo)

contains

  subroutine charge(n,Np,me,i1,iN)
    implicit none
    integer, intent (in)::n,Np,me
    integer, intent (out)::i1,iN
    integer::q
    q=n/Np
    if(me<mod(n,Np))then
       i1=me*(q+1)+1
       iN=(me+1)*(1+q)
    else
       i1=mod(n,Np)+1+me*q
       iN=i1+q-1
    end if
  end subroutine charge

  subroutine prodAx(B,Cx,Cy,X,Nx,Ny,AX)
    implicit none
    integer,intent(in)::Nx,Ny
    real*8,intent(in)::B,Cx,Cy
    real*8,dimension(Nx*Ny),intent(in)::X
    real*8,dimension(Nx*Ny),intent(out)::AX
    integer::i

    AX(1)=B*X(1)+Cx*X(2)+Cy*X(Nx+1)
    do i=2,Nx-1
       AX(i)=B*X(i)+Cx*(X(i-1)+X(i+1))+Cy*X(Nx+i)
    end do
    AX(Nx)=B*X(Nx)+Cx*X(Nx-1)+Cy*X(2*Nx)

    do i=Nx+1,Nx*Ny-Nx
       if(mod(i,Nx)==1)then
          AX(i)=Cy*(X(i-Nx)+X(i+Nx))+Cx*X(i+1)+B*X(i)
       else if(mod(i,Nx)==0)then
          AX(i)=Cy*(X(i-Nx)+X(i+Nx))+Cx*X(i-1)+B*X(i)
       else
          AX(i)=Cy*(X(i-Nx)+X(i+Nx))+Cx*(X(i+1)+X(i-1))+B*X(i)
       end if
    end do

    AX(Nx*Ny-Nx+1)=B*X(Nx*Ny-Nx+1)+Cx*X(Nx*Ny-Nx+2)+Cy*X(Nx*Ny-2*Nx+1)
    do i=Nx*Ny-Nx+2,Nx*Ny-1
       AX(i)=B*X(i)+Cx*(X(i-1)+X(i+1))+Cy*X(i-Nx)
    end do
    AX(Nx*Ny)=B*X(Nx*Ny)+Cx*X(Nx*Ny-1)+Cy*X(Nx*Ny-Nx)

  end subroutine prodAx

!  subroutine MatMulAParallele(AX,X0,B,Cx,Cy,Nx,Ny,i1,iN)
!    implicit none
!    integer, intent (in)::Nx,Ny,i1,iN
!    real*8, intent(in)::B,Cx,Cy
!    real*8, dimension(:),intent(in)::X0
!    real*8, dimension(iN-i1+1),intent(out)::AX
!    integer::j
!    if(i1<Nx)then
!       Do j=1,iN-i1+1
!          !Premier Bloc
!          if((j+i1-2)/(Nx)==0)then
!             if(mod(j+i1-1,Nx)==1)then
!                AX(j)=B*X0(j+i1)+Cx*X0(j+i1+1)+Cy*X0(j+i1+Nx)
!             elseif(mod(j+i1-1,Nx)==0)then
!                AX(j)=B*X0(j+i1)+Cx*X0(j+i1-1)+Cy*X0(j+i1+Nx)
!             else
!                AX(j)=B*X0(j+i1)+Cx*(X0(j+i1-1)+X0(j+i1+1))+Cy*X0(j+i1+Nx)
!             end if
!             !Dernier Bloc
!          elseif((j+i1-2)/Nx==Ny-1)then
!             if(mod(j+i1-1,Nx)==1)then
!                AX(j)=B*X0(j+i1)+Cx*X0(j+i1+1)+Cy*X0(j+i1-Nx)
!             elseif(mod(j+i1-1,Nx)==0)then
!                AX(j)=B*X0(j+i1)+Cx*X0(j+i1-1)+Cy*X0(j+i1-Nx)
!             else
!                AX(j)=B*X0(j+i1)+Cx*(X0(j+i1-1)+X0(j+i1+1))+Cy*X0(j+i1-Nx)
!             end if
!             !Bloc du Milieu
!          else
!             if(mod(j+i1-1,Nx)==1)then
!                AX(j)=B*X0(j+i1)+Cx*X0(j+i1+1)+Cy*(X0(j+i1-Nx)+X0(j+i1+Nx))
!             elseif(mod(j+i1-1,Nx)==0)then
!                AX(j)=B*X0(j+i1)+Cx*X0(j+i1-1)+Cy*(X0(j+i1-Nx)+X0(j+i1+Nx))
!             else
!                AX(j)=B*X0(j+i1)+Cx*(X0(j+i1-1)+X0(j+i1+1))+Cy*(X0(j+i1-Nx)+X0(j+i1+Nx))
!             end if
!          end if
!       end do
!    else
!       Do j=1,iN-i1+1
!          !Premier Bloc
!          if((j+i1-2)/(Nx)==0)then
!             if(mod(j+i1-1,Nx)==1)then
!                AX(j)=B*X0(j+Nx)+Cx*X0(j+Nx+1)+Cy*X0(j+Nx+Nx)
!             elseif(mod(j+i1-1,Nx)==0)then
!                AX(j)=B*X0(j+Nx)+Cx*X0(j+Nx-1)+Cy*X0(j+Nx+Nx)
!             else
!                AX(j)=B*X0(j+Nx)+Cx*(X0(j+Nx-1)+X0(j+Nx+1))+Cy*X0(j+Nx+Nx)
!             end if
!             !Dernier Bloc
!          elseif((j+i1-2)/Nx==Ny-1)then
!             if(mod(j+i1-1,Nx)==1)then
!                AX(j)=B*X0(j+Nx)+Cx*X0(j+Nx+1)+Cy*X0(j+Nx-Nx)
!             elseif(mod(j+i1-1,Nx)==0)then
!                AX(j)=B*X0(j+Nx)+Cx*X0(j+Nx-1)+Cy*X0(j+Nx-Nx)
!             else
!                AX(j)=B*X0(j+Nx)+Cx*(X0(j+Nx-1)+X0(j+Nx+1))+Cy*X0(j+Nx-Nx)
!             end if
!             !Bloc du Milieu
!          else
!             if(mod(j+i1-1,Nx)==1)then
!                AX(j)=B*X0(j+Nx)+Cx*X0(j+Nx+1)+Cy*(X0(j+Nx-Nx)+X0(j+Nx+Nx))
!             elseif(mod(j+i1-1,Nx)==0)then
!                AX(j)=B*X0(j+Nx)+Cx*X0(j+Nx-1)+Cy*(X0(j+Nx-Nx)+X0(j+Nx+Nx))
!             else
!                AX(j)=B*X0(j+Nx)+Cx*(X0(j+Nx-1)+X0(j+Nx+1))+Cy*(X0(j+Nx-Nx)+X0(j+Nx+Nx))
!             end if
!          end if
!       end do
!    end if
!  end subroutine MatMulAParallele

!  subroutine MatMulAParallele(AX,X0,B,Cx,Cy,Nx,Ny,i1,iN)
!    implicit none
!    integer, intent (in)::Nx,Ny,i1,iN
!    real*8, intent(in)::B,Cx,Cy
!    real*8, dimension(:),intent(in)::X0
!    real*8, dimension(iN-i1+1),intent(out)::AX
!    integer::j,proc
!
!    proc=i1/Nx+1
!    Do j=1,iN-i1+1
!       !Premier Bloc
!       if((j+i1-2)/Nx==0)then
!          if(mod(j+i1-1,Nx)==1)then
!             AX(j)=B*X0(j)+Cx*X0(j+1)+Cy*X0(j+Nx)
!          elseif(mod(j+i1-1,Nx)==0)then
!             AX(j)=B*X0(j)+Cx*X0(j-1)+Cy*X0(j+Nx)
!          else
!             AX(j)=B*X0(j)+Cx*(X0(j-1)+X0(j+1))+Cy*X0(j+Nx)
!          end if
!          !Dernier Bloc
!       elseif((j+i1-2)/Nx==Ny-1)then
!          if(mod(j+i1-1,Nx)==1)then !1ere ligne
!             AX(j)=B*X0(j)+Cx*X0(j+1)+Cy*X0(j-Nx)
!          elseif(mod(j+i1-1,Nx)==0)then !dernière ligne
!             AX(j)=B*X0(j)+Cx*X0(j-1)+Cy*X0(j-Nx)
!          else
!             AX(j)=B*X0(j)+Cx*(X0(j-1)+X0(j+1))+Cy*X0(j-Nx)
!          end if
!          !Bloc du Milieu
!       else
!          if(mod(j+i1-1,Nx)==1)then ! 1ere ligne
!             AX(j)=B*X0(j)+Cx*X0(j+1)+Cy*(X0(j-Nx)+X0(j+Nx))
!          elseif(mod(j+i1-1,Nx)==0)then !derniere ligne
!             AX(j)=B*X0(j)+Cx*X0(j-1)+Cy*(X0(j-Nx)+X0(j+Nx))
!          else
!             AX(j)=B*X0(j)+Cx*(X0(j-1)+X0(j+1))+Cy*(X0(j-Nx)+X0(j+Nx))
!          end if
!       end if
!       print*,AX(j),j+i1-1
!    end do
!  end subroutine MatMulAParallele

subroutine MatMulP(AX,X,Xp,Xn,B,Cx,Cy,Nx,Ny,i1,iN)
    integer,intent(in)::i1,iN,Nx,Ny
    real*8,dimension(iN-i1+1),intent(in)::X
    real*8,dimension(Nx),intent(in)::Xp,Xn
    real*8,dimension(iN-i1+1),intent(out)::AX
    real*8,intent(in)::B,Cx,Cy

    do j=1,iN-i1+1

        !Premier Bloc
        if((j+i1-2)/Nx==0)then
            if(mod(j+i1-1,Nx)==1)then !Première ligne
                AX(j)=B*X(j)
                if(j+1>iN)then
                   AX(j)=AX(j)+Cx*Xn(1)
                else
                   AX(j)=AX(j)+Cx*X(j+1)
                end if
                if(j+Nx>iN)then
                    AX(j)=AX(j)+Cy*Xn(j)
                else
                    AX(j)=AX(j)+Cy*X(j+Nx)
                end if

            elseif(mod(j+i1-1,Nx)==0)then !Dernière ligne
                AX(j)=B*X(j)
                if(j-1<i1)then
                    AX(j)=AX(j)+Cx*Xp(Nx)
                else
                    AX(j)=AX(j)+Cx*X(j-1)
                end if
                if(j+Nx>iN)then
                    AX(j)=AX(j)+Cy*Xn(j)
                else
                    AX(j)=AX(j)+Cy*X(j+Nx)
                end if

            else
                AX(j)=B*X(j)
                if(j-1<i1)then
                    AX(j)=AX(j)+Cx*Xp(Nx)
                else
                    AX(j)=AX(j)+Cx*X(j-1)
                end if
                if(j+1>iN)then
                    AX(j)=AX(j)+Cx*Xn(1)
                else
                    AX(j)=AX(j)+Cx*X(j+1)
                end if
                if(j+Nx>iN)then
                    AX(j)=AX(j)+Cy*Xn(j)
                else
                    AX(j)=AX(j)+Cy*X(j+Nx)
                end if
            end if

        !Dernier Bloc
        elseif((j+i1-2)/Nx==Ny-1)then

            if(mod(j+i1-1,Nx)==1)then !Première ligne
                AX(j)=B*X(j)
                if(j-Nx<i1)then
                    AX(j)=AX(j)+Cy*Xp(j)
                else
                    AX(j)=AX(j)+Cy*X(j-Nx)
                end if
                if(j+1>iN)then
                    AX(j)=AX(j)+Cx*Xn(1)
                else
                    AX(j)=AX(j)+Cx*X(j+1)
                end if

            elseif(mod(j+i1-1,Nx)==0)then !Dernière ligne
                AX(j)=B*X(j)
                if(j-Nx<i1)then
                    AX(j)=AX(j)+Cy*Xp(j)
                else
                    AX(j)=AX(j)+Cy*X(j-Nx)
                end if
                if(j-1<i1)then
                    AX(j)=AX(j)+Cx*Xp(Nx)
                else
                    AX(j)=AX(j)+Cx*X(j-1)
                end if

            else
                AX(j)=B*X(j)
                if(j-Nx<i1)then
                    AX(j)=AX(j)+Cy*Xp(j)
                else
                    AX(j)=AX(j)+Cy*X(j-Nx)
                end if
                if(j-1<i1)then
                    AX(j)=AX(j)+Cx*Xp(Nx)
                else
                    AX(j)=AX(j)+Cx*X(j-1)
                end if
                if(j+1>iN)then
                    AX(j)=AX(j)+Cx*Xn(1)
                else
                    AX(j)=AX(j)+Cx*X(j+1)
                end if
            end if

          !Bloc du Milieu
        else

          if(mod(j+i1-1,Nx)==1)then ! Première ligne
            AX(j)=B*X(j)
            if(j-Nx<i1)then
                AX(j)=AX(j)+Cy*Xp(j)
            else
                AX(j)=AX(j)+Cy*X(j-Nx)
            end if
            if(j+1>iN)then
                AX(j)=AX(j)+Cx*Xn(1)
            else
                AX(j)=AX(j)+Cx*X(j+1)
            end if
            if(j+Nx>iN)then
                AX(j)=AX(j)+Cy*Xn(j)
            else
                AX(j)=AX(j)+Cy*X(j+Nx)
            end if

          elseif(mod(j+i1-1,Nx)==0)then !Dernière ligne
            AX(j)=B*X(j)
            if(j-Nx<i1)then
                AX(j)=AX(j)+Cy*Xp(j)
            else
                AX(j)=AX(j)+Cy*X(j-Nx)
            end if
            if(j-1<i1)then
                AX(j)=AX(j)+Cx*Xp(Nx)
            else
                AX(j)=AX(j)+Cx*X(j-1)
            end if
            if(j+Nx>iN)then
                AX(j)=AX(j)+Cy*Xn(j)
            else
                AX(j)=AX(j)+Cy*X(j+Nx)
            end if

          else
            AX(j)=B*X(j)
            if(j-Nx<i1)then
                AX(j)=AX(j)+Cy*Xp(j)
            else
                AX(j)=AX(j)+Cy*X(j-Nx)
            end if
            if(j-1<i1)then
                AX(j)=AX(j)+Cx*Xp(Nx)
            else
                AX(j)=AX(j)+Cx*X(j-1)
            end if
            if(j+1>iN)then
                AX(j)=AX(j)+Cx*Xn(1)
            else
                AX(j)=AX(j)+Cx*X(j+1)
            end if
            if(j+Nx>iN)then
                AX(j)=AX(j)+Cy*Xn(j)
            else
                AX(j)=AX(j)+Cy*X(j+Nx)
            end if

          end if
       end if

    end do

end subroutine



end program parallele
