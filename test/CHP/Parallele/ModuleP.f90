module ModuleP
use fonctions
!use MPI
implicit none

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

  subroutine GradConj(X0,X,B,Cx,Cy,rhs,Nx,Ny,i1,iN,Np,me,statinfo,status)
    implicit none
    include "mpif.h"
    integer,intent(in)::Nx,Ny,i1,iN,Np,me,statinfo
    integer,dimension(:),intent(in):: status
    real*8,intent(in)::B,Cx,Cy
    real*8,dimension(iN-i1+1),intent(in)::rhs,X0
    real*8,dimension(size(X0)),intent(out)::X
    real*8,dimension(iN-i1+1)::d,W,R,R_next
    real*8,dimension(Nx)::Xp,Xn
    real*8::alpha,beta,e,f,Ep,Fp,NormeR
    real*8,parameter::eps=1.e-5
    integer,parameter::N=10000
    integer::k,j,tag


    !Initialisation
    tag=100
    Xp=0.
    Xn=0.
    R=0.
    W=0.
    !R=matmul(A,X0)-rhs
    !calculer Xp et Xn
    if((me/=0) .and. (me/=Np-1))then !on envoie Xp à me-1 et Xn à me+1
       call MPI_SEND(X0(iN-i1+1-Nx+1:iN-i1+1),Nx,MPI_REAL8,me+1,tag,MPI_COMM_WORLD,statinfo)
       call MPI_SEND(X0(1:Nx),Nx,MPI_REAL8,me-1,tag,MPI_COMM_WORLD,statinfo)
    else if(me==0)then !on envoie Xn à me=1
       call MPI_SEND(X0(iN-i1+1-Nx+1:iN-i1+1),Nx,MPI_REAL8,me+1,tag,MPI_COMM_WORLD,statinfo)
    else if(me==Np-1)then !on envoie Xp à me=Np-2
       call MPI_SEND(X0(1:Nx),Nx,MPI_REAL8,me-1,tag,MPI_COMM_WORLD,statinfo)
    end if

    if((me/=0) .and. (me/=Np-1))then !on reçoit Xp de me-1 et Xn de me+1
       call MPI_RECV(Xn,Nx,MPI_REAL8,me+1,tag,MPI_COMM_WORLD,status,statinfo)
       call MPI_RECV(Xp,Nx,MPI_REAL8,me-1,tag,MPI_COMM_WORLD,status,statinfo)
    else if(me==0)then !on reçoit Xn de me=2
       call MPI_RECV(Xn,Nx,MPI_REAL8,me+1,tag,MPI_COMM_WORLD,status,statinfo)
       Xp=0
    else if(me==Np-1)then !on reçoit Xp de me=Nx-2
       call MPI_RECV(Xp,Nx,MPI_REAL8,me-1,tag,MPI_COMM_WORLD,status,statinfo)
       Xn=0
    end if

    call MatMulAParallele(R,X0,Xp,Xn,B,Cx,Cy,Nx,Ny,i1,iN) !R=A*X0
    R=R-rhs
    d=R
    X=X0
    k=0
    !Calculer Norme^2 de R 
    call DotProdParallele(e,R,R,Nx,Ny,i1,iN)
    call MPI_ALLREDUCE(e,NormeR,1,MPi_REAL8,MPi_SUM,MPi_COMM_WORLD,statinfo)



!!!!! Début de la boucle
    Do While ((k<=N) .and. (sqrt(NormeR)>eps/Np))
       k=k+1
       !!W=matmul(A,d)  
       !calcule Xp et Xn pour max(iN-i1+1)> Nx

       if((me/=0) .and. (me/=Np-1))then
          call MPI_SEND(d(iN-i1+1-Nx+1:iN-i1+1),Nx,MPI_REAL8,me+1,tag,MPI_COMM_WORLD,statinfo)
          call MPI_SEND(d(1:Nx),Nx,MPI_REAL8,me-1,tag,MPI_COMM_WORLD,statinfo)
       else if(me==0)then
          call MPI_SEND(d(iN-i1+1-Nx+1:iN-i1+1),Nx,MPI_REAL8,me+1,tag,MPI_COMM_WORLD,statinfo)
       else if(me==Np-1)then
          call MPI_SEND(d(1:Nx),Nx,MPI_REAL8,me-1,tag,MPI_COMM_WORLD,statinfo)
       end if

       if((me/=0) .and. (me/=Np-1))then
          call MPI_RECV(Xn,Nx,MPI_REAL8,me+1,tag,MPI_COMM_WORLD,status,statinfo)
          call MPI_RECV(Xp,Nx,MPI_REAL8,me-1,tag,MPI_COMM_WORLD,status,statinfo)
       else if(me==0)then
          call MPI_RECV(Xn,Nx,MPI_REAL8,me+1,tag,MPI_COMM_WORLD,status,statinfo)
          Xp=0
       else if(me==Np-1)then
          call MPI_RECV(Xp,Nx,MPI_REAL8,me-1,tag,MPI_COMM_WORLD,status,statinfo)
          Xn=0
       end if

       call MatMulAParallele(W,d,Xp,Xn,B,Cx,Cy,Nx,Ny,i1,iN)

       !alpha=dot_product(d,R)/dot_product(d,W)
       call DotProdParallele(e,d,R,Nx,Ny,i1,iN) !e=produit de d et R
       call DotProdParallele(f,d,W,Nx,Ny,i1,iN) !f=produit de d et W
       call MPi_REDUCE(e,Ep,1,MPi_REAL8,MPi_SUM,0,MPi_COMM_WORLD,statinfo)
       call MPi_REDUCE(f,Fp,1,MPi_REAL8,MPi_SUM,0,MPi_COMM_WORLD,statinfo)
       if(me==0)then
          alpha=E/F
          call MPI_BCAST(alpha,1,MPI_REAL8,0,MPI_COMM_WORLD,statinfo)
       end if

       !X=X-alpha*d
       X=X-alpha*d

       !R_next=R-alpha*W
       R_next=R-alpha*W

       !beta=Norme(R_next)**2/Norme(R)**2
       call DotProdParallele(e,R_next,R_next,Nx,Ny,i1,iN)
       call DotProdParallele(f,R,R,Nx,Ny,i1,iN)
       call MPi_REDUCE(e,Ep,1,MPi_REAL8,MPi_SUM,0,MPi_COMM_WORLD,statinfo)
       call MPi_REDUCE(f,Fp,1,MPi_REAL8,MPi_SUM,0,MPi_COMM_WORLD,statinfo)     
       if(me==0)then
          alpha=E/F
          call MPI_BCAST(beta,1,MPI_REAL8,0,MPI_COMM_WORLD,statinfo)
       end if

       !R=R_next
       R=R_next

       !d=R+beta*d  
       d=beta*d+R

       !!Norme R
       call DotProdParallele(e,R,R,Nx,Ny,i1,iN)
       call MPi_ALLREDUCE(e,NormeR,1,MPi_REAL8,MPi_SUM,MPi_COMM_WORLD,statinfo)
       print*,k,NormeR
    end do
!!!rassemblement de X dans le main
    print*,X,me
  end subroutine GradConj
  
    subroutine MatMulAParallele(AX,X,Xp,Xn,B,Cx,Cy,Nx,Ny,i1,iN)
    integer,intent(in)::i1,iN,Nx,Ny
    real*8,dimension(iN-i1+1),intent(in)::X
    real*8,dimension(Nx),intent(in)::Xp,Xn
    real*8,dimension(iN-i1+1),intent(out)::AX
    real*8,intent(in)::B,Cx,Cy
    integer::j
    do j=1,iN-i1+1

       !Premier Bloc
       if((j+i1-2)/Nx==0)then
          if(mod(j+i1-1,Nx)==1)then !Première ligne

             AX(j)=B*X(j)
             if(j+1>iN-i1+1)then
                AX(j)=AX(j)+Cx*Xn(1)
             else
                AX(j)=AX(j)+Cx*X(j+1)
             end if
             if(j+Nx>iN-i1+1)then
                AX(j)=AX(j)+Cy*Xn(j+Nx-(iN-i1+1))
             else
                AX(j)=AX(j)+Cy*X(j+Nx)
             end if

          elseif(mod(j+i1-1,Nx)==0)then !Dernière ligne

             AX(j)=B*X(j)
             if(j-1<1)then
                AX(j)=AX(j)+Cx*Xp(Nx)
             else
                AX(j)=AX(j)+Cx*X(j-1)
             end if
             if(j+Nx>iN-i1+1)then
                AX(j)=AX(j)+Cy*Xn(j+Nx-iN)
             else
                AX(j)=AX(j)+Cy*X(j+Nx)
             end if

          else

             AX(j)=B*X(j)
             if(j-1<1)then
                AX(j)=AX(j)+Cx*Xp(Nx)
             else
                AX(j)=AX(j)+Cx*X(j-1)
             end if
             if(j+1>iN-i1+1)then
                AX(j)=AX(j)+Cx*Xn(1)
             else
                AX(j)=AX(j)+Cx*X(j+1)
             end if
             if(j+Nx>iN-i1+1)then
                AX(j)=AX(j)+Cy*Xn(j+Nx-iN)
             else
                AX(j)=AX(j)+Cy*X(j+Nx)
             end if

          end if

          !Dernier Bloc
       elseif((j+i1-2)/Nx==Ny-1)then

          if(mod(j+i1-1,Nx)==1)then !Première ligne

             AX(j)=B*X(j)
             if(j-Nx<1)then
                AX(j)=AX(j)+Cy*Xp(j)
             else
                AX(j)=AX(j)+Cy*X(j-Nx)
             end if
             if(j+1>iN-i1+1)then
                AX(j)=AX(j)+Cx*Xn(1)
             else
                AX(j)=AX(j)+Cx*X(j+1)
             end if

          elseif(mod(j+i1-1,Nx)==0)then !Dernière ligne

             AX(j)=B*X(j)
             if(j-Nx<1)then
                AX(j)=AX(j)+Cy*Xp(j)
             else
                AX(j)=AX(j)+Cy*X(j-Nx)
             end if
             if(j-1<1)then
                AX(j)=AX(j)+Cx*Xp(Nx)
             else
                AX(j)=AX(j)+Cx*X(j-1)
             end if

          else

             AX(j)=B*X(j)
             if(j-Nx<1)then
                AX(j)=AX(j)+Cy*Xp(j)
             else
                AX(j)=AX(j)+Cy*X(j-Nx)
             end if
             if(j-1<1)then
                AX(j)=AX(j)+Cx*Xp(Nx)
             else
                AX(j)=AX(j)+Cx*X(j-1)
             end if
             if(j+1>iN-i1+1)then
                AX(j)=AX(j)+Cx*Xn(1)
             else
                AX(j)=AX(j)+Cx*X(j+1)
             end if

          end if

          !Bloc du Milieu
       else

          if(mod(j+i1-1,Nx)==1)then ! Première ligne

             AX(j)=B*X(j)
             if(j-Nx<1)then
                AX(j)=AX(j)+Cy*Xp(j)
             else
                AX(j)=AX(j)+Cy*X(j-Nx)
             end if
             if(j+1>iN-i1+1)then
                AX(j)=AX(j)+Cx*Xn(1)
             else
                AX(j)=AX(j)+Cx*X(j+1)
             end if
             if(j+Nx>iN-i1+1)then
                AX(j)=AX(j)+Cy*Xn(j+Nx-(iN-i1+1))
             else
                AX(j)=AX(j)+Cy*X(j+Nx)
             end if

          elseif(mod(j+i1-1,Nx)==0)then !Dernière ligne

             AX(j)=B*X(j)
             if(j-Nx<1)then
                AX(j)=AX(j)+Cy*Xp(j)
             else
                AX(j)=AX(j)+Cy*X(j-Nx)
             end if
             if(j-1<1)then
                AX(j)=AX(j)+Cx*Xp(Nx)
             else
                AX(j)=AX(j)+Cx*X(j-1)
             end if
             if(j+Nx>(iN-i1+1))then
                AX(j)=AX(j)+Cy*Xn(j+Nx-(iN-i1+1))
             else
                AX(j)=AX(j)+Cy*X(j+Nx)
             end if

          else

             AX(j)=B*X(j)
             if(j-Nx<1)then
                AX(j)=AX(j)+Cy*Xp(j)
             else
                AX(j)=AX(j)+Cy*X(j-Nx)
             end if
             if(j-1<1)then
                AX(j)=AX(j)+Cx*Xp(Nx)
             else
                AX(j)=AX(j)+Cx*X(j-1)
             end if
             if(j+1>iN-i1+1)then
                AX(j)=AX(j)+Cx*Xn(1)
             else
                AX(j)=AX(j)+Cx*X(j+1)
             end if
             if(j+Nx>(iN-i1+1))then
                AX(j)=AX(j)+Cy*Xn(j+Nx-(iN-i1+1))
             else
                AX(j)=AX(j)+Cy*X(j+Nx)
             end if

          end if

       end if

    end do

  end subroutine MatMulAParallele


!! DotProdParallele censé OK
  subroutine DotProdParallele(Resultat,A,B,Nx,Ny,i1,iN)
    implicit none
    integer, intent (in)::Nx,Ny,i1,iN
    real*8, dimension(:),intent(in)::A,B
    real*8,intent(out)::Resultat
    integer::j
    Resultat=0
    Do j=1,iN-i1+1
       Resultat=A(j)*B(j)
    End do
 
  end subroutine DotProdParallele

  subroutine CalculRHS(Cx,Cy,Nx,Ny,h,g,RHS,me,i1,iN)
    implicit none
    real*8,intent(in)::Cx,Cy
    real*8,dimension(:),intent(in)::h,g
    real*8,dimension(iN-i1+1),intent(inout)::RHS
    integer,intent(in)::i1,iN,me,Nx,Ny
    integer::i   
    do i=1, iN-i1+1
       if (i<=Nx) then
          RHS(i) = RHS(i) - g(i)*Cy
       elseif (i>(Nx*Ny-Nx)) then
          RHS(i) = RHS(i) - g(i-(Nx*Ny-2*Nx))*Cy
       end if

       if (mod(i,Nx) == 1) then
          RHS(i) = RHS(i) - h(1+(i/Nx))*Cx
       elseif (mod(i, Nx) == 0) then
          RHS(i) = RHS(i) - h(Ny +(i/Nx))*Cx
       end if

    enddo
  end subroutine CalculRHS



end module ModuleP



