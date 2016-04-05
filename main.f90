program projet
use proj_1
implicit none
real :: Lx, Ly, dx, dy
integer :: Nx, Ny, i, j
real, dimension(:), allocatable :: U, F, X, Y;
real*8,dimension(4,4)::A
real*8,dimension(4)::b
real*8,dimension(4)::X2
real*8,dimension(:),allocatable::X1,C
!X2=0.
!print*,X2
!b=1.
!A(1,:)=(/10.,-1.,-1.,-1./)
!A(2,:)=(/-1.,10.,-1.,-1./)
!A(3,:)=(/-1.,-1.,10.,-1./)
!A(4,:)=(/-1.,-1.,-1.,10./)

!call grad_conj(4,A,X2,b)
!print*,matmul(A,X2)

!coucou jafar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! TEST DE LA BOUCLE MULT !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Lx=1
Ly=1

write(*,*)'combien de mailles selon x ? '
read*,Nx
write(*,*)'combien de mailles selon y ? '
read*,Ny

dx = Lx/Nx
dy = Ly/Ny

allocate(U(1:Nx*Ny))
allocate(F(1:Nx*Ny))
allocate(X1(Nx*Ny))
allocate(C(Nx*Ny))

allocate(X(1:Nx))
allocate(Y(1:Ny))

do i=1, Nx
    X(i) = i*dx
enddo

do i=1,Ny
    Y(i) = i*dy
enddo

do i=1,Ny
 do j=1,Nx
    F((i-1)*Nx+j) = Y(i) - Y(i)**2 + X(j) - X(j)**2
 end do
end do

open( unit=10, &
file = "f.txt", &
action = "write", &
status = "unknown")
do i=1,Ny
    do j=1,Nx
    write(10,*) Y(i),X(j),F((i-1)*Nx+j)
    end do
end do

!X1=1.
!call mult(1.0D0,1.0D0,1.0D0,X1,C,Nx,Ny)

!write(*,*) 'Le produit vaut '
!do i=1,NX*Ny
    !print*,C(i)
!end do


end program
