program projet
use proj_1
implicit none
real :: Lx, Ly, dx, dy, dt, D
real*8::alpha, beta, gama
integer :: Nx, Ny, i, j
real*8, dimension(:), allocatable :: U, Test, Uexacte, F, X, Y;
real*8,dimension(4,4)::A
real*8,dimension(4)::b
real*8,dimension(4)::X2
real*8,dimension(:),allocatable::X1, C


Lx=1
Ly=1

write(*,*)'combien de mailles selon x ? '
read*,Nx
write(*,*)'combien de mailles selon y ? '
read*,Ny

dx = Lx/Nx
dy = Ly/Ny

allocate(U(1:Nx*Ny))
allocate(Test(1:Nx*Ny))
allocate(Uexacte(1:Nx*Ny))
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
    Uexacte((i-1)*Nx+j)=X(j)*(1-X(j))*Y(i)*(1-Y(i))
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
    write(10,*) " "
end do
close(10)

U=1.
D=1.
dt=1.
alpha=1./dt + 2.*D/dx**2 + 2.*D/dy**2
beta=-D/dx**2.
gama=-D/dy**2.
call grad_conj(Nx, Ny, alpha, beta, gama, U, F)  
!!$call mult(1.D0, 1.D0, 1.D0, U, Test, Nx, Ny)
!!$print*, Test

open( unit=12, &
file = "sol.txt", &
action = "write", &
status = "unknown")
do i=1,Ny
    do j=1,Nx
    write(12,*) Y(i),X(j),U((i-1)*Nx+j)
    end do
    write(12,*) " "
end do
close(12)

!!$open( unit=13, &
!!$file = "solexacte.txt", &
!!$action = "write", &
!!$status = "unknown")
!!$do i=1,Ny
!!$    do j=1,Nx       
!!$    write(13,*) Y(i),X(j),Uexacte((i-1)*Nx+j)
!!$    end do
!!$    write(13,*) " "
!!$end do
!!$close(13)


end program
