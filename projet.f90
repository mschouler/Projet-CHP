program projet

implicit none
real :: Lx, Ly, dx, dy
integer :: Nx, Ny, i, j
real, dimension(:), allocatable :: U, F, X, Y;

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

allocate(X(1:Nx))
allocate(Y(1:Ny))

do i=1, Nx
    X(i) = i*dx
enddo

do i=1,Ny
    Y(i) = i*dy
enddo

do i=1, Nx*Ny
    F(i) = Y(i) - Y(i)**2 + X(i) - X(i)**2
enddo

write(*,*)F

contains


subroutine mult(alpha, beta, gama, X, B, n)
implicit none
real, intent(in) :: alpha, beta, gama
integer, intent(in) :: n
integer :: i, k, l
real, dimension(n), intent(in) :: X
real, dimension(n), intent(out) :: B


do i=1,5
    B(i) = alpha*X(i) + beta*X(i+1) + gama*X(i+5)
enddo

do i=5, Nx*Ny-6
    if (MOD(i, 5) == 1) then
        k=0
        l=1

    else if (MOD(i, 5) == 0) then
        k=1
        l=0
    else
        k=1
        l=1
    endif

B(i) = gama*(X(i-5) + X(i+5)) + beta*(k*X(i-1) + l*X(i+1)) + gama*X(i)
enddo

do i=Nx*Ny-5, Nx*Ny

        B(i) = gama*X(i-5) + beta*X(i-1) + alpha*X(i)
enddo

end subroutine





end program