program projet
use proj_1
implicit none
real*8:: Lx, Ly, dx, dy, dt, D, Tf, t
real*8::alpha, beta, gama
integer :: Nx, Ny, i, j, Niter, nb_iter, mode, ini
real*8, dimension(:), allocatable :: U, U_0, Test, Uexacte, F, X, Y;
real*8,dimension(4,4)::A
real*8,dimension(4)::b
real*8,dimension(4)::X2
real*8,dimension(:),allocatable::X1, C

! Initialisation des paramètres
Niter=10000
call lec(Nx,Ny,Lx,Ly,D)
!Lx=10.
!Ly=10.

!write(*,*) 'combien de mailles selon x et y ? '
!read*,Nx, Ny
write(*,*) 'Choisir la modelisation en 1_stationnaire ou 2_instationnaire'
read*,mode

dx = Lx/Nx
dy = Ly/Ny

allocate(U(1:Nx*Ny))
allocate(U_0(1:Nx*Ny))
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

!ecriture de la solution exacte
!open( unit=13, &
!file = "solexacte.txt", &
!action = "write", &
!status = "unknown")
!do i=1,Ny
!    do j=1,Nx
!    write(13,*) Y(i),X(j),Uexacte((i-1)*Nx+j)
!    end do
!    write(13,*) " "
!end do
!close(13)


if (mode==1) then
    !cas stationnaire
    write(*,*) 'choisir la condition initiale 1 ou 2'
    read*, ini
    U=1.
    D=1.
    dt=1.
    alpha=1./dt + 2.*D/dx**2 + 2.*D/dy**2
    beta=-D/dx**2.
    gama=-D/dy**2.
    F=ini_f(Nx,Ny,X,Y,t,mode,ini)
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
    call grad_conj(Nx, Ny, alpha, beta, gama, U, F)
    !ecriture de la solution
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

end if
if (mode==2) then
    !cas instationnaire
    write(*,*) 'Choisir le temps final Tf'
    read*,Tf

    dt=Tf/Niter
    U_0=0.
    D=1.
    dt=1.
    alpha=1./dt + 2.*D/dx**2 + 2.*D/dy**2
    beta=-D/dx**2.
    gama=-D/dy**2.
    nb_iter=1

    do while(nb_iter<Niter)

    t=dt*nb_iter !on declare t
    F=ini_f(Nx,Ny,X,Y,t,mode,0) !on initialise le second membre

    do i=1,Nx*Ny
        F(i)=F(i)+U_0(i)/dt
    end do

    call grad_conj(Nx, Ny, alpha, beta, gama, U, F)
    U_0=U

    nb_iter=nb_iter+1
    end do
    !ecriture de la solution
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

end if

!open( unit=10, &
!file = "f.txt", &
!action = "write", &
!status = "unknown")
!do i=1,Ny
!    do j=1,Nx
!    write(10,*) Y(i),X(j),F((i-1)*Nx+j)
!    end do
!    write(10,*) " "
!end do
!close(10)

end program
