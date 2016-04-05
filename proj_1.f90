module proj_1
implicit none


contains
  subroutine grad_conj(Nx, Ny, alpha, beta, gama, X, b)
    implicit none
    integer,intent(in)::Nx, Ny
    real*8,intent(inout),dimension(Nx*Ny)::X
    real*8,intent(in),dimension(Nx*Ny)::b
    real*8,intent(in)::alpha,beta,gama
    real*8,dimension(Nx*Ny)::w,d,R,R1
    integer::l=0, n
    integer,parameter::niter=10000
    real*8,parameter::eps=1e-5
    real*8::beta_,alpha_
    n=Nx*Ny
   
    !initialisation 
    call mult(alpha, beta, gama, X, R, Nx, Ny)
    R = R - b
    d=R
   
    !boucle 
    do while(norme(R,n)>eps .and. l<=niter)
       call mult(alpha, beta, gama, d, w, Nx, Ny)
       alpha_=scal(d,R,n)/scal(d,w,n)
       X=X-alpha_*d
       R1=R-alpha_*w
       beta_=norme(R1,n)**2/norme(R,n)**2
       d=R1+beta_*d
       R=R1
       l=l+1
    end do
    print*,'iteration',l
    end subroutine

    function norme(X,n)
      implicit none
      real*8,dimension(n)::X
      real*8::norme
      integer::i,n  
      norme=0.
      do i=1,n
         norme=norme+X(i)**2
      end do
      norme=sqrt(norme)

      end function

      function scal(X,Y,n)
        implicit none
        real*8,dimension(n)::X,Y
        real*8::scal
        integer::i,n
        scal=0.
        do i=1,n
           scal=scal+X(i)*Y(i)
        end do
        end function

        subroutine mult(alpha, beta, gama, X, B, Nx,Ny)
          implicit none
          real*8, intent(in) :: alpha, beta, gama
          integer, intent(in) :: Nx,Ny
          integer :: i, k, l,n
          real*8, dimension(Nx*Ny), intent(in) :: X
          real*8, dimension(Nx*Ny), intent(out) :: B
          n=Nx*Ny

        B(1)=alpha*X(1) + beta*X(2) + gama*X(1+Nx)

          do i=2,Nx-1
             B(i) = alpha*X(i) + beta*X(i-1)+beta*X(i+1) + gama*X(i+Nx)
          enddo

        B(Nx)=alpha*X(Nx) + beta*X(Nx-1) + gama*X(2*Nx)

          do i=Nx+1, n-(Nx+1)
             if (MOD(i, Nx) == 1) then
                k=0
                l=1
                
             else if (MOD(i, Nx) == 0) then
                k=1
                l=0
             else
                k=1
                l=1
             endif
             
             B(i) = gama*(X(i-Nx) + X(i+Nx)) + beta*(k*X(i-1) + l*X(i+1)) + gama*X(i)
          enddo

        B(n-Nx)=gama*(X(n-2*Nx) + X(n)) + alpha*X(n-Nx) + beta*X(n-Nx+1)

          do i=n-Nx+1, n-1
             
             B(i) = gama*X(i-Nx) + beta*X(i-1) + beta*X(i+1) + alpha*X(i)
          enddo

        B(n)=gama*X(n-Nx) + beta*X(n-1) + alpha*X(n)

        end subroutine mult
      
end module






