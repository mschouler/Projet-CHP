module proj_1
implicit none


contains
  subroutine grad_conj(n,A,X,b)
    implicit none
    integer,intent(in)::n
    real*8,intent(in),dimension(n,n)::A
    real*8,intent(inout),dimension(n)::X
    real*8,intent(in),dimension(n)::b
    real*8,dimension(n)::w,d,R,R1
    integer::l=0
    integer,parameter::niter=10000
    real*8,parameter::eps=1e-5
    real*8::beta,alpha
   
    !initialisation 
    R=matmul(A,X)-b
    d=R
   
    !boucle 
    do while(norme(R,n)>eps .and. l<=niter)
       print*,'iteration',l
       w=matmul(A,d)
       print*,w
       alpha=scal(d,R,n)/scal(d,w,n)
       X=X-alpha*d
       R1=R-alpha*w
       beta=norme(R1,n)**2/norme(R,n)**2
       d=R1+beta*d
       R=R1
       l=l+1
    end do
    
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

          do i=1,Nx
             B(i) = alpha*X(i) + beta*X(i+1) + gama*X(i+Nx)
          enddo
          
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

          do i=n-Nx, n
             
             B(i) = gama*X(i-Nx) + beta*X(i-1) + alpha*X(i)
          enddo

        end subroutine mult
      
end module






