program proj_1
implicit none





contains
  subroutine grad_conj(A,X,b,n)
    implicit none
    real*8,intent(in),dimension(n,n)::A
    real*8,dimension(n)::X,b
    real*8,dimension(n)::w,alpha,d,R,R1
    integer::i,j
    integer,parameter::niter=10000
    real*8,parameter::eps=1e-5
    
    !initialisation 
    R=matmul(A,X)-b
    d=R
    
    do while(norme(R)>eps .and. l<=niter)
       w=matmul(A,d)
       alpha=matmul(d,R)/matmul(d,w)
       x=x-alpha*d
       R1=R-alpha*w
       beta=norme(R1,n)**2/norme(R,n)**2
       d=R1+beta*d
       R=R1
    end do
    
    end subroutine

    function norme(X,n)
      real*8,dimension(n)::X
      real
    
    








