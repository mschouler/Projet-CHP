program proj_1
implicit none
integer :: i
real*8,dimension(4,4)::A
real*8,dimension(4)::b
real*8,dimension(4)::X
X=0.
print*,X
b=1.
A(1,:)=(/10.,-1.,-1.,-1./)
A(2,:)=(/-1.,10.,-1.,-1./)
A(3,:)=(/-1.,-1.,10.,-1./)
A(4,:)=(/-1.,-1.,-1.,10./)

do i=1,4
   write(*,*) A(i,1:4)
enddo
call grad_conj(4,A,X,b)
print*,X

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
        
      
end program






