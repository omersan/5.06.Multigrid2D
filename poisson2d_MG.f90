!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> V-cycle Multigrid method for solving 2D Poisson equation
!     d2u/dx2 + d2u/dy2 = f(x,y)
!     Drichlet b.c.

!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012)
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Oct. 29, 2015
!-----------------------------------------------------------------------------!

program poisson2d
implicit none
integer::i,j,nx,ny,isolver
real*8,dimension(:,:),allocatable ::u,f,ue,e
real*8,dimension(:),allocatable ::x,y
real*8 ::dx,dy,tol,rms,x0,xL,y0,yL

!Domain
x0 =-1.0d0 !left
xL = 1.0d0 !right

y0 =-1.0d0 !bottom
yL = 1.0d0 !up

!number of points
nx = 128 !number of grid points in x (i.e., should be power of 2)
ny = nx  !number of grid points in y

!grid spacing (spatial)
dx = (xL-x0)/dfloat(nx)
dy = (yL-y0)/dfloat(ny)

!spatial coordinates 
allocate(x(0:nx))
do i=0,nx
x(i) = x0 + dfloat(i)*dx
end do

allocate(y(0:ny))
do j=0,ny
y(j) = y0 + dfloat(j)*dy
end do





!Tolerance
tol= 1.0d-4


allocate(u(0:nx,0:ny))
allocate(f(0:nx,0:ny))
allocate(e(0:nx,0:ny))
allocate(ue(0:nx,0:ny))

!---------------------------------------------!
!Exact solution (test case from Moin's textbook):
!---------------------------------------------!
do j=0,ny
do i=0,nx
f(i,j) =-2.0d0*(2.0d0-x(i)*x(i)-y(j)*y(j))
ue(i,j)= (x(i)*x(i)-1.0d0)*(y(j)*y(j)-1.0d0)
end do
end do


!Numerical solution:
do i=0,nx
do j=0,ny
u(i,j)=0.0d0
end do
end do

!Boundary conditions has to satisfy exact solution
do i=0,nx
u(i,0)  = ue(i,0)	
u(i,ny) = ue(i,ny)					  					  	
end do

do j=0,ny
u(0,j)  = ue(0,j)		
u(nx,j) = ue(nx,j)						  	
end do


open(19,file='output.txt')

!----------------------!
!Solver:
!----------------------!
isolver = 1
if (isolver.eq.0) then !Gauss-Seidel scheme
	call GS(nx,ny,dx,dy,f,u,tol)
else
	call MG5(nx,ny,dx,dy,f,u,tol)
end if


!----------------------!
!Error analysis:
!----------------------!
do i=0,nx
do j=0,ny
e(i,j) = dabs(u(i,j)-ue(i,j))
end do 
end do


!L-2 Norm:
call l2norm(nx,ny,e,rms)

write(*,*)"L2-norm =",rms
write(19,*)"L2-norm =",rms

!maximum norm
write(*,*)"Max-norm =",maxval(e)
write(19,*)"Max-norm =",maxval(e)
close(19)

!Plot field
open(10,file='field.plt')
write(10,*) 'variables ="x","y","f","u","ue"'
write(10,*)'zone f=point i=',nx+1,',j=',ny+1
do j=0,ny
do i=0,nx
write(10,*) x(i),y(j),f(i,j),u(i,j),ue(i,j)
end do
end do
close(10)


end




!---------------------------------------------------------------------------!
!Relaxation formula for Poisson equation
!Uses GS relaxation
!Works for Drichlet boundary conditions (Boundary points never updated)
!---------------------------------------------------------------------------!
SUBROUTINE relax(nx,ny,dx,dy,f,u)
implicit none
integer::nx,ny
real*8 ::dx,dy
real*8, dimension(0:nx,0:ny)::u,f
real*8 ::a
integer::i,j

a = -2.0d0/(dx*dx) - 2.0d0/(dy*dy)
  
do i=1,nx-1
do j=1,ny-1
u(i,j) = (1.0d0/a)*(f(i,j) &
                   - (u(i+1,j)+u(i-1,j))/(dx*dx) &
                   - (u(i,j+1)+u(i,j-1))/(dy*dy) )
end do
end do

return
end



!---------------------------------------------------------------------------!
!Residual formula for Poisson equation
!Works for Drichlet boundary conditions (Boundary points never updated)
!---------------------------------------------------------------------------!
SUBROUTINE resid(nx,ny,dx,dy,f,u,r)
implicit none
integer::nx,ny
real*8 ::dx,dy
real*8, dimension(0:nx,0:ny)::u,f,r
integer::i,j

do i=1,nx-1
do j=1,ny-1
r(i,j) = f(i,j) - (u(i+1,j) - 2.0d0*u(i,j) + u(i-1,j))/(dx*dx) &
		        - (u(i,j+1) - 2.0d0*u(i,j) + u(i,j-1))/(dy*dy) 
end do
end do

!Boundary conditions for residuals
do i=0,nx
r(i,0)  = 0.0d0	
r(i,ny) = 0.0d0					  					  	
end do

do j=0,ny
r(0,j)  = 0.0d0		
r(nx,j) = 0.0d0						  	
end do

return
end

!---------------------------------------------------------------------------!
!Compute L2-norm for an array
!---------------------------------------------------------------------------!
SUBROUTINE l2norm(nx,ny,r,rms)
implicit none
integer::Nx,Ny
real*8, dimension(0:nx,0:ny)::r
integer::i,j
real*8 ::rms

rms=0.0d0
do i=1,nx-1
do j=1,ny-1
rms = rms + r(i,j)*r(i,j)
end do 
end do
rms= dsqrt(rms/dfloat((nx-1)*(ny-1)))

return
end


!---------------------------------------------------------------------------!
!Restriction operators
!---------------------------------------------------------------------------!
SUBROUTINE rest(nxf,nyf,nxh,nyh,r,f)
implicit none
integer::nxf,nyf,nxh,nyh
real*8, dimension(0:nxf,0:nyf)::r	!on higher grid
real*8, dimension(0:nxh,0:nyh)::f	!on lower grid
integer::i,j
integer::ireo

ireo = 3

if (ireo.eq.1) then !simply injection

do i=1,nxh-1
do j=1,nyh-1
f(i,j) = r(2*i,2*j) 							  	
end do
end do

else if (ireo.eq.2) then !half-weight

do i=1,nxh-1
do j=1,nyh-1
f(i,j) = 1.0d0/8.0d0*( 4.0d0*r(2*i,2*j) &
	     + 1.0d0*(r(2*i+1,2*j)+r(2*i-1,2*j)+r(2*i,2*j+1)+r(2*i,2*j-1)) )							  	
end do
end do


else !full-weight (trapezoidal)

do i=1,nxh-1
do j=1,nyh-1
f(i,j) = 1.0d0/16.0d0*( 4.0d0*r(2*i,2*j) &
	     + 2.0d0*(r(2*i+1,2*j)+r(2*i-1,2*j)+r(2*i,2*j+1)+r(2*i,2*j-1)) &
	     + 1.0d0*(r(2*i+1,2*j+1)+r(2*i-1,2*j-1)+r(2*i-1,2*j+1)+r(2*i+1,2*j-1)))							  	
end do
end do

end if


!update boundaries
do i=0,nxh
f(i,0)   = r(2*i,0) 	
f(i,nyh) = r(2*i,nyf) 					  					  	
end do

do j=0,nyh
f(0,j)   = r(0,2*j)		
f(nxh,j) = r(nxf,2*j)						  	
end do

  
return
end


!---------------------------------------------------------------------------!
!Prolongation operator
!bilinear interpolation
!---------------------------------------------------------------------------!
SUBROUTINE prol(nxh,nyh,nxf,nyf,u,p)
implicit none
integer::nxf,nyf,nxh,nyh
real*8, dimension(0:nxf,0:nyf)::p	!on higher grid
real*8, dimension(0:nxh,0:nyh)::u	!on lower grid
integer::i,j


do i=0,nxh-1
do j=0,nyh-1
p(2*i,2*j)    = u(i,j)
p(2*i+1,2*j)  = 1.0d0/2.0d0*(u(i,j)+u(i+1,j))
p(2*i,2*j+1)  = 1.0d0/2.0d0*(u(i,j)+u(i,j+1))
p(2*i+1,2*j+1)= 1.0d0/4.0d0*(u(i,j)+u(i,j+1)+u(i+1,j)+u(i+1,j+1))
end do
end do

do j=0,nyh
p(nxf,2*j)    = u(nxh,j)
end do
do i=0,nxh
p(2*i,nyf)    = u(i,nyh)
end do

return
end



!---------------------------------------------------------------------------!
!Gauss Seidel scheme (1 level)
!---------------------------------------------------------------------------!
SUBROUTINE GS(nx,ny,dx,dy,f,u,tol)
implicit none
integer::nx,ny
real*8 ::dx,dy
real*8 ::tol
real*8,dimension(0:nx,0:ny)::u,f
real*8,dimension(:,:),allocatable:: r
real*8 ::rms0,rms
integer::k,ke,wl,nI

nI = 100000  !maximum number of iteration

allocate(r(0:nx,0:ny))

!Compute initial resitual:
call resid(nx,ny,dx,dy,f,u,r)

!and its l2 norm:
call l2norm(nx,ny,r,rms0)

open(66,file='residual.plt')
write(66,*) 'variables ="k","rms","rms/rms0"'

do k=1,nI

	call relax(nx,ny,dx,dy,f,u)
	call resid(nx,ny,dx,dy,f,u,r)
    
	! Check for convergence on smallest grid	
	call l2norm(nx,ny,r,rms)
	if (rms/rms0.le.tol) goto 10

    ! Write residual history
	write(66,*) k,rms,rms/rms0
    write(*,*) k,rms,rms/rms0     
end do

10 continue
close(66)

deallocate(r)

ke=k

!work load (total number of operations)
wl = ke*(nx*ny)

write(19,*)"outer number of iteration = ",ke
write(19,*)"normalized workload       = ",dfloat(wl)/dfloat(nx*ny)
write(*,*)"outer number of iteration = ",ke
write(*,*)"normalized workload       = ",dfloat(wl)/dfloat(nx*ny)

return  
end  




!---------------------------------------------------------------------------!
!Multigrid scheme (5 level)
!Full-weighting is used as restriction operator
!Bilinear interpolation procedure is used as prolongation operator
!---------------------------------------------------------------------------!
SUBROUTINE MG5(nx,ny,dx,dy,f,u,tol)
implicit none
integer::nx,ny
real*8 ::dx,dy
real*8 ::tol
integer::nI,v1,v2,v3
real*8,dimension(0:nx,0:ny)::u,f
real*8,dimension(:,:),allocatable:: r,r2,r3,r4,r5
real*8,dimension(:,:),allocatable:: p,p2,p3,p4
real*8,dimension(:,:),allocatable:: u2,u3,u4,u5
real*8,dimension(:,:),allocatable:: f2,f3,f4,f5
real*8 ::dx2,dy2,dx3,dy3,dx4,dy4,dx5,dy5,rms0,rms,rmsc
integer::i,j,k,ke,me,wl,m,nx2,ny2,nx3,ny3,nx4,ny4,nx5,ny5

nI = 100000 !maximum number of outer iteration
v1 = 2   	!number of relaxation for restriction in V-cycle
v2 = 2   	!number of relaxation for prolongation in V-cycle
v3 = 100 	!number of relaxation at coarsest level


dx2=dx*2.0d0
dy2=dy*2.0d0

dx3=dx*4.0d0
dy3=dy*4.0d0

dx4=dx*8.0d0
dy4=dy*8.0d0

dx5=dx*16.0d0
dy5=dy*16.0d0

nx2=nx/2
ny2=ny/2

nx3=nx/4
ny3=ny/4

nx4=nx/8
ny4=ny/8

nx5=nx/16
ny5=ny/16

me = 0

if (nx5.lt.2.or.ny5.lt.2) then
write(*,*)"5 level is high for this grid.."
stop
end if

allocate(r (0:nx ,0:ny))
allocate(p (0:nx ,0:ny))

allocate(u2(0:nx2,0:ny2))
allocate(f2(0:nx2,0:ny2))
allocate(r2(0:nx2,0:ny2))
allocate(p2(0:nx2,0:ny2))

allocate(u3(0:nx3,0:ny3))
allocate(f3(0:nx3,0:ny3))
allocate(r3(0:nx3,0:ny3))
allocate(p3(0:nx3,0:ny3))

allocate(u4(0:nx4,0:ny4))
allocate(f4(0:nx4,0:ny4))
allocate(r4(0:nx4,0:ny4))
allocate(p4(0:nx4,0:ny4))

allocate(u5(0:nx5,0:ny5))
allocate(f5(0:nx5,0:ny5))
allocate(r5(0:nx5,0:ny5))

!Compute initial resitual:
call resid(nx,ny,dx,dy,f,u,r)
!and its l2 norm:
call l2norm(nx,ny,r,rms0)

open(66,file='residual.plt')
write(66,*) 'variables ="k","rms","rms/rms0"'


do k=1,nI

!1.Relax v1 times
do m=1,v1
call relax(nx,ny,dx,dy,f,u)			
end do

! Compute residual
call resid(nx,ny,dx,dy,f,u,r)

! Check for convergence on finest grid	
call l2norm(nx,ny,r,rms)
write(66,*) k,rms,rms/rms0  
write(*,*) k,rms,rms/rms0 
if (rms/rms0.le.tol) goto 10

!1r.Restriction	
call rest(nx,ny,nx2,ny2,r,f2)

!Set zero
do i=0,nx2
do j=0,ny2
u2(i,j)=0.0d0
end do
end do


!2.Relax v1 times
do m=1,v1
call relax(nx2,ny2,dx2,dy2,f2,u2)
end do

! Compute residual
call resid(nx2,ny2,dx2,dy2,f2,u2,r2)

!2r.Restriction
call rest(nx2,ny2,nx3,ny3,r2,f3)


!Set zero
do i=0,nx3
do j=0,ny3
u3(i,j)=0.0d0
end do
end do


!3.Relax v1 times
do m=1,v1
call relax(nx3,ny3,dx3,dy3,f3,u3)
end do

! Compute residual
call resid(nx3,ny3,dx3,dy3,f3,u3,r3)


!3r.Restriction
call rest(nx3,ny3,nx4,ny4,r3,f4)


!Set zero
do i=0,nx4
do j=0,ny4
u4(i,j)=0.0d0
end do
end do

!4.Relax v1 times
do m=1,v1
call relax(nx4,ny4,dx4,dy4,f4,u4)
end do


! Compute residual
call resid(nx4,ny4,dx4,dy4,f4,u4,r4)


!4r.Restriction
call rest(nx4,ny4,nx5,ny5,r4,f5)


!Set zero
do i=0,nx5
do j=0,ny5
u5(i,j)=0.0d0
end do
end do

!5.Relax v3 times (or it can be solved exactly)
!call initial residual:
	call resid(nx5,ny5,dx5,dy5,f5,u5,r5)
	call l2norm(nx5,ny5,r5,rmsc)
do m=1,v3
	call relax(nx5,ny5,dx5,dy5,f5,u5)
	call resid(nx5,ny5,dx5,dy5,f5,u5,r5)
	! Check for convergence on smallest grid	
	call l2norm(nx5,ny5,r5,rms)
	if (rms/rmsc.le.tol) goto 11
end do
	11 continue

me = me + m

!4p.Prolongation
call prol(nx5,ny5,nx4,ny4,u5,p4)


!Correct
do i=1,nx4-1
do j=1,ny4-1
u4(i,j) = u4(i,j) + p4(i,j)
end do
end do

!4.Relax v2 times
do m=1,v2
call relax(nx4,ny4,dx4,dy4,f4,u4)
end do

!3p.Prolongation
call prol(nx4,ny4,nx3,ny3,u4,p3)

!Correct
do i=1,nx3-1
do j=1,ny3-1
u3(i,j) = u3(i,j) + p3(i,j)
end do
end do

!3.Relax v2 times
do m=1,v2
call relax(nx3,ny3,dx3,dy3,f3,u3)
end do

!2p.Prolongation
call prol(nx3,ny3,nx2,ny2,u3,p2)


!Correct
do i=1,nx2-1
do j=1,ny2-1
u2(i,j) = u2(i,j) + p2(i,j)
end do
end do

!2.Relax v2 times
do m=1,v2
call relax(nx2,ny2,dx2,dy2,f2,u2)
end do

!1p.Prolongation
call prol(nx2,ny2,nx,ny,u2,p)


!Correct
do i=1,nx-1
do j=1,ny-1
u(i,j) = u(i,j) + p(i,j)
end do
end do


!1.Relax v2 times
do m=1,v2
call relax(nx,ny,dx,dy,f,u)
end do

end do	! Outer iteration loop

10 continue

close(66)


deallocate(r,p,u2,f2,r2,p2,u3,f3,r3,p3,u4,f4,r4,p4,u5,f5,r5)

ke=k

!work load (total number of operations)
wl = ke*(v1+v2)*(nx*ny + nx2*ny2 + nx3*ny3 + nx4*ny4) + me*(nx5*ny5) 

write(19,*)"outer number of iteration = ",ke
write(19,*)"normalized workload       = ",dfloat(wl)/dfloat(nx*ny)
write(*,*)"outer number of iteration = ",ke
write(*,*)"normalized workload       = ",dfloat(wl)/dfloat(nx*ny)

return  
end


