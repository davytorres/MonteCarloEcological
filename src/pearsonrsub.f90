subroutine pearsonrsub(r,slope,x,y,n,sdx,sdy,yint,xbar,ybar)
! Computes Pearson R coefficient and linear regression line given
! a set of original scores
 implicit none

! Input
! Length of set of original scores
 integer n
! Original scores 
 double precision x(n),y(n)

! Output

! Standard deviation of x- and y-scores
 double precision sdx,sdy
! Mean of x- and y-scores
 double precision xbar,ybar
! y-intercept of linear regression line
 double precision yint
! Slope of linear regression line
 double precision slope

! Local variables
 integer i
 double precision sumx,sumy,sumx2,sumy2,sumxy
 double precision r,rn,rtop,rbot
 double precision sigmax,sigmay

 sumx = 0.d0
 sumy = 0.d0
 sumx2 = 0.d0
 sumy2 = 0.d0
 sumxy = 0.d0
 do i = 1,n
    sumx = sumx + x(i)
    sumy = sumy + y(i)
    sumx2 = sumx2 + x(i)**2
    sumy2 = sumy2 + y(i)**2
    sumxy = sumxy + x(i)*y(i)
 end do

 rn = 1.d0/dble(n)

 ybar = sumy*rn
 xbar = sumx*rn

 rtop = sumxy - sumx*sumy/n
 sigmax = sumx2 - rn*(sumx**2)
 sigmay = sumy2 - rn*(sumy**2)
      
 rbot = dsqrt(sigmax*sigmay)

 r = rtop/rbot

 sdx = dsqrt(sigmax*rn)
 sdy = dsqrt(sigmay*rn)

 slope = r*dsqrt(sigmay/sigmax)

 yint = ybar - slope*xbar


 return
end
