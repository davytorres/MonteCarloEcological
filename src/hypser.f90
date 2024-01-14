      SUBROUTINE hypser(a,b,c,z,series,deriv)
      implicit none
      INTEGER n
      double precision a,b,c,z,series,deriv,aa,bb,cc,fac,temp
!     Returns the hypergeometric series 2F1 and its derivative, iterating to machine accuracy
!     For cabs(z) ? 1/2 convergence is quite rapid.                              
                                                                                  
      deriv=0.d0
      fac=1.d0
      temp=fac
      aa=a
      bb=b
      cc=c
      do n=1,1000
         fac=((aa*bb)/cc)*fac
         deriv=deriv+fac
         fac=fac*z/n
                                                                                  
         series=temp+fac
         if (series.eq.temp) return
         temp=series
         aa=aa+1.
         bb=bb+1.
         cc=cc+1.
      end do
      print*,'convergence failure in hypser'
      stop
      END
