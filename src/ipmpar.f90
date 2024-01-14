      INTEGER FUNCTION ipmpar(i)
!-----------------------------------------------------------------------                             
!                                                                                                    
!     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER                                 
!     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER                                  
!     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...                                     
!                                                                                                    
!  INTEGERS.                                                                                         
!                                                                                                    
!     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM                                    
!                                                                                                    
!               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )                                       
!                                                                                                    
!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.                                            
!                                                                                                    
!     IPMPAR(1) = A, THE BASE.                                                                       
!                                                                                                    
!     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS.                                                    
!                                                                                                    
!     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE.                                                   
!                                                                                                    
!  FLOATING-POINT NUMBERS.                                                                           
!                                                                                                    
!     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING                                    
!     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE                                      
!     NONZERO NUMBERS ARE REPRESENTED IN THE FORM                                                    
!                                                                                                    
!               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)                                             
!                                                                                                    
!               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,                                              
!               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.                                              
!                                                                                                    
!     IPMPAR(4) = B, THE BASE.                                                                       
!                                                                                                    
!  SINGLE-PRECISION                                                                                  
!                                                                                                    
!     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.                                                    
!                                                                                                    
!     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.                                                     
!                                                                                                    
!     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.                                     

!  DOUBLE-PRECISION                                                                                  
!                                                                                                    
!     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.                                                    
!                                                                                                    
!     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.                                                     
!                                                                                                    
!     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.                                                     
!                                                                                                    
!-----------------------------------------------------------------------                             
!                                                                                                    
!     TO DEFINE THIS FUNCTION FOR THE COMPUTER BEING USED, ACTIVATE                                  
!     THE DATA STATMENTS FOR THE COMPUTER BY REMOVING THE C FROM                                     
!     COLUMN 1. (ALL THE OTHER DATA STATEMENTS SHOULD HAVE C IN                                      
!     COLUMN 1.)                                                                                     
!                                                                                                    
!-----------------------------------------------------------------------                             
!                                                                                                    
!     IPMPAR IS AN ADAPTATION OF THE FUNCTION I1MACH, WRITTEN BY                                     
!     P.A. FOX, A.D. HALL, AND N.L. SCHRYER (BELL LABORATORIES).                                     
!     IPMPAR WAS FORMED BY A.H. MORRIS (NSWC). THE CONSTANTS ARE                                     
!     FROM BELL LABORATORIES, NSWC, AND OTHER SOURCES.                                               
!                                                                                                    
!-----------------------------------------------------------------------                             
!     .. Scalar Arguments ..
      implicit none  
      INTEGER i
!     ..                                                                                             
!     .. Local Arrays ..                                                                             
      INTEGER imach(10)
!     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T                               
!     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T                                  
!     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).                                   
!                                                                                                    
      DATA IMACH( 1) /     2 /
      DATA IMACH( 2) /    31 /
      DATA IMACH( 3) / 2147483647 /
      DATA IMACH( 4) /     2 /
      DATA IMACH( 5) /    24 /
      DATA IMACH( 6) /  -125 /
      DATA IMACH( 7) /   128 /
      DATA IMACH( 8) /    53 /
      DATA IMACH( 9) / -1021 /
      DATA IMACH(10) /  1024 /

      ipmpar = imach(i)
      RETURN
      END
