! Copyright (c) 2007 Free Software Foundation, Inc.
! This program creates sampling distributions of Pearson R and slope for random samples selected
! without replacement from a finite population
program correlation
 implicit none
 include 'mpif.h'

! Input variables
 character*4 strz       
 integer n,np,nk,mp
! rho from bivariate distribution
 double precision rho
 double precision rminz,rmaxz 

 integer mp_end
 
! Logical variable to determine if a distribution should be created
 logical create_distribution
! Logical variable to determine if the sampling should be done with replacement
 logical nonreplacement

! Number of populations 
 integer num_populations

! Number of intervals in analytical and simulated distributions 
 integer m, num_intervals
 parameter(m=10000,num_intervals=60)
! x- and y-coordinates of analytical distribution Pearson R distribution
 double precision ra(m+1)
 double precision analytical(m+1)
! Arrays and variables to build simulated distributions
 double precision smid(num_intervals),fdist(num_intervals)      

 logical failed
      
 integer ijk,i,iii,jjj,ii,nitr,iijj
 integer j,itimes,irs,ija
 integer ij,jj,iiata,itimesnum
 integer iiat,jjk,jjk_save
 
! Variables used to determine if sample is unique 
 integer maxindex,minindex,sumindex
 integer maxisum,minisum,minnp

! Analytical expectation and variance of analytical Pearson R distribution 
 double precision expectation,var
! Analytical expectation and variance of analytical slope distribution 
 double precision exp_slope,varb2
! Analytical expectation and variance of analytical Fisher distribution 
 double precision zexprho,varz_a

! Variables used to compute the area difference in simulated and normal Fisher distribution 
 double precision zminz,zmaxz
 double precision area1,area2,area,zvalue
 double precision areadif
 double precision dzr,sumdist,cumarea
 
! Variables used to compute the area difference in simulated and slope distribution
 double precision zminzb,zmaxzb,dzrb
 double precision cumareab
 double precision areadifb 
 double precision nu,facgnu,expgnu,top,bot,fgnu
 double precision facnu,dbconv
 double precision dz,v 

! Variables for normal and slope simulated distribution 
 double precision sminz,smaxz,sminzb,smaxzb 

! Parallel arrays and variables
 integer status
 integer num_processors,image,ierror
 double precision amax(4),amin(4)
 double precision a(4),aa(8),asum(8)

! Maximum iterations to achieve R_{individual} = rho
 integer nitdiff
! Variables used to sample from normal distriution
 real rrr
 double precision diff 
 double precision meanz,stz 
 double precision pnor,qnor,meannor,sdnor
 double precision boundnor,meannor2,sdnor2
 double precision sdnort,sdnor2t
 double precision xmsum,xms,ymsum,yms
 double precision eps
 double precision xmm,ymm
 
      
! Standard deviation of R, z, and slope for distribution
 double precision sigmar,sigmaz,sigma_slope
! Variance of R, z, and slope for distribution
 double precision varr,varz,var_slope

 ! Variable used to compute sampling distribution statistics
 integer nnn,imultiple
 double precision rn     
 double precision ravg,zavg,slopeavg 
 double precision rpairsum2,zpairsum2,slopeavg2
 double precision rminza,rmaxza
 double precision slopemin,slopemax

! Pearson R distribution variables 
 double precision rpairavg,rpair2avg

! Slope distribution variables
 integer nnnm1,nnnm1s,nnns  
 double precision sstop,ssbot,expr2
 double precision facq,hyp,dhyp,ssz,fac,s
 double precision rnnns

! Analytical slope distribution variables 
 double precision sumb
 double precision smida

! Analytical distribution of Pearson R variables 
 double precision pi,facdist,facdist1,series,deriv,dr
 double precision meanu,sigmau

! Array that stores sample sizes 
 integer, allocatable :: npa(:)

!*********************************************************************      
! Arrays that are written to files      
 double precision, allocatable :: pearsonr(:)
 double precision, allocatable :: pearsonr_max(:)
 double precision, allocatable :: pearsonr_min(:)
 double precision, allocatable :: pearsonr_var(:)
 double precision, allocatable :: pearsonr_var_min(:)
 double precision, allocatable :: pearsonr_var_max(:)
     
 double precision, allocatable :: least_squares_slope(:)
 double precision, allocatable :: least_squares_slope_max(:)
 double precision, allocatable :: least_squares_slope_min(:)
 double precision, allocatable :: least_squares_slope_var(:)      
 double precision, allocatable :: least_squares_slope_var_max(:)
 double precision, allocatable :: least_squares_slope_var_min(:)            
 double precision, allocatable :: areab_array(:)

 double precision, allocatable :: fisherz(:),fisherz_var(:)
 double precision, allocatable :: fisherz_area(:)
!*********************************************************************      
! Arrays to determine if sample is unique         
 double precision, allocatable :: arrayr(:)
 double precision, allocatable :: arrayi(:,:)
! Original scores
 double precision, allocatable :: xm(:),ym(:)
 double precision, allocatable :: xma(:,:),yma(:,:)


! R and linear regression parameters for original scores
 double precision r,slope
 double precision, allocatable :: slopeu(:)
 double precision yint,xbar,ybar

! standard deviation of x- and y-scores
 double precision sdx,sdy,rr,sumrr

 integer, allocatable :: iia(:),indexs(:)
 integer, allocatable :: index_save(:)
      
 double precision, allocatable :: xavg(:),yavg(:)
 double precision, allocatable :: rs(:),zs(:)

 double precision, allocatable :: rsave(:),sdnorx(:),sdnory(:)
      
 call MPI_INIT(ierror)
 call MPI_COMM_SIZE(MPI_COMM_WORLD, num_processors, ierror)
 call MPI_COMM_RANK(MPI_COMM_WORLD, image, ierror)

 open(unit=9,file="input_swr")
! Output files that generated
 if (image .eq. 0) then
    open(unit=85,file="pearsonr")
    open(unit=86,file="slope")
    open(unit=87,file="fisher")
    open(unit=96,file="rdistsim")
    open(unit=98,file="sdistsim")
    open(unit=97,file="rdistanalytical")
    open(unit=99,file="sdistanalytical")
 end if

 print*,'image = ',image
! n: length of original scores
! np: size of the group
! nnn : sample size - now calculated
! nk: number of iterations to generate sampling distribution
! mp: number of different sample sizes used      

! Read input file "input_swr"      
 read(9,*) strz,n
 read(9,*) strz,np
 read(9,*) strz,nk
 read(9,*) strz,rho
 read(9,*) strz,rminz
 read(9,*) strz,rmaxz
 read(9,*) strz,mp

! Currently 112 populations are used
 num_populations = 112
 itimes = num_populations/num_processors

 create_distribution = .false.
 create_distribution = .true. 
 if (create_distribution) then
    mp = 1
    itimes = 1
 end if

 itimesnum = itimes*num_processors
      
 nonreplacement = .true.

 allocate(arrayr(-1000002:1000002))
 allocate(arrayi(-1000002:1000002,3))

 eps = 1.d-14

 if (image .eq. 0) then
    print*,'n = ',n,'np = ',np,'nk = ',nk, &
    'rho = ',rho,'rmin = ',rminz,'rmax = ',rmaxz,'mp = ',mp
 end if


 allocate(rsave(itimes),sdnorx(itimes),sdnory(itimes))      
 allocate(indexs(n),index_save(n))

 do jjj = 1,n
    index_save(jjj) = jjj
 end do

 allocate(npa(mp))

 allocate(pearsonr(mp))
 allocate(pearsonr_max(mp))
 allocate(pearsonr_min(mp))
 allocate(pearsonr_var(mp))      
 allocate(pearsonr_var_max(mp))      
 allocate(pearsonr_var_min(mp))      
      
 allocate(least_squares_slope(mp))
 allocate(least_squares_slope_max(mp))
 allocate(least_squares_slope_min(mp))
 allocate(least_squares_slope_var(mp))
 allocate(least_squares_slope_var_max(mp))
 allocate(least_squares_slope_var_min(mp))            
 allocate(areab_array(mp))
      
 allocate(fisherz(mp),fisherz_var(mp))      
 allocate(fisherz_area(mp))
 
 allocate(xm(n),ym(n))
 allocate(xma(n,itimes),yma(n,itimes))

 meanz = 0.d0
 stz = 1.d0

! Bivariate parameters      
 meannor = 0.d0
 sdnor = 1.d0

 meannor2 = 0.d0
 sdnor2 = 1.0d0

! These lines are executed so each image selects different populations      
 sumrr = 0.d0
 do ii = 1,image+20
    rr = rand()
 end do
 if (sumrr > 200.d0) then
    stop
 end if

!**********************************************************************************            
!**********************************************************************************      
! Create itimes populations of size n on each image      
  do ii = 1,itimes
              
     diff = 1.d0
     nitdiff = 0

!    Iterate until rho and R agree to within .0001
     do while (diff .gt. .0001 .and. nitdiff .lt. 100000)
        nitdiff = nitdiff + 1

        xmsum = 0.d0
        xms = 0.d0
        ymsum = 0.d0
        yms = 0.d0
        do j = 1,n
!          Create a random number between 0 and 1
           pnor = min(max(rand(),eps),1.d0-eps)
           qnor = 1.d0 - pnor

!          Sample x from a standard normal distribution               
           call cdfnor(2,pnor,qnor,xmm,meanz,stz,status,boundnor)
                       
           xm(j) = sdnor*xmm + meannor

           xmsum = xmsum + xm(j)
           xms = xms + xm(j)**2

           pnor = min(max(rand(),eps),1.d0-eps)
           qnor = 1.d0 - pnor

!          Sample y from a standard normal distribution               
           call cdfnor(2,pnor,qnor,ymm,meanz,stz,status,boundnor)

!          Correlate xm(j) and ym(j) at level of correlation rho
           ym(j) = sdnor2*(rho*xmm + dsqrt(1.d0-rho**2)*ymm) + meannor2

           ymsum = ymsum + ym(j)
           yms = yms + ym(j)**2

        end do

!       Standard deviation of x            
        sdnort = dsqrt((xms - xmsum**2/dble(n))/dble(n))
!       Standard deviation of y            
        sdnor2t = dsqrt((yms - ymsum**2/dble(n))/dble(n))
            
        do i = 1,n
           xma(i,ii) = xm(i)
           yma(i,ii) = ym(i)
        end do

!       Compute the Pearson R coefficient of the sample            
        call pearsonrsub(r,slope,xm,ym,n,sdx,sdy,yint,xbar,ybar)

        diff = abs(r-rho)

     end do
!    Store the sample Pearson R value and standard deviation in the x- and y-variables
     rsave(ii) = r

     sdnorx(ii) = sdnort
     sdnory(ii) = sdnor2t
         
     if (nitdiff .ge. 100000) then
        print*,'STOP ERROR'
        stop
     end if
         
  end do
!**********************************************************************************
!**********************************************************************************            


 imultiple = n/(mp*np)
! Create mp different sample sizes so the sample percent of populationranges between 0 and 100      
 do iijj = 1,mp
    npa(iijj) = imultiple*iijj
 end do
 if (create_distribution) npa(1) = 50

 
 allocate(rs(nk),slopeu(nk),zs(nk))      
 allocate(iia(np))

 print*,'mp = ',mp
! Loop over different sample sizes      
  do iijj = 1,mp

!    nnn is the sample size         
     nnn = npa(iijj)

     if (image .eq. 0) then
        print*,'iijj = ',iijj,npa(iijj),n
     end if
         
     allocate(xavg(nnn),yavg(nnn))

     pearsonr(iijj) = 0.d0
     pearsonr_max(iijj) = -1.d+15
     pearsonr_min(iijj) = 1.d+15         
     pearsonr_var(iijj) = 0.d0         
     pearsonr_var_max(iijj) = -1.d+15
     pearsonr_var_min(iijj) = 1.d+15
         
     least_squares_slope(iijj) = 0.d0
     least_squares_slope_max(iijj) = -1.d+15
     least_squares_slope_min(iijj) = 1.d+15                  
     least_squares_slope_var(iijj) = 0.d0
     least_squares_slope_var_max(iijj) = -1.d+15
     least_squares_slope_var_min(iijj) = 1.d+15
     areab_array(iijj) = 0.d0

     fisherz(iijj) = 0.d0
     fisherz_var(iijj) = 0.d0         
     fisherz_area(iijj) = 0.d0

         
!    Loop over all the populations on an image
     do ija = 1,itimes

        arrayr = 3.d0
        arrayi = -1
            
        ravg = 0.d0
        zavg = 0.d0
        slopeavg = 0.d0
        slopeavg2 = 0.d0
        rpairsum2 = 0.d0
        zpairsum2 = 0.d0
            
        rminza = 1.d+20
        rmaxza = -1.d+20
        slopemin = 1.d+20
        slopemax = -1.d+20

!       Select nk different samples from the population for the sampling distribution            
        do ij = 1,nk

           if (nonreplacement) then
              indexs = index_save
              jjk = 0
              jjk_save = jjk
           end if
               
           failed = .true.
           nitr = 0
           do while (failed .and. nitr .le. 1000)

              nitr = nitr + 1

              maxindex = -1
              minindex = 100000000
              sumindex = 0
                  
              if (.not. nonreplacement) then
!                Choose nsamples of length nnn with replacement                     
                 do i = 1,nnn

                    xavg(i) = 0.d0
                    yavg(i) = 0.d0

                    do ii = 1,np
                       pnor = min(max(rand(),eps),1.d0-eps)         
                       iii = pnor*dble(n)            
                       iii = iii + 1

                       xavg(i) = xavg(i) + xma(iii,ija)
                       yavg(i) = yavg(i) + yma(iii,ija)
                    end do

                    xavg(i) = xavg(i)/dble(np)
                    yavg(i) = yavg(i)/dble(np)

                 end do

              else
!                Choose nsamples of length nnn without replacement
!                Each element is composed of of an average of np elements
                 do j = 1,nnn

                    maxisum = 0
                    minisum = 0
                    minnp = 100000000
!                   Choose np random unique elements from population (np is the group size)
                    do jj = 1,np
                     
                       rrr = rand()
                      
                       rn = dble(n-jjk)
                       iiat = int(min(max(rrr*rn,eps),rn-eps))+1

                       iiata = iiat
                       iiat = indexs(iiat)

                       indexs(iiata) = indexs(n-jjk)
                     
                       jjk = jjk + 1
                       iia(jj) = iiat
                       minnp = min(minnp,iia(jj))
                       maxisum = maxisum + iia(jj)
                       minisum = minisum + iia(jj)

                    end do

                    maxindex = max(maxisum,maxindex)
                    minindex = min(minisum,minindex)
                    sumindex = sumindex + minnp

                    xavg(j) = 0.d0
                    yavg(j) = 0.d0
                    do jj = 1,np
                       xavg(j) = xavg(j) + xma(iia(jj),ija)
                       yavg(j) = yavg(j) + yma(iia(jj),ija)                  
                    end do

                    xavg(j) = xavg(j)/dble(np)
                    yavg(j) = yavg(j)/dble(np)

                 end do
                     
              end if

!             (xavg,yavg, i = 1,nnn) are the x- and y-values              
!             rs(ij) and slopeu(ij) are the Pearson R coefficient and slope of the sample
              call pearsonrsub(rs(ij),slopeu(ij),xavg,yavg,nnn,sdx,sdy,yint,xbar,ybar)

!             Make sure the choice of the sample is unique                  
              if (nonreplacement) then
                 irs = 1000000*rs(ij)
            
                 if (arrayr(irs) .lt. 2.d0) then
                    diff = abs(rs(ij)-arrayr(irs))
                    if (diff .lt. 1.d-12) then
                       if (maxindex-arrayi(irs,1) .eq. 0 .and. &
                           minindex-arrayi(irs,2) .eq. 0 .and. &
                           sumindex-arrayi(irs,3) .eq. 0) then
!                          print*,'failed = ',failed
!                          print*,'rs = ',rs(ij),arrayr(irs)
                           failed = .true.
                           indexs = index_save 
                           jjk = jjk_save
                        else
                           failed = .false.
                        end if
                    else
                        failed = .false.
                    end if
                 else
                    failed = .false.
                 end if
              else
                 failed = .false.
              end if

           end do  ! End do While loop

           if (nitr .ge. 1000) then
              print*,'Failed to choose a sample'
              stop
           end if
         
           if (nonreplacement) then
!             Save characteristics for each sample
              arrayr(irs) = rs(ij)
              arrayi(irs,1) = maxindex
              arrayi(irs,2) = minindex
              arrayi(irs,3) = sumindex
           end if

!          Calculate the Fisher transformed R value 
           zs(ij) = .5d0*log((1.d0 + rs(ij))/(1.d0 - rs(ij)))
               
           ravg = ravg + rs(ij)
           zavg = zavg + zs(ij)
           rpairsum2 = rpairsum2 + rs(ij)**2
           zpairsum2 = zpairsum2 + zs(ij)**2
           slopeavg = slopeavg + slopeu(ij)
           slopeavg2 = slopeavg2 + slopeu(ij)**2
         
           rminza = min(rminza,rs(ij))
           rmaxza = max(rmaxza,rs(ij))
           slopemin = min(slopemin,slopeu(ij))
           slopemax = max(slopemax,slopeu(ij))

!       nk loop    
        end do

!       Compute standard deviations of Pearson R, Z, and slope for the sample distribution            
        sigmar = dsqrt( (rpairsum2 - ravg**2/dble(nk))/dble(nk) )
        sigmaz = dsqrt((zpairsum2 - zavg**2/dble(nk))/dble(nk))            
        sigma_slope = dsqrt((slopeavg2 - slopeavg**2/dble(nk))/dble(nk))

!       Compute variances            
        varr = sigmar**2
        var_slope = sigma_slope**2
        varz = sigmaz**2
            
!       average of r average scores
        ravg = ravg/dble(nk)

!       average of z scores
        zavg = zavg/dble(nk)

!       average of slope scores            
        slopeavg = slopeavg/dble(nk)            


!**********************************************************************************
!**********************************************************************************                        
!       Compute the area difference from a normal distribution            
        zminz = -3.5d0
        zmaxz =  3.5d0

        dzr = (zmaxz-zminz)/dble(num_intervals)
            
        fdist = 0.d0
        cumarea = 0.d0
        areadif = 0.d0
!       Use num_intervals intervals in distribution
        do i = 1,num_intervals
           sminz = zminz + dble(i-1)*dzr
           smaxz = zminz + dble(i)*dzr
           smid(i) = .5d0*(sminz+smaxz)

           call cdfnor(1,area1,qnor,sminz,0.d0,1.d0,status,boundnor)
           call cdfnor(1,area2,qnor,smaxz,0.d0,1.d0,status,boundnor)
           area = area2 - area1
           cumarea = cumarea + area
           do ijk = 1,nk
              zvalue = (zs(ijk)-zavg)/sigmaz
              if (zvalue >= sminz .and. zvalue < smaxz) then                  
                 fdist(i) = fdist(i) + 1.d0
              end if
           end do
           areadif = areadif + abs(fdist(i)/dble(nk)-area)
        end do
!**********************************************************************************
!**********************************************************************************            


!**********************************************************************************
!**********************************************************************************                        
!       Find area difference between sample and analytical b distribution
!       The variance is adjusted to account for the sample percent of the population      
        ii = 2
        fac = 1.d0

        facq = dble(nnn)/dble(n)
        facq = 1.d0 - facq
        rnnns = dble(nnn-3)/facq + 3.d0
        nnns = nint(rnnns)            

        s = rnnns
        ssz = s/2.d0
        
        do while (ssz .gt. 3)
           fac = fac*(s-dble(ii))/(s-dble(ii+1))
           ssz = (s-dble(ii))/2.d0
           ii = ii + 2
        end do

        call GAMMA(ssz,sstop)
        call GAMMA(ssz-.5d0,ssbot)
        fac = fac*sstop/ssbot
        facnu = fac

        nnnm1s = nnns - 1

        nu = s - 1.d0
        facgnu = facnu/dsqrt(nu*acos(-1.d0))
        expgnu = -(nu + 1.d0)/2.d0

        bot = dsqrt(1.d0 - rsave(ija)**2)

        dbconv = dsqrt(nu)*sdnorx(ija)/(bot*sdnory(ija))

        zminzb = slopemin
        zmaxzb = slopemax

        dzrb = (zmaxzb-zminzb)/dble(num_intervals)

        fdist = 0.d0
        cumareab = 0.d0
        areadifb = 0.d0
!       Use num_intervals intervals
        do i = 1,num_intervals
           sminzb = zminzb + dble(i-1)*dzrb
           smaxzb = zminzb + dble(i)*dzrb

           top = (sminzb*sdnorx(ija)/sdnory(ija) - rsave(ija))
           fgnu = top*dsqrt(nu)/bot
           call hypser(0.5d0,(nu+1.d0)/2.d0,1.5d0,-fgnu**2/nu,hyp,dhyp)
           area1 = .5d0 + fgnu*facgnu*hyp

           top = (smaxzb*sdnorx(ija)/sdnory(ija) - rsave(ija))
           fgnu = top*dsqrt(nu)/bot
           call hypser(0.5d0,(nu+1.d0)/2.d0,1.5d0,-fgnu**2/nu,hyp,dhyp)
           area2 = .5d0 + fgnu*facgnu*hyp               
               
           area = area2 - area1
           cumareab = cumareab + area
           do ijk = 1,nk
              if (slopeu(ijk) >= sminzb .and. slopeu(ijk) < smaxzb) then                 
                 fdist(i) = fdist(i) + 1.d0
              end if
           end do
           areadifb = areadifb + abs(fdist(i)/dble(nk)-area)
        end do
!**********************************************************************************
!**********************************************************************************            

!       Calculate expectation (expectation) and  variance (var)
!       of analytical Pearson R distribution
        call mct(nnn,rsave(ija),sdnorx(ija),sdnory(ija),expectation,var)            
            

!       Analytical expectation (exp_slope) of analytical slope distribution            
        exp_slope = rsave(ija)*sdnory(ija)/sdnorx(ija)            

!       Variance of analytical slope distribution            
        varb2 = (1.d0 - rsave(ija)**2)/(dble(nnn) - 3.d0)
        varb2 = varb2*sdnory(ija)**2/sdnorx(ija)**2            

!       Mean of Fisher distribution
        zexprho = .5d0*log((1.d0+rsave(ija))/(1.d0-rsave(ija)))
!       Variance of Fisher distribution
        varz_a = 1.d0/dble(nnn-3) 

!       Sum, Max, Min each metric for each population
            
        pearsonr(iijj) = pearsonr(iijj) + (ravg-expectation)
        pearsonr_max(iijj) = max(pearsonr_max(iijj),(ravg-expectation))
        pearsonr_min(iijj) = min(pearsonr_min(iijj),(ravg-expectation))            

        pearsonr_var(iijj) = pearsonr_var(iijj) + (varr-var)/var
        pearsonr_var_min(iijj) = min(pearsonr_var_min(iijj),(varr-var)/var)           
        pearsonr_var_max(iijj) = max(pearsonr_var_max(iijj),(varr-var)/var)           

        least_squares_slope(iijj) = least_squares_slope(iijj) + (slopeavg-exp_slope)        
        least_squares_slope_max(iijj) = max(least_squares_slope_max(iijj),(slopeavg-exp_slope))
        least_squares_slope_min(iijj) = min(least_squares_slope_min(iijj),(slopeavg-exp_slope))                        
            
        least_squares_slope_var(iijj) = least_squares_slope_var(iijj) + ((var_slope-varb2)/varb2)            
        least_squares_slope_var_max(iijj) = max(least_squares_slope_var_max(iijj),((var_slope-varb2)/varb2))            
        least_squares_slope_var_min(iijj) = min(least_squares_slope_var_min(iijj),((var_slope-varb2)/varb2))            

        areab_array(iijj) = areab_array(iijj) + areadifb            
            

        fisherz(iijj) = fisherz(iijj) + (zavg - zexprho)            
        fisherz_var(iijj) = fisherz_var(iijj) + (varz - varz_a)/varz_a
        fisherz_area(iijj) = fisherz_area(iijj) + areadif
            
!    itimes loop      
     end do

!    Average sum of population values     
     pearsonr(iijj) = pearsonr(iijj)/dble(itimes)         
     pearsonr_var(iijj) = pearsonr_var(iijj)/dble(itimes)

     least_squares_slope(iijj) = least_squares_slope(iijj)/dble(itimes)
     least_squares_slope_var(iijj) = least_squares_slope_var(iijj)/dble(itimes)         

     areab_array(iijj) = areab_array(iijj)/dble(itimes)

     fisherz(iijj) = fisherz(iijj)/dble(itimes)
     fisherz_var(iijj) = fisherz_var(iijj)/dble(itimes)
     fisherz_area(iijj) = fisherz_area(iijj)/dble(itimes)
            
     aa(1) = pearsonr(iijj)         
     aa(2) = pearsonr_var(iijj)
     aa(3) = least_squares_slope(iijj)
     aa(4) = least_squares_slope_var(iijj)
     aa(5) = areab_array(iijj)                  
     aa(6) = fisherz(iijj)
     aa(7) = fisherz_var(iijj)
     aa(8) = fisherz_area(iijj)


     call MPI_REDUCE(aa,asum,8,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)

     a(1) = pearsonr_max(iijj)         
     a(2) = pearsonr_var_max(iijj)
     a(3) = least_squares_slope_max(iijj)
     a(4) = least_squares_slope_var_max(iijj)

     call MPI_REDUCE(a,amax,4,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierror)

     a(1) = pearsonr_min(iijj)         
     a(2) = pearsonr_var_min(iijj)
     a(3) = least_squares_slope_min(iijj)
     a(4) = least_squares_slope_var_min(iijj)

     call MPI_REDUCE(a,amin,4,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierror)         

     if (image .eq. 0) then
         
        do i = 1,8
           asum(i) = asum(i)/dble(num_processors)
        end do

        pearsonr(iijj) = asum(1)
        pearsonr_var(iijj) = asum(2)
        least_squares_slope(iijj) = asum(3)
        least_squares_slope_var(iijj) = asum(4)
        areab_array(iijj) = asum(5)            
        fisherz(iijj) = asum(6)
        fisherz_var(iijj) = asum(7)
        fisherz_area(iijj) = asum(8)
            
        pearsonr_max(iijj) = amax(1)
        pearsonr_var_max(iijj) = amax(2)
        least_squares_slope_max(iijj) = amax(3)            
        least_squares_slope_var_max(iijj) = amax(4)

        pearsonr_min(iijj) = amin(1)
        pearsonr_var_min(iijj) = amin(2)
        least_squares_slope_min(iijj) = amin(3)            
        least_squares_slope_var_min(iijj) = amin(4)
            
     end if

     deallocate(xavg,yavg)
      
  end do ! mp loop
  
  if (image .eq. 0) then
     
     do iijj = 1,mp
        write(85,*) dble(np*npa(iijj)/dble(n)), &
                    pearsonr(iijj), &
                    pearsonr_max(iijj), &
                    pearsonr_min(iijj), &
                    pearsonr_var(iijj), &
                    pearsonr_var_max(iijj), &
                    pearsonr_var_min(iijj)
        write(86,*) dble(np*npa(iijj)/dble(n)), &
                    least_squares_slope(iijj), &
                    least_squares_slope_max(iijj), &
                    least_squares_slope_min(iijj), &
                    least_squares_slope_var(iijj), &
                    least_squares_slope_var_max(iijj), &
                    least_squares_slope_var_min(iijj), &
                    areab_array(iijj)
        write(87,*) dble(np*npa(iijj)/dble(n)), &
                    fisherz(iijj), &
                    fisherz_var(iijj), &
                    fisherz_area(iijj)
     end do   
            

     if (create_distribution) then

!**********************************************************************************
!**********************************************************************************                        
!       Sample distribution of simulated Pearson R            
        fdist = 0.d0
        sumdist = 0.d0
        dzr = (rmaxz - rminz)/dble(num_intervals)
            
!       Use num_intervals intervals
        do i = 1,num_intervals
           sminz = rminz + dble(i-1)*dzr
           smaxz = rminz + dble(i)*dzr
           smid(i) = .5d0*(sminz+smaxz)
           do ijk = 1,nk
              if (rs(ijk) >= sminz .and. rs(ijk) < smaxz) then                  
                 fdist(i) = fdist(i) + 1.d0
              end if
           end do
           sumdist = sumdist + fdist(i)/dble(nk)
           write(96,*) smid(i),fdist(i)/(dble(dzr*nk))
        end do
!**********************************************************************************
!**********************************************************************************                        


!**********************************************************************************
!**********************************************************************************                        
!       Sample distribution of linear regression slope            
        slopemax = slopeavg + 3.5d0*sigma_slope
        slopemin = slopeavg - 3.5d0*sigma_slope
            
        dz = (slopemax - slopemin)/dble(num_intervals)

        fdist = 0.d0
        sumdist = 0.d0
        do i = 1,num_intervals
           sminz = slopemin + dble(i-1)*dz
           smaxz = slopemin + dble(i)*dz
           smid(i) = .5d0*(sminz+smaxz)
           do ijk = 1,nk
              if (slopeu(ijk) >= sminz .and. slopeu(ijk) < smaxz) then                  
                 fdist(i) = fdist(i) + 1.d0
              end if
           end do
           sumdist = sumdist + fdist(i)/dble(nk)
           write(98,*) smid(i),fdist(i)/(dble(dz*nk))
        end do
!**********************************************************************************
!**********************************************************************************                        


!**********************************************************************************
!**********************************************************************************                        
!       Analytical distribution of Pearson R
        pi = acos(-1.d0)
        facdist = 1.d0
        do i = nnn-2,1,-1
           facdist = facdist*dble(i)/(dble(i)+.5d0)
        end do
        facdist = facdist/(.5d0*dsqrt(pi))
        facdist = facdist*dble(nnn-2)/dsqrt(2.d0*pi)
        facdist1 = facdist*(1.d0 - rsave(1)**2)**(dble(nnn-1)/2.d0)
        dr = 1.999d0/dble(m)
        rpairavg = 0.d0
        rpair2avg = 0.d0
        do i = 1,m+1
           ra(i) = -1.d0 + dr*dble(i-1) + 1.d-8
           facdist = facdist1*(1.d0 - rsave(1)*ra(i))**(-dble(nnn)+1.5d0)
           facdist =facdist*(1.d0 - ra(i)**2)**(dble(nnn)/2.d0-2.d0)
           call hypser(0.5d0,0.5d0,dble(nnn)-0.5d0,0.5d0*(1.d0+rsave(1)*ra(i)),series,deriv)
           analytical(i) = facdist*series
           if (ra(i) .gt. rminz .and. ra(i) .lt. rmaxz) then
              write(97,*) ra(i),analytical(i)
           end if
           rpairavg = rpairavg + ra(i)*analytical(i)
           rpair2avg = rpair2avg + (ra(i)**2)*analytical(i)
        end do
        meanu = rpairavg*dr
        rpair2avg = rpair2avg*dr
        sigmau = dsqrt(rpair2avg - meanu**2)
!**********************************************************************************
!**********************************************************************************                        


!**********************************************************************************
!**********************************************************************************                        
!       Analytical distribution of linear regression slope

        dz = (slopemax - slopemin)/dble(m)
        print*,'nnn = ',nnn,rsave(1)
        s = dble(nnn)
        ssz = s/2.d0

        ii = 2
        fac = 1.d0
        do while (ssz .gt. 3)
           fac = fac*(s-dble(ii))/(s-dble(ii+1))
           ssz = (s-dble(ii))/2.d0
           ii = ii + 2
        end do

        call GAMMA(ssz,sstop)
        call GAMMA(ssz-.5d0,ssbot)
        fac = fac*sstop/ssbot
        fac = fac*((1.d0-rho**2)**((s-1.d0)/2.d0))
        fac = fac/sqrt(acos(-1.d0))
        
        sumb = 0.d0
        do i = 1,m
           sminz = slopemin + dble(i-1)*dz
           smaxz = slopemin + dble(i)*dz
           smida = .5d0*(sminz+smaxz)
           v = fac*(1.d0 - rsave(1)**2 + (rsave(1) - (sdnorx(1)/sdnory(1))*smida)**2)**(-s/2.d0)
           v = v*sdnorx(1)/sdnory(1)
           write(99,*) smida,v
           sumb = sumb + v*dz
        end do
!**********************************************************************************
!**********************************************************************************            
     end if   
      
  end if

  call MPI_FINALIZE(ierror)         

end
