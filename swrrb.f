c     Copyright (c) 2007 Free Software Foundation, Inc.
      program correlation
      implicit none
      include 'mpif.h'

      integer nnnm1,ndnp
      logical create_distribution
      double precision s,fac,ssz
      double precision rpairavg,rpair2avg,meanu
      double precision sstop,ssbot,expr2
      double precision hyp,dhyp

      integer nnns,nnnm1s
      double precision facq,rnnns
      
      integer m,num_intervals,ijk
      double precision dzr,cumfdist,sumdist,cumarea
      double precision sminz,smaxz,sminzb,smaxzb
      double precision dz,sigmau,v
      
      parameter(m=10000,num_intervals=60)
c     x- and y-coordinates of analytical distribution
      double precision ra(m+1)
      double precision analytical(m+1)
c     Arrays and variables to build partition distribution
      double precision smid(num_intervals),fdist(num_intervals)      
      integer num_processors,image,ierror
      
c     Number of points in analytical distribution
      integer nnn,mp,imultiple
      
      character*4 strz      
      logical failed,nonreplacement
      
      integer i,nk,iii,jjj,ii,nitr,iijj
      integer n,np,j,itimes,irs,ija
      integer ij,jj,iiata,itimesnum
      integer iiat,jjk,jjk_save
      integer maxindex,minindex,sumindex
      integer maxisum,minisum,minnp
      
      real rrr      
      double precision expectation,var,zexp,zexprho
      double precision rminza,rmaxza,zminz,zmaxz
      double precision area1,area2,area,zvalue
      double precision areadif
      double precision slopemin,slopemax
      double precision meanz,stz
      double precision xmsum,xms,ymsum,yms
      double precision eps,rn,dist_ks
      double precision xmm,ymm
      double precision svarb2,svar
      double precision varr,var_slope
      double precision varz,varz_a
c      double precision a(8),asum(8),
      double precision amax(10),amin(10)
      double precision a(10),aa(21),asum(21)

      double precision standard_dev_r
      double precision standard_dev_rs
      double precision standard_dev_s
      double precision standard_dev_ss
      
c     arrays to sample from standard normal distriution
      integer status
      double precision pnor,qnor,meannor,sdnor
      double precision boundnor,meannor2,sdnor2
      double precision sdnort,sdnor2t
      
c     standard deviation of R average values
      double precision sigmar,sigmaz

      double precision sumr_error,sums_error
      double precision rpairsum2,zpairsum2
      double precision diff
      double precision slopeavg2
      double precision ravg,zavg,slopeavg
      double precision varb2,sigma_slope,exp_slope
      double precision rminz,rmaxz
c     rho from bivariate distribution
      double precision rho
      double precision sumrs_error,sumss_error

      double precision nu,facgnu,expgnu,top,bot,fgnu
      double precision facnu,dbconv,sumb

      double precision zminzb,zmaxzb,dzrb
      double precision cumfdistb,cumareab
      double precision areadifb
            
      integer num_populations

c      double precision slope_sd

      double precision pi,facdist,facdist1,series,deriv,dr

      
      integer, allocatable :: npa(:)

      double precision, allocatable :: ravg_array_abs(:)
      double precision, allocatable :: ravg_array_abs_max(:)
      double precision, allocatable :: ravg_array_abs_min(:)
      double precision, allocatable :: ravg_array_dif(:)
      double precision, allocatable :: ravg_array_value(:)

      double precision, allocatable :: slope_array_abs(:)
      double precision, allocatable :: slope_array_abs_max(:)
      double precision, allocatable :: slope_array_abs_min(:)
      double precision, allocatable :: slope_array_dif(:)
      double precision, allocatable :: slope_array_value(:)      
      
      double precision, allocatable :: sigmar_array(:),ravg_array(:)
      double precision, allocatable :: sigmar_arrayd(:),ravg_arrayd(:)      
      double precision, allocatable :: slope_array(:),d_array(:)
      double precision, allocatable :: zrho_array(:),varz_array(:)
      double precision, allocatable :: z_array(:),zq_array(:)
      double precision, allocatable :: area_array(:),areab_array(:)
      double precision, allocatable :: sigma_slope_array(:)
      double precision, allocatable :: slope_arrayd(:)
      double precision, allocatable :: sigma_slope_arrayd(:)

      double precision, allocatable :: sigmar_array_max(:)
      double precision, allocatable :: ravg_array_max(:)
      double precision, allocatable :: sigmar_arrayd_max(:)
      double precision, allocatable :: ravg_arrayd_max(:)      
      double precision, allocatable :: slope_array_max(:)
      double precision, allocatable :: sigma_slope_array_max(:)
      double precision, allocatable :: slope_arrayd_max(:)
      double precision, allocatable :: sigma_slope_arrayd_max(:)

      double precision, allocatable :: sigmar_array_min(:)
      double precision, allocatable :: ravg_array_min(:)
      double precision, allocatable :: sigmar_arrayd_min(:)
      double precision, allocatable :: ravg_arrayd_min(:)      
      double precision, allocatable :: slope_array_min(:)
      double precision, allocatable :: sigma_slope_array_min(:)
      double precision, allocatable :: slope_arrayd_min(:)
      double precision, allocatable :: sigma_slope_arrayd_min(:)            
      

         
      double precision, allocatable :: arrayr(:)
      double precision, allocatable :: arrayi(:,:)
c     Original scores
      double precision, allocatable :: xm(:),ym(:)
      double precision, allocatable :: xma(:,:),yma(:,:)

c     Maximum iterations to achieve R_{individual} = rho
      integer nitdiff

c     R and linear regression parameters for original scores
      double precision r,slope
      double precision, allocatable :: slopeu(:)
      double precision yint,xbar,ybar

c     standard deviation of x- and y-scores
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
c     Output files that generated
      if (image .eq. 0) then
         open(unit=83,file="rnp")
         open(unit=84,file="snp")
         open(unit=85,file="re")
         open(unit=86,file="se")
         open(unit=96,file="rdistsim")
         open(unit=98,file="sdistsim")
         open(unit=97,file="rdista")
         open(unit=99,file="sdista")
      end if

      print*,'image = ',image
c     n: length of original scores
c     np: size of the group
c     nnn : sample size - now calculated
c     nk: number of iterations to generate sampling distribution
c     mp: number of different sample sizes used      

c     Read input file "input_swr"      
      read(9,*) strz,n
      read(9,*) strz,np
      read(9,*) strz,nk
      read(9,*) strz,rho
      read(9,*) strz,rminz
      read(9,*) strz,rmaxz
      read(9,*) strz,mp

c     Currently 112 populations are used
      num_populations = 112
      itimes = num_populations/num_processors

      itimesnum = itimes*num_processors
      
      nonreplacement = .true.
      create_distribution = .false.
      allocate(arrayr(-1000002:1000002))
      allocate(arrayi(-1000002:1000002,3))

      eps = 1.d-14

      standard_dev_r = 0.d0
      standard_dev_rs = 0.d0
      standard_dev_s = 0.d0
      standard_dev_ss = 0.d0

      
      if (image .eq. 0) then
         print*,'n = ',n,'np = ',np,'nk = ',nk,
     &   'rho = ',rho,'rmin = ',rminz,'rmax = ',rmaxz,'mp = ',mp
      end if

      allocate(rsave(itimes),sdnorx(itimes),sdnory(itimes))      
      allocate(indexs(n),index_save(n))

      do jjj = 1,n
         index_save(jjj) = jjj
      end do

      allocate(npa(mp))

      allocate(ravg_array_abs(mp))
      allocate(ravg_array_abs_max(mp))
      allocate(ravg_array_abs_min(mp))
      allocate(ravg_array_dif(mp))
      allocate(ravg_array_value(mp))

      allocate(slope_array_abs(mp))
      allocate(slope_array_abs_max(mp))
      allocate(slope_array_abs_min(mp))
      allocate(slope_array_dif(mp))
      allocate(slope_array_value(mp))                        
      
      allocate(sigmar_array(mp),ravg_array(mp))
      allocate(sigmar_arrayd(mp),ravg_arrayd(mp))      
      allocate(slope_array(mp),sigma_slope_array(mp))
      allocate(slope_arrayd(mp),sigma_slope_arrayd(mp))
      allocate(d_array(mp),area_array(mp),z_array(mp))
      allocate(areab_array(mp))
      allocate(zrho_array(mp),varz_array(mp),zq_array(mp))
      allocate(sigmar_array_max(mp),ravg_array_max(mp))
      allocate(sigmar_arrayd_max(mp),ravg_arrayd_max(mp))      
      allocate(slope_array_max(mp),sigma_slope_array_max(mp))
      allocate(slope_arrayd_max(mp),sigma_slope_arrayd_max(mp))

      allocate(sigmar_array_min(mp),ravg_array_min(mp))
      allocate(sigmar_arrayd_min(mp),ravg_arrayd_min(mp))      
      allocate(slope_array_min(mp),sigma_slope_array_min(mp))
      allocate(slope_arrayd_min(mp),sigma_slope_arrayd_min(mp))            
      
      allocate(xm(n),ym(n))
      allocate(xma(n,itimes),yma(n,itimes))

      meanz = 0.d0
      stz = 1.d0

c     Bivariate parameters      
      meannor = 0.d0
      sdnor = 1.d0

      meannor2 = 0.d0
      sdnor2 = 1.0d0

c     These lines are executed so each image selects different populations      
      sumrr = 0.d0
      do ii = 1,image+20
         rr = rand()
      end do
      if (sumrr > 200.d0) then
         stop
      end if

c     Create itimes populations of size n on each image      
      do ii = 1,itimes
              
         diff = 1.d0
         nitdiff = 0

c        Iterate until rho and R agree to within .0001
         do while (diff .gt. .0001 .and. nitdiff .lt. 100000)
            nitdiff = nitdiff + 1

            xmsum = 0.d0
            xms = 0.d0
            ymsum = 0.d0
            yms = 0.d0
            do j = 1,n
c              Create a random number between 0 and 1
               pnor = min(max(rand(),eps),1.d0-eps)
               qnor = 1.d0 - pnor

c              Sample x from a standard normal distribution               
               call cdfnor(2,pnor,qnor,xmm,meanz,stz,
     &                     status,boundnor)
                       
               xm(j) = sdnor*xmm + meannor

               xmsum = xmsum + xm(j)
               xms = xms + xm(j)**2

               pnor = min(max(rand(),eps),1.d0-eps)
               qnor = 1.d0 - pnor

c              Sample y from a standard normal distribution               
               call cdfnor(2,pnor,qnor,ymm,meanz,stz,
     &                     status,boundnor)

c              correlate xm(j) and ym(j) at level of correlation rho
c              ym(j) = rho*xm(j) + dsqrt(1.d0 - rho**2)*ym(j)

               ym(j) = sdnor2*(rho*xmm + dsqrt(1.d0-rho**2)*ymm)
     &               + meannor2

               ymsum = ymsum + ym(j)
               yms = yms + ym(j)**2

            end do

c           Standard deviation of x            
            sdnort = dsqrt((xms - xmsum**2/dble(n))/dble(n))
c           Standard deviation of y            
            sdnor2t = dsqrt((yms - ymsum**2/dble(n))/dble(n))
            
            do i = 1,n
               xma(i,ii) = xm(i)
               yma(i,ii) = ym(i)
            end do

c           Compute the Pearson R coefficient of the sample            
            call pearsonr(r,slope,
     &                    xm,ym,n,sdx,sdy,yint,xbar,ybar)

            diff = abs(r-rho)

         end do
c        Store the sample Pearson R value and standard deviation in the x- and y-variables
         rsave(ii) = r

         sdnorx(ii) = sdnort
         sdnory(ii) = sdnor2t
         
c         rsave(ii) = rho

         if (nitdiff .ge. 100000) then
            print*,'STOP ERROR'
            stop
         end if
         
      end do


      imultiple = n/(mp*np)
c     Create mp different sample sizes so the sample percent of populationranges between 0 and 100      
      do iijj = 1,mp
         npa(iijj) = imultiple*iijj
      end do

      allocate(rs(nk),slopeu(nk),zs(nk))      
      allocate(iia(np))

c     Loop over different sample sizes      
      do iijj = 1,mp

c        nnn is the sample size         
         nnn = npa(iijj)

         if (image .eq. 0) then
            print*,'iijj = ',iijj,npa(iijj),n
         end if
         
         allocate(xavg(nnn),yavg(nnn))

         ravg_array_abs(iijj) = 0.d0
         ravg_array_dif(iijj) = 0.d0
         ravg_array_value(iijj) = 0.d0

         slope_array_abs(iijj) = 0.d0
         slope_array_dif(iijj) = 0.d0
         slope_array_value(iijj) = 0.d0         

         area_array(iijj) = 0.d0
         areab_array(iijj) = 0.d0
         
         sigmar_array(iijj) = 0.d0
         ravg_array(iijj) = 0.d0
         sigmar_arrayd(iijj) = 0.d0
         ravg_arrayd(iijj) = 0.d0
         slope_array(iijj) = 0.d0
         sigma_slope_array(iijj) = 0.d0
         slope_arrayd(iijj) = 0.d0
         sigma_slope_arrayd(iijj) = 0.d0
         d_array(iijj) = 0.d0
         z_array(iijj) = 0.d0
         zq_array(iijj) = 0.d0
         zrho_array(iijj) = 0.d0
         varz_array(iijj) = 0.d0

         ravg_array_abs_max(iijj) = -1.d+15
         ravg_array_abs_min(iijj) = 1.d+15

         slope_array_abs_max(iijj) = -1.d+15
         slope_array_abs_min(iijj) = 1.d+15         
         
         sigmar_array_max(iijj) = -1.d+15
         ravg_array_max(iijj) = -1.d+15
         sigmar_arrayd_max(iijj) = -1.d+15
         ravg_arrayd_max(iijj) = -1.d+15
         slope_array_max(iijj) = -1.d+15
         sigma_slope_array_max(iijj) = -1.d+15
         slope_arrayd_max(iijj) = -1.d+15
         sigma_slope_arrayd_max(iijj) = -1.d+15

         sigmar_array_min(iijj) = 1.d+15
         ravg_array_min(iijj) = 1.d+15
         sigmar_arrayd_min(iijj) = 1.d+15
         ravg_arrayd_min(iijj) = 1.d+15
         slope_array_min(iijj) = 1.d+15
         sigma_slope_array_min(iijj) = 1.d+15
         slope_arrayd_min(iijj) = 1.d+15
         sigma_slope_arrayd_min(iijj) = 1.d+15
         
c        Loop over all the populations on an image
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

c           Select nk different samples from population for sampling distribution            
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
c                    Choose nsamples of length nnn with replacement                     
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
c                    Choose nsamples of length nnn without replacement
                     do j = 1,nnn

                        maxisum = 0
                        minisum = 0
                        minnp = 100000000
c                       Choose np random unique elements from population (np is the group size)
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

c                 rs(ij) and slopeu(ij) are the Pearson R coefficient and slope of the sample                  
                  call pearsonr(rs(ij),slopeu(ij),
     &                 xavg,yavg,nnn,sdx,sdy,yint,xbar,ybar)

c                 Make sure the choice of the sample is unique                  
                  if (nonreplacement) then
                     irs = 1000000*rs(ij)
            
                     if (arrayr(irs) .lt. 2.d0) then
                        diff = abs(rs(ij)-arrayr(irs))
                        if (diff .lt. 1.d-12) then
                           if (maxindex-arrayi(irs,1) .eq. 0 .and.
     &                         minindex-arrayi(irs,2) .eq. 0 .and.
     &                         sumindex-arrayi(irs,3) .eq. 0) then
c                              print*,'failed = ',failed
c                              print*,'rs = ',rs(ij),arrayr(irs)
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

               end do

               if (nitr .ge. 1000) then
                  print*,'Failed to choose a sample'
                  stop
               end if
         
               if (nonreplacement) then
                  arrayr(irs) = rs(ij)
                  arrayi(irs,1) = maxindex
                  arrayi(irs,2) = minindex
                  arrayi(irs,3) = sumindex
               end if

c              Calculate the Fisher transformed R value 
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

c           nk loop    
            end do

c           Compute standard deviations of Pearson R, Z, and slope for the sample distribution            
            sigmar =
     &      dsqrt( (rpairsum2 - ravg**2/dble(nk))/dble(nk) )
            sigmaz =
     &      dsqrt((zpairsum2 - zavg**2/dble(nk))/dble(nk))            
            sigma_slope =
     &      dsqrt((slopeavg2 - slopeavg**2/dble(nk))/dble(nk))

c           Compute variances            
            varr = sigmar**2
            var_slope = sigma_slope**2
            varz = sigmaz**2
            
c           average of r average scores
            ravg = ravg/dble(nk)

c           average of z scores
            zavg = zavg/dble(nk)

c           average of slope scores            
            slopeavg = slopeavg/dble(nk)            

c           Compute the area difference from a normal distribution            
            zminz = -3.5d0
            zmaxz =  3.5d0

            dzr = (zmaxz-zminz)/dble(num_intervals)
            
            fdist = 0.d0
            cumfdist = 0.d0
            cumarea = 0.d0
            areadif = 0.d0
c           Use num_intervals intervals
            do i = 1,num_intervals
               sminz = zminz + dble(i-1)*dzr
               smaxz = zminz + dble(i)*dzr
               smid(i) = .5d0*(sminz+smaxz)

               call cdfnor(1,area1,qnor,sminz,0.d0,1.d0,
     &                     status,boundnor)
               call cdfnor(1,area2,qnor,smaxz,0.d0,1.d0,
     &                     status,boundnor)
               area = area2 - area1
               cumarea = cumarea + area
               do ijk = 1,nk
                  zvalue = (zs(ijk)-zavg)/sigmaz
                  if (zvalue >= sminz .and.
     &                zvalue < smaxz) then                  
                     fdist(i) = fdist(i) + 1.d0
                     cumfdist = cumfdist + 1.d0
                  end if
               end do
               areadif = areadif + abs(fdist(i)/dble(nk)-area)
            end do


c           Find area difference between sample and analytical b distribution
c           The variance is adjusted to account for the sample percent of the population      
            ii = 2
            fac = 1.d0

            facq = dble(nnn)/dble(n)
            facq = 1.d0 - facq
            rnnns = dble(nnn-3)/facq + 3.d0
            nnns = nint(rnnns)            
            
c            s = dble(nnns)
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

c            nu = dble(nnnm1s)
            nu = s - 1.d0
            facgnu = facnu/dsqrt(nu*acos(-1.d0))
            expgnu = -(nu + 1.d0)/2.d0

            bot = dsqrt(1.d0 - rsave(ija)**2)

            dbconv = dsqrt(nu)*sdnorx(ija)/(bot*sdnory(ija))

            zminzb = slopemin
            zmaxzb = slopemax

            dzrb = (zmaxzb-zminzb)/dble(num_intervals)

            fdist = 0.d0
            cumfdistb = 0.d0
            cumareab = 0.d0
            areadifb = 0.d0
c           Use num_intervals intervals
            do i = 1,num_intervals
               sminzb = zminzb + dble(i-1)*dzrb
               smaxzb = zminzb + dble(i)*dzrb

               top = (sminzb*sdnorx(ija)/sdnory(ija) - rsave(ija))
               fgnu = top*dsqrt(nu)/bot
               call hypser(0.5d0,(nu+1.d0)/2.d0,1.5d0,-fgnu**2/nu,
     &              hyp,dhyp)
               area1 = .5d0 + fgnu*facgnu*hyp

               top = (smaxzb*sdnorx(ija)/sdnory(ija) - rsave(ija))
               fgnu = top*dsqrt(nu)/bot
               call hypser(0.5d0,(nu+1.d0)/2.d0,1.5d0,-fgnu**2/nu,
     &              hyp,dhyp)
               area2 = .5d0 + fgnu*facgnu*hyp               
               
               area = area2 - area1
               cumareab = cumareab + area
               do ijk = 1,nk
                  if (slopeu(ijk) >= sminzb .and.
     &                slopeu(ijk) < smaxzb) then                  
                     fdist(i) = fdist(i) + 1.d0
                     cumfdistb = cumfdistb + 1.d0
                  end if
               end do
               areadifb = areadifb + abs(fdist(i)/dble(nk)-area)
            end do

c           The Kolmogorov-Smirnov test is currently not used
            dist_ks = 0.d0
c            call ksone(zs,nk,dist_ks,prob,zavg,sigmaz)
            

c           Analytical expectation of analytical slope distribution            
            exp_slope = rsave(ija)*sdnory(ija)/sdnorx(ija)            

c           Variance of analytical slope distribution            
            varb2 = (1.d0 - rsave(ija)**2)/(dble(nnn) - 3.d0)
c            varb2 = varb2*sdnor2**2/sdnor**2
            varb2 = varb2*sdnory(ija)**2/sdnorx(ija)**2            

c           Calculate expectation (expectation) and  variance (var)
c           of analytical Pearson R distribution
            call mct(nnn,rsave(ija),
     &      sdnorx(ija),sdnory(ija),expectation,var)            

c           Mean of Fisher distribution
            zexp = .5d0*log((1.d0+expectation)/(1.d0-expectation))
            zexprho = .5d0*log((1.d0+rsave(ija))/(1.d0-rsave(ija)))
c           Variance of Fisher distribution
            varz_a = 1.d0/dble(nnn-3) 

            
            svarb2 = varb2
            svar = var

c           Sum each metric for each population
            
            d_array(iijj) = d_array(iijj) + dist_ks
            area_array(iijj) = area_array(iijj) + areadif
            areab_array(iijj) = areab_array(iijj) + areadifb            
            z_array(iijj) = z_array(iijj) + (zavg - zexp)
            zrho_array(iijj) = z_array(iijj) + (zavg - zexprho)
            
            varz_array(iijj) = varz_array(iijj) + (varz - varz_a)/varz_a
c            varz_array(iijj) = varz_array(iijj) + varz            
            zq_array(iijj) = zq_array(iijj) + varz

            slope_array_abs(iijj) = slope_array_abs(iijj) +
     &      abs(slopeavg-exp_slope)

            slope_array_abs_max(iijj) = max(slope_array_abs_max(iijj),
     &      (slopeavg-exp_slope))

            slope_array_abs_min(iijj) = min(slope_array_abs_min(iijj),
     &      (slopeavg-exp_slope))                        

            slope_array_dif(iijj) = slope_array_dif(iijj) +
     &      (slopeavg-exp_slope)

            slope_array_value(iijj) = slope_array_value(iijj) +
     &      slopeavg                        
            

            slope_array(iijj) = slope_array(iijj) +
     &      (abs(slopeavg-exp_slope)/exp_slope)

            sigma_slope_array(iijj) = sigma_slope_array(iijj) +
     &      (abs(var_slope-svarb2)/svarb2)

            slope_arrayd(iijj) = slope_arrayd(iijj) +
     &      ((slopeavg-exp_slope)/exp_slope)

            sigma_slope_arrayd(iijj) = sigma_slope_arrayd(iijj) +
     &      ((var_slope-svarb2)/svarb2)            

            
            ravg_array_abs(iijj) = ravg_array_abs(iijj) +
     &      abs(ravg-expectation)

            ravg_array_abs_max(iijj) = max(ravg_array_abs_max(iijj),
     &      (ravg-expectation))

            ravg_array_abs_min(iijj) = min(ravg_array_abs_min(iijj),
     &      (ravg-expectation))            

            ravg_array_dif(iijj) = ravg_array_dif(iijj) +
     &      (ravg-expectation)

            ravg_array_value(iijj) = ravg_array_value(iijj) +
     &      ravg                        
            

            ravg_array(iijj) = ravg_array(iijj) +
     &      abs((ravg-expectation)/expectation)

            sigmar_array(iijj) = sigmar_array(iijj) +
     &      abs((varr-svar)/svar)

            ravg_arrayd(iijj) = ravg_arrayd(iijj) +
     &      (ravg-expectation)/expectation
            
            sigmar_arrayd(iijj) = sigmar_arrayd(iijj) +
     &           (varr-svar)/svar


            slope_array_max(iijj) = max(slope_array_max(iijj),
     &      (abs(slopeavg-exp_slope)/exp_slope))

            sigma_slope_array_max(iijj)=max(sigma_slope_array_max(iijj),
     &      (abs(var_slope-svarb2)/svarb2))

            slope_arrayd_max(iijj) = max(slope_arrayd_max(iijj),
     &      ((slopeavg-exp_slope)/exp_slope))

            sigma_slope_arrayd_max(iijj) =
     &      max(sigma_slope_arrayd_max(iijj),
     &      ((var_slope-svarb2)/svarb2))            

            ravg_array_max(iijj) = max(ravg_array_max(iijj),
     &      abs((ravg-expectation)/expectation))            

            sigmar_array_max(iijj) = max(sigmar_array_max(iijj),
     &      abs((varr-svar)/svar))

            ravg_arrayd_max(iijj) = max(ravg_arrayd_max(iijj),
     &      (ravg-expectation)/expectation)
            
            sigmar_arrayd_max(iijj) = max(sigmar_arrayd_max(iijj),
     &      (varr-svar)/svar)           


            slope_array_min(iijj) = min(slope_array_min(iijj),
     &      (abs(slopeavg-exp_slope)/exp_slope))

            sigma_slope_array_min(iijj)=min(sigma_slope_array_min(iijj),
     &      (abs(var_slope-svarb2)/svarb2))

            slope_arrayd_min(iijj) = min(slope_arrayd_min(iijj),
     &      ((slopeavg-exp_slope)/exp_slope))

            sigma_slope_arrayd_min(iijj) =
     &      min(sigma_slope_arrayd_min(iijj),
     &      ((var_slope-svarb2)/svarb2))            

            ravg_array_min(iijj) = min(ravg_array_min(iijj),
     &      abs((ravg-expectation)/expectation))            

            sigmar_array_min(iijj) = min(sigmar_array_min(iijj),
     &      abs((varr-svar)/svar))

            ravg_arrayd_min(iijj) = min(ravg_arrayd_min(iijj),
     &      (ravg-expectation)/expectation)
            
            sigmar_arrayd_min(iijj) = min(sigmar_arrayd_min(iijj),
     &      (varr-svar)/svar)           
            
c        itimes loop      
         end do
         
         ravg_array(iijj) = ravg_array(iijj)/dble(itimes)         
         sigmar_array(iijj) = sigmar_array(iijj)/dble(itimes)
         slope_array(iijj) = slope_array(iijj)/dble(itimes)
         sigma_slope_array(iijj) = sigma_slope_array(iijj)/dble(itimes)
         slope_arrayd(iijj) = slope_arrayd(iijj)/dble(itimes)
         sigma_slope_arrayd(iijj) =sigma_slope_arrayd(iijj)/dble(itimes)         
         ravg_arrayd(iijj) = ravg_arrayd(iijj)/dble(itimes)
         sigmar_arrayd(iijj) = sigmar_arrayd(iijj)/dble(itimes)
         d_array(iijj) = d_array(iijj)/dble(itimes)
         area_array(iijj) = area_array(iijj)/dble(itimes)
         areab_array(iijj) = areab_array(iijj)/dble(itimes)
         z_array(iijj) = z_array(iijj)/dble(itimes)
         zrho_array(iijj) = zrho_array(iijj)/dble(itimes)
         varz_array(iijj) = varz_array(iijj)/dble(itimes)
         zq_array(iijj) = zq_array(iijj)/dble(itimes)

         ravg_array_abs(iijj) = ravg_array_abs(iijj)/dble(itimes)
         ravg_array_dif(iijj) = ravg_array_dif(iijj)/dble(itimes)
         ravg_array_value(iijj) = ravg_array_value(iijj)/dble(itimes)

         slope_array_abs(iijj) = slope_array_abs(iijj)/dble(itimes)
         slope_array_dif(iijj) = slope_array_dif(iijj)/dble(itimes)
         slope_array_value(iijj) = slope_array_value(iijj)/dble(itimes)
            
         
         aa(1) = ravg_array(iijj)
         aa(2) = sigmar_array(iijj)
         aa(3) = slope_array(iijj)
         aa(4) = sigma_slope_array(iijj)
         aa(5) = ravg_arrayd(iijj)
         aa(6) = sigmar_arrayd(iijj)
         aa(7) = slope_arrayd(iijj)
         aa(8) = sigma_slope_arrayd(iijj)
         aa(9) = d_array(iijj)
         aa(10) = area_array(iijj)
         aa(11) = z_array(iijj)
         aa(12) = zrho_array(iijj)
         aa(13) = varz_array(iijj)
         aa(14) = zq_array(iijj)
         aa(15) = ravg_array_abs(iijj)
         aa(16) = ravg_array_dif(iijj)
         aa(17) = ravg_array_value(iijj)
         aa(18) = slope_array_abs(iijj)
         aa(19) = slope_array_dif(iijj)
         aa(20) = slope_array_value(iijj)
         aa(21) = areab_array(iijj)

         
         call MPI_REDUCE(aa,asum,21,MPI_DOUBLE_PRECISION,
     &                   MPI_SUM,0,MPI_COMM_WORLD,ierror)

         a(1) = ravg_array_max(iijj)
         a(2) = sigmar_array_max(iijj)
         a(3) = slope_array_max(iijj)
         a(4) = sigma_slope_array_max(iijj)
         a(5) = ravg_arrayd_max(iijj)
         a(6) = sigmar_arrayd_max(iijj)
         a(7) = slope_arrayd_max(iijj)
         a(8) = sigma_slope_arrayd_max(iijj)
         a(9) = slope_array_abs_max(iijj)
         a(10) = ravg_array_abs_max(iijj)

         call MPI_REDUCE(a,amax,10,MPI_DOUBLE_PRECISION,
     &        MPI_MAX,0,MPI_COMM_WORLD,ierror)

         a(1) = ravg_array_min(iijj)
         a(2) = sigmar_array_min(iijj)
         a(3) = slope_array_min(iijj)
         a(4) = sigma_slope_array_min(iijj)
         a(5) = ravg_arrayd_min(iijj)
         a(6) = sigmar_arrayd_min(iijj)
         a(7) = slope_arrayd_min(iijj)
         a(8) = sigma_slope_arrayd_min(iijj)
         a(9) = slope_array_abs_min(iijj)
         a(10) = ravg_array_abs_min(iijj)         

         call MPI_REDUCE(a,amin,10,MPI_DOUBLE_PRECISION,
     &                   MPI_MIN,0,MPI_COMM_WORLD,ierror)         

         if (image .eq. 0) then
         
            do i = 1,21
               asum(i) = asum(i)/dble(num_processors)
            end do

            ravg_array(iijj) = asum(1)   
            sigmar_array(iijj) = asum(2)
            slope_array(iijj) = asum(3)
            sigma_slope_array(iijj) = asum(4)
            ravg_arrayd(iijj) = asum(5)
            sigmar_arrayd(iijj) = asum(6)
            slope_arrayd(iijj) = asum(7)
            sigma_slope_arrayd(iijj) = asum(8)
            d_array(iijj) = asum(9)
            area_array(iijj) = asum(10)
            z_array(iijj) = asum(11)
            zrho_array(iijj) = asum(12)
            varz_array(iijj) = asum(13)
            zq_array(iijj) = asum(14)
            ravg_array_abs(iijj) = asum(15)
            ravg_array_dif(iijj) = asum(16)
            ravg_array_value(iijj) = asum(17)
            slope_array_abs(iijj) = asum(18)
            slope_array_dif(iijj) = asum(19)
            slope_array_value(iijj) = asum(20)
            areab_array(iijj) = asum(21)
            
            ravg_array_max(iijj) = amax(1)   
            sigmar_array_max(iijj) = amax(2)
            slope_array_max(iijj) = amax(3)
            sigma_slope_array_max(iijj) = amax(4)
            ravg_arrayd_max(iijj) = amax(5)
            sigmar_arrayd_max(iijj) = amax(6)
            slope_arrayd_max(iijj) = amax(7)
            sigma_slope_arrayd_max(iijj) = amax(8)
            slope_array_abs_max(iijj) = amax(9)
            ravg_array_abs_max(iijj) = amax(10)

            ravg_array_min(iijj) = amin(1)   
            sigmar_array_min(iijj) = amin(2)
            slope_array_min(iijj) = amin(3)
            sigma_slope_array_min(iijj) = amin(4)
            ravg_arrayd_min(iijj) = amin(5)
            sigmar_arrayd_min(iijj) = amin(6)
            slope_arrayd_min(iijj) = amin(7)
            sigma_slope_arrayd_min(iijj) = amin(8)
            slope_array_abs_min(iijj) = amin(9)
            ravg_array_abs_min(iijj) = amin(10)            

            
         end if

         deallocate(xavg,yavg)
      
c     mp loop      
      end do

      if (image .eq. 0) then

         sumr_error = 0.d0
         sumrs_error = 0.d0
         sums_error = 0.d0
         sumss_error = 0.d0
         do iijj = 1,mp
            sumr_error = sumr_error + ravg_array(iijj)
            sumrs_error = sumrs_error + sigmar_array(iijj)
            sums_error = sums_error + sigmar_array(iijj)
            sumss_error = sumss_error + sigma_slope_array(iijj)

            write(83,*) npa(iijj),
     &           ravg_array(iijj),sigmar_array(iijj),
     &           ravg_arrayd(iijj),sigmar_arrayd(iijj),
     &           ravg_array_max(iijj),sigmar_array_max(iijj),
     &           ravg_array_min(iijj),sigmar_array_min(iijj),                        
     &           ravg_arrayd_max(iijj),sigmar_arrayd_max(iijj),
     &           ravg_arrayd_min(iijj),sigmar_arrayd_min(iijj),
     &           standard_dev_r,standard_dev_rs,
     &           d_array(iijj),area_array(iijj),z_array(iijj),
     &           zrho_array(iijj),varz_array(iijj),zq_array(iijj)
            write(84,*) npa(iijj),
     &           slope_array(iijj),sigma_slope_array(iijj),
     &           slope_arrayd(iijj),sigma_slope_arrayd(iijj),
     &           slope_array_max(iijj),sigma_slope_array_max(iijj),
     &           slope_array_min(iijj),sigma_slope_array_min(iijj),                        
     &           slope_arrayd_max(iijj),sigma_slope_arrayd_max(iijj),
     &           slope_arrayd_min(iijj),sigma_slope_arrayd_min(iijj),
     &           standard_dev_s,standard_dev_ss
            write(85,*) dble(np*npa(iijj))/dble(n),
     &           ravg_array(iijj),sigmar_array(iijj),
     &           ravg_arrayd(iijj),sigmar_arrayd(iijj),
     &           ravg_array_max(iijj),sigmar_array_max(iijj),
     &           ravg_array_min(iijj),sigmar_array_min(iijj),                        
     &           ravg_arrayd_max(iijj),sigmar_arrayd_max(iijj),
     &           ravg_arrayd_min(iijj),sigmar_arrayd_min(iijj),
     &           standard_dev_r,standard_dev_rs,
     &           d_array(iijj),area_array(iijj),z_array(iijj),
     &           zrho_array(iijj),varz_array(iijj),zq_array(iijj),
     &           ravg_array_abs(iijj),ravg_array_abs_max(iijj),
     &           ravg_array_abs_min(iijj),ravg_array_dif(iijj),
     &           ravg_array_value(iijj)
            write(86,*) dble(np*npa(iijj))/dble(n),
     &           slope_array(iijj),sigma_slope_array(iijj),
     &           slope_arrayd(iijj),sigma_slope_arrayd(iijj),
     &           slope_array_max(iijj),sigma_slope_array_max(iijj),
     &           slope_array_min(iijj),sigma_slope_array_min(iijj),                        
     &           slope_arrayd_max(iijj),sigma_slope_arrayd_max(iijj),
     &           slope_arrayd_min(iijj),sigma_slope_arrayd_min(iijj),
     &           standard_dev_s,standard_dev_ss,
     &           slope_array_abs(iijj),slope_array_abs_max(iijj),
     &           slope_array_abs_min(iijj),slope_array_dif(iijj),
     &           slope_array_value(iijj),areab_array(iijj)            

         end do   
            

         print*,'r error = ',sumr_error/dble(mp),sumrs_error/dble(mp)
         print*,'s error = ',sums_error/dble(mp),sumss_error/dble(mp)

         if (create_distribution) then
            dzr = (rmaxz - rminz)/dble(num_intervals)

c           Sample distribution of simulated Pearson R            
            fdist = 0.d0
            cumfdist = 0.d0
            sumdist = 0.d0
c           Use num_intervals intervals
            do i = 1,num_intervals
               sminz = rminz + dble(i-1)*dzr
               smaxz = rminz + dble(i)*dzr
               smid(i) = .5d0*(sminz+smaxz)
               do ijk = 1,nk
                  if (rs(ijk) >= sminz .and.
     &                rs(ijk) < smaxz) then                  
                     fdist(i) = fdist(i) + 1.d0
                     cumfdist = cumfdist + 1.d0
                  end if
               end do
               sumdist = sumdist + fdist(i)/dble(nk)
               write(96,*) smid(i),fdist(i)/(dble(dzr*nk)),         
     &                     fdist(i),dble(nk)
            end do

c           Sample distribution of linear regression slope            
            slopemax = slopeavg + 3.5d0*sigma_slope
            slopemin = slopeavg - 3.5d0*sigma_slope
            
            dz = (slopemax - slopemin)/dble(num_intervals)

            fdist = 0.d0
            cumfdist = 0.d0
            sumdist = 0.d0
c           Use num_intervals intervals
            do i = 1,num_intervals
               sminz = slopemin + dble(i-1)*dz
               smaxz = slopemin + dble(i)*dz
               smid(i) = .5d0*(sminz+smaxz)
               do ijk = 1,nk
                  if (slopeu(ijk) >= sminz .and.
     &                slopeu(ijk) < smaxz) then                  
                     fdist(i) = fdist(i) + 1.d0
                     cumfdist = cumfdist + 1.d0
                  end if
               end do
               sumdist = sumdist + fdist(i)/dble(nk)
               write(98,*) smid(i),fdist(i)/(dble(dz*nk)),         
     &                     fdist(i),dble(nk)
            end do

c           Analytical distribution of Pearson R
            pi = acos(-1.d0)
            facdist = 1.d0
            do i = nnn-2,1,-1
               facdist = facdist*dble(i)/(dble(i)+.5d0)
            end do
            facdist = facdist/(.5d0*dsqrt(pi))
            facdist = facdist*dble(nnn-2)/dsqrt(2.d0*pi)
            facdist1 = facdist*(1.d0 - rho**2)**(dble(nnn-1)/2.d0)
            dr = 1.999d0/dble(m)
            rpairavg = 0.d0
            rpair2avg = 0.d0
            do i = 1,m+1
               ra(i) = -1.d0 + dr*dble(i-1) + 1.d-8
               facdist = facdist1*(1.d0 - rho*ra(i))**(-dble(nnn)+1.5d0)
               facdist =facdist*(1.d0 - ra(i)**2)**(dble(nnn)/2.d0-2.d0)
               call hypser(0.5d0,0.5d0,dble(nnn)-0.5d0,
     &                     0.5d0*(1.d0+rho*ra(i)),series,deriv)
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


            ndnp = n/np
c            print*,rho - rho*(1.d0 - rho**2)/(2.d0*dble(ndnp-1))

c           Analytical distribution of linear regression slope
            ii = 2
            fac = 1.d0

            s = dble(nnn)
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

            expectation = 2.d0*fac**2/dble(nnn-1)

            fac = fac*((1.d0-rho**2)**((s-1.d0)/2.d0))
            fac = fac/sqrt(acos(-1.d0))

            expectation = expectation*rho

            call hypser(0.5d0,0.5d0,dble(nnn+1)/2.d0,rho**2,hyp,dhyp)
            expectation = expectation*hyp

            nnnm1 = nnn - 1
            call hypser(1.d0,1.d0,dble(nnnm1)/2.d0+1.d0,rho**2,hyp,dhyp)
            expr2 = 1.d0 - hyp*(1.d0-rho**2)*dble(nnnm1-1)/dble(nnnm1)

            var = expr2 - expectation**2
      
            ija = 1
            sdnorx(1) = 1.d0
            sdnory(1) = 1.d0

            sumb = 0.d0
            do i = 1,num_intervals
               sminz = slopemin + dble(i-1)*dz
               smaxz = slopemin + dble(i)*dz
               smid(i) = .5d0*(sminz+smaxz)
               v = fac*(1.d0 - rho**2 + (rho -
     &              (sdnorx(ija)/sdnory(ija))*smid(i))**2)**(-s/2.d0)
               v = v*sdnorx(1)/sdnory(1)
               write(99,*) smid(i),v
               sumb = sumb + v*dz
            end do

         end if   
      
      end if

      call MPI_FINALIZE(ierror)         

      end

      subroutine mct(nnn,rho,sdnor,sdnor2,expectation,var)
      implicit none
c     Input
      integer nnn
      double precision rho
      double precision sdnor,sdnor2
c     Output
      double precision expectation,var
c     Local
      integer ii,nnnm1
      double precision hyp,dhyp,expr2
      double precision fac,s,ssz,sstop,ssbot
      
      ii = 2
      fac = 1.d0
      s = dble(nnn)
      ssz = s/2.d0

      do while (ssz .gt. 3)
         fac = fac*(s-dble(ii))/(s-dble(ii+1))
         ssz = (s-dble(ii))/2.d0
         ii = ii + 2
      end do

      call GAMMA(ssz,sstop)
      call GAMMA(ssz-.5d0,ssbot)      
      fac = fac*sstop/ssbot
      
      expectation = 2.d0*(fac**2)/dble(nnn-1)      

      fac = fac*((1.d0-rho**2)**((s-1.d0)/2.d0))
      fac = fac/sqrt(acos(-1.d0))
      fac = fac*sdnor/sdnor2

      expectation = expectation*rho

      call hypser(0.5d0,0.5d0,dble(nnn+1)/2.d0,rho**2,hyp,dhyp)
      expectation = expectation*hyp

      nnnm1 = nnn - 1
      call hypser(1.d0,1.d0,dble(nnnm1)/2.d0+1.d0,rho**2,hyp,dhyp)
      expr2 = 1.d0 - hyp*(1.d0-rho**2)*dble(nnnm1-1)/dble(nnnm1)

      var = expr2 - expectation**2

      return
      end
      

      SUBROUTINE hypser(a,b,c,z,series,deriv)
      INTEGER n
      double precision a,b,c,z,series,deriv,aa,bb,cc,fac,temp
c     Returns the hypergeometric series 2F1 and its derivative, iterating to machine accuracy
c     For cabs(z) ? 1/2 convergence is quite rapid.                              
                                                                                  
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
      
      subroutine pearsonr(r,slope,x,y,n,sdx,sdy,yint,xbar,ybar)
c     Computes Pearson R coefficient and linear regression line given
c     a set of original scores
      implicit none

c     Input
c     Length of set of original scores
      integer n
c     Original scores 
      double precision x(n),y(n)

c     Output

c     Standard deviation of x- and y-scores
      double precision sdx,sdy
c     Mean of x- and y-scores
      double precision xbar,ybar
c     y-intercept of linear regression line
      double precision yint
c     slope of linear regression line
      double precision slope

c     Local variables
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

      SUBROUTINE GAMMA(X,GA)

!       ==================================================
!       Purpose: Compute gamma function â(x)
!       Input :  x  --- Argument of â(x)
!                       ( x is not equal to 0,-1,-2,úúú)
!       Output:  GA --- â(x)
!       ==================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO K=2,M1
                 GA=GA*K
              END DO
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO K=1,M
                 R=R*(Z-K)
              END DO
              Z=Z-M
           ELSE
              Z=X
           END IF
           DATA G/1.0D0,0.5772156649015329D0,  
     &           -0.6558780715202538D0, -0.420026350340952D-1, 
     &           0.1665386113822915D0,-.421977345555443D-1, 
     &           -.96219715278770D-2, .72189432466630D-2, 
     &           -.11651675918591D-2, -.2152416741149D-3, 
     &           .1280502823882D-3, -.201348547807D-4, 
     &           -.12504934821D-5, .11330272320D-5, 
     &           -.2056338417D-6, .61160950D-8, 
     &           .50020075D-8, -.11812746D-8, 
     &           .1043427D-9, .77823D-11, 
     &           -.36968D-11, .51D-12, 
     &           -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO K=25,1,-1
              GR=GR*Z+G(K)
           END DO
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
      RETURN
      END
      

      subroutine cumbinomial(x,xx)
c     Computes the cumulative bimodal distribution      
      implicit none
      integer n,i
      double precision aa,bb,dx,summ,a,b
      double precision x,xx,bimodal,eps 

      n = 1000
      aa = -8.d0
      bb = 8.d0
      dx = (bb - aa)/dble(n)
      summ = 0.d0
      i = 1
      xx = aa
      eps = 1.d-12

      do while (summ .lt. x-eps)
         a = aa + dble(i-1)*dx
         b = aa + dble(i)*dx
         summ = summ + .5d0*(bimodal(a)+bimodal(b))*dx
         xx = .5d0*(a+b)
         i = i + 1
      end do

      return
      end

      double precision function bimodal(x)
c     Computes a bimodal normal distribution      
      implicit none
      double precision x
      double precision sigma1,mu1
      double precision sigma2,mu2
      double precision pi
      double precision aa1,aa2
      double precision bi1,bi2

      pi = acos(-1.d0)
        
      mu1 = -1d0
      sigma1 = 0.5d0

      mu2 = 1.d0
      sigma2 = 0.5d0
        
      aa1 = (x - mu1)/sigma1
      aa2 = (x - mu2)/sigma2
        
      bi1 = exp(-.5d0*aa1**2)/(sigma1*sqrt(2.d0*pi))
      bi2 = exp(-.5d0*aa2**2)/(sigma2*sqrt(2.d0*pi))        
        
      bimodal = 0.5d0*(bi1 + bi2)
        
      return
      end


      

      SUBROUTINE cdfnor(which,p,q,x,mean,sd,status,bound)
C**********************************************************************
C
C      SUBROUTINE CDFNOR( WHICH, P, Q, X, MEAN, SD, STATUS, BOUND )
C               Cumulative Distribution Function
C               NORmal distribution
C
C
C                              Function
C
C
C     Calculates any one parameter of the normal
C     distribution given values for the others.
C
C
C                              Arguments
C
C
C     WHICH  --> Integer indicating  which of the  next  parameter
C     values is to be calculated using values  of the others.
C     Legal range: 1..4
C               iwhich = 1 : Calculate P and Q from X,MEAN and SD
C               iwhich = 2 : Calculate X from P,Q,MEAN and SD
C               iwhich = 3 : Calculate MEAN from P,Q,X and SD
C               iwhich = 4 : Calculate SD from P,Q,X and MEAN
C                    INTEGER WHICH
C
C     P <--> The integral from -infinity to X of the normal density.
C            Input range: (0,1].
C                    DOUBLE PRECISION P
C
C     Q <--> 1-P.
C            Input range: (0, 1].
C            P + Q = 1.0.
C                    DOUBLE PRECISION Q
C
C     X < --> Upper limit of integration of the normal-density.
C             Input range: ( -infinity, +infinity)
C                    DOUBLE PRECISION X
C
C     MEAN <--> The mean of the normal density.
C               Input range: (-infinity, +infinity)
C                    DOUBLE PRECISION MEAN
C
C     SD <--> Standard Deviation of the normal density.
C             Input range: (0, +infinity).
C                    DOUBLE PRECISION SD
C
C     STATUS <-- 0 if calculation completed correctly
C               -I if input parameter number I is out of range
C                1 if answer appears to be lower than lowest
C                  search bound
C                2 if answer appears to be higher than greatest
C                  search bound
C                3 if P + Q .ne. 1
C                    INTEGER STATUS
C
C     BOUND <-- Undefined if STATUS is 0
C
C               Bound exceeded by parameter number I if STATUS
C               is negative.
C
C               Lower search bound if STATUS is 1.
C
C               Upper search bound if STATUS is 2.
C
C
C                              Method
C
C
C
C
C     A slightly modified version of ANORM from
C
C     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
C     Package of Special Function Routines and Test Drivers"
C     acm Transactions on Mathematical Software. 19, 22-32.
C
C     is used to calulate the  cumulative standard normal distribution.
C
C     The rational functions from pages  90-95  of Kennedy and Gentle,
C     Statistical  Computing,  Marcel  Dekker, NY,  1980 are  used  as
C     starting values to Newton's Iterations which compute the inverse
C     standard normal.  Therefore no  searches  are necessary for  any
C     parameter.
C
C     For X < -15, the asymptotic expansion for the normal is used  as
C     the starting value in finding the inverse standard normal.
C     This is formula 26.2.12 of Abramowitz and Stegun.
C
C
C                              Note
C
C
C      The normal density is proportional to
C      exp( - 0.5 * (( X - MEAN)/SD)**2)
C
C
C**********************************************************************
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION bound,mean,p,q,sd,x
      INTEGER status,which
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION z,pq
C     ..
C     .. External Functions ..

      DOUBLE PRECISION dinvnr,spmpar
      EXTERNAL dinvnr,spmpar
C     ..
C     .. External Subroutines ..
      EXTERNAL cumnor
C     ..
C     .. Executable Statements ..
C
C     Check arguments
C
      status = 0
      IF (.NOT. ((which.LT.1).OR. (which.GT.4))) GO TO 30
      IF (.NOT. (which.LT.1)) GO TO 10
      bound = 1.0D0
      GO TO 20

   10 bound = 4.0D0
   20 status = -1
      RETURN

   30 IF (which.EQ.1) GO TO 70
C
C     P
C
      IF (.NOT. ((p.LE.0.0D0).OR. (p.GT.1.0D0))) GO TO 60
      IF (.NOT. (p.LE.0.0D0)) GO TO 40
      bound = 0.0D0
      GO TO 50

   40 bound = 1.0D0
   50 status = -2
      RETURN

   60 CONTINUE
   70 IF (which.EQ.1) GO TO 110
C
C     Q
C
      IF (.NOT. ((q.LE.0.0D0).OR. (q.GT.1.0D0))) GO TO 100
      IF (.NOT. (q.LE.0.0D0)) GO TO 80
      bound = 0.0D0
      GO TO 90

   80 bound = 1.0D0
   90 status = -3
      RETURN

  100 CONTINUE
  110 IF (which.EQ.1) GO TO 150
C
C     P + Q
C
      pq = p + q
      IF (.NOT. (abs(((pq)-0.5D0)-0.5D0).GT.
     +    (3.0D0*spmpar(1)))) GO TO 140
      IF (.NOT. (pq.LT.0.0D0)) GO TO 120
      bound = 0.0D0
      GO TO 130

  120 bound = 1.0D0
  130 status = 3
      RETURN

  140 CONTINUE
  150 IF (which.EQ.4) GO TO 170
C
C     SD
C
      IF (.NOT. (sd.LE.0.0D0)) GO TO 160
      bound = 0.0D0
      status = -6
      RETURN

  160 CONTINUE
C
C     Calculate ANSWERS
C
  170 IF ((1).EQ. (which)) THEN
C
C     Computing P
C
          z = (x-mean)/sd
          CALL cumnor(z,p,q)

      ELSE IF ((2).EQ. (which)) THEN
C
C     Computing X
C
          z = dinvnr(p,q)
          x = sd*z + mean

      ELSE IF ((3).EQ. (which)) THEN
C
C     Computing the MEAN
C
          z = dinvnr(p,q)
          mean = x - sd*z

      ELSE IF ((4).EQ. (which)) THEN
C
C     Computing SD
C
          z = dinvnr(p,q)
          sd = (x-mean)/z
      END IF

      RETURN

      END

      SUBROUTINE cumnor(arg,result,ccum)
C**********************************************************************
C
C     SUBROUINE CUMNOR(X,RESULT,CCUM)
C
C
C                              Function
C
C
C     Computes the cumulative  of    the  normal   distribution,   i.e.,
C     the integral from -infinity to x of
C          (1/sqrt(2*pi)) exp(-u*u/2) du
C
C     X --> Upper limit of integration.
C                                        X is DOUBLE PRECISION
C
C     RESULT <-- Cumulative normal distribution.
C                                        RESULT is DOUBLE PRECISION
C
C     CCUM <-- Compliment of Cumulative normal distribution.
C                                        CCUM is DOUBLE PRECISION
C
C
C     Renaming of function ANORM from:
C
C     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
C     Package of Special Function Routines and Test Drivers"
C     acm Transactions on Mathematical Software. 19, 22-32.
C
C     with slight modifications to return ccum and to deal with
C     machine constants.
C
C**********************************************************************
C
C
C Original Comments:
C------------------------------------------------------------------
C
C This function evaluates the normal distribution function:
C
C                              / x
C                     1       |       -t*t/2
C          P(x) = ----------- |      e       dt
C                 sqrt(2 pi)  |
C                             /-oo
C
C   The main computation evaluates near-minimax approximations
C   derived from those in "Rational Chebyshev approximations for
C   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
C   This transportable program uses rational functions that
C   theoretically approximate the normal distribution function to
C   at least 18 significant decimal digits.  The accuracy achieved
C   depends on the arithmetic system, the compiler, the intrinsic
C   functions, and proper selection of the machine-dependent
C   constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants.
C
C   MIN   = smallest machine representable number.
C
C   EPS   = argument below which anorm(x) may be represented by
C           0.5  and above which  x*x  will not underflow.
C           A conservative value is the largest machine number X
C           such that   1.0 + X = 1.0   to machine precision.
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns  ANORM = 0     for  ARG .LE. XLOW.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, EXP
C
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C  Latest modification: March 15, 1992
C
C------------------------------------------------------------------
      INTEGER i
      DOUBLE PRECISION a,arg,b,c,d,del,eps,half,p,one,q,result,sixten,
     +                 temp,sqrpi,thrsh,root32,x,xden,xnum,y,xsq,zero,
     +                 min,ccum
      DIMENSION a(5),b(4),c(9),d(8),p(6),q(5)
C------------------------------------------------------------------
C  External Function
C------------------------------------------------------------------
      DOUBLE PRECISION spmpar
      EXTERNAL spmpar
C------------------------------------------------------------------
C  Mathematical constants
C
C  SQRPI = 1 / sqrt(2*pi), ROOT32 = sqrt(32), and
C  THRSH is the argument for which anorm = 0.75.
C------------------------------------------------------------------
      DATA one,half,zero,sixten/1.0D0,0.5D0,0.0D0,1.60D0/,
     +     sqrpi/3.9894228040143267794D-1/,thrsh/0.66291D0/,
     +     root32/5.656854248D0/
C------------------------------------------------------------------
C  Coefficients for approximation in first interval
C------------------------------------------------------------------
      DATA a/2.2352520354606839287D00,1.6102823106855587881D02,
     +     1.0676894854603709582D03,1.8154981253343561249D04,
     +     6.5682337918207449113D-2/
      DATA b/4.7202581904688241870D01,9.7609855173777669322D02,
     +     1.0260932208618978205D04,4.5507789335026729956D04/
C------------------------------------------------------------------
C  Coefficients for approximation in second interval
C------------------------------------------------------------------
      DATA c/3.9894151208813466764D-1,8.8831497943883759412D00,
     +     9.3506656132177855979D01,5.9727027639480026226D02,
     +     2.4945375852903726711D03,6.8481904505362823326D03,
     +     1.1602651437647350124D04,9.8427148383839780218D03,
     +     1.0765576773720192317D-8/
      DATA d/2.2266688044328115691D01,2.3538790178262499861D02,
     +     1.5193775994075548050D03,6.4855582982667607550D03,
     +     1.8615571640885098091D04,3.4900952721145977266D04,
     +     3.8912003286093271411D04,1.9685429676859990727D04/
C------------------------------------------------------------------
C  Coefficients for approximation in third interval
C------------------------------------------------------------------
      DATA p/2.1589853405795699D-1,1.274011611602473639D-1,
     +     2.2235277870649807D-2,1.421619193227893466D-3,
     +     2.9112874951168792D-5,2.307344176494017303D-2/
      DATA q/1.28426009614491121D00,4.68238212480865118D-1,
     +     6.59881378689285515D-2,3.78239633202758244D-3,
     +     7.29751555083966205D-5/
C------------------------------------------------------------------
C  Machine dependent constants
C------------------------------------------------------------------
      eps = spmpar(1)*0.5D0
      min = spmpar(2)
C------------------------------------------------------------------
      x = arg
      y = abs(x)
      IF (y.LE.thrsh) THEN
C------------------------------------------------------------------
C  Evaluate  anorm  for  |X| <= 0.66291
C------------------------------------------------------------------
          xsq = zero
          IF (y.GT.eps) xsq = x*x
          xnum = a(5)*xsq
          xden = xsq
          DO 10 i = 1,3
              xnum = (xnum+a(i))*xsq
              xden = (xden+b(i))*xsq
   10     CONTINUE
          result = x* (xnum+a(4))/ (xden+b(4))
          temp = result
          result = half + temp
          ccum = half - temp
C------------------------------------------------------------------
C  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
C------------------------------------------------------------------
      ELSE IF (y.LE.root32) THEN
          xnum = c(9)*y
          xden = y
          DO 20 i = 1,7
              xnum = (xnum+c(i))*y
              xden = (xden+d(i))*y
   20     CONTINUE
          result = (xnum+c(8))/ (xden+d(8))
          xsq = aint(y*sixten)/sixten
          del = (y-xsq)* (y+xsq)
          result = exp(-xsq*xsq*half)*exp(-del*half)*result
          ccum = one - result
          IF (x.GT.zero) THEN
              temp = result
              result = ccum
              ccum = temp
          END IF
C------------------------------------------------------------------
C  Evaluate  anorm  for |X| > sqrt(32)
C------------------------------------------------------------------
      ELSE
          result = zero
          xsq = one/ (x*x)
          xnum = p(6)*xsq
          xden = xsq
          DO 30 i = 1,4
              xnum = (xnum+p(i))*xsq
              xden = (xden+q(i))*xsq
   30     CONTINUE
          result = xsq* (xnum+p(5))/ (xden+q(5))
          result = (sqrpi-result)/y
          xsq = aint(x*sixten)/sixten
          del = (x-xsq)* (x+xsq)
          result = exp(-xsq*xsq*half)*exp(-del*half)*result
          ccum = one - result
          IF (x.GT.zero) THEN
              temp = result
              result = ccum
              ccum = temp
          END IF

      END IF

      IF (result.LT.min) result = 0.0D0
      IF (ccum.LT.min) ccum = 0.0D0
C------------------------------------------------------------------
C  Fix up for negative argument, erf, etc.
C------------------------------------------------------------------
C----------Last card of ANORM ----------
      END
      

      DOUBLE PRECISION FUNCTION spmpar(i)
C-----------------------------------------------------------------------                                   
C                                                                                                          
C     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR                                           
C     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT                                             
C     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE                                          
C     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND                                       
C     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN                                           
C                                                                                                          
C        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,                                                    
C                                                                                                          
C        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,                                                
C                                                                                                          
C        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.                                         
C                                                                                                          
C-----------------------------------------------------------------------                                   
C     WRITTEN BY                                                                                           
C        ALFRED H. MORRIS, JR.                                                                             
C        NAVAL SURFACE WARFARE CENTER                                                                      
C        DAHLGREN VIRGINIA                                                                                 
C-----------------------------------------------------------------------                                   
C-----------------------------------------------------------------------                                   
C     MODIFIED BY BARRY W. BROWN TO RETURN DOUBLE PRECISION MACHINE                                        
C     CONSTANTS FOR THE COMPUTER BEING USED.  THIS MODIFICATION WAS                                        
C     MADE AS PART OF CONVERTING BRATIO TO DOUBLE PRECISION                                                
C-----------------------------------------------------------------------                                   
C     .. Scalar Arguments ..                                                                               
      INTEGER i
C     ..                                                                                                   
C     .. Local Scalars ..                                                                                  
      DOUBLE PRECISION b,binv,bm1,one,w,z
      INTEGER emax,emin,ibeta,m
C     ..                                                                                                   
C     .. External Functions ..                                                                             
      INTEGER ipmpar
      EXTERNAL ipmpar
C     ..                                                                                                   
C     .. Intrinsic Functions ..                                                                            
      INTRINSIC dble
C     .. Executable Statements ..                                                                          
C                                                                                                          
      IF (i.GT.1) GO TO 10
      b = ipmpar(4)
      m = ipmpar(8)
      spmpar = b** (1-m)
      RETURN
C                                                                                                          
 10    IF (i.GT.2) GO TO 20
      b = ipmpar(4)
      emin = ipmpar(9)
      one = dble(1)
      binv = one/b
      w = b** (emin+2)
      spmpar = ((w*binv)*binv)*binv
      RETURN
C                                                                                                          
 20    ibeta = ipmpar(4)
      m = ipmpar(8)
      emax = ipmpar(10)
C                                                                                                          
      b = ibeta
      bm1 = ibeta - 1
      one = dble(1)
      z = b** (m-1)
      w = ((z-one)*b+bm1)/ (b*z)
C                                                                                                          
      z = b** (emax-2)
      spmpar = ((w*z)*b)*b
      RETURN

      END
      

      DOUBLE PRECISION FUNCTION dinvnr(p,q)
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION DINVNR(P,Q)
C     Double precision NoRmal distribution INVerse
C
C
C                              Function
C
C
C     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
C     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
C
C
C                              Arguments
C
C
C     P --> The probability whose normal deviate is sought.
C                    P is DOUBLE PRECISION
C
C     Q --> 1-P
C                    P is DOUBLE PRECISION
C
C
C                              Method
C
C
C     The  rational   function   on  page 95    of Kennedy  and  Gentle,
C     Statistical Computing, Marcel Dekker, NY , 1980 is used as a start
C     value for the Newton method of finding roots.
C
C
C                              Note
C
C
C     If P or Q .lt. machine EPS returns +/- DINVNR(EPS)
C
C**********************************************************************
C     .. Parameters ..
      INTEGER maxit
      PARAMETER (maxit=100)
      DOUBLE PRECISION eps
      PARAMETER (eps=1.0D-13)
      DOUBLE PRECISION r2pi
      PARAMETER (r2pi=0.3989422804014326D0)
      DOUBLE PRECISION nhalf
      PARAMETER (nhalf=-0.5D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION p,q
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION strtx,xcur,cum,ccum,pp,dx
      INTEGER i
      LOGICAL qporq
C     ..
C     .. External Functions ..
      DOUBLE PRECISION stvaln
      EXTERNAL stvaln
C     ..
C     .. External Subroutines ..
      EXTERNAL cumnor
C     ..
C     .. Statement Functions ..
      DOUBLE PRECISION dennor,x

      dennor(x) = r2pi*exp(nhalf*x*x)
C     ..
C     .. Executable Statements ..
C
C     FIND MINIMUM OF P AND Q
C
      qporq = p .LE. q
      IF (.NOT. (qporq)) GO TO 10
      pp = p
      GO TO 20

   10 pp = q
C
C     INITIALIZATION STEP
C
   20 strtx = stvaln(pp)
      xcur = strtx
C
C     NEWTON INTERATIONS
C
      DO 30,i = 1,maxit
          CALL cumnor(xcur,cum,ccum)
          dx = (cum-pp)/dennor(xcur)
          xcur = xcur - dx
          IF (abs(dx/xcur).LT.eps) GO TO 40
   30 CONTINUE
      dinvnr = strtx
C
C     IF WE GET HERE, NEWTON HAS FAILED
C
      IF (.NOT.qporq) dinvnr = -dinvnr
      RETURN
C
C     IF WE GET HERE, NEWTON HAS SUCCEDED
C
   40 dinvnr = xcur
      IF (.NOT.qporq) dinvnr = -dinvnr
      RETURN

      END
      
      INTEGER FUNCTION ipmpar(i)
C-----------------------------------------------------------------------                             
C                                                                                                    
C     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER                                 
C     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER                                  
C     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...                                     
C                                                                                                    
C  INTEGERS.                                                                                         
C                                                                                                    
C     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM                                    
C                                                                                                    
C               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )                                       
C                                                                                                    
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.                                            
C                                                                                                    
C     IPMPAR(1) = A, THE BASE.                                                                       
C                                                                                                    
C     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS.                                                    
C                                                                                                    
C     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE.                                                   
C                                                                                                    
C  FLOATING-POINT NUMBERS.                                                                           
C                                                                                                    
C     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING                                    
C     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE                                      
C     NONZERO NUMBERS ARE REPRESENTED IN THE FORM                                                    
C                                                                                                    
C               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)                                             
C                                                                                                    
C               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,                                              
C               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.                                              
C                                                                                                    
C     IPMPAR(4) = B, THE BASE.                                                                       
C                                                                                                    
C  SINGLE-PRECISION                                                                                  
C                                                                                                    
C     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.                                                    
C                                                                                                    
C     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.                                                     
C                                                                                                    
C     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.                                     

C  DOUBLE-PRECISION                                                                                  
C                                                                                                    
C     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.                                                    
C                                                                                                    
C     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.                                                     
C                                                                                                    
C     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.                                                     
C                                                                                                    
C-----------------------------------------------------------------------                             
C                                                                                                    
C     TO DEFINE THIS FUNCTION FOR THE COMPUTER BEING USED, ACTIVATE                                  
C     THE DATA STATMENTS FOR THE COMPUTER BY REMOVING THE C FROM                                     
C     COLUMN 1. (ALL THE OTHER DATA STATEMENTS SHOULD HAVE C IN                                      
C     COLUMN 1.)                                                                                     
C                                                                                                    
C-----------------------------------------------------------------------                             
C                                                                                                    
C     IPMPAR IS AN ADAPTATION OF THE FUNCTION I1MACH, WRITTEN BY                                     
C     P.A. FOX, A.D. HALL, AND N.L. SCHRYER (BELL LABORATORIES).                                     
C     IPMPAR WAS FORMED BY A.H. MORRIS (NSWC). THE CONSTANTS ARE                                     
C     FROM BELL LABORATORIES, NSWC, AND OTHER SOURCES.                                               
C                                                                                                    
C-----------------------------------------------------------------------                             
C     .. Scalar Arguments ..                                                                         
      INTEGER i
C     ..                                                                                             
C     .. Local Arrays ..                                                                             
      INTEGER imach(10)
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T                               
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T                                  
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).                                   
C                                                                                                    
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

      DOUBLE PRECISION FUNCTION stvaln(p)
C                                                                                                          
C**********************************************************************                                    
C                                                                                                          
C     DOUBLE PRECISION FUNCTION STVALN(P)                                                                  
C                    STarting VALue for Neton-Raphon                                                       
C                calculation of Normal distribution Inverse                                                
C                                                                                                          
C                                                                                                          
C                              Function                                                                    
C                                                                                                          
C                                                                                                          
C     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -                                   
C     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P                                                  
C                                                                                                          
C                                                                                                          
C                              Arguments                                                                   
C                                                                                                          
C                                                                                                          
C     P --> The probability whose normal deviate is sought.                                                
C                    P is DOUBLE PRECISION                                                                 
C                                                                                                          
C                                                                                                          
C                              Method                                                                      
C                                                                                                          
C                                                                                                          
C     The  rational   function   on  page 95    of Kennedy  and  Gentle,                                   
C     Statistical Computing, Marcel Dekker, NY , 1980.                                                     
C                                                                                                          
C**********************************************************************                                    
C                                                                                                          
C     .. Scalar Arguments ..                                                                               
      DOUBLE PRECISION p
C     ..                                                                                                   
C     .. Local Scalars ..                                                                                  
      DOUBLE PRECISION sign,y,z
C     ..                                                                                                   
C     .. Local Arrays ..                                                                                   
      DOUBLE PRECISION xden(5),xnum(5)
C     ..                                                                                                   
C     .. External Functions ..                                       
      DOUBLE PRECISION devlpl
      EXTERNAL devlpl
C     ..                                                                                                   
C     .. Intrinsic Functions ..                                                                            
      INTRINSIC dble,log,sqrt
C     ..                                                                                                   
C     .. Data statements ..                                                                                
      DATA xnum/-0.322232431088D0,-1.000000000000D0,-0.342242088547D0,
     +     -0.204231210245D-1,-0.453642210148D-4/
      DATA xden/0.993484626060D-1,0.588581570495D0,0.531103462366D0,
     +     0.103537752850D0,0.38560700634D-2/
C     ..                                                                                                   
C     .. Executable Statements ..                                                                          
      IF (.NOT. (p.LE.0.5D0)) GO TO 10
      sign = -1.0D0
      z = p
      GO TO 20

 10    sign = 1.0D0
      z = 1.0D0 - p
 20    y = sqrt(-2.0D0*log(z))
      stvaln = y + devlpl(xnum,5,y)/devlpl(xden,5,y)
      stvaln = sign*stvaln
      RETURN

      END
      
      DOUBLE PRECISION FUNCTION devlpl(a,n,x)
C**********************************************************************                              
C                                                                                                    
C     DOUBLE PRECISION FUNCTION DEVLPL(A,N,X)                                                        
C              Double precision EVALuate a PoLynomial at X                                           
C                                                                                                    
C                                                                                                    
C                              Function                                                              
C                                                                                                    
C                                                                                                    
C     returns                                                                                        
C          A(1) + A(2)*X + ... + A(N)*X**(N-1)                                                       
C                                                                                                    
C                                                                                                    
C                              Arguments                                                             
C                                                                                                    
C                                                                                                    
C     A --> Array of coefficients of the polynomial.                                                 
C                                        A is DOUBLE PRECISION(N)                                    
C                                                                                                    
C     N --> Length of A, also degree of polynomial - 1.                                              
C                                        N is INTEGER                                                
C                                                                                                    
C     X --> Point at which the polynomial is to be evaluated.                                        
C                                        X is DOUBLE PRECISION                                       
C                                                                                                    
C**********************************************************************                              
C                                                                                                    
C     .. Scalar Arguments ..                                                                         
      DOUBLE PRECISION x
      INTEGER n
C     ..                                                                                             
C     .. Array Arguments ..                                                                          
      DOUBLE PRECISION a(n)
C     ..                                                                                             
C     .. Local Scalars ..                                                                            
      DOUBLE PRECISION term
      INTEGER i
C     ..                                                                                             
C     .. Executable Statements ..                                                                    
      term = a(n)
      DO 10,i = n - 1,1,-1
          term = a(i) + term*x
 10   CONTINUE
      devlpl = term
      RETURN

      END
