      program correlation
      implicit none
      include 'mpif.h'

      integer ss,nnnm1,ndnp
      logical create_distribution
      double precision varhyper,dzp,s,fac,ssz
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
      double precision cdf(num_intervals),trap(num_intervals)
      integer num_processors,image,ierror
      
c     Number of points in analytical distribution
      integer nnn,ndim,mp,imultiple
      
      character*4 strz      
      logical failed,partition
      
      integer i,k,nk,iii,jjj,ii,nitr,iijj,iijja
      integer n,nit,np,j,p,itimes,irs,ija
      integer ij,jj,iiata,itimesnum
      integer iiat,jjk,jjk_save,ijb
      integer maxindex,minindex,sumindex
      integer maxisum,minisum,minnp
      
      real rrr      
      double precision expectation,var,zexp,zexprho
      double precision rminza,rmaxza,zminz,zmaxz
      double precision area1,area2,area,zvalue
      double precision areadif,areadif_sum,areadifb_sum
      double precision slopemin,slopemax
      double precision meanz,stz
      double precision xmsum,xms,ymsum,yms
      double precision eps,rn,prob,dist_ks
      double precision dist_ks_sum
      double precision xmm,ymm
      double precision svarb2,svar
      double precision varr,var_slope
      double precision varz,varz_sum,varz_a
c      double precision a(8),asum(8),
      double precision amax(10),amin(10)
      double precision a(10),aa(21),asum(21)
      double precision sumintegral
      
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
      double precision rminz,rmaxz,r2
c     rho from bivariate distribution
      double precision rho
      double precision sumrs_error,sumss_error

      double precision nu,facgnu,expgnu,top,bot,fgnu
      double precision facnu,dbconv,sumb

      double precision zminzb,zmaxzb,dzrb
      double precision fdistb,cumfdistb,cumareab
      double precision areadifb
            
      integer kja,niterations
      double precision sigmar_sum,sigma_slope_sum
      double precision sigma_slope_sum2,zavg_sum
      double precision ravg_sum,slopeavg_sum
      double precision sigmaz_sum,sigmaz_sum2
c      double precision slope_sd
      
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
      
      double precision, allocatable :: rerror(:,:),rserror(:,:)
      double precision, allocatable :: slopea(:,:),sslope(:,:)
      double precision, allocatable :: rsum(:),rsum2(:)
      double precision, allocatable :: rssum(:),rssum2(:)
      double precision, allocatable :: slopesum(:),slopesum2(:)
      double precision, allocatable :: sslopesum(:),sslopesum2(:)
      double precision, allocatable :: standard_dev_r(:)
      double precision, allocatable :: standard_dev_rs(:)
      double precision, allocatable :: standard_dev_s(:)
      double precision, allocatable :: standard_dev_ss(:)
         
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

      double precision, allocatable :: xfull(:,:),yfull(:),solfull(:)
      double precision, allocatable :: yhatfull(:)

      double precision, allocatable :: rsave(:),sdnorx(:),sdnory(:)
      double precision, allocatable :: all(:),allp(:)
      
      call MPI_INIT(ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, num_processors, ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, image, ierror)

      open(unit=9,file="input_swr")
      if (image .eq. 0) then
         open(unit=83,file="rnp")
         open(unit=84,file="snp")
         open(unit=85,file="re")
         open(unit=86,file="se")
         open(unit=87,file="all")
         open(unit=88,file="input_process")
         open(unit=96,file="rdistsim")
         open(unit=97,file="rdista")
         open(unit=98,file="sdistsim")
         open(unit=99,file="sdista")
         open(unit=107,file="bdist")
      end if

      print*,'image = ',image
c     n: length of original scores
c     np: size of the group
c     nnn : sample size
c     nk: number of iterations
c     nit: number of samples

      read(9,*) strz,n
      read(9,*) strz,np
      read(9,*) strz,nk
      read(9,*) strz,nit
      read(9,*) strz,rho
      read(9,*) strz,rminz
      read(9,*) strz,rmaxz
      read(9,*) strz,mp

      niterations = nit
      itimes = 112/num_processors

      itimesnum = itimes*num_processors

      
      allocate(all(itimes*4),allp(itimes*4*num_processors))
      partition = .true.
      create_distribution = .false.
      allocate(arrayr(-1000002:1000002))
      allocate(arrayi(-1000002:1000002,3))

      eps = 1.d-14

      write(88,*) itimes
      write(88,*) mp
      write(88,*) num_processors
      
      if (image .eq. 0) then
         print*,'n = ',n,'np = ',np,'nk = ',nk,'nit = ',nit,
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
      
      allocate(rerror(itimes*num_processors,mp))
      allocate(rserror(itimes*num_processors,mp))
      allocate(slopea(itimes*num_processors,mp))
      allocate(sslope(itimes*num_processors,mp))
      allocate(rsum(mp))
      allocate(rsum2(mp))
      allocate(rssum(mp))
      allocate(rssum2(mp))
      allocate(slopesum(mp))
      allocate(slopesum2(mp))
      allocate(sslopesum(mp))
      allocate(sslopesum2(mp))
      allocate(standard_dev_r(mp))
      allocate(standard_dev_rs(mp))
      allocate(standard_dev_s(mp))
      allocate(standard_dev_ss(mp))
      
      ndim = 1
      allocate(xfull(n,ndim),yfull(n),solfull(ndim),yhatfull(n))

      allocate(xm(n),ym(n))
      allocate(xma(n,itimes),yma(n,itimes))

      meanz = 0.d0
      stz = 1.d0

      meannor = 0.d0
      sdnor = 1.d0

      meannor2 = 0.d0
      sdnor2 = 1.0d0

      sumrr = 0.d0
      do ii = 1,image+20
         rr = rand()
      end do
      if (sumrr > 200.d0) then
         stop
      end if
      
      do ii = 1,itimes
              
         diff = 1.d0
         nitdiff = 0

c        Iterate until rho and R_individual agree to within .0005
         do while (diff .gt. .01 .and. nitdiff .lt. 100000)
            nitdiff = nitdiff + 1

            xmsum = 0.d0
            xms = 0.d0
            ymsum = 0.d0
            yms = 0.d0
            do j = 1,n
c              create a random number xm(i) sampled from a standard normal distribution      
               pnor = min(max(rand(),eps),1.d0-eps)
               qnor = 1.d0 - pnor

c               xmm = pnor
c               xm(j) = pnor
c               call cdfnor(2,pnor,qnor,xmm,meanz,stz,
c     &                     status,boundnor)
               call cumbinomial(pnor,xmm)
                       
c               xm(j) = sdnor*xmm + meannor

c               xmm = -log(1.d0 - pnor)
               xm(j) = xmm
            
               xmsum = xmsum + xm(j)
               xms = xms + xm(j)**2

               pnor = min(max(rand(),eps),1.d0-eps)
               qnor = 1.d0 - pnor

c               ymm = -log(1.d0 - pnor)

               call cdfnor(2,pnor,qnor,ymm,meanz,stz,
     &                     status,boundnor)
c               call cumbinomial(pnor,ymm)
                       
c              correlate xm(j) and ym(j) at level of correlation rho
c              ym(j) = rho*xm(j) + dsqrt(1.d0 - rho**2)*ym(j)

               ym(j) = sdnor2*(rho*xmm + dsqrt(1.d0-rho**2)*ymm)
     &               + meannor2

               ymsum = ymsum + ym(j)
               yms = yms + ym(j)**2

            end do

            sdnort = dsqrt((xms - xmsum**2/dble(n))/dble(n))
            sdnor2t = dsqrt((yms - ymsum**2/dble(n))/dble(n))
            
            do i = 1,n
               do p = 1,ndim
                  xfull(i,p) = xm(i)
               end do
            end do

            do i = 1,n
               yfull(i) = ym(i)
            end do

            do i = 1,n
               xma(i,ii) = xm(i)
               yma(i,ii) = ym(i)
            end do

            call multiple(xfull,yfull,solfull,yhatfull,r2,n,ndim) 

            call pearsonr(r,slope,
     &                    xm,ym,n,sdx,sdy,yint,xbar,ybar)

            diff = abs(r**2-rho**2)
            diff = abs(r**2 - .7d0**2)

         end do
         rsave(ii) = r
         print*,'rsave(',ii,') = ',rsave(ii),image
         sdnorx(ii) = sdnort
         sdnory(ii) = sdnor2t
         
c         rsave(ii) = rho

         if (nitdiff .ge. 100000) then
            print*,'STOP ERROR'
            stop
         end if
         
c         print*,'rsave(ii) = ',ii,rsave(ii),rsave(ii)-rho,
c     &   sdnorx(ii)-sdnor,sdnory(ii)-sdnor2

      end do


      imultiple = n/(mp*np)
      do iijj = 1,mp
         npa(iijj) = imultiple*iijj
      end do

      allocate(rs(nk),slopeu(nk),zs(nk))      
      allocate(iia(np))
      
c     do iijj = 1,mp
      do iijj = 1,mp

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
         

         do ija = 1,itimes

            varz_sum = 0.d0
            sigmar_sum = 0.d0
            sigmaz_sum = 0.d0
            sigmaz_sum2 = 0.d0
            sigma_slope_sum = 0.d0
            ravg_sum = 0.d0
            zavg_sum = 0.d0
            slopeavg_sum = 0.d0
            dist_ks_sum = 0.d0
            areadif_sum = 0.d0
            areadifb_sum = 0.d0

            arrayr = 3.d0
            arrayi = -1
            
            do kja = 1,niterations
            
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
         
            do ij = 1,nk

               if (partition) then
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
                  
                  if (.not. partition) then
c                    Choose nsamples of length nnn from vector                     
                     do i = 1,nnn

                        xavg(i) = 0.d0
                        yavg(i) = 0.d0

c                        do ii = 1,np
c
c                           pnor = min(max(rand(),eps),1.d0-eps)
c                           qnor = 1.d0 - pnor
c
c                           call cdfnor(2,pnor,qnor,xmm,meanz,stz,
c     &                                 status,boundnor)
c                       
c                           pnor = min(max(rand(),eps),1.d0-eps)
c                           qnor = 1.d0 - pnor
c
c                           call cdfnor(2,pnor,qnor,ymm,meanz,stz,
c     &                                 status,boundnor)
c
c                           ymm =
c     &                     sdnor2*(rho*xmm + dsqrt(1.d0-rho**2)*ymm)
c     &                     + meannor2
c
c                           xavg(i) = xavg(i) + xmm
c                           yavg(i) = yavg(i) + ymm
c                        end do


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
c                    Choose nsamples of length nnn from vector
                     do j = 1,nnn

                        maxisum = 0
                        minisum = 0
                        minnp = 100000000
c                       Choose np random unique elements from vector
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

                  call pearsonr(rs(ij),slopeu(ij),
     &                 xavg,yavg,nnn,sdx,sdy,yint,xbar,ybar)
                  
                  if (partition) then
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
         
               if (partition) then
                  arrayr(irs) = rs(ij)
                  arrayi(irs,1) = maxindex
                  arrayi(irs,2) = minindex
                  arrayi(irs,3) = sumindex
               end if

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

c            if (image .eq. 0) then
c               print*,'nitr itimes = ',nitr,ija
c            end if
            
            sigmar =
     &      dsqrt( (rpairsum2 - ravg**2/dble(nk))/dble(nk) )
            sigmaz =
     &      dsqrt((zpairsum2 - zavg**2/dble(nk))/dble(nk))            
            sigma_slope =
     &      dsqrt((slopeavg2 - slopeavg**2/dble(nk))/dble(nk))

c            sigmar = sigmar**2
c            sigma_slope = sigma_slope**2

            varr = sigmar**2
            var_slope = sigma_slope**2
            varz = sigmaz**2
            
c           average of r average scores
            ravg = ravg/dble(nk)

c           average of z scores
            zavg = zavg/dble(nk)

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

c            sumintegral = 0.d0
c            do i = 1,num_intervals
c               sminz = slopemin + dble(i-1)*dz
c               smaxz = slopemin + dble(i)*dz
c               smid(i) = .5d0*(sminz+smaxz)
c               top = (smid(i)*sdnorx(ija)/sdnory(ija) - rsave(ija))
c               fgnu = top*dsqrt(nu)/bot
c               v = dbconv*facgnu*(1.d0 + fgnu**2/nu)**expgnu
c
c               sumintegral = sumintegral + v*dz
c               call hypser(0.5d0,(nu+1.d0)/2.d0,1.5d0,-fgnu**2/nu,
c     &                     hyp,dhyp)
c               trap(i) = sumintegral
c               cdf(i) = .5d0 + fgnu*facgnu*hyp
cc               v = facgnu*(1.d0 + smid(i)**2/dsqrt(nu))**expgnu
c               write(107,*) smid(i),v
c            end do
c            do i = 1,num_intervals
c               write(108,*) smid(i),trap(i),cdf(i)
c            end do

            
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
c            print*,'cumareab cumfdistb areadifb = ',
c     &      cumareab,cumfdistb/dble(nk),areadifb
            
            call ksone(zs,nk,dist_ks,prob,zavg,sigmaz)
            
            slopeavg = slopeavg/dble(nk)
            
            sigmar_sum = sigmar_sum + varr
            sigmaz_sum = sigmaz_sum + sigmaz
            varz_sum = varz_sum + varz
            sigma_slope_sum = sigma_slope_sum + var_slope
            ravg_sum = ravg_sum + ravg
            zavg_sum = zavg_sum + zavg
            slopeavg_sum = slopeavg_sum + slopeavg
            dist_ks_sum = dist_ks_sum + dist_ks
            areadif_sum = areadif_sum + areadif
            areadifb_sum = areadifb_sum + areadifb
            
            end do

            varr = sigmar_sum/dble(niterations)
            sigmaz = sigmaz_sum/dble(niterations)
            varz = varz_sum/dble(niterations)
            var_slope = sigma_slope_sum/dble(niterations)
            ravg = ravg_sum/dble(niterations)
            zavg = zavg_sum/dble(niterations)
            slopeavg = slopeavg_sum/dble(niterations)
            dist_ks = dist_ks_sum/dble(niterations)
            areadif = areadif_sum/dble(niterations)
            areadifb = areadifb_sum/dble(niterations)
            
            exp_slope = rsave(ija)*sdnory(ija)/sdnorx(ija)            
            
            varb2 = (1.d0 - rsave(ija)**2)/(dble(nnn) - 3.d0)
c            varb2 = varb2*sdnor2**2/sdnor**2
            varb2 = varb2*sdnory(ija)**2/sdnorx(ija)**2            
            
c            call mct(nnn,rsave(ija),
c     &           sdnor,sdnor2,expectation,var)

            call mct(nnn,rsave(ija),
     &      sdnorx(ija),sdnory(ija),expectation,var)            


            zexp = .5d0*log((1.d0+expectation)/(1.d0-expectation))
            zexprho = .5d0*log((1.d0+rsave(ija))/(1.d0-rsave(ija)))
c            svarb2 = sqrt(varb2)
c            svar = sqrt(var)

            svarb2 = varb2
            svar = var
            varz_a = 1.d0/dble(nnn-3) 

            all(1 + 4*(ija-1)) = (ravg - expectation)/expectation
            all(2 + 4*(ija-1)) = (varr - svar)/svar
            all(3 + 4*(ija-1)) = (slopeavg - exp_slope)/exp_slope
            all(4 + 4*(ija-1)) = (var_slope - svarb2)/svarb2
c            all(5 + 6*(ija-1)) = (zavg - zexp)/zexp
c            all(6 + 6*(ija-1)) = sigmaz_sd

            d_array(iijj) = d_array(iijj) + dist_ks
            area_array(iijj) = area_array(iijj) + areadif
            areab_array(iijj) = areab_array(iijj) + areadifb            
c            z_array(iijj) = z_array(iijj) + (zavg - zexp)/zexp
c            zrho_array(iijj) = z_array(iijj) + (zavg - zexprho)/zexprho
            z_array(iijj) = z_array(iijj) + (zavg - zexp)
            zrho_array(iijj) = z_array(iijj) + (zavg - zexprho)
            
            varz_array(iijj) = varz_array(iijj) + (varz - varz_a)/varz_a
c            varz_array(iijj) = varz_array(iijj) + varz            
            zq_array(iijj) = zq_array(iijj) + varz
c           (zavg - zexp)/zexp

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
         
         call MPI_GATHER(all,itimes*4,MPI_DOUBLE_PRECISION,allp,
     &   itimes*4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)

         if (image .eq. 0) then
            write(87,*) npa(iijj),dble(np*npa(iijj))/dble(n)
            do iijja = 1,num_processors*itimes*4
               write(87,*) allp(iijja)
            end do

            rsum(iijj) = 0.d0
            rsum2(iijj) = 0.d0
            rssum(iijj) = 0.d0
            rssum2(iijj) = 0.d0
            slopesum(iijj) = 0.d0
            slopesum2(iijj) = 0.d0
            sslopesum(iijj) = 0.d0
            sslopesum2(iijj) = 0.d0
         
            ij = 0
            ija = 0
            do k = 1,num_processors
               do i = 1,itimes
                  ija = ija + 1
                  ij = ij + 1
                  rerror(ija,iijj) = allp(ij)
                  ij = ij + 1
                  rserror(ija,iijj) = allp(ij)
                  ij = ij + 1
                  slopea(ija,iijj) = allp(ij)
                  ij = ij + 1
                  sslope(ija,iijj) = allp(ij)

                  rsum(iijj) = rsum(iijj) + rerror(ija,iijj)
                  rsum2(iijj) = rsum2(iijj) + rerror(ija,iijj)**2
                  rssum(iijj) = rssum(iijj) + rserror(ija,iijj)
                  rssum2(iijj) = rssum2(iijj) + rserror(ija,iijj)**2
                  slopesum(iijj) = slopesum(iijj) + slopea(ija,iijj)
                  slopesum2(iijj) = slopesum2(iijj)+ slopea(ija,iijj)**2
                  sslopesum(iijj) = sslopesum(iijj) + sslope(ija,iijj)
                  sslopesum2(iijj)=sslopesum2(iijj)+ sslope(ija,iijj)**2
               end do
            end do
            
            standard_dev_r(iijj) =
     &      rsum2(iijj) - rsum(iijj)**2/dble(itimesnum)
            standard_dev_r(iijj) =
     &      dsqrt(standard_dev_r(iijj)/dble(itimesnum))

            standard_dev_rs(iijj) =
     &      rssum2(iijj) - rssum(iijj)**2/dble(itimesnum)
            standard_dev_rs(iijj)=
     &      dsqrt(standard_dev_rs(iijj)/dble(itimesnum))

            standard_dev_s(iijj) =
     &      slopesum2(iijj) - slopesum(iijj)**2/dble(itimesnum)
            standard_dev_s(iijj) =
     &      dsqrt(standard_dev_s(iijj)/dble(itimesnum))

            standard_dev_ss(iijj) =
     &      sslopesum2(iijj) - sslopesum(iijj)**2/dble(itimesnum)
            standard_dev_ss(iijj) =
     &      dsqrt(standard_dev_ss(iijj)/dble(itimesnum))
         
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
         
c         deallocate(iia)
         deallocate(xavg,yavg)
c         deallocate(rs,zs,slopeu)
      
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
c         write(83,*) dble(np*npa(iijj))/dble(n),
c     &        ravg_array(iijj),sigmar_array(iijj),
c     &        ravg_arrayd(iijj),sigmar_arrayd(iijj)
c         write(84,*) dble(np*npa(iijj))/dble(n),
c     &        slope_array(iijj),sigma_slope_array(iijj)
            write(83,*) npa(iijj),
     &           ravg_array(iijj),sigmar_array(iijj),
     &           ravg_arrayd(iijj),sigmar_arrayd(iijj),
     &           ravg_array_max(iijj),sigmar_array_max(iijj),
     &           ravg_array_min(iijj),sigmar_array_min(iijj),                        
     &           ravg_arrayd_max(iijj),sigmar_arrayd_max(iijj),
     &           ravg_arrayd_min(iijj),sigmar_arrayd_min(iijj),
     &           standard_dev_r(iijj),standard_dev_rs(iijj),
     &           d_array(iijj),area_array(iijj),z_array(iijj),
     &           zrho_array(iijj),varz_array(iijj),zq_array(iijj)
            write(84,*) npa(iijj),
     &           slope_array(iijj),sigma_slope_array(iijj),
     &           slope_arrayd(iijj),sigma_slope_arrayd(iijj),
     &           slope_array_max(iijj),sigma_slope_array_max(iijj),
     &           slope_array_min(iijj),sigma_slope_array_min(iijj),                        
     &           slope_arrayd_max(iijj),sigma_slope_arrayd_max(iijj),
     &           slope_arrayd_min(iijj),sigma_slope_arrayd_min(iijj),
     &           standard_dev_s(iijj),standard_dev_ss(iijj)
            write(85,*) dble(np*npa(iijj))/dble(n),
     &           ravg_array(iijj),sigmar_array(iijj),
     &           ravg_arrayd(iijj),sigmar_arrayd(iijj),
     &           ravg_array_max(iijj),sigmar_array_max(iijj),
     &           ravg_array_min(iijj),sigmar_array_min(iijj),                        
     &           ravg_arrayd_max(iijj),sigmar_arrayd_max(iijj),
     &           ravg_arrayd_min(iijj),sigmar_arrayd_min(iijj),
     &           standard_dev_r(iijj),standard_dev_rs(iijj),
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
     &           standard_dev_s(iijj),standard_dev_ss(iijj),
     &           slope_array_abs(iijj),slope_array_abs_max(iijj),
     &           slope_array_abs_min(iijj),slope_array_dif(iijj),
     &           slope_array_value(iijj),areab_array(iijj)            

         end do   
            

         print*,'r error = ',sumr_error/dble(mp),sumrs_error/dble(mp)
         print*,'s error = ',sums_error/dble(mp),sumss_error/dble(mp)

         if (create_distribution) then
            dzr = (rmaxz - rminz)/dble(num_intervals)

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

            slopemax = slopeavg + 3.5d0*sigma_slope
            slopemin = slopeavg - 3.5d0*sigma_slope
            slopemin = .65d0
            slopemax = 1.15d0
            
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

c           Create analytical distribution based on general hypergeometric function
            print*,'rho analytical = ',rho
            print*,'start rdistribution'
            ss = nnn
            call rdistribution(rho,ss,varhyper,ra,analytical,m)
            print*,'end rdistribution'
      
            rpairavg = 0.d0
            rpair2avg = 0.d0
            dzp = 1.999d0/dble(m)
            do i = 1,m+1
               if (ra(i) .ge. rminz .and. ra(i) .le. rmaxz) then
                  write(97,*) ra(i),analytical(i)
               end if
               rpairavg = rpairavg + ra(i)*analytical(i)
               rpair2avg = rpair2avg + (ra(i)**2)*analytical(i)
            end do
            meanu = rpairavg*dzp

            print*,'mean analytical r = ',meanu
            rpair2avg = rpair2avg*dzp
            sigmau = dsqrt(rpair2avg - meanu**2)
            print*,'standard dev analytical r = ',sigmau**2
            print*,'rpairavg = ',rho-meanu

            ndnp = n/np
            print*,rho - rho*(1.d0 - rho**2)/(2.d0*dble(ndnp-1))


            ii = 2
            fac = 1.d0

c            nnn = 30
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
c            fac = fac*sdnor/sdnor2

            expectation = expectation*rho

            call hypser(0.5d0,0.5d0,dble(nnn+1)/2.d0,rho**2,hyp,dhyp)
            expectation = expectation*hyp

            nnnm1 = nnn - 1
            call hypser(1.d0,1.d0,dble(nnnm1)/2.d0+1.d0,rho**2,hyp,dhyp)
            expr2 = 1.d0 - hyp*(1.d0-rho**2)*dble(nnnm1-1)/dble(nnnm1)

            var = expr2 - expectation**2
      
            print*,'expectation var = ',expectation,var

            ija = 1
            sdnorx(1) = 1.d0
            sdnory(1) = 1.3d0

            sumb = 0.d0
            do i = 1,num_intervals
               sminz = slopemin + dble(i-1)*dz
               smaxz = slopemin + dble(i)*dz
               smid(i) = .5d0*(sminz+smaxz)
c               v = fac*(1.d0 - rho**2 + (rho -
c     &              (sdnor/sdnor2)*smid(i))**2)**(-s/2.d0)

               v = fac*(1.d0 - rho**2 + (rho -
     &              (sdnorx(ija)/sdnory(ija))*smid(i))**2)**(-s/2.d0)
               v = v*sdnorx(1)/sdnory(1)
c               v = fac*(1.d0 - rho**2 + (rho -
c     &         (1.d0/1.d0)*smid(i))**2)**(-s/2.d0)                              
               write(99,*) smid(i),v
               sumb = sumb + v*dz
            end do
            print*,'sumb = ',sumb
c            stop


            
            nu = dble(nnnm1)
            facgnu = facnu/dsqrt(nu*acos(-1.d0))
            expgnu = -(nu + 1.d0)/2.d0
c           print*,'sigma = ',sdnorx(ija),sdnory(ija),rsave(ija)
c            print*,'ratio = ',sdnorx(ija)/sdnory(ija)
c            sdnorx(1) = 1.d0
c            sdnory(1) = 1.d0

c            slopemin = -1.d0
c            dz = 2.d0/dble(num_intervals)
            bot = dsqrt(1.d0 - rsave(ija)**2)

c            print*,'sdn = ',sdnory(ija),sdnorx(ija)
            dbconv = dsqrt(nu)*sdnorx(ija)/(bot*sdnory(ija))
c            dbconv = dsqrt(nu)/bot
c            dbconv = 0.5d0*dsqrt(nu)*sdnory(ija)/(bot*sdnorx(ija))

            sumintegral = 0.d0
            do i = 1,num_intervals
               sminz = slopemin + dble(i-1)*dz
               smaxz = slopemin + dble(i)*dz
               smid(i) = .5d0*(sminz+smaxz)
               top = (smid(i)*sdnorx(ija)/sdnory(ija) - rsave(ija))
c               top = (smid(i)*1.d0/1.d0 - rsave(ija))
               
               fgnu = top*dsqrt(nu)/bot
c               print*,'fgnu = ',facgnu,fgnu,
c     &         (1.d0 + fgnu**2/dsqrt(nu))**expgnu
               v = dbconv*facgnu*(1.d0 + fgnu**2/nu)**expgnu

               sumintegral = sumintegral + v*dz
               call hypser(0.5d0,(nu+1.d0)/2.d0,1.5d0,-fgnu**2/nu,
     &                     hyp,dhyp)
               trap(i) = sumintegral
               cdf(i) = .5d0 + fgnu*facgnu*hyp
c               v = facgnu*(1.d0 + smid(i)**2/dsqrt(nu))**expgnu
               write(107,*) smid(i),v
            end do
            do i = 1,num_intervals
               write(108,*) smid(i),trap(i),cdf(i)
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
c         print*,'fac = ',fac                                                    \
                                                                                  
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
c      print*,'xbar ybar = ',xbar,ybar

      rtop = sumxy - sumx*sumy/n
      sigmax = sumx2 - rn*(sumx**2)
      sigmay = sumy2 - rn*(sumy**2)

c      print*,'sigmax sigmay = ',sigmax,sigmay
      
      rbot = dsqrt(sigmax*sigmay)

      r = rtop/rbot

      sdx = dsqrt(sigmax*rn)
      sdy = dsqrt(sigmay*rn)

c      print*,'sdx sdy = ',sdx,sdy

c      print*,'sigmay sigmax = ',sigmay,sigmax

c      print*,'r = ',r

c      print*,'sqrt = ',dsqrt(sigmay/sigmax)
      
      slope = r*dsqrt(sigmay/sigmax)

c      print*,'slope = ',slope

      yint = ybar - slope*xbar

c      print*,'slope yint = ',slope,yint

      return
      end

      subroutine aggreg(p,n,m,ndm,newpairs)
c     aggreg creates a partition given n original scores
      implicit none

c     Inputs
c     Number of groups
      integer ndm
c     Number of original scores
      integer n
c     Size of the groups
      integer m

c     Outputs
c     p(i,j) stores the original score index in element j of group i
      integer p(ndm,m)
c     pairs(i) stores the group the original index i is stored in
      integer pairs(n)

c     Local variables
      logical found,assigned(ndm)
      integer k,irn,kk,ij,i,j
      integer igroup,igroup_count
      integer map(n),newpairs(n)
      double precision r,eps


      eps = 1.d-14

      do i = 1,n
         map(i) = i
      end do

      p = 0

c     Index i refers to the group number
      do i = 1,ndm

         j = 0
         do while (j .lt. m) 

c           Create a random number (irn) from 1 to n - (i-1)*m
            r = rand(0)
            r = max(eps,r)
            r = min(r,1.d0-eps)
            irn = r*dble(n - (i-1)*m)
            irn = irn + 1


c           Ensure no duplicate members in group i
            found = .false.
            if (j .gt. 0) then
               do k = 1,j
                  if (p(i,k) .eq. map(irn)) then
                     found = .true.
                  end if
               end do
            end if

            if (.not. found) then
               j = j + 1
               p(i,j) = map(irn)
               pairs(p(i,j)) = i
            end if

         end do

         ij = 0

c        Collect all indices that have not been selected in array map
         do k = 1,n

            found = .false.
            do j = 1,i
               do kk = 1,m
                  if (k .eq. p(j,kk)) then
                     found = .true.
                  end if
               end do
            end do

            if (.not. found) then
               ij = ij + 1
               map(ij) = k
            end if

         end do

      end do

c     Reorder group assignments
      igroup_count = 0

      do k = 1,ndm
         assigned(k) = .false.
      end do

      do k = 1,n
         igroup = pairs(k)
         if (.not. assigned(igroup)) then
            assigned(igroup) = .true.
            igroup_count = igroup_count + 1
            do j = 1,n
               if (pairs(j) .eq. igroup) then
                  newpairs(j) = igroup_count
               end if
            end do
         end if
      end do

      end

      subroutine rdistribution(rho,ss,varhyper,ra,analytical,m)
      implicit none
      integer ss,n,nterms,k,i,m,ijk
      integer irn,irn32,ijka
      double precision daa,daa32,top2a,top3a
      double precision aa,bb,er,er32
      double precision ric,rick,prod
      double precision fachg,fachgk,factop2,factop3
      double precision rn,r,rho
      double precision coef,sum,sumr,pi,dr
      double precision fac3,addition
      double precision nsum,rsum,rsum2
      double precision varhyper
      double precision ff(m+1),ra(m+1)
      double precision analytical(m+1)

      pi = acos(-1.d0)


      n = ss - 1
      rn = dble(n)

      dr = 1.999d0/dble(m)
      nsum = 0.d0
      rsum = 0.d0
      rsum2 = 0.d0

      aa = 1.d0 - rho**2

      er = rn/2.d0
      er32 = (rn - 3.d0)/2.d0

      irn = er
      irn32 = er32

      daa = er - dble(irn)
      daa32 = er32 - dble(irn32)

      top2a = aa**daa

      do i = 1,m+1

         r = -1.d0 + dr*dble(i-1) + 1.d-8

         bb = 1.d0 - r**2
         top3a = bb**daa32
         coef = top2a*top3a/pi

         sum = 0.d0
         addition = 1.d0
         nterms = 0
         k = 0

         do while ((abs(addition) .gt. 1.d-14
     &              .and. nterms .lt. 1000) .or.
     &             nterms .lt. 200)

            ric = dble(n+k)/2.d0
            rick = dble(k+1)
            prod = 1.d0
            ijka = 1
            do ijk = n-2,1,-1

               if (ijka .le. irn) then
                  factop2 = aa
               else
                  factop2 = 1.d0
               end if

               if (ijka .le. irn32) then
                  factop3 = bb
               else
                  factop3 = 1.d0
               end if

               ijka = ijka + 1

               ric = ric - 1.d0
               if (ric .gt. 1.1d0) then
                  fachg = ric
               else
                  if (dabs(ric-.5d0) .lt. 1.d-12) then
                     fachg = .5d0*dsqrt(pi)
                  else
                     fachg = 1.d0
                  end if
               end if

               rick = rick - 1.d0
               if (rick .gt. 1.1d0) then
                  fachgk = rick
               else
                  fachgk = 1.d0
               end if

               prod = prod*fachg*fachg*2.d0/(fachgk*dble(ijk))
               prod = prod*factop2*factop3

            end do

            do while (ric .gt. 1.1d0 .or. rick .gt. 1.1d0)
               ric = ric - 1.d0
               if (ric .gt. 1.1d0) then
                  fachg = ric
               else
                  if (dabs(ric-.5d0) .lt. 1.d-12) then
                     fachg = .5d0*dsqrt(pi)
                  else
                     fachg = 1.d0
                  end if
               end if

               rick = rick - 1.d0
               if (rick .gt. 1.1d0) then
                  fachgk = rick
               else
                  fachgk = 1.d0
               end if
               prod = prod*fachg*fachg/fachgk

            end do

            fac3 = (2.d0*r*rho)**k
            addition = coef*fac3*prod

            nterms = nterms + 1
            sum = sum + addition

            k = k + 1
         end do

         ra(i) = r
         ff(i) = sum

         rsum = rsum + r*sum
         rsum2 = rsum2 + (r**2)*sum
         nsum = nsum + sum
      end do

      sumr = 0.d0
      do i = 1,m
         sumr = sumr + .5d0*dr*(ff(i)+ff(i+1))
      end do

      varhyper = ( rsum2 - rsum**2/nsum )/(nsum - 1.d0)

      do i = 1,m+1
         analytical(i) = ff(i)/sumr
      end do

c      print*,'sumr = ',sumr
c      print*,'nsum = ',nsum
c      print*,'varhyper = ',varhyper

      end
      
      subroutine sorte(a,n)
      integer n,a(n)
      integer i,j,amin,atemp

      do i = 1,n-1
         amin = a(i)
         jmin = 0
         do j = i+1,n
            if (a(j) .lt. amin) then
               jmin = j
               amin = a(j)
            end if
         end do
c         print*,'jmin amin = ',jmin,amin,i
         if (jmin .gt. 0) then
            atemp = a(i)
            a(i) = a(jmin)
            a(jmin) = atemp
         end if
c         do j = 1,n
c            print*,'a(',j,') = ',a(j)
c         end do
      end do

      return
      end

      
      subroutine freqdist(x,y,n,ymin,dy,freq)
      implicit none
      integer i,n,ii
      double precision ymin,ymax,dy
      double precision x(n),y(n)
      double precision freq(10),sumfreq
         
      ymin = 1.d+10
      ymax = -1.d+10
      do i = 1,n
         if (x(i) .ge.-.1d0 .and.
     &       x(i) .le. .1d0) then
            ymin = min(ymin,y(i))
            ymax = max(ymax,y(i))
         end if
      end do

      dy = (ymax - ymin)/10.d0

      freq = 0.d0
      do ii = 1,n
         do i = 1,10
            if (x(ii) .gt. -.1d0 .and.
     &          x(ii) .le. .1d0) then
               if (y(ii) .gt. ymin + dble(i-1)*dy .and.
     &             y(ii) .le. ymin + dble(i)*dy) then
                  freq(i) = freq(i) + 1.d0
               end if
            end if
         end do
      end do

      sumfreq = 0.d0
      do i = 1,10
         sumfreq = sumfreq + freq(i)
      end do

      do i = 1,10
         freq(i) = freq(i)/sumfreq
      end do

         
      return
      end

      
      subroutine permutation(p,n)
      implicit none
      logical alreadyused
      integer int50,is1,is2
      integer p(n),n,n2,i,j
      double precision eps,pv

      n2 = n/2
      eps = 1.d-10

      is1 = 0
      do while (is1 .lt. n2)
         pv = max(min(rand(),1.d0-eps),eps)
         int50 = dble(n)*pv+1
         alreadyused = .false.
         do i = 1,is1
            if (int50 .eq. p(i)) alreadyused = .true.
         end do
         if (.not. alreadyused) then
            is1 = is1 + 1
            p(is1) = int50
         end if
      end do

      is2 = 0
      do i = 1,n
         alreadyused = .false.
         do j = 1,is1
            if (i .eq. p(j)) alreadyused = .true.
         end do
         if (.not. alreadyused) then
            is2 = is2 +1
            p(is2+n2) = i
         end if
      end do

      return
      end

      subroutine tvalue(t,g1,g2,n)
      implicit none
      integer n,i
      double precision g1(n),g2(n),rn,t
      double precision sumg1,sumg2
      double precision sumg1_2,sumg2_2
      double precision g1a,g2a,g1ss,g2ss
      double precision sst

      rn = 1.d0/dble(n)
      sumg1 = 0.d0
      sumg2 = 0.d0
      sumg1_2 = 0.d0
      sumg2_2 = 0.d0
      do i = 1,n
         sumg1 = sumg1 + g1(i)
         sumg1_2 = sumg1_2 + g1(i)**2
         sumg2 = sumg2 + g2(i)
         sumg2_2 = sumg2_2 + g2(i)**2
      end do
      g1a = rn*sumg1
      g2a = rn*sumg2
      g1ss = sumg1_2 - rn*(sumg1**2)
      g2ss = sumg2_2 - rn*(sumg2**2)
      sst = (g1ss + g2ss)/dble(2*n-2)
      t = (g1a - g2a)/dsqrt(sst*2.d0*rn)
      t = g1a - g2a

      return
      end

      subroutine combin(n,m,nchoosem,comb)
      implicit none
      integer n,m
      integer nchoosem,ij,nit,istart,ik(m)
      integer comb(nchoosem,m)            
      external combinations

      ij = 0
      nit = 1
      istart = 1
c123456789012345678901234567890123456789012345678901234567890123456789012
      call combinations(istart,ij,ik,n,m,nit,nchoosem,comb,combinations)

      end

      subroutine combinations(istart,ij,ik,n,m,nit,nchoosem,comb,dsub)
      implicit none
      integer i,k,n,m,ij,nchoosem,nit,istart
      integer comb(nchoosem,m),ik(m)
      external dsub

c      print*,'nit = ',nit,istart,n
      do i = istart,n-m+nit
c         print*,'i = ',i
         if (nit .le. m-1) then
            ik(nit) = i
c            print*,'ik(',nit,') = ',ik(nit)
            nit = nit + 1      
            call dsub(i+1,ij,ik,n,m,nit,nchoosem,comb,dsub)
c            call dsub(1,ij,ik,n,m,nit,nchoosem,comb,dsub)            
         else
            ij = ij + 1
            ik(nit) = i
c            print*,'ik inside(',nit,') = ',ik(nit)
            do k = 1,m
               comb(ij,k) = ik(k)
c               print*,'comb(',ij,k,') = ',ik(k)
            end do
         end if
      end do

      nit = nit - 1

      end

      function fact(n)
      implicit none
      integer fact,n,i
      fact = 1
      do i = 2,n
         fact = fact*i
c         print*,'fact = ',fact
      end do
      end

      subroutine multiple(x,y,sol,yhat,r2,n,k)
      implicit none
      integer n,k,i,j,p
      double precision x(n,k),xmean(k)
      double precision y(n),ymean,yhat(n)
      double precision xa(n,k),xat(k,n)
      double precision xtx(k,k),rhs(k),sol(k)
      double precision sumx
      double precision r2,r2top,r2bot,rn

      rn = 1.d0/dble(n)
      
      ymean = 0.d0
      do i = 1,n
         ymean = ymean + y(i)
      end do
      ymean = ymean/dble(n)
      
      do j = 1,k
         sumx = 0.d0
         do i = 1,n
            sumx = sumx + x(i,j)
         end do
         xmean(j) = sumx/dble(n)
      end do

c      sumy = 0.d0
c      sumxx = 0.d0
c      sumxy = 0.d0
c      do i = 1,n
c         sumxy = sumxy + x(i,1)*y(i)
c         sumxx = sumxx + x(i,1)*x(i,1)
c         sumy = sumy + y(i)
c      end do

c      slope = sumxy - sumx*sumy/dble(n)
c      slope = slope/(sumxx - sumx*sumx/dble(n))
c      print*,'slope = ',slope
      
      do j = 1,k
         do i = 1,n
            xa(i,j) = x(i,j) - xmean(j)
         end do
      end do


      do j = 1,k
         do i = 1,n
            xat(j,i) = xa(i,j)
         end do
      end do

      do j = 1,k
         do p = 1,k
            sumx = 0.d0
            do i = 1,n
               sumx = sumx + xat(j,i)*xa(i,p)
            end do
            xtx(j,p) = sumx
         end do
      end do

      do j = 1,k
         sumx = 0.d0
         do i = 1,n
            sumx = sumx + xat(j,i)*(y(i) - ymean)
         end do
         rhs(j) = sumx
      end do

      call gauss_2(xtx,rhs,sol,k)

c      do p = 1,k
c         print*,'sol(',k,') = ',sol(k),rhs(k),xtx(k,k),rhs(k)/xtx(k,k)
c      end do
c      stop
      
      do i = 1,n
         yhat(i) = 0.d0
         do j = 1,k
            yhat(i) = yhat(i) + sol(j)*(x(i,j) - xmean(j))
         end do
         yhat(i) = yhat(i) + ymean
      end do

      r2top = 0.d0
      r2bot = 0.d0
      do i = 1,n
         r2top = r2top + (yhat(i) - ymean)**2
         r2bot = r2bot + (y(i) - ymean)**2
      end do
      r2 = r2top/r2bot

c      print*,'r2 = ',r2

c      x2 = 0.0
c      y2 = 0.0
c      z2 = 0.0
c      xy = 0.0
c      xz = 0.0
c      yz = 0.0
c      xmean(1) = 0.d0
c      xmean(2) = 0.d0
c      ymean = 0.d0

c      do i = 1,n
c         xmean(1) = xmean(1) + x(i,1)
c         xmean(2) = xmean(2) + x(i,2)
c         ymean = ymean + y(i)
c         x2 = x2 + x(i,1)**2
c         y2 = y2 + x(i,2)**2
c         z2 = z2 + y(i)**2
c         xy = xy + x(i,1)*x(i,2)
c         xz = xz + x(i,1)*y(i)
c         yz = yz + x(i,2)*y(i)
c      end do

c      rxy = (xy - xmean(1)*xmean(2)*rn)/
c     &      dsqrt((x2 - rn*xmean(1)**2)*(y2 - rn*xmean(2)**2))
c      rxz = (xz - xmean(1)*ymean*rn)/
c     &      dsqrt((x2 - rn*xmean(1)**2)*(z2 - rn*ymean**2))
c      ryz = (yz - xmean(2)*ymean*rn)/
c     &      dsqrt((y2 - rn*xmean(2)**2)*(z2 - rn*ymean**2))
c
c      r2 = (rxz**2 - 2.d0*rxy*ryz*rxz + ryz**2)/(1.d0 - rxy**2)
c
c      print*,'r2 = ',r2
      
      end

      subroutine gauss_2(a,b,x,n)
!===========================================================
! Solutions to a system of linear equations A*x=b
! Method: Gauss elimination (with scaling and pivoting)
! Alex G. (November 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! b(n)   - array of the right hand coefficients b
! n      - number of equations (size of matrix A)
! output ...
! x(n)   - solutions
! coments ...
! the original arrays a(n,n) and b(n) will be destroyed 
! during the calculation
!===========================================================
       implicit none 
       integer n
       double precision a(n,n), b(n), x(n)
       double precision s(n)
       double precision c, pivot, store
       integer i, j, k, l

c step 1: begin forward elimination
       do k=1, n-1

c step 2: "scaling"
c s(i) will have the largest element from row i 
          do i= k,n                       ! loop over rows
             s(i) = 0.0
             do j= k,n                    ! loop over elements of row i
                s(i) = max(s(i),abs(a(i,j)))
             end do
          end do

! step 3: "pivoting 1" 
! find a row with the largest pivoting element
          pivot = abs(a(k,k)/s(k))
          l = k
          do j=k+1,n
             if (abs(a(j,k)/s(j)) > pivot) then
                pivot = abs(a(j,k)/s(j))
                l = j
             end if
          end do

c         Check if the system has a sigular matrix
          if (pivot == 0.0) then
             write(*,*) ' The matrix is sigular '
             return
          end if

! step 4: "pivoting 2" interchange rows k and l (if needed)
          if (l /= k) then
             do j=k,n
                store = a(k,j)
                a(k,j) = a(l,j)
                a(l,j) = store
             end do
             store = b(k)
             b(k) = b(l)
             b(l) = store
          end if

! step 5: the elimination (after scaling and pivoting)
          do i=k+1,n
             c=a(i,k)/a(k,k)
             a(i,k) = 0.0
             b(i)=b(i)- c*b(k)
             do j=k+1,n
                a(i,j) = a(i,j)-c*a(k,j)
             end do
          end do
       end do

! step 6: back substiturion 
       x(n) = b(n)/a(n,n)
       do i=n-1,1,-1
          c=0.0
          do j=i+1,n
             c= c + a(i,j)*x(j)
          end do 
          x(i) = (b(i)- c)/a(i,i)
       end do

      end subroutine gauss_2 


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
      

        double precision function bimodal(x)
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
c        bimodal = bi1
        
        return
        end

        subroutine cumbinomial(x,xx)
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

      
c      subroutine ks
c      integer status
c      double precision pnor,qnor,xmm,meannor
c      double precision sdnor,boundnor,stz,meanz
c      double precision d,prob
c      double precision, allocatable :: data(:)
c
c      
c      meanz = 0.d0
c      stz = 1.d0
c     
c      meannor = 0.d0
c      sdnor = 1.d0
c
c      n = 10000
c      allocate(data(n))
c
c      do i =  1,n
c         pnor = min(max(rand(),eps),1.d0-eps)
c         qnor = 1.d0 - pnor
c
c         call cdfnor(2,pnor,qnor,xmm,meanz,stz,
c     &               status,boundnor)
c      
c         xm = sdnor*xmm + meannor
c         data(i) = xm
c      end do
c
c      call ksone(data,n,d,prob)
c      print*,'d prob = ',d,prob
c      
c      end

      SUBROUTINE ksone(data,n,d,prob,zavg,sigmaz)
      implicit none
      integer n
      double precision d,data(n),func,prob
      double precision zavg,sigmaz
      external func
c USES probks,sort
c Given an array data(1:n), and given a user-supplied function of a single variable func
c which is a cumulative distribution function ranging from 0 (for smallest values of its argument) to 1
c (for largest values of its argument), this routine returns the K–S statistic d, and
c the significance level prob. Small values of prob show that the cumulative distribution
c function of data is significantly different from func. The array data is modified by being
c sorted into ascending order.

      integer j
      double precision dt,en,ff,fn,fo,probks
      call sort(n,data) ! If the data are already sorted into ascending oren=n der, then this call can be omitted.
c      do j = 1,n
c         print*,'data(',j,') = ',data(j)
c      end do
      en=n
      d=0.
      fo=0. ! Data’s c.d.f. before the next step.
      do j=1,n ! Loop over the sorted data points.
           fn=j/en                ! Data’s c.d.f. after this step.
           ff=func(data(j),zavg,sigmaz)     ! Compare to the user-supplied function.
           dt=max(abs(fo-ff),abs(fn-ff)) ! Maximum distance.
           if (dt.gt.d) d=dt
           fo=fn
      end do 
      en=sqrt(en)
c      print*,(en+0.12d0+0.11d0/en),d,(en+0.12d0+0.11d0/en)*d
      prob=probks((en+0.12d0+0.11d0/en)*d) ! Compute significance.
      return
      END


      
      function probks(alam)
      implicit none
      double precision probks,alam,EPS1,EPS2
      parameter (EPS1=0.001d0, EPS2=1.d-8)
c     Kolmogorov-Smirnov probability function.
      integer j
      double precision a2,fac,term,termbf
      a2=-2.d0*alam**2
      fac=2.d0
      probks=0.d0
      termbf=0.d0 ! Previous term in sum.
      do j=1,100
         term=fac*exp(a2*j**2)
         probks=probks+term
         if(abs(term).le.EPS1*termbf.or.abs(term).le.EPS2*probks) return
         fac=-fac ! Alternating signs in sum.
         termbf=abs(term)
      end do 
      probks=1.d0 ! Get here only by failing to converge.
      return
      end

      function func(rr,zavg,sigmaz)
      implicit none
      double precision func,result,rr,ccum
      double precision zavg,sigmaz,zz
      zz = (rr - zavg)/sigmaz
      call cumnor(zz,result,ccum)
      func = result
      end 
      
      subroutine sort(n,arr)
      implicit none
      integer n,M,NSTACK
      double precision arr(n)
      parameter (M=7,NSTACK=50)
c     Sorts an array arr(1:n) into ascending numerical order using the Quicksort algorithm.
c     n is input; arr is replaced on output by its sorted rearrangement.
c     Parameters: M is the size of subarrays sorted by straight insertion and NSTACK is the required
c     auxiliary storage.
      integer i,ir,j,jstack,k,l,istack(NSTACK)
      double precision a,temp
      jstack=0
      l=1
      ir=n
 1    if (ir-l.lt.M) then
c       Insertion sort when subarray small enough.
        do j=l+1,ir
          a=arr(j)
          do i=j-1,l,-1
            if (arr(i).le.a) goto 2
              arr(i+1)=arr(i)
          end do 
          i=l-1
 2        arr(i+1)=a
        end do 

        if (jstack.eq.0) return
        ir=istack(jstack) ! Pop stack and begin a new round of partitioning.
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2 
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if (arr(l).gt.arr(ir)) then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        end if
        if (arr(l+1).gt.arr(ir)) then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        end if
        if (arr(l).gt.arr(l+1)) then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        end if
        i=l+1 ! Initialize pointers for partitioning.
        j=ir
        a=arr(l+1) ! Partitioning element.
 3      continue ! Beginning of innermost loop.
           i=i+1 ! Scan up to find element > a.
        if (arr(i) .lt. a) goto 3
 4      continue
           j=j-1 ! Scan down to find element < a.
        if (arr(j).gt.a) goto 4
        if (j.lt.i) goto 5 ! Pointers crossed. Exit with partitioning complete.
        temp=arr(i) ! Exchange elements.
        arr(i) = arr(j)
        arr(j) = temp
        goto 3 ! End of innermost loop.
 5      arr(l+1)=arr(j) ! Insert partitioning element.
        arr(j)=a
        jstack=jstack+2
c        Push pointers to larger subarray on stack, process smaller subarray immediately.
        if (jstack.gt.NSTACK) then
           print*,"NSTACK too small in sort"
           stop
        end if
        if (ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        end if
      end if
      goto 1
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
