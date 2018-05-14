c
c ... maximal number of particles, times, and weight-function bins
       PARAMETER( NKMAX=250001, NT=201, NSMAX=100 ) 
c                                                   -> tav9sym250k
c      PARAMETER( NKMAX=37012, NT=2800, NSMAX=100 ) 
c                                                   -> tav9sym37k
c      PARAMETER( NKMAX=11200, NT=800, NSMAX=100 )
c
c ... maximal number of bins in x, y, and z direction
c                            or r, phi, and z direction
c
       include 'tav10_para.h'
c
c ... charcter constants for filenames
       character cfstat*80, cfc3d*80
       character  c1*1
       character  c4*4
c
c ... averaging fields
       integer ixzC(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       integer ixzN(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)

       real*8  icc0(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       real*8  icc3(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       real*8  icc4(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       real*8  icc6(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       real*8  icc9(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       real*8  iccx(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       real*8  jcc3(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       real*8  jcc4(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       real*8  jcc6(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       real*8  jcc9(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       real*8  jccx(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)

       double precision  xznu(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  xznu_4(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  xCnu(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  vxCnu(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  vyCnu(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  vzCnu(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  xznur(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  xznu1(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  xznu2(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
c
       double precision  pnum1(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  cnum1(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  amom1(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  amom2(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  amom3(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  amom4(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  amom5(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  apnum1(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  aamom1(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  aamom2(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  aamom3(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  aamom4(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  aamom5(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
c
       double precision  f(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  f1(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)

       double precision  e(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  e1(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  et(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  et1(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)

       double precision  sd(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  sd1(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  st(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  st1(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)

       double precision  s(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  s1(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  s2(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)

       double precision Cn(3,3,3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision Cn1(3,3,3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision Cn13(21)
       double precision Ct(3,3,3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision Ct1(3,3,3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision Ct13(45)
c
       integer NS, ishape
       integer ixb,iyb,izb,ixp,iyp,izp
       double precision  cfn,cft,cfl,cfr,cdelta,ctheta
       double precision  cnx,cny,cnz,ccnx,ccny,ccnz
       double precision  ctx,cty,ctz,cctx,ccty,cctz
       double precision  cmass
c
       double precision   px(-1:NKMAX),  py(-1:NKMAX),  pz(-1:NKMAX)
       double precision  pvx(-1:NKMAX), pvy(-1:NKMAX), pvz(-1:NKMAX)
       double precision omex(-1:NKMAX),omey(-1:NKMAX),omez(-1:NKMAX)
       double precision phix(-1:NKMAX),phiy(-1:NKMAX),phiz(-1:NKMAX)
       double precision prad(-1:NKMAX),pxi(-1:NKMAX)

       integer inc(-1:NKMAX), inc0(0:IPCMAX), inc1(0:IPCMAX)
c
       double precision  time0,time,time_old,ctime,ctt0,ttime(NT)
       double precision  xj(-1:NKMAX),yj(-1:NKMAX),zj(-1:NKMAX)
       double precision  dr(-1:NKMAX),dz(-1:NKMAX),dphi(-1:NKMAX)
       double precision  vvrN,vvpN,vvzN,vvx,vvy,vvz
       double precision  wwrN,wwpN,wwzN,pwrN,pwpN,pwzN
       double precision  wwx,wwy,wwz,pwwx,pwwy,pwwz
       double precision  x00(-1:NKMAX),y00(-1:NKMAX),z00(-1:NKMAX)
       integer ijc(-1:NKMAX)

       double precision   ek(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision   epn(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision   ept(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision   de(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision   avfn(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision   avft(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision   avrfn(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision   avrft(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision   avdn(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision   avdn2(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision   avdt(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision   avdt2(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 

       double precision   vr(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision   vp(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision   vz(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision  vr2(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision  vp2(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  vz2(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)

       double precision   wr(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision   wp(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision   wz(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 

       double precision   pwr(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision   pwp(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision   pwz(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 

       double precision  drz_r(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX) 
       double precision  drz_p(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision  drz_z(-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)

       integer   icdr2(-1:IXBINMAX,-1:IZBINMAX,NT)
       double precision     dr1r(-1:IXBINMAX,-1:IZBINMAX,NT)
       double precision     dr1z(-1:IXBINMAX,-1:IZBINMAX,NT)
       double precision    dr1rp(-1:IXBINMAX,-1:IZBINMAX,NT)
       double precision     dr2r(-1:IXBINMAX,-1:IZBINMAX,NT)
       double precision     dr2z(-1:IXBINMAX,-1:IZBINMAX,NT)
       double precision    dr2rp(-1:IXBINMAX,-1:IZBINMAX,NT)
       double precision      dr2(-1:IXBINMAX,-1:IZBINMAX,NT)
       double precision    dr2xy(-1:IXBINMAX,-1:IZBINMAX,NT)

       double precision wt(-1:NSMAX)
c...   define the vector for the contact calculation
       double precision  hao_n1
       double precision  hao_n2
       double precision  hao_n3
       double precision  hao_delta


       PI=2.d0*ASIN(1.d0)
       ACCURA=1.d-7
c
c       print*, 'test'
c
       open(81,file='tav10.data')
       write(*,'(a)') '# tav10.1 - data averaging tool'
       open(92,file='dr2z.data')
       open(93,file='displacement.xb7')
c
c ... cylindrical (11|12) or cartesian (0|1)
c     with (1) or without (0) periodic boundaries
       read(*,*) icoord
c
c ... output control: (0) no stiffness tensor or (1)
c     with stiffness tensor  
       read(*,*) ioutput
       if(ioutput.gt.1)then
          icomout=0
          ioutput=1
       else
          icomout=1
       endif
c
       if(icomout.eq.1)then
          write(81,'(a)') '# tav10.1 - data averaging tool'
       endif
c
c ... read in contact information filename
       read(*,'(a)') cfc3d
       read(*,'(a)') cfstat
c
c ... get length of filenames
       lc10=80
       do 35 ilc=1,80
          if( cfc3d(ilc:ilc) .eq. ' ' )then
             lc10=ilc
             goto 36
          endif
35     continue
36     continue
       lc20=80
       do 37 ilc=1,80
          if( cfstat(ilc:ilc) .eq. ' ' )then
             lc20=ilc
             goto 38
          endif
37     continue
38     continue
c
c ... read in 4-digit begin and end numbers
       read(*,*) ifcnum0,ifcnum1,ifcstep
       if(ifcstep.lt.1) ifcstep=1
c
c ... read in time window
       read(*,*) tmin,tmax
c
c ... read in geometry info for averaging
       read(*,*) ixbin,xmin0,xmax0
       read(*,*) iybin,ymin0,ymax0
       read(*,*) izbin,zmin0,zmax0
c
       write(*,*) '# bins: ',ixbin,iybin,izbin
       write(*,*) '# max:  ',IXBINMAX,IYBINMAX,IZBINMAX
c
c ... read in smoothin method and shape function mode
       read(*,*) NS, ishape
c
       if(NS.gt.NSMAX) NS=NSMAX
c
c ... ishape > 6000 or ishape < -6000 (center weighted)
c ... ishape > 5000 or ishape < -5000 (contact weighted)
c ... ishape > 4000 or ishape < -4000 (center weighted)
c ... ishape > 3000 or ishape < -3000 (contact weighted)
c ... ishape > 2000 
c ... ishape > 1000 or ishape < -1000
c     the first digit is used as ishape, the rest 
c            sscale:=(ishape-3000)/100
c     is used as scaling factor, such that 3100 corresponds to ishape=3
c     but ishape=3215 stretches the line by a factor of 2.15
c     
       sscale=1.0
       if(ishape.lt.-1000)then
          if((ishape.le.-1000).and.(ishape.gt.-2000))then
             sscale=dble(-ishape-1000)/100.
             ishape=-1
          endif
c ... not implemented
c         if((ishape.lt.-2000).and.(ishape.gt.-3000))then
c            sscale=dble(-ishape-2000)/100.
c            ishape=-2
c         endif
          if((ishape.le.-3000).and.(ishape.gt.-4000))then
             sscale=dble(-ishape-3000)/100.
             ishape=-3
          endif
          if((ishape.le.-4000).and.(ishape.gt.-5000))then
             sscale=dble(-ishape-4000)/100.
             ishape=-4
          endif
          if((ishape.le.-5000).and.(ishape.gt.-6000))then
             sscale=dble(-ishape-5000)/100.
             ishape=-5
          endif
          if((ishape.le.-6000).and.(ishape.gt.-7000))then
             sscale=dble(-ishape-6000)/100.
             ishape=-6
          endif
       endif
       if(ishape.gt.1000)then
          if((ishape.ge.1000).and.(ishape.lt.2000))then
             sscale=dble(ishape-1000)/100.
             ishape=1
          endif
          if((ishape.ge.2000).and.(ishape.lt.3000))then
             sscale=dble(ishape-2000)/100.
             ishape=2
          endif
          if((ishape.ge.3000).and.(ishape.lt.4000))then
             sscale=dble(ishape-3000)/100.
             ishape=3
          endif
c ... not implemented
c         if((ishape.ge.5000).and.(ishape.lt.6000))then
c            sscale=dble(ishape-5000)/100.
c            ishape=5
c         endif
       endif
c
c ... ishape==-1,-3,-4,-5,-6 puts equidistant points on the center-contact 
c     line but weighs the contact higher according to the sphere-volume (-3,-4)
c     or according to a Gaussian (-5,-6) and the field wt() contains the weights
c     where (-3,-5) are contact weighted and (-4,-6) are center weighted.
       if(NS.ge.1)then
c
c ... weight according to sphere volume (increasing to contact)
          if((ishape.eq.-3).or.(ishape.eq.-4))then
             write(*,'(2a)') '# weight-function according to ',
     1                         'sphere-mass-distribution ...'
             wwt=1./(2*NS)**3*(NS+1)
             wt(0)=wwt
             wtsum=wt(0)
c
             do 1501 is=1,NS-1
                wt(is)=wwt*((24*is-16)*is+6)
                wtsum=wtsum+wt(is)
1501         continue
             wt(NS)=wwt*((20*NS-18)*NS+5)
             wtsum=wtsum+wt(NS)
c
c ... weight increasing towards contact (Gaussian)
          elseif((ishape.eq.-5).or.(ishape.eq.-6))then
             write(*,'(2a)') '# weight-function according to ',
     1                         'Gaussian with 1/2 scale width ...'
             wt(0)=1.
             wtsum=wt(0)
c            write(*,*) 0,wt(0),wtsum
             do 1505 is=1,NS-1
                wt(is)=exp(-(2.0*is/NS)**2)
                wtsum=wtsum+wt(is)
c               write(*,*) is,wt(is),wtsum
1505         continue
             wt(NS)=0.
             wtsum=wtsum+wt(NS)
             wwt=wtsum/(NS+1)
c            write(*,*) is,wt(is),wtsum
c
c ... normalise
             wt(0)=wt(0)/wwt
             wtsum=wt(0)
             write(*,*) 0,wt(0)
             do  1506 is=1,NS-1
                wt(is)=wt(is)/wwt
                wtsum=wtsum+wt(is)
                write(*,*) is,wt(is)
1506         continue
             wt(NS)=0.
             wtsum=wtsum+wt(NS)
c
c ... weight is constant along branch vector
          else
             write(*,'(a)') '# weight-function constant ... '
             wwt=1.
             wt(0)=wwt
             wtsum=wt(0)
c            write(*,*) 0,wt(0)
             do 1502 is=1,NS-1
                wt(is)=wwt
                wtsum=wtsum+wt(is)
c               write(*,*) is,wt(is)
1502         continue
             wt(NS)=wwt
             wtsum=wtsum+wt(NS)
          endif
          wtsum=wtsum/(NS+1)
       else
          wtsum=1.0
          wt(-1)=1.0
          wt(0)=1.0
          wt(1)=1.0
       endif
c
       write(*,'(a,i4,a,i4,a)') 
     1       '# weight function ',ishape,' with ',NS,' points'
       write(*,'(a,g16.7,a,g16.7)') 
     2       '# checksum=',wtsum,' stretch-scaling=',sscale
c
c ... read in particle density 
       read(*,*) pdensity
c
c ... IF ifcnum1 > ifcnum0
c     read in from ifcnum0 to ifcnum1 - and average over all
c ... ELSE
c     read in from ifcnum0 to -ifcnum1 - and average over each snapshot
c
c ... define logical for empty-line output in case of all-average
       iout_el=1
c
       ifcnum00=ifcnum0
       if(ifcnum1.lt.0)then
          ifcnum1=-ifcnum1
          ifcnum10=ifcnum1
          ifcnumstep=ifcstep
          iout_el=0
       else
c
c ... average over all files from ifcnum0 - ifcnum1
          if(ifcnum1.gt.ifcnum0)then
             ifcnumstep=ifcnum1-ifcnum0
c
c ... average over snapshots 
          else
             ifcnum1=ifcnum0+1
             ifcnumstep=1
          endif
          ifcnum10=ifcnum1
       endif
c
       itd=0
c
c ... switch off empty lines in case of one-d bins
       if(iybin.eq.0)then
          if(ixbin.eq.0)then
             iout_el=-3
          endif
          if(izbin.eq.0)then
             iout_el=-1
          endif
       else
         if(izbin.eq.0)then
            if(ixbin.eq.0)then
               iout_el=-2
            endif
         endif
       endif
c
c ... loop for average-over-all or time-step-average
       do 9000 ifcnum0=ifcnum00,ifcnum10-1,ifcnumstep

c      write(*,*) '... loop ... ',ifcnum0,ifcnum00,ifcnum10-1,ifcnumstep
c
c ... first reset all fields
       include 'tav10_reset.h' 
c
       ifcnum1=ifcnum0+ifcnumstep-1
       write(*,'(a,4i8)') '# start averaging sequence: ',
     1            ifcnum0,ifcnum00,ifcnum10,ifcnumstep
c       
       irnum1=-1
       irnum2=-1
c
c ... use the wall positions from c3d-files if x|y|z_min >= _max
       iwx=0
       iwy=0
       iwz=0
       if(xmin0.ge.xmax0)then 
          iwx=1
       endif
       if(ymin0.ge.ymax0)then
          if(icoord.eq.12)then
             ymin0=0.
             ymax0=2.*PI
          else
             iwy=1
          endif
       endif
       if(zmin0.ge.zmax0)then
          iwz=1
       endif
c
       if(ixbin.lt.0) ixbin=-1
       if(iybin.lt.0) iybin=-1
       if(izbin.lt.0) izbin=-1
       if(ixbin.gt.IXBINMAX) ixbin=IXBINMAX
       if(iybin.gt.IYBINMAX) iybin=IYBINMAX
       if(izbin.gt.IZBINMAX) izbin=IZBINMAX
       dxbin=(xmax0-xmin0)/(ixbin+1)
       dybin=(ymax0-ymin0)/(iybin+1)
       dzbin=(zmax0-zmin0)/(izbin+1)
c
       pcontrol=0.0
       p0control=0.0
c
c ... loop for multiple files input
       ifcc=0
       ifcN=0
       ifcF=0
       ifcc0=0
c
       time=-1.d10
       time_old=-1.d10
       itc=0
       jj=0
c
       do 8000 ifcnum=ifcnum0,ifcnum1,ifcstep
c
       write(*,'(a,4i6)') 
     1    '# ... loop 8000 ... ',ifcnum,ifcnum0,ifcnum1,ifcstep
c
       if((ifcnum0.eq.0).and.(ifcnum10.eq.0))then
          write(*,'(a)') '# do not use numbered files !'
       else
          write(c4,'(i4.4)') ifcnum
          cfc3d=cfc3d(1:lc10-1)//c4
          cfstat=cfstat(1:lc20-1)//c4
          lc1=lc10+4
          lc2=lc20+4
       endif
c      write(*,*) cfc3d,cfstat,len(cfc3d),lc1,lc2
c
c ... read rotation-format info data (7 - no rot-info
c     that means 8 columns, 14 - with rot info that
c     means 14 columns)
       open(90,file='c3d_rot.ini')
       irotinfo=14
       read(90,*,ERR=1090,END=1090) irotinfo
       if(irotinfo.ne.14) irotinfo=8
c      write(*,'(a,i4)') '# irotinfo = ',irotinfo
1090   continue
       close(90)
c
c ... reset averaging variables, which are used per loop, always
       do 1001 ixp=-1,IXBINMAX
          do 1002 iyp=-1,IYBINMAX
             do 1003 izp=-1,IZBINMAX
                xznu(ixp,iyp,izp)=0.d0
                xznu_4(ixp,iyp,izp)=0.d0
1003   continue
1002   continue
1001   continue
c
c ... read from xb7 coordinates file
       time0=time
c
c      write(*,'(2a)') '# Open file: ',cfc3d
       open(90,file=cfc3d)
       read(90,*) NK,time,wx10,wy10,wz10,wx11,wy11,wz11
       time_old=time0
       if(NK.gt.NKMAX) STOP '# NKMAX - Particle number too large !!!'
c
c ... strain calculations
       if(itd.eq.0)then
          sys_Lx00=wx11-wx10
          sys_Ly00=wy11-wy10
          sys_Lz00=wz11-wz10
          sys_Lx=wx11-wx10
          sys_Ly=wy11-wy10
          sys_Lz=wz11-wz10
       endif
       sys_Lx0=sys_Lx
       sys_Ly0=sys_Ly
       sys_Lz0=sys_Lz
       sys_Lx=wx11-wx10
       sys_Ly=wy11-wy10
       sys_Lz=wz11-wz10
c
c      if(itd.gt.0)then
          sys_Sx=(sys_Lx-sys_Lx0)/sys_Lx0
          sys_Sy=(sys_Ly-sys_Ly0)/sys_Ly0
          sys_Sz=(sys_Lz-sys_Lz0)/sys_Lz0
          sys_Sx0=(sys_Lx-sys_Lx00)/sys_Lx00
          sys_Sy0=(sys_Ly-sys_Ly00)/sys_Ly00
          sys_Sz0=(sys_Lz-sys_Lz00)/sys_Lz00
c      else
c         sys_Sx=0.
c         sys_Sy=0.
c         sys_Sz=0.
c         sys_Sx0=0.
c         sys_Sy0=0.
c         sys_Sz0=0.
c      endif

c      write(*,*) itd,sys_Lx00,sys_Lx0,sys_Lx,sys_Sx,sys_Sx0

       write(*,'(a,2g16.6)') '# read time: ',time,time_old
c
c ... do averaging only if ( tmax > time > tmin )
       if((time.ge.tmin).and.(time.le.tmax))then
c
       itc=itc+1
       itd=itd+1
       if(irnum1.eq.-1) irnum1=ifcnum
c
c ... set the min-max cell coordinates to wall-values
       if(iwx.gt.0)then
          xmin=wx10
          xmax=wx11
       else
          xmin=xmin0
          xmax=xmax0
       endif
       dxbin=(xmax-xmin)/(ixbin+1)
       if(iwy.gt.0)then
          ymin=wy10
          ymax=wy11
       else
          ymin=ymin0
          ymax=ymax0
       endif
       dybin=(ymax-ymin)/(iybin+1)
       if(iwz.gt.0)then
          zmin=wz10
          zmax=wz11
       else
          zmin=zmin0
          zmax=zmax0
       endif
       dzbin=(zmax-zmin)/(izbin+1)
c
c ... reset some fields for this loop
       include 'tav10_reloop.h'
c
c ... read from contact info-file for the first time (to count contacts)
       open(91,file=cfstat)
       read(91,*,ERR=9011) c1,ctime,xnu
9011   read(91,*,ERR=9012) c1,wx0,wy0,wz0,wx1,wy1,wz1
c ... 9012   read(91,*,ERR=9013) c1,rmin,rmax,r1,r2,r3,r4,r5
9012   read(91,*,ERR=9013) c1,rmin,rmax,r1,r2,r3,r4

9013   continue
c
       do 1009 i=0,NK-1
          inc(i)=0
1009   continue
c
c ... first read in of the contact file
c     in order to assign particle coordination number
       inumc=0
       do 1091 i=1,ICMAX
          read(91,*,ERR=1097,END=1191)
     1       ctt0,ic,jc,xcc,ycc,zcc,cdelta,ctheta,
     2       cfn,cft,cnx,cny,cnz,ctx,cty,ctz
c ... calculate the contact point from positions
c      hao_n1 = (px(ic)-px(jc))/sqrt((px(ic)-px(jc))**2
c     &+(py(ic)-py(jc))**2+(pz(ic)-pz(jc))**2)
c      hao_n2 = (py(ic)-py(jc))/sqrt((px(ic)-px(jc))**2
c     &+(py(ic)-py(jc))**2+(pz(ic)-pz(jc))**2)
c      hao_n3 = (pz(ic)-pz(jc))/sqrt((px(ic)-px(jc))**2
c     &+(py(ic)-py(jc))**2+(pz(ic)-pz(jc))**2)
c      hao_delta = (prad(ic)+prad(jc)) - ((px(ic)-px(jc))*hao_n1
c     &+(py(ic)-py(jc))*hao_n2+(pz(ic)-pz(jc))*hao_n3)
c      xcc = px(ic)-(hao_n1*(prad(ic)-0.5*hao_delta))
c      ycc = py(ic)-(hao_n2*(prad(ic)-0.5*hao_delta))
c      zcc = pz(ic)-(hao_n3*(prad(ic)-0.5*hao_delta))

          inumc=inumc+1
          inc(ic)=inc(ic)+1
c         write(*,*) ic,jc,inc(ic)
1097      continue
1091   continue
1191   continue
       close(91)
       write(*,'(a,i7,a)') '# read in ',inumc,' contacts ... first ...'
c
       if(abs(ctime-time).gt.ACCURA)then
          write(*,'(2a,2g18.9)') 
     1      '# WARNING: time-error ',cfc3d(1:lc10+4),time,time_old
          write(*,'(2a,g18.9)') 
     1      '# WARNING: time-error ',cfstat(1:lc20+4),ctime
       endif
c
c ... read in particles
       do 10 i=0,NK-1
          if(irotinfo.ne.14)then
            read(90,*,ERR=11,END=11) 
     1        px(i),py(i),pz(i),pvx(i),pvy(i),pvz(i),prad(i),pxi(i)
          else
            read(90,*,ERR=11,END=11) 
     1        px(i),py(i),pz(i),pvx(i),pvy(i),pvz(i),prad(i),
     2        omex(i),omey(i),omez(i),phix(i),phiy(i),phiz(i),pxi(i)
          endif
c
c ... compute cylindrical coordinates (11) and
c     project r,z results into the x,z plane
          zp=pz(i)
          if(icoord.ge.11)then
c ... radial distance and normal
             rip=sqrt(px(i)**2+py(i)**2)
             pnx=px(i)/rip
             pny=py(i)/rip
c            pnz=pz(i)/rip
             xp=rip
c ... tangential distance and normal
             prphi=atan(py(i)/px(i))
c ... full circle
             if(icoord.eq.12)then
                if((px(i).lt.0).and.(py(i).ge.0))then
                   prphi=prphi+PI
                endif
                if((px(i).lt.0).and.(py(i).lt.0))then
                   prphi=prphi+PI
                endif
                if((px(i).ge.0).and.(py(i).lt.0))then
                   prphi=prphi+PI*2.
                endif
c               write(*,'(a,i9,5g16.6)') 
c    1           '# mode-12 ',i,px(i),py(i),xp,prphi,pz(i)
             endif
             yp=prphi
c            write(*,'(a,i9,3g16.6)') 
c    1            '# mode-11 ',i,px(i),py(i),pz(i)
c            write(*,'(a,2i9,5g16.6)') 
c    1            '# mode-11 ',i,icoord,xp,yp,zp,rip,prphi
c
             vvx=pnx*pvx(i)+pny*pvy(i)
             vvy=-pnx*pvy(i)+pny*pvx(i)
c            write(*,'(a,i9,6g16.6)') 
c    1            '# mode-11 ',i,pnx,pny,pvx(i),pvy(i),vvx,vvy
             wwx=pnx*omex(i)+pny*omey(i)
             wwy=-pnx*omey(i)+pny*omex(i)
             pwx=pnx*phix(i)+pny*phiy(i)
             pwy=-pnx*phiy(i)+pny*phix(i)
          else
             xp=px(i)
             yp=py(i)
c            write(*,'(a,2i9,3g16.6)') 
c    1            '# mode-0|1 ',i,icoord,xp,yp,zp
             vvx=pvx(i)
             vvy=pvy(i)
             wwx=omex(i)
             wwy=omey(i)
             pwx=phix(i)
             pwy=phiy(i)
          endif
          vvz=pvz(i)
          wwz=omez(i)
          pwwz=phiz(i)
c
c ... get the cell indices
          ixp=int((xp-xmin)/dxbin)
          iyp=int((yp-ymin)/dybin)
          izp=int((zp-zmin)/dzbin)
c ... non-periodic
          if((icoord.eq.0).or.(icoord.eq.12))then
c            write(*,'(a,3i5)') '# 12: ',ixp,iyp,izp
             if(xp-xmin.lt.0) ixp=-1
             if(yp-ymin.lt.0) iyp=-1
             if(zp-zmin.lt.0) izp=-1
             if(ixp.ge.IXBINMAX) ixp=IXBINMAX
             if(iyp.ge.IYBINMAX) iyp=IYBINMAX
             if(izp.ge.IZBINMAX) izp=IZBINMAX
c ... periodic
          else
             if(xp-xmin.lt.0)then
                ixp=int((xp-xmin+xmax-xmin)/dxbin)
             endif
             if(yp-ymin.lt.0)then
                iyp=int((yp-ymin+ymax-ymin)/dybin)
             endif
             if(zp-zmin.lt.0)then
                izp=int((zp-zmin+zmax-zmin)/dzbin)
             endif
c
             if(xp.ge.xmax)then
                ixp=int((xp-(xmax-xmin))/dxbin)
             endif
             if(yp.ge.ymax)then
                iyp=int((yp-(ymax-ymin))/dybin)
             endif
             if(zp.ge.zmax)then
                izp=int((zp-(zmax-zmin))/dzbin)
             endif
c
             if(ixp.lt.0)then
                write(*,*) '... WARNING ixp = ',ixp
                ixp=-1
             endif
             if(iyp.lt.0)then
                write(*,*) '... WARNING iyp = ',iyp
                iyp=-1
             endif
             if(izp.lt.0)then
                write(*,*) '... WARNING izp = ',izp
                izp=-1
             endif
c
             if(ixp.ge.IXBINMAX)then
                write(*,*) '... WARNING ixp = ',ixp
                ixp=IXBINMAX
             endif
             if(iyp.ge.IYBINMAX)then
                write(*,*) '... WARNING iyp = ',iyp
                iyp=IYBINMAX
             endif
             if(izp.ge.IZBINMAX)then
                write(*,*) '... WARNING izp = ',izp
                izp=IZBINMAX
             endif
c
          endif
c
          ifcN=ifcN+1
          ixzN(ixp,iyp,izp)=ixzN(ixp,iyp,izp)+1
          xznu(ixp,iyp,izp)=xznu(ixp,iyp,izp)+4./3.*PI*prad(i)**3
          if(inc(i).lt.4)
     1     xznu_4(ixp,iyp,izp)=xznu_4(ixp,iyp,izp)+4./3.*PI*prad(i)**3
c
c ... take care of C=0 particles here ...
          iNS=NS+1
          if(inc(i).eq.0)then
c
             pcontrol=pcontrol+1.0
             p0control=p0control+1.0
             icc0(ixp,iyp,izp)=icc0(ixp,iyp,izp)+1+NS
c
             if(NS.le.0) iNS=1
             pViwt=4./3.*PI*prad(i)**3*iNS
             xCnu(ixp,iyp,izp)=xCnu(ixp,iyp,izp)+pViwt
             vxCnu(ixp,iyp,izp)=vxCnu(ixp,iyp,izp) + pViwt*vvx
             vyCnu(ixp,iyp,izp)=vyCnu(ixp,iyp,izp) + pViwt*vvy
             vzCnu(ixp,iyp,izp)=vzCnu(ixp,iyp,izp) + pViwt*vvz
c
             apnum1(ixp,iyp,izp)=apnum1(ixp,iyp,izp)+1
             aamom1(ixp,iyp,izp)=aamom1(ixp,iyp,izp)+prad(i)
             aamom2(ixp,iyp,izp)=aamom2(ixp,iyp,izp)+prad(i)**2
             aamom3(ixp,iyp,izp)=aamom3(ixp,iyp,izp)+prad(i)**3
             aamom4(ixp,iyp,izp)=aamom4(ixp,iyp,izp)+prad(i)**4
             aamom5(ixp,iyp,izp)=aamom5(ixp,iyp,izp)+prad(i)**5
c
c ... average the dynamic stress tensor = (1/V) \sum_i m_i v_i v_i 
             cmass=4./3.*PI*pdensity*prad(i)**3*iNS
c OLD        call sdaverage(ixp,iyp,izp,pvx(i),pvy(i),pvz(i),cmass,sd)
             call sdaverage(ixp,iyp,izp,vvx,vvy,vvz,cmass,sd)
c            write(*,*) 1,1,ixp,iyp,izp,cmass,pViwt,
c    1                  sd(1,1,ixp,iyp,izp),vxCnu(ixp,iyp,izp)
c
          else
c
c ... average the dynamic stress tensor = (1/V) \sum_i m_i v_i v_i
c     contact averaging has a small error - therefore average here ...
            cmass=4./3.*PI*pdensity*prad(i)**3*iNS
c OLD       call sdaverage(ixp,iyp,izp,pvx(i),pvy(i),pvz(i),cmass,sd)
            call sdaverage(ixp,iyp,izp,vvx,vvy,vvz,cmass,sd)
c           write(*,*) 1,1,ixp,iyp,izp,sd(1,1,ixp,iyp,izp)

          endif
c
c ... kinetic energy
          ek(ixp,iyp,izp)=ek(ixp,iyp,izp)+0.5*cmass*
     1                      (vvx**2+vvy**2+vvz**2)
c           write(*,*) ixp,iyp,izp,ek(ixp,iyp,izp)
c
c ... average particle velocity and square thereof
          vr(ixp,iyp,izp)=vr(ixp,iyp,izp)+vvx
          vp(ixp,iyp,izp)=vp(ixp,iyp,izp)+vvy
          vz(ixp,iyp,izp)=vz(ixp,iyp,izp)+vvz
          vr2(ixp,iyp,izp)=vr2(ixp,iyp,izp)+vvx**2
          vp2(ixp,iyp,izp)=vp2(ixp,iyp,izp)+vvy**2
          vz2(ixp,iyp,izp)=vz2(ixp,iyp,izp)+vvz**2
          wr(ixp,iyp,izp)=wr(ixp,iyp,izp)+wwx
          wp(ixp,iyp,izp)=wp(ixp,iyp,izp)+wwy
          wz(ixp,iyp,izp)=wz(ixp,iyp,izp)+wwz
          pwr(ixp,iyp,izp)=pwr(ixp,iyp,izp)+pwx
          pwp(ixp,iyp,izp)=pwp(ixp,iyp,izp)+pwy
          pwz(ixp,iyp,izp)=pwz(ixp,iyp,izp)+pwwz
c
c ... average mean square displacements
c     relative to the first data file
          if(itc.eq.1)then
             ttime(itc)=time
             x00(i)=px(i)
             y00(i)=py(i)
             z00(i)=pz(i)
          else
             ttime(itc)=time
             icdr2(ixp,izp,itc)=icdr2(ixp,izp,itc)+1
c
             dr1z(ixp,izp,itc)=dr1z(ixp,izp,itc)+(pz(i)-z00(i))
             dr2z(ixp,izp,itc)=dr2z(ixp,izp,itc)+(pz(i)-z00(i))**2
c
             if(icoord.ge.11)then
                rip00=sqrt(x00(i)**2+y00(i)**2)
                dr1r(ixp,izp,itc)=dr1r(ixp,izp,itc)
     1                           +(rip-rip00)
                dr2r(ixp,izp,itc)=dr2r(ixp,izp,itc)
     1                           +(rip-rip00)**2
c OLD
c                dphii=(prphi-atan(y00(i)/x00(i)))
c                if(dphii.lt.-PI/4) dphii=dphii+PI/2
c                if(dphii.gt.PI/4)  dphii=dphii-PI/2
                phi00=atan(y00(i)/x00(i))
                if(icoord.eq.12)then
                   if((x00(i).lt.0).and.(y00(i).ge.0))then
                      phi00=phi00+PI
                   endif
                   if((x00(i).lt.0).and.(y00(i).lt.0))then
                      phi00=phi00+PI
                   endif
                   if((x00(i).ge.0).and.(y00(i).lt.0))then
                      phi00=phi00+PI*2.
                   endif
                   dphii=prphi-phi00
                   if(dphii.lt.-PI/4) dphii=dphii+PI*2
                   if(dphii.gt.PI/4)  dphii=dphii-PI*2
                else
                   dphii=prphi-phi00
                   if(dphii.lt.-PI/4) dphii=dphii+PI/2
                   if(dphii.gt.PI/4)  dphii=dphii-PI/2
                endif

                dr1rp(ixp,izp,itc)=dr1rp(ixp,izp,itc)+rip*dphii
                dr2rp(ixp,izp,itc)=dr2rp(ixp,izp,itc)+(rip*dphii)**2
             else
                dr1r(ixp,izp,itc)=dr1r(ixp,izp,itc)+(px(i)-x00(i))
                dr2r(ixp,izp,itc)=dr2r(ixp,izp,itc)+(px(i)-x00(i))**2
                dr1rp(ixp,izp,itc)=0.
                dr1r(ixp,izp,itc)=dr1r(ixp,izp,itc)+(px(i)-x00(i))
                dr2r(ixp,izp,itc)=dr2r(ixp,izp,itc)+(px(i)-x00(i))**2
                dr1rp(ixp,izp,itc)=0.
                dr2rp(ixp,izp,itc)=0.
             endif
c
             dr2xy(ixp,izp,itc)=dr2xy(ixp,izp,itc)
     1        +(px(i)-x00(i))**2+(py(i)-y00(i))**2
             dr2(ixp,izp,itc)=dr2(ixp,izp,itc)
     1        +(px(i)-x00(i))**2+(py(i)-y00(i))**2+(pz(i)-z00(i))**2
          endif
c
c ... average over relative displacement since last output
c         write(*,*) '# jj,itd ', jj,itd
c
c ... use this line for long-time averages
c         if((jj.gt.0))then
c ... use this line for single-time-averages
          if((jj.gt.0).or.(itd.gt.1))then
c
             if(icoord.ge.11)then
                rjp=sqrt(xj(i)*xj(i)+yj(i)*yj(i))

                dr(i)=dr(i)+rip-rjp
                drz_r(ixp,iyp,izp)=drz_r(ixp,iyp,izp)+rip-rjp
 
c               phij=atan(yj(i)/xj(i))
c               dphii=prphi-phij
c               if(dphii.lt.-PI/4) dphii=dphii+PI/2
c               if(dphii.gt.PI/4)  dphii=dphii-PI/2

                phij=atan(yj(i)/xj(i))
                if(icoord.eq.12)then
                   if((xj(i).lt.0).and.(yj(i).ge.0))then
                      phij=phij+PI
                   endif
                   if((xj(i).lt.0).and.(yj(i).lt.0))then
                      phij=phij+PI
                   endif
                   if((xj(i).ge.0).and.(yj(i).lt.0))then
                      phij=phij+PI*2.
                   endif
                   dphii=prphi-phij
                   if(dphii.lt.-PI/4) dphii=dphii+PI*2
                   if(dphii.gt.PI/4)  dphii=dphii-PI*2
             write(*,*) 'mode12: ', xj(i),yj(i),zj(i),irjp,dphii
                else
                   dphii=prphi-phij
                   if(dphii.lt.-PI/4) dphii=dphii+PI/2
                   if(dphii.gt.PI/4)  dphii=dphii-PI/2
                endif

                dphi(i)=dphi(i)+(rip+rjp)/2.*dphii
                drz_p(ixp,iyp,izp)=drz_p(ixp,iyp,izp)
     1                             +(rip+rjp)/2.*dphii
              write(*,*) 'dphii: ', ixp,iyp,izp,drz_p(ixp,iyp,izp),dphii
             else
                ddx=px(i)-xj(i)
                if(icoord.gt.0)then
                   dwx=(xmax-xmin)/2.
                   if(ddx.gt. dwx) ddx=ddx-2.*dwx
                   if(ddx.lt.-dwx) ddx=ddx+2.*dwx
                endif
                dr(i)=dr(i)+ddx
                drz_r(ixp,iyp,izp)=drz_r(ixp,iyp,izp)+ddx
c               write(*,*) ixp,iyp,izp,drz_r(ixp,iyp,izp),ddx
c
                ddy=py(i)-yj(i)
                if(icoord.gt.0)then
                   wwy=(ymax-ymin)/2.
                   if(ddy.gt. wwy) ddy=ddy-2.*wwy
                   if(ddy.lt.-wwy) ddy=ddy+2.*wwy
                endif
                dphi(i)=dphi(i)+ddy
                drz_p(ixp,iyp,izp)=drz_p(ixp,iyp,izp)+ddy
c            write(*,*) 'ddy: ', ixp,iyp,izp,drz_p(ixp,iyp,izp),ddy
             endif
c
             ddz=pz(i)-zj(i)
             if(icoord.gt.0)then
                wwz=(zmax-zmin)/2.
                if(ddz.gt. wwz) ddz=ddz-2.*wwz
                if(ddz.lt.-wwz) ddz=ddz+2.*wwz
             endif
             drz_z(ixp,iyp,izp)=drz_z(ixp,iyp,izp)+ddz
             dz(i)=dz(i)+ddz
c            write(*,*) i,pz(i),zj(i),ddz
c
             ijc(i)=ijc(i)+1
          endif
c
          xj(i)=px(i)
          yj(i)=py(i)
          zj(i)=pz(i)
c
10     continue
       write(*,'(a,2i8,a,i8)') 
     1   '# read in ',i,NK,' particles ... ',ifcnum
       goto 12
11     continue
       write(*,'(a,a,i8)') '# ERROR reading: ',cfc3d(1:lc1),i
12     continue
       close(90)
c
c ... compute the mean and the square of the density
       do 1011 ixp=-1,IXBINMAX
          do 1012 iyp=-1,IYBINMAX
             do 1013 izp=-1,IZBINMAX
                xznur(ixp,iyp,izp)=xznur(ixp,iyp,izp)
     1                            +xznu_4(ixp,iyp,izp)
                xznu1(ixp,iyp,izp)=xznu1(ixp,iyp,izp)
     1                            +xznu(ixp,iyp,izp)
                xznu2(ixp,iyp,izp)=xznu2(ixp,iyp,izp)
     1                            +xznu(ixp,iyp,izp)**2
1013   continue
1012   continue
1011   continue
       ifcF=ifcF+1
c
c ... reset =0 the averaging fields which are used per loop
c     this was formerly placed here (tav7-version)
c       include 'tav10_reloop.h'
c
c ... read from contact info-file for the first time (to count contacts)
c
c       open(91,file=cfstat)
c       read(91,*) c1,ctime,xnu
c       read(91,*) c1,wx0,wy0,wz0,wx1,wy1,wz1
c       read(91,*) c1,rmin,rmax,r1,r2,r3,r4
c
c       if(abs(ctime-time).gt.ACCURA)then
c          write(*,'(2a,2g18.9)') 
c     1      '# WARNING: time-error ',cfc3d,time,time_old
c          write(*,'(2a,g18.9)') 
c     1      '# WARNING: time-error ',cfstat,ctime
c       endif
c
c ... first read in of the contact file
c     in order to assign particle cooridnation number
c       inumc=0
c       do 1091 i=1,ICMAX
c          read(91,*,ERR=1097,END=1191)
c     1       ctime,ic,jc,xcc,ycc,zcc,cdelta,ctheta,
c     2       cfn,cft,cnx,cny,cnz,ctx,cty,ctz
c          inumc=inumc+1
c          inc(ic)=inc(ic)+1
c          write(*,*) ic,jc,inc(ic)
c1097      continue
c1091   continue
c1191   continue
c       close(91)
c
c ... assign coordination numbers to particle-counter
       do 1092 i=0,NK-1
          incic=inc(i)
          if(incic.gt.IPCMAX) incic=IPCMAX
          inc0(incic)=inc0(incic)+1
          inc1(incic)=inc1(incic)+incic
1092   continue
c
c ... read from contact info-file (again)
       open(91,file=cfstat)
       read(91,*,ERR=8011) c1,ctime,xnu
8011   read(91,*,ERR=8012) c1,wx0,wy0,wz0,wx1,wy1,wz1
8012   read(91,*,ERR=8013) c1,rmin,rmax,r1,r2,r3,r4
c8012   read(91,*,ERR=8013) c1,rmin,rmax,r1,r2,r3,r4,r5
8013   continue

c
       if(abs(ctime-time).gt.ACCURA)then
          write(*,'(2a,2g18.9)') 
     1      '# WARNING: time-error ',cfc3d,time,time_old
          write(*,'(2a,g18.9)') 
     1      '# WARNING: time-error ',cfstat,ctime
       endif
c
       do 100 i=1,inumc
          read(91,*,ERR=97,END=101)
     1       ctime,ic,jc,xcc,ycc,zcc,cdelta,ctheta,
     2       cfn,cft,cnx,cny,cnz,ctx,cty,ctz
c ... calculate the contact point from positions
      hao_n1 = (px(ic)-px(jc))/sqrt((px(ic)-px(jc))**2
     &+(py(ic)-py(jc))**2+(pz(ic)-pz(jc))**2)
      hao_n2 = (py(ic)-py(jc))/sqrt((px(ic)-px(jc))**2
     &+(py(ic)-py(jc))**2+(pz(ic)-pz(jc))**2)
      hao_n3 = (pz(ic)-pz(jc))/sqrt((px(ic)-px(jc))**2
     &+(py(ic)-py(jc))**2+(pz(ic)-pz(jc))**2)
      hao_delta = (prad(ic)+prad(jc)) - ((px(ic)-px(jc))*hao_n1
     &+(py(ic)-py(jc))*hao_n2+(pz(ic)-pz(jc))*hao_n3)
      xcc = px(ic)-(hao_n1*(prad(ic)-0.5*hao_delta))
      ycc = py(ic)-(hao_n2*(prad(ic)-0.5*hao_delta))
      zcc = pz(ic)-(hao_n3*(prad(ic)-0.5*hao_delta))
          goto 98
 97       continue
          abscn=sqrt(cnx*cnx+cny*cny+cnz*cnz)
          if((abscn.gt.1.+ACCURA).or.(abscn.lt.1.d0-ACCURA))then
             write(*,*) '# Error reading n ',i,ic,jc,abscn
             goto 99
          endif
          absct=sqrt(ctx*ctx+cty*cty+ctz*ctz)
          if((absct.gt.1.+ACCURA).or.(absct.lt.1.d0-ACCURA))then
             if((absct.gt.ACCURA).or.(absct.lt.-ACCURA))then
                write(*,*) '# Error reading t ',i,ic,jc,absct
                goto 99
             endif
          endif
 98       continue
c         write(*,*) 'read ',i,' of ',inumc,' contacts ...'
c
c ... Loop for smoothing along the branch vector
c
c ... NS is the averaging method:
c     NS=-1 means: average if pq-contact in av. volume
c     NS= 0 means: average if  p-center  in av. volume
c     NS >0 means: average along center-contact line with NS+1 steps
c                  and count if step-point lies in the volume
c
c ... ishape=0 is step shape function
c     ishape>0 or <0 allows for (future) other smoothing 
c                    shape functions, e.g. Gaussian
c     ishape=   -6 weigh according to a Gaussian (contact weighted)
c     ishape=   -5 weigh according to a Gaussian (particle center weighted)
c     ishape=+3/-3 weigh according to a sphere (the contacts are
c                  weighed higher than the center.  ishape positive
c                  shifts the NS points (each point same weight),
c                  while ishape negative keeps the av. points
c                  equidistant but puts according weight on them (xwt)
c     ishape=+2    same as +3, but weight according to a disk particle
c
c ... NEW - go along the contact line - start at center
c
          NS0=0
          if(NS.eq.-1) NS0=-1
c
          do 202 is=NS0,NS
c         write(*,*) ' loop ',is,NS0,NS
c
c ... weighting factor (for spheres ishape==3)
             xwt=1.
c
             if(is.eq.0)then
                xc=px(ic)
                yc=py(ic)
                zc=pz(ic)
                xwt=wt(0)
             else
c
c ... default step on the center-contact line (constant) 
                xis=dble(is)/NS
c
c ... selection of smoothing method
                if(ishape.gt.0)then
c
c ... this makes the steps smaller at the contact-point - wider at center
                   if(ishape.ge.2)then
                      xis=xis**(1./ishape)
                   endif
                else
c ... this uses a weighting factor for spheres with constant steps
                   if(ishape.eq.-6)then
                      xwt=wt(is)
                   elseif(ishape.eq.-5)then
                      xwt=wt(NS-is)
                   elseif(ishape.eq.-4)then
                      xwt=wt(NS-is)
                   elseif(ishape.eq.-3)then
                      xwt=wt(is)
                   elseif(ishape.eq.-2)then
c ... not implemented yet
                      xwt=1.0
                   elseif(ishape.eq.-1)then
                      xwt=wt(is)
                   else
c ... not implemented yet
                      ishape=0
                   endif
                endif
c
                if(sscale.eq.1.0)then
                   xc=px(ic)-xis*(prad(ic)-cdelta/2.)*cnx
                   yc=py(ic)-xis*(prad(ic)-cdelta/2.)*cny
                   zc=pz(ic)-xis*(prad(ic)-cdelta/2.)*cnz
                else
c
c ... contact weighted
                   if((ishape.eq.-1).or.(ishape.eq.-3)
     1                              .or.(ishape.eq.-5))then
                      xc=px(ic)+prad(ic)*(sscale-1.0)*cnx
     1                         -xis*(prad(ic)*(sscale)-cdelta/2.)*cnx
                      yc=py(ic)+prad(ic)*(sscale-1.0)*cny
     1                         -xis*(prad(ic)*(sscale)-cdelta/2.)*cny
                      zc=pz(ic)+prad(ic)*(sscale-1.0)*cnz
     1                         -xis*(prad(ic)*(sscale)-cdelta/2.)*cnz
c
c ... particle center weighted
                   elseif((ishape.eq.-4).or.(ishape.eq.-6))then
                      xc=px(ic)-xis*prad(ic)*(sscale)*cnx
                      yc=py(ic)-xis*prad(ic)*(sscale)*cny
                      zc=pz(ic)-xis*prad(ic)*(sscale)*cnz
                   endif
                endif
c
             endif
c
c          write(*,*) is,ic,jc,xcc,xc,px(jc),px(ic),prad(ic),cdelta
c          write(*,*) " Y.   ",ycc,yc,py(jc),py(ic),cnx,cny,cnz
c          write(*,*) " Z.   ",zcc,zc,pz(jc),pz(ic),xwt

c
c ... compute cylindrical coordinates (11|12) and
c     project r,z results into the x,z plane
c     or use cartesian coordinates for mode (0|1)
c
          if(icoord.ge.11)then
c ... radial distance and normal
             r11=sqrt(xc*xc+yc*yc)
             rnx=xc/r11
             rny=yc/r11
c ... tangential distance and normal
             rphi=atan(yc/xc)
             if(icoord.eq.12)then
                if((xc.lt.0).and.(yc.ge.0))then
                   rphi=rphi+PI
                endif
                if((xc.lt.0).and.(yc.lt.0))then
                   rphi=rphi+PI
                endif
                if((xc.ge.0).and.(yc.lt.0))then
                   rphi=rphi+PI*2.
                endif
c               write(*,'(a,5g14.6)') '# m12: ',xc,yc,zc,r11,rphi
             endif
             rtx=-rny
             rty=rnx
c
             x=r11
             y=rphi
c ... unit normal vector in cylindrical coordinates
             ccnr=rnx*cnx+rny*cny
             ccnt=rtx*cnx+rty*cny
             ccnx=ccnr
             ccny=ccnt
c ... unit tangential vector in cylindrical coordinates
             cctr=rnx*ctx+rny*cty
             cctt=rtx*ctx+rty*cty
             cctx=cctr
             ccty=cctt
          else
c
c ... Cartesian coordinates
             x=xc
             y=yc
             ccnx=cnx
             ccny=cny
             cctx=ctx
             ccty=cty
          endif
          z=zc
          ccnz=cnz
          cctz=ctz
c         write(*,*) '# cm12: ',ixb,iyb,izb,prad(ic),ccnx,ccny,ccnz
c
c ... compute bin number
c
          ixb=int((x-xmin)/dxbin)
          iyb=int((y-ymin)/dybin)
          izb=int((z-zmin)/dzbin)
c         write(*,'(3i4,3g15.5)') ixb,iyb,izb,x,y,z
c
c ... non-periodic
          if((icoord.eq.0).or.(icoord.eq.12))then
             if(x-xmin.lt.0) ixb=-1
             if(y-ymin.lt.0) iyb=-1
             if(z-zmin.lt.0) izb=-1
             if(ixb.ge.IXBINMAX) ixb=IXBINMAX
             if(iyb.ge.IYBINMAX) iyb=IYBINMAX
             if(izb.ge.IZBINMAX) izb=IZBINMAX
c ... periodic
          else
c
             if(x-xmin.lt.0)then
                ixb=int((x-xmin+xmax-xmin)/dxbin)
c               write(*,*) x,xmin,xmax,dxbin,ixb
             endif
             if(y-ymin.lt.0)then
                iyb=int((y-ymin+ymax-ymin)/dybin)
             endif
             if(z-zmin.lt.0)then
                izb=int((z-zmin+zmax-zmin)/dzbin)
             endif
c
             if(x.ge.xmax)then
                ixb=int((x-(xmax-xmin))/dxbin)
             endif
             if(y.ge.ymax)then
                iyb=int((y-(ymax-ymin))/dybin)
             endif
             if(z.ge.zmax)then
                izb=int((z-(zmax-zmin))/dzbin)
             endif
c
             if(ixb.lt.0)then
                write(*,*) i,is,'... WARNING ixb = ',ixb
                ixb=-1
             endif
             if(iyb.lt.0)then
                write(*,*) i,is,'... WARNING iyb = ',iyb
                iyb=-1
             endif
             if(izb.lt.0)then
                write(*,*) i,is,'... WARNING izb = ',izb
                izb=-1
             endif
             if(ixb.gt.IXBINMAX)then
                write(*,*) i,is,'... WARNING ixb = ',ixb
                ixb=IXBINMAX
             endif
             if(iyb.gt.IYBINMAX)then
                write(*,*) i,is,'... WARNING iyb = ',iyb
                iyb=IYBINMAX
             endif
             if(izb.gt.IZBINMAX)then
                write(*,*) i,is,'... WARNING izb = ',izb
                izb=IZBINMAX
             endif
          endif
c         write(*,'(4i4,3g15.5)') i,IXBINMAX,IYBINMAX,IZBINMAX
c         write(*,'(4i4,3g15.5)') icoord,ixb,iyb,izb,x,y,z
c
c ... update fields
          ixzC(ixb,iyb,izb)=ixzC(ixb,iyb,izb)+1
          de(ixb,iyb,izb)=de(ixb,iyb,izb)+cdelta
          avfn(ixb,iyb,izb)=avfn(ixb,iyb,izb)+cfn
          avft(ixb,iyb,izb)=avft(ixb,iyb,izb)+cft
          avrfn(ixb,iyb,izb)=avrfn(ixb,iyb,izb)+cfn*prad(ic)
          avrft(ixb,iyb,izb)=avrft(ixb,iyb,izb)+cft*prad(ic)
          avdn(ixb,iyb,izb)=avdn(ixb,iyb,izb)+cdelta
          avdn2(ixb,iyb,izb)=avdn2(ixb,iyb,izb)+cdelta**2
          avdt(ixb,iyb,izb)=avdt(ixb,iyb,izb)+ctheta
          avdt2(ixb,iyb,izb)=avdt2(ixb,iyb,izb)+ctheta**2
          epn(ixb,iyb,izb)=epn(ixb,iyb,izb)+0.5*cdelta*cfn
          ept(ixb,iyb,izb)=ept(ixb,iyb,izb)+0.5*ctheta*cft
c
c         write(*,*) ic,jc,inc(ic)
c
c ... compute cylindrical coordinates (11|12) and
c     project r,z results into the x,z plane
          zp=pz(ic)
          if(icoord.ge.11)then
c
c ... radial distance and normal
            rip=sqrt(px(ic)**2+py(ic)**2)
            pnx=px(ic)/rip
            pny=py(ic)/rip
c           pnz=pz(ic)/rip
            xp=rip
c
c ... tangential distance and normal
            prphi=atan(py(ic)/px(ic))
            if(icoord.eq.12)then
               if((px(ic).lt.0).and.(py(ic).ge.0))then
                  prphi=prphi+PI
               endif
               if((px(ic).lt.0).and.(py(ic).lt.0))then
                  prphi=prphi+PI
               endif
               if((px(ic).ge.0).and.(py(ic).lt.0))then
                  prphi=prphi+PI*2.
               endif
            endif
            yp=prphi
c           write(*,'(a,2i9,5g16.6)')
c    1           '# mode-11 ',i,icoord,xp,yp,zp,rip,prphi

            vvx=pnx*pvx(ic)+pny*pvy(ic)
            vvy=-pnx*pvy(ic)+pny*pvx(ic)
            wwx=pnx*omex(ic)+pny*omey(ic)
            wwy=-pnx*omey(ic)+pny*omex(ic)
            pwx=pnx*phix(i)+pny*phiy(i)
            pwy=-pnx*phiy(i)+pny*phix(i)
c
c ... Cartesian distances and velocities
          else
            xp=px(ic)
            yp=py(ic)
c           write(*,'(a,2i9,3g16.6)')
c    1           '# mode-0|1 ',i,icoord,xp,yp,zp
            vvx=pvx(ic)
            vvy=pvy(ic)
            wwx=omex(ic)
            wwy=omey(ic)
            pwx=phix(ic)
            pwy=phiy(ic)
          endif
          vvz=pvz(ic)
          wwz=omez(ic)
          pwwz=phiz(ic)
c
c ... smooth density calculation - note: does NOT go over C=0 particles
          if(inc(ic).gt.0)then
             pViwt=4./3.*PI*prad(ic)**3 /inc(ic) * xwt
             xCnu(ixb,iyb,izb)=xCnu(ixb,iyb,izb) + pViwt
c            write(*,*)'$',ic,ixb,iyb,izb,inc(ic),pViwt,
c    1             xCnu(ixb,iyb,izb),vxCnu(ixb,iyb,izb),xwt
c
             pcontrol=pcontrol+1.0/inc(ic)
c
c ... compute average velocity smoothly
             vxCnu(ixb,iyb,izb)=vxCnu(ixb,iyb,izb) + pViwt*vvx
             vyCnu(ixb,iyb,izb)=vyCnu(ixb,iyb,izb) + pViwt*vvy
             vzCnu(ixb,iyb,izb)=vzCnu(ixb,iyb,izb) + pViwt*vvz
c
c ... this case takes care of C=0 particles
c     NO this is taken care of before ... see above ...
c         else
c            pViwt=4./3.*PI*prad(ic)**3
c            xCnu(ixb,iyb,izb)=xCnu(ixb,iyb,izb) + pViwt
c
c            p0control=p0control+1.0
c            pcontrol=pcontrol+1.0
c
c            vxCnu(ixb,iyb,izb)=vxCnu(ixb,iyb,izb) + pViwt*vvx
c            vyCnu(ixb,iyb,izb)=vyCnu(ixb,iyb,izb) + pViwt*vvy
c            vzCnu(ixb,iyb,izb)=vzCnu(ixb,iyb,izb) + pViwt*vvz
c
          endif
c
c         write(*,'(g15.5,i7,a,3i4,a,$)') pcontrol,
c    1          ixzC(ixb,iyb,izb),'(',ixb,iyb,izb,')'
c
c ... count "rattlers" - that is particles with less-equal 3 contacts
          if(inc(ic).le.3)then
             if(inc(ic).gt.0)then
                icc3(ixb,iyb,izb)=icc3(ixb,iyb,izb)+1./inc(ic)
                jcc3(ixb,iyb,izb)=jcc3(ixb,iyb,izb)+1.
c ... the C=0 counting is done somewhere else
             endif
          else
c ... count particles with more ( <=4, <=6, <=9, >9) contacts
            if(inc(ic).le.4)then
              icc4(ixb,iyb,izb)=icc4(ixb,iyb,izb)+1./inc(ic)
              jcc4(ixb,iyb,izb)=jcc4(ixb,iyb,izb)+1.
            else
              if(inc(ic).le.6)then
                icc6(ixb,iyb,izb)=icc6(ixb,iyb,izb)+1./inc(ic)
                jcc6(ixb,iyb,izb)=jcc6(ixb,iyb,izb)+1.
              else
                if(inc(ic).le.9)then
                   icc9(ixb,iyb,izb)=icc9(ixb,iyb,izb)+1./inc(ic)
                   jcc9(ixb,iyb,izb)=jcc9(ixb,iyb,izb)+1.
                else
                   iccx(ixb,iyb,izb)=iccx(ixb,iyb,izb)+1./inc(ic)
                   jccx(ixb,iyb,izb)=jccx(ixb,iyb,izb)+1.
                endif
              endif
            endif
          endif
c
c         write(*,*) ic,jc,inc(ic),icc3(ixb,iyb,izb),
c    1       icc4(ixb,iyb,izb),icc6(ixb,iyb,izb),icc9(ixb,iyb,izb)
c
c ... stress average - compute branch vector
          if(jc.ge.0)then
             cfl=sqrt((xcc-px(ic))**2+
     1                (ycc-py(ic))**2+
     2                (zcc-pz(ic))**2)
             cfr=prad(ic)-cdelta*0.5
c ... check for periodic boundaries !!!
             if(cfl.gt.cfr) cfl=cfr
          else
             cfl=prad(ic)-cdelta*0.5
          endif
c         if(jc.ge.0)then
c            cfl=sqrt((px(jc)-px(ic))**2+
c    1                (py(jc)-py(ic))**2+
c    2                (pz(jc)-pz(ic))**2)
c            cfr=prad(ic)+prad(jc)
c ... check for periodic boundaries !!!
c            if(cfl.gt.cfr) cfl=cfr
c         else
c            cfl=prad(ic)
c         endif
c
c ... fabric average
c         write(*,*) ixb,iyb,izb,prad(ic),ccnx,ccny,ccnz
          call faverage(ixb,iyb,izb,prad(ic),ccnx,ccny,ccnz,xwt,f)
c
c ... elastic strain average
          call eaverage(ixb,iyb,izb,cdelta,cfl,ccnx,ccny,ccnz,xwt,e)
          call eaverage(ixb,iyb,izb,ctheta,cfl,cctx,ccty,cctz,xwt,et)
c
c ... stress average (normal and tangential) contact contribution
          call saverage(ixb,iyb,izb,cfn,cfl,ccnx,ccny,ccnz,xwt,s)
          call staverage(ixb,iyb,izb,cft,cfl,
     1                   ccnx,ccny,ccnz,cctx,ccty,cctz,xwt,st)
c
c ... average the dynamic stress tensor = (1/V) \sum_i m_i v_i v_i 
c     WARNING - here is a small error O(1/200) ??? average in p-loop
c         cmass=4./3.*PI*pdensity*prad(ic)**3
c         call sdaverage(ixb,iyb,izb,pvx(ic),pvy(ic),pvz(ic),cmass,sd)
c
c ... stiffness tensor average (normal) contact contribution
          call Caverage(ixb,iyb,izb,cfl,ccnx,ccny,ccnz,xwt,Cn)
c
c ... stiffness tensor average (tangential) contact contribution
          call Ctaverage(ixb,iyb,izb,cfl,
     1                  ccnx,ccny,ccnz,cctx,ccty,cctz,xwt,Ct)
c
c ... weighted averages of radii-distribution-moments
         cincic=1.d10
         if(inc(ic).gt.0)then
           cincic=inc(ic)
           pnum1(ixb,iyb,izb)=pnum1(ixb,iyb,izb)+xwt/cincic
           apnum1(ixb,iyb,izb)=pnum1(ixb,iyb,izb)+xwt/cincic
           cnum1(ixb,iyb,izb)=cnum1(ixb,iyb,izb)+xwt
c         
           amom1(ixb,iyb,izb)=amom1(ixb,iyb,izb)+prad(ic)*xwt/cincic
           amom2(ixb,iyb,izb)=amom2(ixb,iyb,izb)+prad(ic)**2*xwt/cincic
           amom3(ixb,iyb,izb)=amom3(ixb,iyb,izb)+prad(ic)**3*xwt/cincic
           amom4(ixb,iyb,izb)=amom4(ixb,iyb,izb)+prad(ic)**4*xwt/cincic
           amom5(ixb,iyb,izb)=amom5(ixb,iyb,izb)+prad(ic)**5*xwt/cincic
          aamom1(ixb,iyb,izb)=aamom1(ixb,iyb,izb)+prad(ic)*xwt/cincic
          aamom2(ixb,iyb,izb)=aamom2(ixb,iyb,izb)+prad(ic)**2*xwt/cincic
          aamom3(ixb,iyb,izb)=aamom3(ixb,iyb,izb)+prad(ic)**3*xwt/cincic
          aamom4(ixb,iyb,izb)=aamom4(ixb,iyb,izb)+prad(ic)**4*xwt/cincic
          aamom5(ixb,iyb,izb)=aamom5(ixb,iyb,izb)+prad(ic)**5*xwt/cincic
         endif
c
c ... end contact smoothing loop
202    continue
c
 99    continue
100    continue
       write(*,'(a,i8,a)') '# success: read in ',inumc,' contacts.'
       goto 102
101    continue
       write(*,'(a,i8,a)') '# end-of-file: read in ',inumc,' contacts.'
102    continue
       close(91)
c
c ... second loop over all particles to account for particles with C=0
c       do 1801 ic=0,NK-1
c          if(inc(ic).eq.0)then
c             xc=px(ic)
c             yc=py(ic)
c             zc=pz(ic)
cc
cc ... compute cylindrical coordinates (11) and
cc     project r,z results into the x,z plane
cc     or use cartesian coordinates for mode (0|1)
cc
c             if(icoord.eq.11)then
cc ... radial distance and normal
c                r11=sqrt(xc*xc+yc*yc)
cc ... tangential distance and normal
c                rphi=atan(yc/xc)
cc
c                x=r11
c                y=rphi
c             else
cc
cc ... Cartesian coordinates
c                x=xc
cc                y=yc
c             endif
c             z=zc
c
c ... compute bin number
c             ixb=int((x-xmin)/dxbin)
c             iyb=int((y-ymin)/dybin)
c             izb=int((z-zmin)/dzbin)
c             if(ixb.lt.0) ixb=-1
c             if(iyb.lt.0) iyb=-1
c             if(izb.lt.0) izb=-1
c             if(x-xmin.lt.0) ixb=-1
c             if(y-ymin.lt.0) iyb=-1
c             if(z-zmin.lt.0) izb=-1
c             if(ixb.ge.IXBINMAX) ixb=IXBINMAX
c             if(iyb.ge.IYBINMAX) iyb=IYBINMAX
c             if(izb.ge.IZBINMAX) izb=IZBINMAX
c
c ... account for C=0 particles
c             icc0(ixb,iyb,izb)=icc0(ixb,iyb,izb)+1+NS
c             xCnu(ixb,iyb,izb)=xCnu(ixb,iyb,izb)
c     1                        +4./3.*PI*prad(ic)**3*(1+NS)
c             p0control=p0control+1.0+NS
c
c            write(*,*) ic,inc(ic),ixb,iyb,izb,pcontrol,p0control,ifcc
c          endif
c1801   continue
c
c ... compute the mean and the square of the stresses
       do 2011 ixp=-1,IXBINMAX
         do 2012 iyp=-1,IYBINMAX
           do 2013 izp=-1,IZBINMAX
             do 2014 i1=1,3
               do 2015 i2=1,3
                 f1(i1,i2,ixp,iyp,izp)=f1(i1,i2,ixp,iyp,izp)
     1                                +f(i1,i2,ixp,iyp,izp)

                 e1(i1,i2,ixp,iyp,izp)=e1(i1,i2,ixp,iyp,izp)
     1                                +e(i1,i2,ixp,iyp,izp)
                 et1(i1,i2,ixp,iyp,izp)=et1(i1,i2,ixp,iyp,izp)
     1                                 +et(i1,i2,ixp,iyp,izp)

                 sd1(i1,i2,ixp,iyp,izp)=sd1(i1,i2,ixp,iyp,izp)
     1                                 +sd(i1,i2,ixp,iyp,izp)
c                write(*,*) i1,i2,ixp,iyp,izp, sd(i1,i2,ixp,iyp,izp)

                 s1(i1,i2,ixp,iyp,izp)=s1(i1,i2,ixp,iyp,izp)
     1                                +s(i1,i2,ixp,iyp,izp)
                 s2(i1,i2,ixp,iyp,izp)=s2(i1,i2,ixp,iyp,izp)
     1                                +s(i1,i2,ixp,iyp,izp)**2

                 st1(i1,i2,ixp,iyp,izp)=st1(i1,i2,ixp,iyp,izp)
     1                                 +st(i1,i2,ixp,iyp,izp)

             do 2016 i3=1,3
               do 2017 i4=1,3
                 Cn1(i1,i2,i3,i4,ixp,iyp,izp)
     1              =Cn1(i1,i2,i3,i4,ixp,iyp,izp)
     2              +Cn(i1,i2,i3,i4,ixp,iyp,izp)
                 Ct1(i1,i2,i3,i4,ixp,iyp,izp)
     1              =Ct1(i1,i2,i3,i4,ixp,iyp,izp)
     2              +Ct(i1,i2,i3,i4,ixp,iyp,izp)
2017   continue
2016   continue
2015   continue
2014   continue
2013   continue
2012   continue
2011   continue
c
       ifcc=ifcc+1
       jj=jj+1
c
c ... if time larger than max-time stop 
       else
          close(90)
          if(time.gt.tmax) goto 8099
          jj=jj+1
          if(irnum2.eq.-1) irnum2=ifcnum
c
c ... end tmin-tmax if-then-else-endif
       endif
c
c ... end loop for multi-file input
8000   continue
8099   continue
c
c .. create data on output file
c
       if(irnum2.eq.-1) irnum2=ifcnum-1
c
       if(icomout.eq.1)then
          write(81,'(a,a)') '# read until file: ',cfc3d(1:lc1)
          write(81,'(a,a)') '# read until file: ',cfstat(1:lc2)
          write(81,'(a,2i6)') '# number-range: ',irnum1,irnum2
          write(81,'(a,i6)') '# icoord-mode: ',icoord
       endif
c
c ... smoothing correction
       ifcc0=ifcc
       if(NS.gt.0) ifcc=ifcc*(NS+1)
c
       if(icomout.eq.1)then
          write(81,'(a,4i8)') '# number of read files = ',
     1                             ifcc,ifcc0,ifcN,ifcF
          write(81,'(a,i4)') '# number of points on branch-vec. = 1+',NS
          write(81,'(a,i4)') '# smoothing+shape method = ',ishape
          write(81,'(a,g16.6)') '#        stretch-scaling = ',sscale
          write(81,'(a,i4)') '# irotinfo (from c3d_rot.ini) = ',irotinfo
          write(81,'(a,g16.6)') '# particle mat. density = ',pdensity
c
          write(81,'(a,g16.6,22i8,i11,21i8)') '# CN-stat: ',
     1       xnu,NK,(inc0(i),i=0,20),inumc,(inc1(i),i=0,20)
       endif
c
       if(ioutput.eq.1)then  
       if(icoord.ge.11)then
         if(icomout.eq.1)then
          write(81,'(a,i7,3g16.6)') '# r-bins:',
     1                               ixbin,xmin,xmax,dxbin
          write(81,'(a,i7,3g16.6)') '# p-bins:',
     1                               iybin,ymin,ymax,dybin
          write(81,'(a,i7,3g16.6)') '# z-bins:',
     1                               izbin,zmin,zmax,dzbin
          write(81,'(a)')
     1     '# averages in the r-p-z cubes'
          write(81,'(33a,$)')
     1     '#  1:ir   2:ip   3:iz   4:r    5:phi  6:z    7:Vol  ',
     2       ' 8:Nc   9:Cc  10:nu  11:nuS 12:nuC 13:C_0 14:C_3 ',
     3       '15:C_4 16:C_6 17:C_9 18:C_x 19:p_n 20:pS  21:Srr ',
     4       '22:Sry 23:Srz 24:Syr 25:Syy 26:Syz 27:Szr 28:Szy ',
     5       '29:Szz 30:p_t 31: rr 32: ry 33: rz 34: yr 35: yy ',
     6       '36: yz 37: zr 38: zy 39: zz 40:p_d 41: rr 42: ry ',
     7       '43: rz 44: yr 45: yy 46: yz 47: zr 48: zy 49: zz ',
     8       '50: t  51:Vr  52:Vy  53:Vz  54:V2r 55:V2y 56:V2z ',
     9       '57:dr  58:dy  59:dz  60:Dab 61: rr 62: ry 63: rz ',
     1       '64: yr 65: yy 66: yz 67: zr 68: zy 69: zz 70:trF ',
     2       '71:Frr 72:Fry 73:Frz 74:Fyr 75:Fyy 76:Fyz 77:Fzr ',
     3       '78:Fzy 79:Fzz 80:p_E 81:Err 82:Ery 83:Erz 84:Eyr ',
     4       '85:Eyy 86:Eyz 87:Ezr 88:Ezy 89:Ezz 90:pEt 91: rr ',
     5       '92: ry 93: rz 94: yr 95: yy 96: yz 97: zr 98: zy ',
     6       '99: zz 100: pw 101: cw 102: a1 103: a2 104: a3 ',
     7      '105: a4 106: a5 107: Ex 108: Ey 109: Ez 110:    ',
     8      '111: wr 112: wy 113: wz 114:<Vx>115:<Vy>',
     9      '116:<Vz>117:    118:    119:    120:    ',
     1      '121:C11 122:C12 123:C13 124:C14 125:C15 ',
     2      '126:C16 127:C22 128:C23 129:C24 130:C25 ',
     3      '131:C26 132:C33 133:C34 134:C35 135:C36 ',
     4      '136:C44 137:C45 138:C46 139:C55 140:C56 ',
     5      '141:C66 142:    143:    144:    145:    ',
     6      '146:    147:    148:    149:    150:    ', 
     7      '151:Ct11 152:Ct12 153:Ct13 154:Ct14 155:Ct15 ',
     8      '156:Ct16 157:Ct17 158:Ct18 159:Ct19 160:Ct22 ',
     9      '161:Ct23 162:Ct24 163:Ct25 164:Ct26 165:Ct27 ',
     1      '166:Ct28 167:Ct29 168:Ct33 169:Ct34 170:Ct35 ',
     2      '171:Ct36 172:Ct37 173:Ct38 174:Ct39 175:Ct44 ',
     3      '176:Ct45 177:Ct46 178:Ct47 179:Ct48 180:Ct49 ',
     4      '181:Ct55 182:Ct56 183:Ct57 184:Ct58 185:Ct59 ',
     5      '186:Ct66 187:Ct67 188:Ct68 189:Ct69 190:Ct77 ',
     6      '191:Ct78 192:Ct79 193:Ct88 194:Ct89 195:Ct99 '
        write(81,'(20a)')
     7      '196:     197:nu_r 198:Ntot 199:Ctot 200: n0  ',
     8      '201: n1  202: n2  203: n3  204: n4  205: n5  ',
     9      '206: n6  207: n7  208: n8  209: n9  210: n10 ',
     1      '211: n11 212: n12 213: n13 214: n14 215: n15 ',
     2      '216: n16 217: n17 218: n18 219: n19 ',
     3      '220: c0  221: c1  222: c2  223: c3  224: c4  ',
     4      '225: c5  226: c6  227: c7  228: c8  229: c9  ',
     5      '230: c10 231: c11 232: c12 233: c13 234: c14 ',
     6      '235: c15 236: c16 237: c17 238: c18 239: c19 ',
     7      '240:Nc0  241:     242:     243:Nc3  244:Nc4  ',
     8      '245:     246:Nc6  247:Nc9  248:Ncx  249:     250:     ',
     9      '251:     252:     253:Cc3  254:Cc4  255:     ',
     1      '256:Cc6  257:Cc9  258:Ccx  259:     260: Z0  ',
     2      '261:     262:     263: Z3  264: Z4  265:     ',
     3      '266: Z6  267:     268:     269: Z00 270:     ',
     4      '271: aa1 272: aa2 273: aa3 274: aa4 275: aa5 ',
     5      '276: deN 277: ek  278: epn 279: ept 280:     ',
     6      '281:<fn> 282:<ft> 283: <d> 284: <t> 285:<d^2>',
     7      '286:<t^2>287:<fnr>288:<ftr>289:     290:     ',
     8      '291:<pr> 292:<pp> 293:<pz> <EOL>'
        endif
c     
       else
        if(icomout.eq.1)then
          write(81,'(a,i7,3g16.6)') '# x-bins:',
     1                               ixbin,xmin,xmax,dxbin
          write(81,'(a,i7,3g16.6)') '# y-bins:',
     1                               iybin,ymin,ymax,dybin
          write(81,'(a,i7,3g16.6)') '# z-bins:',
     1                               izbin,zmin,zmax,dzbin
          write(81,'(a)')
     1     '# averages in the x-y-z cubes'
          write(81,'(33a,$)')
     1     '#  1:ix   2:iy   3:iz   4:x    5:y    6:z    7:Vol  ',
     2       ' 8:Nc   9:Cc  10:nu  11:nuS 12:nuC 13:C_0 14:C_3 ',
     3       '15:C_4 16:C_6 17:C_9 18:C_x 19:p   20:pS  21:Sxx ',
     4       '22:Sxy 23:Sxz 24:Syx 25:Syy 26:Syz 27:Szx 28:Szy ',
     5       '29:Szz 30:p_t 31: xx 32: xy 33: xz 34: yx 35: yy ',
     6       '36: yz 37: zx 38: zy 39: zz 40:p_d 41: xx 42: xy ',
     7       '43: xz 44: yx 45: yy 46: yz 47: zx 48: zy 49: zz ',
     8       '50: t  51:Vx  52:Vy  53:Vz  54:V2x 55:V2y 56:V2z ',
     9       '57:dx  58:dy  59:dz  60:Dab 61: xx 62: xy 63: xz ',
     1       '64: yx 65: yy 66: yz 67: zx 68: zy 69: zz 70:trF ',
     2       '71:Fxx 72:Fxy 73:Fxz 74:Fyx 75:Fyy 76:Fyz 77:Fzx ',
     3       '78:Fzy 79:Fzz 80:p_E 81:Exx 82:Exy 83:Exz 84:Eyx ',
     4       '85:Eyy 86:Eyz 87:Ezx 88:Ezy 89:Ezz 90:pEt 91: xx ',
     5       '92: xy 93: xz 94: yx 95: yy 96: yz 97: zx 98: zy ',
     6       '99: zz 100: pw 101: cw 102: a1 103: a2 104: a3 ',
     7      '105: a4 106: a5 107:dEx 108:dEy 109:dEz 110:    ',
     8      '111: wx 112: wy 113: wz 114:<Vx>115:<Vy>',
     9      '116:<Vz>117: Ex 118: Ey 119: Ez 120:    ',
     1      '121:C11 122:C12 123:C13 124:C14 125:C15 ',
     2      '126:C16 127:C22 128:C23 129:C24 130:C25 ',
     3      '131:C26 132:C33 133:C34 134:C35 135:C36 ',
     4      '136:C44 137:C45 138:C46 139:C55 140:C56 ',
     5      '141:C66 142:    143:    144:    145:    ',
     6      '146:    147:    148:    149:    150:    ', 
     7      '151:Ct11 152:Ct12 153:Ct13 154:Ct14 155:Ct15 ',
     8      '156:Ct16 157:Ct17 158:Ct18 159:Ct19 160:Ct22 ',
     9      '161:Ct23 162:Ct24 163:Ct25 164:Ct26 165:Ct27 ',
     1      '166:Ct28 167:Ct29 168:Ct33 169:Ct34 170:Ct35 ',
     2      '171:Ct36 172:Ct37 173:Ct38 174:Ct39 175:Ct44 ',
     3      '176:Ct45 177:Ct46 178:Ct47 179:Ct48 180:Ct49 ',
     4      '181:Ct55 182:Ct56 183:Ct57 184:Ct58 185:Ct59 ',
     5      '186:Ct66 187:Ct67 188:Ct68 189:Ct69 190:Ct77 ',
     6      '191:Ct78 192:Ct79 193:Ct88 194:Ct89 195:Ct99  '
          write(81,'(20a)')
     7      '196:     197:nu_r 198:Ntot 199:Ctot 200: n0  ',
     8      '201: n1  202: n2  203: n3  204: n4  205: n5  ',
     9      '206: n6  207: n7  208: n8  209: n9  210: n10 ',
     1      '211: n11 212: n12 213: n13 214: n14 215: n15 ',
     2      '216: n16 217: n17 218: n18 219: n19 ',
     3      '220: c0  221: c1  222: c2  223: c3  224: c4  ',
     4      '225: c5  226: c6  227: c7  228: c8  229: c9  ',
     5      '230: c10 231: c11 232: c12 233: c13 234: c14 ',
     6      '235: c15 236: c16 237: c17 238: c18 239: c19 ',
     7      '240:Nc0  241:     242:     243:Nc3  244:Nc4  ',
     8      '245:     246:Nc6  247:Nc9  248:Ncx  249:     250:     ',
     9      '251:     252:     253:Cc3  254:Cc4  255:     ',
     1      '256:Cc6  257:Cc9  258:Ccx  259:     260: Z0  ',
     2      '261:     262:     263: Z3  264: Z4  265:     ',
     3      '266: Z6  267:     268:     269: Z00 270:     ',
     4      '271: aa1 272: aa2 273: aa3 274: aa4 275: aa5 ',
     5      '276: deN 277: ek  278: epn 279: ept 280:     ',
     6      '281:<fn> 282:<ft> 283: <d> 284: <t> 285:<d^2>',
     7      '286:<t^2>287:<fnr>288:<ftr>289:     290:     ',
     8      '291:<px> 292:<py> 293:<pz> <EOL>'
        endif
c
       endif
c
       else
       if(icoord.ge.11)then
        if(icomout.eq.1)then
          write(81,'(a,i7,3g16.6)') '# r-bins:',
     1                               ixbin,xmin,xmax,dxbin
          write(81,'(a,i7,3g16.6)') '# p-bins:',
     1                               iybin,ymin,ymax,dybin
          write(81,'(a,i7,3g16.6)') '# z-bins:',
     1                               izbin,zmin,zmax,dzbin
          write(81,'(a)')
     1     '# averages in the r-p-z cubes'
          write(81,'(22a)')
     1     '#  1:ir   2:ip   3:iz   4:r    5:phi  6:z    7:Vol  ',
     2       ' 8:Nc   9:Cc  10:nu  11:nuS 12:nuC 13:C_0 14:C_3 ',
     3       '15:C_4 16:C_6 17:C_9 18:C_x 19:p_n 20:pS  21:Srr ',
     4       '22:Sry 23:Srz 24:Syr 25:Syy 26:Syz 27:Szr 28:Szy ',
     5       '29:Szz 30:p_t 31: rr 32: ry 33: rz 34: yr 35: yy ',
     6       '36: yz 37: zr 38: zy 39: zz 40:p_d 41: rr 42: ry ',
     7       '43: rz 44: yr 45: yy 46: yz 47: zr 48: zy 49: zz ',
     8       '50: t  51:Vr  52:Vy  53:Vz  54:V2r 55:V2y 56:V2z ',
     9       '57:dr  58:dy  59:dz  60:Dab 61: rr 62: ry 63: rz ',
     1       '64: yr 65: yy 66: yz 67: zr 68: zy 69: zz 70:trF ',
     2       '71:Frr 72:Fry 73:Frz 74:Fyr 75:Fyy 76:Fyz 77:Fzr ',
     3       '78:Fzy 79:Fzz 80:p_E 81:Err 82:Ery 83:Erz 84:Eyr ',
     4       '85:Eyy 86:Eyz 87:Ezr 88:Ezy 89:Ezz 90:pEt 91: rr ',
     5       '92: ry 93: rz 94: yr 95: yy 96: yz 97: zr 98: zy ',
     6       '99: zz 100: pw 101: cw 102: a1 103: a2 104: a3 ',
     7      '105: a4 106: a5 107: Er 108: Ep 109: Ez 110:    ',
     8      '111: wr 112: wp 113: wz 114:<Vx>115:<Vy>',
     9      '116:<Vz>117:    118:    119:    120:    ',
     5      '121: de 122: ek 123:epn 124:ept 125:     ',
     6      '126:<fn> 127:<ft> 128: <d> 129: <t> 130:<d^2>',
     7      '131:<t^2>132:<fnr>133:<ftr>134:     135:     ',
     8      '136:<pr> 137:<pp> 138:<pz> <EOL>'
        endif
c     
       else
        if(icomout.eq.1)then
          write(81,'(a,i7,3g16.6)') '# x-bins:',
     1                               ixbin,xmin,xmax,dxbin
          write(81,'(a,i7,3g16.6)') '# y-bins:',
     1                               iybin,ymin,ymax,dybin
          write(81,'(a,i7,3g16.6)') '# z-bins:',
     1                               izbin,zmin,zmax,dzbin
          write(81,'(a)')
     1     '# averages in the x-y-z cubes'
          write(81,'(22a)')
     1     '#  1:ix   2:iy   3:iz   4:x    5:y    6:z    7:Vol  ',
     2       ' 8:Nc   9:Cc  10:nu  11:nuS 12:nuC 13:C_0 14:C_3 ',
     3       '15:C_4 16:C_6 17:C_9 18:C_x 19:p   20:pS  21:Sxx ',
     4       '22:Sxy 23:Sxz 24:Syx 25:Syy 26:Syz 27:Szx 28:Szy ',
     5       '29:Szz 30:p_t 31: xx 32: xy 33: xz 34: yx 35: yy ',
     6       '36: yz 37: zx 38: zy 39: zz 40:p_d 41: xx 42: xy ',
     7       '43: xz 44: yx 45: yy 46: yz 47: zx 48: zy 49: zz ',
     8       '50: t  51:Vx  52:Vy  53:Vz  54:V2x 55:V2y 56:V2z ',
     9       '57:dx  58:dy  59:dz  60:Dab 61: xx 62: xy 63: xz ',
     1       '64: yx 65: yy 66: yz 67: zx 68: zy 69: zz 70:trF ',
     2       '71:Fxx 72:Fxy 73:Fxz 74:Fyx 75:Fyy 76:Fyz 77:Fzx ',
     3       '78:Fzy 79:Fzz 80:p_E 81:Exx 82:Exy 83:Exz 84:Eyx ',
     4       '85:Eyy 86:Eyz 87:Ezx 88:Ezy 89:Ezz 90:pEt 91: xx ',
     5       '92: xy 93: xz 94: yx 95: yy 96: yz 97: zx 98: zy ',
     6       '99: zz 100: pw 101: cw 102: a1 103: a2 104: a3 ',
     7      '105: a4 106: a5 107: Ex 108: Ey 109: Ez 110:    ',
     8      '111: wx 112: wy 113: wz 114:<Vx>115:<Vy>',
     9      '116:<Vz>117:    118:    119:    120:    ',
     1      '121: de 122: ek 123:epn 124:ept 125:    ',
     6      '126:<fn> 127:<ft> 128: <d> 129: <t> 130:<d^2>',
     7      '131:<t^2>132:<fnr>133:<ftr>134:     135:     ',
     8      '136:<px> 137:<py> 138:<pz> <EOL>'
        endif
c
       endif
       endif
c
       chk_isum=0
       chk_xnu =0
       chk_xCnu=0
c
       do 1100 iix=0,ixbin
         do 1101 iiy=0,iybin
           do 1102 iiz=0,izbin
c
c ... set the phi- or y-length of the cell
             xpos=xmin+iix*dxbin+0.5*dxbin
             ypos=ymin+iiy*dybin+0.5*dybin

             if(icoord.ge.11)then
c ... OLD ...   dy=0.5*PI*xpos/(iybin+1)
                dy=dybin*xpos
                ypos=ymin+(ymax-ymin)/(iybin+1) * (iiy+0.5)
             else
                dy=dybin
             endif
             zpos=zmin+iiz*dzbin+0.5*dzbin
             Vxyz=dxbin*dzbin*dy
c            write(*,*) iix,iiy,iiz,xpos,ypos,zpos
c
c ... output to file (no correction for ixzN and xznu)
             Vifcc=1.d0/(Vxyz*ifcc)
             Vifcc2=1.d0/(Vxyz**2*ifcc)
             trf=(f1(1,1,iix,iiy,iiz)
     1           +f1(2,2,iix,iiy,iiz)
     2           +f1(3,3,iix,iiy,iiz))*Vifcc
             p1=(s1(1,1,iix,iiy,iiz)
     1          +s1(2,2,iix,iiy,iiz)
     2          +s1(3,3,iix,iiy,iiz))*Vifcc/3
             p2=(s2(1,1,iix,iiy,iiz)
     1          +s2(2,2,iix,iiy,iiz)
     2          +s2(3,3,iix,iiy,iiz))*Vifcc2/3
             pt=(st1(1,1,iix,iiy,iiz)
     1          +st1(2,2,iix,iiy,iiz)
     2          +st1(3,3,iix,iiy,iiz))*Vifcc/3
             pd=(sd1(1,1,iix,iiy,iiz)
     1          +sd1(2,2,iix,iiy,iiz)
     2          +sd1(3,3,iix,iiy,iiz))*Vifcc/3
             xnur=xznur(iix,iiy,iiz)/Vxyz/dble(ifcF)
             xnu =xznu1(iix,iiy,iiz)/Vxyz/dble(ifcF)
             xnu2=xznu2(iix,iiy,iiz)/Vxyz**2/dble(ifcF)
c
             if(ixzN(iix,iiy,iiz).gt.0)then
                deN=de(iix,iiy,iiz)/ixzC(iix,iiy,iiz)
                avfnN=avfn(iix,iiy,iiz)/ixzC(iix,iiy,iiz)
                avftN=avft(iix,iiy,iiz)/ixzC(iix,iiy,iiz)
                avrfnN=avrfn(iix,iiy,iiz)/ixzC(iix,iiy,iiz)
                avrftN=avrft(iix,iiy,iiz)/ixzC(iix,iiy,iiz)
                avdnN=avdn(iix,iiy,iiz)/ixzC(iix,iiy,iiz)
                avdtN=avdt(iix,iiy,iiz)/ixzC(iix,iiy,iiz)
                avdnN2=avdn2(iix,iiy,iiz)/ixzC(iix,iiy,iiz)
                avdtN2=avdt2(iix,iiy,iiz)/ixzC(iix,iiy,iiz)
                CZ00=ixzC(iix,iiy,iiz)/dble(ifcc) /
     1               (ixzN(iix,iiy,iiz)/dble(ifcF)) /2.
                vvrN=vr(iix,iiy,iiz)/ixzN(iix,iiy,iiz)
                vvpN=vp(iix,iiy,iiz)/ixzN(iix,iiy,iiz)
                vvzN=vz(iix,iiy,iiz)/ixzN(iix,iiy,iiz)
                vvr2N=vr2(iix,iiy,iiz)/ixzN(iix,iiy,iiz)
                vvp2N=vp2(iix,iiy,iiz)/ixzN(iix,iiy,iiz)
                vvz2N=vz2(iix,iiy,iiz)/ixzN(iix,iiy,iiz)
                wwrN=wr(iix,iiy,iiz)/ixzN(iix,iiy,iiz)
                wwpN=wp(iix,iiy,iiz)/ixzN(iix,iiy,iiz)
                wwzN=wz(iix,iiy,iiz)/ixzN(iix,iiy,iiz)
                pwrN=pwr(iix,iiy,iiz)/ixzN(iix,iiy,iiz)
                pwpN=pwp(iix,iiy,iiz)/ixzN(iix,iiy,iiz)
                pwzN=pwz(iix,iiy,iiz)/ixzN(iix,iiy,iiz)
                dt=time-time_old
c
                if(dt.ne.0.d0)then
                   drx=drz_r(iix,iiy,iiz)/ixzN(iix,iiy,iiz)/dt
                   dry=drz_p(iix,iiy,iiz)/ixzN(iix,iiy,iiz)/dt
                   drz=drz_z(iix,iiy,iiz)/ixzN(iix,iiy,iiz)/dt
                else
                   drx=0.d0
                   dry=0.d0
                   drz=0.d0
                endif
c               write(*,*) 'dxyz: ',iix, iiy, iiz, dt, drx, dry, drz
c
                if((iix.eq.1).and.(iiy.eq.1).and.(iiz.eq.1))then
                   ddxdx=(drz_r(iix+1,iiy,iiz)-drz_r(iix-1,iiy,iiz))
     1                     /ixzN(iix,iiy,iiz)/dt
                   ddxdy=(drz_r(iix,iiy+1,iiz)-drz_r(iix,iiy-1,iiz))
     1                     /ixzN(iix,iiy,iiz)/dt
                   ddxdz=(drz_r(iix,iiy,iiz+1)-drz_r(iix,iiy,iiz-1))
     1                     /ixzN(iix,iiy,iiz)/dt
                   ddydx=(drz_p(iix+1,iiy,iiz)-drz_r(iix-1,iiy,iiz))
     1                     /ixzN(iix,iiy,iiz)/dt
                   ddydy=(drz_p(iix,iiy+1,iiz)-drz_r(iix,iiy-1,iiz))
     1                     /ixzN(iix,iiy,iiz)/dt
                   ddydz=(drz_p(iix,iiy,iiz+1)-drz_r(iix,iiy,iiz-1))
     1                     /ixzN(iix,iiy,iiz)/dt
                   ddzdx=(drz_z(iix+1,iiy,iiz)-drz_r(iix-1,iiy,iiz))
     1                     /ixzN(iix,iiy,iiz)/dt
                   ddzdy=(drz_z(iix,iiy+1,iiz)-drz_r(iix,iiy-1,iiz))
     1                     /ixzN(iix,iiy,iiz)/dt
                   ddzdz=(drz_z(iix,iiy,iiz+1)-drz_r(iix,iiy,iiz-1))
     1                     /ixzN(iix,iiy,iiz)/dt
                   write(*,*) ddxdx,ddxdy,ddxdz
                   write(*,*) ddydx,ddydy,ddydz
                   write(*,*) ddzdx,ddzdy,ddzdz
                endif
             else
                deN=0.
                avfnN=0.
                avftN=0.
                avrfnN=0.
                avrftN=0.
                avdnN=0.
                avdtN=0.
                avdnN2=0.
                avdtN2=0.
                vvrN=0.
                vvpN=0.
                vvzN=0.
                vvr2N=0.
                vvp2N=0.
                vvz2N=0.
                wwrN=0.
                wwpN=0.
                wwzN=0.
                pwrN=0.
                pwpN=0.
                pwzN=0.
                drx=0.
                dry=0.
                drz=0.
             endif
c            write(*,*) 'dxyz: ',iix, iiy, iiz, dt, drx, dry, drz
c
             if(p2.gt.p1*p1)then
                pstdev=sqrt(p2-p1*p1)
             else
                pstdev=0.
             endif
c
             if(xnu2.gt.xnu*xnu)then
                xnustdev=sqrt(xnu2-xnu*xnu)
             else
                xnustdev=0.
             endif
c
             tek=ek(iix,iiy,iiz)/dble(ifcc)
             tepn=epn(iix,iiy,iiz)/dble(ifcc)
             tept=ept(iix,iiy,iiz)/dble(ifcc)
c            write(*,*) 'de,ek,ep: ',deN,tek,tepn,tept

             xCnu_av=xCnu(iix,iiy,iiz)/dble(ifcc)/Vxyz
             vxCnu_av=vxCnu(iix,iiy,iiz)/dble(ifcc)/Vxyz
             vyCnu_av=vyCnu(iix,iiy,iiz)/dble(ifcc)/Vxyz
             vzCnu_av=vzCnu(iix,iiy,iiz)/dble(ifcc)/Vxyz
c
             cc0n=icc0(iix,iiy,iiz)/dble(ifcc)
             cc3n=icc3(iix,iiy,iiz)/dble(ifcc)
             cc4n=icc4(iix,iiy,iiz)/dble(ifcc)
             cc6n=icc6(iix,iiy,iiz)/dble(ifcc)
             cc9n=icc9(iix,iiy,iiz)/dble(ifcc)
             ccxn=iccx(iix,iiy,iiz)/dble(ifcc)
c
             cc3c=jcc3(iix,iiy,iiz)/dble(ifcc)
             cc4c=jcc4(iix,iiy,iiz)/dble(ifcc)
             cc6c=jcc6(iix,iiy,iiz)/dble(ifcc)
             cc9c=jcc9(iix,iiy,iiz)/dble(ifcc)
             ccxc=jccx(iix,iiy,iiz)/dble(ifcc)
c
c            write(*,*) iix,iiy,iiz,xpos,ypos,zpos,
c    1                  icc3(iix,iiy,iiz),ifcc
c
             pw1=pnum1(iix,iiy,iiz)/dble(ifcc)
             cw1=cnum1(iix,iiy,iiz)/dble(ifcc)
             if((cw1.gt.0.0).and.(ifcc.gt.0))then
                cw1inv=1./cw1
c OLD?          cw1inv=1./cw1/dble(ifcc)
                if(apnum1(iix,iiy,iiz).gt.0)then
                   aam1=aamom1(iix,iiy,iiz)/apnum1(iix,iiy,iiz)
                   aam2=aamom2(iix,iiy,iiz)/apnum1(iix,iiy,iiz)
                   aam3=aamom3(iix,iiy,iiz)/apnum1(iix,iiy,iiz)
                   aam4=aamom4(iix,iiy,iiz)/apnum1(iix,iiy,iiz)
                   aam5=aamom5(iix,iiy,iiz)/apnum1(iix,iiy,iiz)
                else
                   aam1=0
                   aam2=0
                   aam3=0
                   aam4=0
                   aam5=0
                endif
                if(pnum1(iix,iiy,iiz).gt.0)then
                   am1=amom1(iix,iiy,iiz)/pnum1(iix,iiy,iiz)
                   am2=amom2(iix,iiy,iiz)/pnum1(iix,iiy,iiz)
                   am3=amom3(iix,iiy,iiz)/pnum1(iix,iiy,iiz)
                   am4=amom4(iix,iiy,iiz)/pnum1(iix,iiy,iiz)
                   am5=amom5(iix,iiy,iiz)/pnum1(iix,iiy,iiz)
                else
                   am1=0
                   am2=0
                   am3=0
                   am4=0
                   am5=0
                endif
c
                ee1=(e1(1,1,iix,iiy,iiz)
     1           +e1(2,2,iix,iiy,iiz)
     2           +e1(3,3,iix,iiy,iiz))/3.0*cw1inv
                eet1=(et1(1,1,iix,iiy,iiz)
     1           +et1(2,2,iix,iiy,iiz)
     2           +et1(3,3,iix,iiy,iiz))/3.0*cw1inv
             else
                cw1inv=0
                am1=0
                am2=0
                am3=0
                am4=0
                am5=0
                ee1=0
                eet1=0
             endif
c            write(*,*) iix,iiy,iiz,cw1,am1,am2,am3,am4,am5
c
c ... compute strain from the differential displacement data
c     x-derivatives
             if(iix.gt.0)then
                m1=1
             else
                m1=0
             endif
c
             if((dt.ne.0.d0).and.(ixzN(iix-m1,iiy,iiz).gt.0))then
                drx1=vxCnu(iix-m1,iiy,iiz)/xCnu(iix-m1,iiy,iiz)
                dry1=vyCnu(iix-m1,iiy,iiz)/xCnu(iix-m1,iiy,iiz)
                drz1=vzCnu(iix-m1,iiy,iiz)/xCnu(iix-m1,iiy,iiz)
c               drx1=drz_r(iix-m1,iiy,iiz)/ixzN(iix-m1,iiy,iiz)/dt
c               dry1=drz_p(iix-m1,iiy,iiz)/ixzN(iix-m1,iiy,iiz)/dt
c               drz1=drz_z(iix-m1,iiy,iiz)/ixzN(iix-m1,iiy,iiz)/dt
             else
                drx1=0.
                dry1=0.
                drz1=0.
             endif
c
             if(iix.lt.ixbin)then
                m1=1
             else
                m1=0
             endif
c
             if((dt.ne.0.d0).and.(ixzN(iix+m1,iiy,iiz).gt.0))then
                drx2=vxCnu(iix+m1,iiy,iiz)/xCnu(iix+m1,iiy,iiz)
                dry2=vyCnu(iix+m1,iiy,iiz)/xCnu(iix+m1,iiy,iiz)
                drz2=vzCnu(iix+m1,iiy,iiz)/xCnu(iix+m1,iiy,iiz)
c               drx2=drz_r(iix+m1,iiy,iiz)/ixzN(iix+m1,iiy,iiz)/dt
c               dry2=drz_p(iix+m1,iiy,iiz)/ixzN(iix+m1,iiy,iiz)/dt
c               drz2=drz_z(iix+m1,iiy,iiz)/ixzN(iix+m1,iiy,iiz)/dt
             else
                drx2=0.
                dry2=0.
                drz2=0.
             endif
c
c ... WARNING ERROR - CHECK WHY NOT DIVIDED BY 2
c     NEW version - division by 2
             dxbin2=dxbin*2.0
c
             es11=(drx2-drx1)/dxbin2
             es21=(dry2-dry1)/dxbin2
             es31=(drz2-drz1)/dxbin2
c
c     y-derivatives
             if(iiy.gt.0)then
                m1=1
             else
                m1=0
             endif
c
             if((dt.ne.0.d0).and.(ixzN(iix,iiy-m1,iiz).gt.0))then
                drx1=vxCnu(iix,iiy-m1,iiz)/xCnu(iix,iiy-m1,iiz)
                dry1=vyCnu(iix,iiy-m1,iiz)/xCnu(iix,iiy-m1,iiz)
                drz1=vzCnu(iix,iiy-m1,iiz)/xCnu(iix,iiy-m1,iiz)
c               drx1=drz_r(iix,iiy-m1,iiz)/ixzN(iix,iiy-m1,iiz)/dt
c               dry1=drz_p(iix,iiy-m1,iiz)/ixzN(iix,iiy-m1,iiz)/dt
c               drz1=drz_z(iix,iiy-m1,iiz)/ixzN(iix,iiy-m1,iiz)/dt
             else
                drx1=0.
                dry1=0.
                drz1=0.
             endif
c
             if(iiy.lt.iybin)then
                m1=1
             else
                m1=0
             endif
c
             if((dt.ne.0.d0).and.(ixzN(iix,iiy+m1,iiz).gt.0))then
                drx2=vxCnu(iix,iiy+m1,iiz)/xCnu(iix,iiy+m1,iiz)
                dry2=vyCnu(iix,iiy+m1,iiz)/xCnu(iix,iiy+m1,iiz)
                drz2=vzCnu(iix,iiy+m1,iiz)/xCnu(iix,iiy+m1,iiz)
c               drx2=drz_r(iix,iiy+m1,iiz)/ixzN(iix,iiy+m1,iiz)/dt
c               dry2=drz_p(iix,iiy+m1,iiz)/ixzN(iix,iiy+m1,iiz)/dt
c               drz2=drz_z(iix,iiy+m1,iiz)/ixzN(iix,iiy+m1,iiz)/dt
             else
                drx2=0.
                dry2=0.
                drz2=0.
             endif
c
             dybin2=2.0*dybin
             es12=(drx2-drx1)/dybin2
             es22=(dry2-dry1)/dybin2
             es32=(drz2-drz1)/dybin2
c
c     z-derivatives
             if(iiz.gt.0)then
                m1=1
             else
                m1=0
             endif
c
             if((dt.ne.0.d0).and.(ixzN(iix,iiy,iiz-m1).gt.0))then
                drx1=vxCnu(iix,iiy,iiz-m1)/xCnu(iix,iiy,iiz-m1)
                dry1=vyCnu(iix,iiy,iiz-m1)/xCnu(iix,iiy,iiz-m1)
                drz1=vzCnu(iix,iiy,iiz-m1)/xCnu(iix,iiy,iiz-m1)
c               drx1=drz_r(iix,iiy,iiz-m1)/ixzN(iix,iiy,iiz-m1)/dt
c               dry1=drz_p(iix,iiy,iiz-m1)/ixzN(iix,iiy,iiz-m1)/dt
c               drz1=drz_z(iix,iiy,iiz-m1)/ixzN(iix,iiy,iiz-m1)/dt
             else
                drx1=0.
                dry1=0.
                drz1=0.
             endif
c
             if(iiz.lt.izbin)then
                m1=1
             else
                m1=0
             endif
c
             if((dt.ne.0.d0).and.(ixzN(iix,iiy,iiz+m1).gt.0))then
                drx2=vxCnu(iix,iiy,iiz+m1)/xCnu(iix,iiy,iiz+m1)
                dry2=vyCnu(iix,iiy,iiz+m1)/xCnu(iix,iiy,iiz+m1)
                drz2=vzCnu(iix,iiy,iiz+m1)/xCnu(iix,iiy,iiz+m1)
c               drx2=drz_r(iix,iiy,iiz+m1)/ixzN(iix,iiy,iiz+m1)/dt
c               dry2=drz_p(iix,iiy,iiz+m1)/ixzN(iix,iiy,iiz+m1)/dt
c               drz2=drz_z(iix,iiy,iiz+m1)/ixzN(iix,iiy,iiz+m1)/dt
             else
                drx2=0.
                dry2=0.
                drz2=0.
             endif
c
             dzbin2=2.0*dzbin
             es13=(drx2-drx1)/dzbin2
             es23=(dry2-dry1)/dzbin2
             es33=(drz2-drz1)/dzbin2
c
c      trace of derivatives
             esV=(es11+es22+es33)/3.
c      normal stiffness tensor entries that are not necessary independent
c      column numbers in tav.data are 120+index of the Cn(index) field
             Cn13(1) =Cn1(1,1,1,1,iix,iiy,iiz)
             Cn13(2) =Cn1(1,1,2,2,iix,iiy,iiz)
             Cn13(3) =Cn1(1,1,3,3,iix,iiy,iiz)
             Cn13(4) =Cn1(1,1,2,3,iix,iiy,iiz)
             Cn13(5) =Cn1(1,1,1,3,iix,iiy,iiz)
             Cn13(6) =Cn1(1,1,1,2,iix,iiy,iiz)
             Cn13(7) =Cn1(2,2,2,2,iix,iiy,iiz)
             Cn13(8) =Cn1(2,2,3,3,iix,iiy,iiz)
             Cn13(9) =Cn1(2,2,2,3,iix,iiy,iiz)
             Cn13(10)=Cn1(2,2,1,3,iix,iiy,iiz)
             Cn13(11)=Cn1(2,2,1,2,iix,iiy,iiz)
             Cn13(12)=Cn1(3,3,3,3,iix,iiy,iiz)
             Cn13(13)=Cn1(3,3,2,3,iix,iiy,iiz)
             Cn13(14)=Cn1(3,3,1,3,iix,iiy,iiz)
             Cn13(15)=Cn1(3,3,1,2,iix,iiy,iiz)
             Cn13(16)=Cn1(2,3,2,3,iix,iiy,iiz)
             Cn13(17)=Cn1(2,3,1,3,iix,iiy,iiz)
             Cn13(18)=Cn1(2,3,1,2,iix,iiy,iiz)
             Cn13(19)=Cn1(1,3,1,3,iix,iiy,iiz)
             Cn13(20)=Cn1(1,3,1,2,iix,iiy,iiz)
             Cn13(21)=Cn1(1,2,1,2,iix,iiy,iiz) 
c
c      tangential stiffness tensor entries that are not necessary independent
c      column numbers in tav.data are 150+index of the Ct(index) field
             Ct13(1) =Ct1(1,1,1,1,iix,iiy,iiz)
             Ct13(2) =Ct1(1,1,2,2,iix,iiy,iiz)
             Ct13(3) =Ct1(1,1,3,3,iix,iiy,iiz)
             Ct13(4) =Ct1(1,1,2,3,iix,iiy,iiz)
             Ct13(5) =Ct1(1,1,3,2,iix,iiy,iiz)
             Ct13(6) =Ct1(1,1,1,3,iix,iiy,iiz)
             Ct13(7) =Ct1(1,1,3,1,iix,iiy,iiz)
             Ct13(8) =Ct1(1,1,1,2,iix,iiy,iiz)
             Ct13(9) =Ct1(1,1,2,1,iix,iiy,iiz)
             Ct13(10)=Ct1(2,2,2,2,iix,iiy,iiz)
             Ct13(11)=Ct1(2,2,3,3,iix,iiy,iiz)
             Ct13(12)=Ct1(2,2,2,3,iix,iiy,iiz)
             Ct13(13)=Ct1(2,2,3,2,iix,iiy,iiz)
             Ct13(14)=Ct1(2,2,1,3,iix,iiy,iiz)
             Ct13(15)=Ct1(2,2,3,1,iix,iiy,iiz)
             Ct13(16)=Ct1(2,2,1,2,iix,iiy,iiz)
             Ct13(17)=Ct1(2,2,2,1,iix,iiy,iiz)
             Ct13(18)=Ct1(3,3,3,3,iix,iiy,iiz)
             Ct13(19)=Ct1(3,3,2,3,iix,iiy,iiz)
             Ct13(20)=Ct1(3,3,3,2,iix,iiy,iiz)
             Ct13(21)=Ct1(3,3,1,3,iix,iiy,iiz)
             Ct13(22)=Ct1(3,3,3,1,iix,iiy,iiz)
             Ct13(23)=Ct1(3,3,1,2,iix,iiy,iiz)
             Ct13(24)=Ct1(3,3,2,1,iix,iiy,iiz)
             Ct13(25)=Ct1(2,3,2,3,iix,iiy,iiz)
             Ct13(26)=Ct1(2,3,3,2,iix,iiy,iiz)
             Ct13(27)=Ct1(2,3,1,3,iix,iiy,iiz)
             Ct13(28)=Ct1(2,3,3,1,iix,iiy,iiz)
             Ct13(29)=Ct1(2,3,1,2,iix,iiy,iiz)
             Ct13(30)=Ct1(2,3,2,1,iix,iiy,iiz)
             Ct13(31)=Ct1(3,2,3,2,iix,iiy,iiz)
             Ct13(32)=Ct1(3,2,1,3,iix,iiy,iiz)
             Ct13(33)=Ct1(3,2,3,1,iix,iiy,iiz)
             Ct13(34)=Ct1(3,2,1,2,iix,iiy,iiz)
             Ct13(35)=Ct1(3,2,2,1,iix,iiy,iiz)
             Ct13(36)=Ct1(1,3,1,3,iix,iiy,iiz)
             Ct13(37)=Ct1(1,3,3,1,iix,iiy,iiz)
             Ct13(38)=Ct1(1,3,1,2,iix,iiy,iiz)
             Ct13(39)=Ct1(1,3,2,1,iix,iiy,iiz)
             Ct13(40)=Ct1(3,1,3,1,iix,iiy,iiz)
             Ct13(41)=Ct1(3,1,1,2,iix,iiy,iiz)
             Ct13(42)=Ct1(3,1,2,1,iix,iiy,iiz)
             Ct13(43)=Ct1(1,2,1,2,iix,iiy,iiz)
             Ct13(44)=Ct1(1,2,2,1,iix,iiy,iiz)
             Ct13(45)=Ct1(2,1,2,1,iix,iiy,iiz)              
c
         chk_isum=chk_isum+1
         chk_xnu =chk_xnu +xnu
         chk_xCnu=chk_xCnu+xCnu_av
c
         if(ioutput.eq.1)then   
c
           CZ0=-1
           CZ3=-1
           CZ4=-1
           CZ6=-1
c
           if((cc0n+cc3n+cc4n+cc6n+cc9n+ccxn).gt.0)
     1        CZ0=(cc3c+cc4c+cc6c+cc9c+ccxc)/
     2            (cc0n+cc3n+cc4n+cc6n+cc9n+ccxn)
           if((cc3n+cc4n+cc6n+cc9n+ccxn).gt.0)
     1        CZ3=(cc3c+cc4c+cc6c+cc9c+ccxc)/
     2            (cc3n+cc4n+cc6n+cc9n+ccxn)
           if((cc4n+cc6n+cc9n+ccxn).gt.0)
     1        CZ4=(cc4c+cc6c+cc9c+ccxc)/
     2            (cc4n+cc6n+cc9n+ccxn)
           if((cc6n+cc9n+ccxn).gt.0)
     1        CZ6=(cc6c+cc9c+ccxc)/
     2            (cc6n+cc9n+ccxn)
c
           write(81,'(3i6,194g16.8,98g16.8)') 
     1       iix,iiy,iiz,xpos,ypos,zpos,Vxyz,
     2       ixzN(iix,iiy,iiz)/dble(ifcF),
     3       ixzC(iix,iiy,iiz)/dble(ifcc),
     4       xnu,xnustdev,xCnu_av,
     6       cc0n,cc3n,cc4n,cc6n,cc9n,ccxn,p1,pstdev,
     7          (( s1(n,m,iix,iiy,iiz)*Vifcc, m=1,3), n=1,3),
     8       pt,(( st1(n,m,iix,iiy,iiz)*Vifcc, m=1,3), n=1,3),
     9       pd,(( sd1(n,m,iix,iiy,iiz)*Vifcc, m=1,3), n=1,3),
     1       time,vvrN,vvpN,vvzN,vvr2N,vvp2N,vvz2N,drx,dry,drz,
     2       esV,es11,es12,es13,es21,es22,es23,es31,es32,es33,
     3       trf,(( f1(n,m,iix,iiy,iiz)*Vifcc, m=1,3), n=1,3),
     4       ee1,((e1(n,m,iix,iiy,iiz)*cw1inv,m=1,3),n=1,3),
     5       eet1,((et1(n,m,iix,iiy,iiz)*cw1inv,m=1,3),n=1,3),
     6       pw1,cw1,am1,am2,am3,am4,am5,sys_Sx,sys_Sy,sys_Sz,
     7       0.,wwrN,wwpN,wwzN,vxCnu_av,vyCnu_av,vzCnu_av,
     8       sys_Sx0,sys_Sy0,sys_Sz0,dble(itd),
     8                                  (Cn13(n1)*Vifcc,n1=1,21),
     9       0.,0.,0.,0.,0.,0.,0.,0.,0.,(Ct13(n1)*Vifcc,n1=1,45),
     1       0.,xnur,NK,inumc,(inc0(i),i=0,19),
     2                        (inc1(i),i=0,19),
     3       cc0n,0.,0.,cc3n,cc4n,0.,cc6n,cc9n,ccxn,0.,0.,
     4            0.,0.,cc3c,cc4c,0.,cc6c,cc9c,ccxc,0.,
     5       CZ0,0.,0.,CZ3,CZ4,0.,CZ6,0.,0.,CZ00,0.,
     6       aam1,aam2,aam3,aam4,aam5,deN,tek,tepn,tept,0.,
     7       avfnN,avftN,avdnN,avdtN,avdnN2,avdtN2,avrfnN,avrftN,
     8       0.,0.,pwrN,pwpN,pwzN
c
         else
          write(81,'(3i6,135g16.8)') 
     1       iix,iiy,iiz,xpos,ypos,zpos,Vxyz,
     2       ixzN(iix,iiy,iiz)/dble(ifcF),
     3       ixzC(iix,iiy,iiz)/dble(ifcc),
     4       xnu,xnustdev,xCnu_av,
     6       cc0n,cc3n,cc4n,cc6n,cc9n,ccxn,p1,pstdev,
     7          (( s1(n,m,iix,iiy,iiz)*Vifcc, m=1,3), n=1,3),
     8       pt,(( st1(n,m,iix,iiy,iiz)*Vifcc, m=1,3), n=1,3),
     9       pd,(( sd1(n,m,iix,iiy,iiz)*Vifcc, m=1,3), n=1,3),
     1       time,vvrN,vvpN,vvzN,vvr2N,vvp2N,vvz2N,drx,dry,drz,
     2       esV,es11,es12,es13,es21,es22,es23,es31,es32,es33,
     3       trf,(( f1(n,m,iix,iiy,iiz)*Vifcc, m=1,3), n=1,3),
     4       ee1,((e1(n,m,iix,iiy,iiz)*cw1inv,m=1,3),n=1,3),
     5       eet1,((et1(n,m,iix,iiy,iiz)*cw1inv,m=1,3),n=1,3),
     6       pw1,cw1,am1,am2,am3,am4,am5,sys_Sx,sys_Sy,sys_Sz,
     7       0.,wwrN,wwpN,wwzN,vxCnu_av,vyCnu_av,vzCnu_av,
     8       sys_Sx0,sys_Sy0,sys_Sz0,dble(itd),deN,tek,tepn,tept,0,
     9       avfnN,avftN,avdnN,avdtN,avdnN2,avdtN2,avrfnN,avrftN,
     1       0.,0.,pwrN,pwpN,pwzN
c
         endif 
c
1102       continue
c
c ... uncomment if all 3Dblocks should be separate
c       comment if gnuplot desired with contour
c          write(81,'(a,$)') char(10)
1101     continue
c
c ... create an empty line here (if over-all averaging is on)
c     this is appropriate for space-space-plots
         if(iout_el.eq.1)
     1      write(81,'(a,$)') char(10)

1100   continue
c
c ... create an empty line here (if snap-shot averaging is on)
c     this is appropriate for space-time-plots
       if(iout_el.eq.0)
     1    write(81,'(a,$)') char(10)
c
       write(*,'(a,i8,3g16.6)') '# density check: ',
     1      NK,chk_isum,chk_xnu/chk_isum,chk_xCnu/chk_isum
       write(*,'(a,2i8,3g16.6)') '# control number: ',
     1      NK,ifcc,pcontrol/ifcc,p0control/ifcc,
     2              (pcontrol-p0control)/ifcc
c
c .. write dr2.data
       write(92,'(a)')
     1   '# ix iz it-1 x z time ic dr1,2 dz1,2 drp1,2 dxy2 dxyz2'
       do 300 ix=0,IZBINMAX
        do 301 iz=0,IZBINMAX
          do 302 it=2,itc
             x=xmin+(ix+0.5)*(xmax-xmin)/(ixbin+1)
             z=zmin+(iz+0.5)*(zmax-zmin)/(izbin+1)
             if(icdr2(ix,iz,it).gt.0)then
               write(92,'(3i5,3g18.8,i6,8g18.8)')
     1           ix,iz,it-1,x,z,ttime(it)-ttime(1),icdr2(ix,iz,it),
     2  dr1r(ix,iz,it)/icdr2(ix,iz,it),dr2r(ix,iz,it)/icdr2(ix,iz,it),
     3  dr1z(ix,iz,it)/icdr2(ix,iz,it),dr2z(ix,iz,it)/icdr2(ix,iz,it),
     4  dr1rp(ix,iz,it)/icdr2(ix,iz,it),dr2rp(ix,iz,it)/icdr2(ix,iz,it),
     5  dr2xy(ix,iz,it)/icdr2(ix,iz,it),dr2(ix,iz,it)/icdr2(ix,iz,it)
             else
                write(92,'(3i5,3g18.8,i6,8g18.8)')
     1             ix,iz,it-1,x,z,ttime(it)-ttime(1),0,
     2             0,0,0,0,0,0,0,0
             endif
302       continue
          write(92,'(a,$)') char(10)
301    continue
300    continue
c
c ... output of displacement data for xballs visualisation
       dt=time-time_old
c      write(*,*) dt,time,itc
       write(93,'(i9,7g16.7)') 
     1      NK,time,wx10,wy10,wz10,wx11,wy11,wz11
       do 400 i=0,NK-1
          if(dt.ne.0.d0)then
          write(93,'(8g16.7)') xj(i),yj(i),zj(i),
     1      dr(i)/(ijc(i))/dt,dphi(i)/(ijc(i))/dt,
     2      dz(i)/(ijc(i))/dt,prad(i),pxi(i)
          else
          write(93,'(8g16.7)') xj(i),yj(i),zj(i),
     1      0, 0, 0, prad(i),pxi(i)
          endif
400    continue

c
c ... end read in if time > tmax
       if(time.gt.tmax) goto 9001
c
9000   continue
9001   continue
       close(81)
       close(92)
       close(93)

       end

       subroutine faverage(ixb,iyb,izb, radius, cnx,cny,cnz, xw,f)
c
       include 'tav10_para.h'

       integer ixb, iyb, izb
       double precision radius,cnx,cny,cnz
       double precision f(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision cn(3)
       double precision PI
       cn(1)=cnx
       cn(2)=cny
       cn(3)=cnz
       PI=asin(1.d0)*2.d0

       do 10 iii=1,3
          do 11 jjj=1,3
             f(iii,jjj,ixb,iyb,izb)=f(iii,jjj,ixb,iyb,izb)
     1               +cn(iii)*cn(jjj)*(4./3.*PI)*radius**3*xw
11     continue
10     continue

       return
       end

       subroutine saverage(ixb,iyb,izb, cfn, cfl, cnx,cny,cnz, xw,s)
c
       include 'tav10_para.h'

       integer ixb, iyb, izb
       double precision cfn,cfl,cnx,cny,cnz
       double precision s(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision cn(3)
       cn(1)=cnx
       cn(2)=cny
       cn(3)=cnz

       do 10 iii=1,3
          do 11 jjj=1,3
             s(iii,jjj,ixb,iyb,izb)=s(iii,jjj,ixb,iyb,izb)
     1                             +cn(iii)*cn(jjj)*cfn*cfl*xw
c           write(*,*) iii,jjj,ixb,iyb,izb,cnx,cny,cnz,cfn,cfl
c           write(*,*) s(iii,jjj,ixb,iyb,izb),cn(iii)*cn(jjj)*cfn*cfl
11     continue
10     continue

       return
       end

       subroutine staverage(ixb, iyb, izb, cft, cfl, 
     1                      cnx, cny, cnz, ctx, cty, ctz, xw, s)
c
       include 'tav10_para.h'

       integer ixb, iyb, izb
       double precision cft,cfl, cnx,cny,cnz, ctx,cty,ctz
       double precision s(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision cn(3),ct(3)
       cn(1)=cnx
       cn(2)=cny
       cn(3)=cnz
       ct(1)=ctx
       ct(2)=cty
       ct(3)=ctz

       do 10 iii=1,3
          do 11 jjj=1,3
             s(iii,jjj,ixb,iyb,izb)=s(iii,jjj,ixb,iyb,izb)
     1                             +cn(iii)*ct(jjj)*cft*cfl*xw
c           write(*,*) iii,jjj,ixb,iyb,izb,cnx,cny,cnz,cft,cfl
c           write(*,*) s(iii,jjj,ixb,iyb,izb),cn(iii)*cn(jjj)*cft*cfl
11     continue
10     continue

       return
       end

       subroutine eaverage(ixb,iyb,izb, cfd, cfl, cnx,cny,cnz, xw, e)
c
       include 'tav10_para.h'

       integer ixb, iyb, izb
       double precision cfd,cfl,cnx,cny,cnz
       double precision e(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision cn(3)
       cn(1)=cnx
       cn(2)=cny
       cn(3)=cnz

       do 10 iii=1,3
          do 11 jjj=1,3
             e(iii,jjj,ixb,iyb,izb)=e(iii,jjj,ixb,iyb,izb)
     1                             +cn(iii)*cn(jjj)*cfd/cfl*xw
c           write(*,*) iii,jjj,ixb,iyb,izb,cnx,cny,cnz,cfn,cfl
c           write(*,*) s(iii,jjj,ixb,iyb,izb),cn(iii)*cn(jjj)*cfn*cfl
11     continue
10     continue

       return
       end

       subroutine sdaverage(ixb, iyb, izb, cvx, cvy, cvz, cmass, s)
c
       include 'tav10_para.h'

       integer ixb, iyb, izb
       double precision cvx,cvy,cvz,cmass
       double precision s(3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision cv(3)
       cv(1)=cvx
       cv(2)=cvy
       cv(3)=cvz

       do 10 iii=1,3
          do 11 jjj=1,3
             s(iii,jjj,ixb,iyb,izb)=s(iii,jjj,ixb,iyb,izb)
     1                             +cv(iii)*cv(jjj)*cmass
c           write(*,*) iii,jjj,ixb,iyb,izb,cnx,cny,cnz,cfn,cfl
c           write(*,*) s(iii,jjj,ixb,iyb,izb),cn(iii)*cn(jjj)*cfn*cfl
11     continue
10     continue

       return
       end

c
c ... the stiffness tensor has to be multiplied by k^n !
c
       subroutine Caverage(ixb,iyb,izb, cfl, cnx,cny,cnz, xw,C)
c
       include 'tav10_para.h'

       integer iii, jjj, kkk, lll
       integer ixb, iyb, izb
       double precision cfl, cnx,cny,cnz
       double precision C(3,3,3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision cn(3)
       cn(1)=cnx
       cn(2)=cny
       cn(3)=cnz

       do 10 iii=1,3
          do 11 jjj=1,3
             do 12 kkk=1,3
                do 13 lll=1,3
                   C(iii,jjj,kkk,lll,ixb,iyb,izb)
     1                =C(iii,jjj,kkk,lll,ixb,iyb,izb)
     2                +cn(iii)*cn(jjj)*cn(kkk)*cn(lll)*cfl*cfl*xw
c           write(*,*) iii,jjj,ixb,iyb,izb,cnx,cny,cnz,cfl
c           write(*,*) C(iii,jjj,kkk,lll,ixb,iyb,izb)
13     continue
12     continue
11     continue
10     continue

       return
       end

       subroutine Ctaverage(ixb, iyb, izb, cfl,
     1                      cnx, cny, cnz, ctx, cty, ctz, xw, D)
c
       include 'tav10_para.h'

       integer iii, jjj, kkk, lll
       integer ixb, iyb, izb
       double precision cfl, cnx,cny,cnz,ctx,cty,ctz
       double precision D(3,3,3,3,-1:IXBINMAX,-1:IYBINMAX,-1:IZBINMAX)
       double precision cn(3),ct(3)
       cn(1)=cnx
       cn(2)=cny
       cn(3)=cnz
       ct(1)=ctx
       ct(2)=cty
       ct(3)=ctz

       do 10 iii=1,3
          do 11 jjj=1,3
             do 12 kkk=1,3
                do 13 lll=1,3
                   D(iii,jjj,kkk,lll,ixb,iyb,izb)
     1                =D(iii,jjj,kkk,lll,ixb,iyb,izb)
     2                +cn(iii)*ct(jjj)*cn(kkk)*ct(lll)*cfl*cfl*xw
c           write(*,*) iii,jjj,ixb,iyb,izb,cnx,cny,cnz,cfl
c           write(*,*) C(iii,jjj,kkk,lll,ixb,iyb,izb)
13     continue
12     continue
11     continue
10     continue

       return
       end


