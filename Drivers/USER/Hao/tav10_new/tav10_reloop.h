c
c ... reset =0 the averaging fields which are used per loop
c
       do 2001 ixp=-1,IXBINMAX
         do 2002 iyp=-1,IYBINMAX
           do 2003 izp=-1,IZBINMAX
             do 2004 i1=1,3
               do 2005 i2=1,3
                 f(i1,i2,ixp,iyp,izp)=0.d0
                 e(i1,i2,ixp,iyp,izp)=0.d0
                 et(i1,i2,ixp,iyp,izp)=0.d0
                 s(i1,i2,ixp,iyp,izp)=0.d0
                 sd(i1,i2,ixp,iyp,izp)=0.d0
                 st(i1,i2,ixp,iyp,izp)=0.d0
                 do 2006 i3=1,3
                   do 2007 i4=1,3
                     Cn(i1,i2,i3,i4,ixp,iyp,izp)=0.d0
                     Ct(i1,i2,i3,i4,ixp,iyp,izp)=0.d0
2007   continue
2006   continue
2005   continue
2004   continue
2003   continue
2002   continue
2001   continue

