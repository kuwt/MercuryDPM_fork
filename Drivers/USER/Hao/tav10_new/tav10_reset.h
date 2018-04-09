      write(*,'(a,$)') '# reset data-structures ...' 
c
      do 1111 m1=-1,IXBINMAX
         do 1112 m2=-1,IYBINMAX
            do 1113 m3=-1,IZBINMAX
               ixzC(m1,m2,m3)=0
               ixzN(m1,m2,m3)=0
               icc0(m1,m2,m3)=0
               icc3(m1,m2,m3)=0
               icc4(m1,m2,m3)=0
               icc6(m1,m2,m3)=0
               icc9(m1,m2,m3)=0
               iccx(m1,m2,m3)=0
               jcc3(m1,m2,m3)=0
               jcc4(m1,m2,m3)=0
               jcc6(m1,m2,m3)=0
               jcc9(m1,m2,m3)=0
               jccx(m1,m2,m3)=0
               xznu(m1,m2,m3)=0
               xCnu(m1,m2,m3)=0
               vxCnu(m1,m2,m3)=0
               vyCnu(m1,m2,m3)=0
               vzCnu(m1,m2,m3)=0
               xznur(m1,m2,m3)=0
               xznu1(m1,m2,m3)=0
               xznu2(m1,m2,m3)=0
c
               apnum1(m1,m2,m3)=0
               pnum1(m1,m2,m3)=0
               cnum1(m1,m2,m3)=0
               aamom1(m1,m2,m3)=0
               aamom2(m1,m2,m3)=0
               aamom3(m1,m2,m3)=0
               aamom4(m1,m2,m3)=0
               aamom5(m1,m2,m3)=0
               amom1(m1,m2,m3)=0
               amom2(m1,m2,m3)=0
               amom3(m1,m2,m3)=0
               amom4(m1,m2,m3)=0
               amom5(m1,m2,m3)=0
               do 1115 n1=1,3
                  do 1116 n2=1,3
                     f(n1,n2,m1,m2,m3)=0
                     f1(n1,n2,m1,m2,m3)=0
                     e(n1,n2,m1,m2,m3)=0
                     e1(n1,n2,m1,m2,m3)=0
                     et(n1,n2,m1,m2,m3)=0
                     et1(n1,n2,m1,m2,m3)=0
                     s(n1,n2,m1,m2,m3)=0
                     s1(n1,n2,m1,m2,m3)=0
                     s2(n1,n2,m1,m2,m3)=0
                     sd(n1,n2,m1,m2,m3)=0
                     sd1(n1,n2,m1,m2,m3)=0
                     st(n1,n2,m1,m2,m3)=0
                     st1(n1,n2,m1,m2,m3)=0
                     do 1117 n3=1,3
                        do 1118 n4=1,3
                           Cn1(n1,n2,n3,n4,m1,m2,m3)=0
                           Ct1(n1,n2,n3,n4,m1,m2,m3)=0
1118                    continue
1117                 continue
1116              continue
1115           continue
               ek(m1,m2,m3)=0 
               epn(m1,m2,m3)=0 
               ept(m1,m2,m3)=0 
               de(m1,m2,m3)=0 
               avfn(m1,m2,m3)=0 
               avft(m1,m2,m3)=0 
               avrfn(m1,m2,m3)=0 
               avrft(m1,m2,m3)=0 
               avdn(m1,m2,m3)=0 
               avdt(m1,m2,m3)=0 
               avdn2(m1,m2,m3)=0 
               avdt2(m1,m2,m3)=0 
               vr(m1,m2,m3)=0 
               vp(m1,m2,m3)=0 
               vz(m1,m2,m3)=0 
               wr(m1,m2,m3)=0 
               wp(m1,m2,m3)=0 
               wz(m1,m2,m3)=0 
               pwr(m1,m2,m3)=0 
               pwp(m1,m2,m3)=0 
               pwz(m1,m2,m3)=0 
               vr2(m1,m2,m3)=0 
               vp2(m1,m2,m3)=0 
               vz2(m1,m2,m3)=0 
               drz_r(m1,m2,m3)=0 
               drz_p(m1,m2,m3)=0 
               drz_z(m1,m2,m3)=0 
1113        continue
1112     continue
1111  continue

      do 2100 m1=0,IPCMAX
         inc0(m1)=0
         inc1(m1)=0
2100  continue
c
       write(*,'(a)') '. finished.' 

