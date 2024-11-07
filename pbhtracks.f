      character*5 pbh(10)
      dimension x(1000),y(1000),s(12),t(12),l(10),xm(1000,10),start(10)
      dimension kend(10),k(10),im(10),!xc(10),yc(10),
     1 c(2),xx(300),yy(300),tart(10),xart(10)
      data pbh/'m23','m21','m19','m17','m15','m13','m1 ','m9','m7','m5'/
      data start/-23.,-21.,-19.,-17.,-15.,-13.,-1.,-9.,-7.,-5./,j/0/
      data im/10*0/,xm/10000*0./
      b=alog10(8.62)!Boltzmann
      teq=57.*3.14!teq4
      const=10.125/6.67/3.14159!mass loss
      zeq=alog10(3500.)!Planck
      call pgbegin(0,'?',1,1)
      call pgsch(1.6)
      do i=1,10
      tart(i)=start(i)
      start(i)=alog10(6.62e-27*27.e30/6.67e23/8./1.38e-16/3.14159**2)
     1 +b-5.-start(i)!Hawking temperature
      sch=tart(i)+alog10(6.67e-10/9.e20*2.e33)
      t(1)=alog10(8./3.e10)+sch!r=ct/8
      xart(i)=alog10 (1.3*27./13.34) + (19.-tart(i))/2
      print *,pbh(i),xart(i),start(i)
      end do
      s(1)=xart(1)
      t(1)=start(1)
      s(2 )=xart( 7)
      t(2 )=start( 7)
      call pgenv(21.,-5.9,-17.,25.,0,0)!make the axes
      call pgsci(3)
      call pgline(2,s,t)!birth line
      call pgsci(1)
      call labels0(b,c,d,s,t,xt8)!label the axes
      sgns=-1!mass loss on
      do m=1,10 !loop over masses
      call setzeros(x,y,iflag,it8,iz,icmbflag,ian
     1,mmflag,jflag,kflag,mflag,lflag)
      y0=start(m)
      do i=1,1000!loop over PBS temperature
      if(i.eq.1)then!initial conditions
      x(i)=xart(m)
      y(i)=y0!y0 is the initial temp, must stay @Hawking level
      xm(i,m)=tart(m)!these are the original starting masses
      print *,'start',m,'with mass',int(xm(i,m))
      write(11,*)'start',m,'with mass',   (xm(i,m)),x(i),y(i)
      dlmdlt=1.e-35
      else
      x(i)=x(i-1)-0.1
      end if
      if(i.gt.1)then
              call alf(xm,alpha,i,m)
              if(iz.eq.0)dla=-0.1!normal resolution before DE
              expon=(3.*(xm(i-1,m)+20.)+2.*x(i-1))
      if(abs(expon).gt.38.)then
              dlmdlt=0
      else
              dlmdlt=-13.16*alpha/10.**expon
      end if
              dlt=dla/2.!RE
              dm=dlmdlt*dlt
              xm(i,m)=xm(i-1,m)-dm   
              y(i)=y(i-1)+dm
              xm(i,m)=xm(i-1,m)-dm  
      end if
      if(x(i).lt.0..and.iflag.eq.0)then
      iflag=1
      n=i
      end if
      if(i.eq.9)then
      t(1)=alog10(6.62e-27*27.e30/6.67e23/8./1.38e-16/3.14159**2)
     1 +b-5.-xm(i,m)!Hawking temperature change
      end if
      if(i.eq.17)then
      t(2)=alog10(6.62e-27*27.e30/6.67e23/8./1.38e-16/3.14159**2)
     1 +b-5.-xm(i,m)!Hawking temperature change
      end if
      if(x(i).lt.3.-zeq+b.and.iflag.eq.1)then!end of RE
              call alf(xm,alpha,i,m)
              if(iz.eq.0)dla=-0.1
              dlt=4.*dla/3.!MD    take the half and turn it into 2/3
              dm=dlmdlt*dlt
              xm(i,m)=xm(i-1,m)-dm   
              y(i)=y(i-1)+dm
      if(x(i).lt.3.-zeq+b.and.jflag.eq.0)then!end of RE
               jflag=1
      l(m)=i
             print *,'end of RE',i,xm(i-1,m)!              
             write(11,*)'end of RE',i,xm(i-1,m),x(i),y(i)
      end if
      end if
      if(x(i).lt.3.-alog10(10./1.3)-zeq+b.and.kflag.eq.0)then
       kflag=1
c      icmb=i
       end if
       z0=10.**(x(i)+3.633)
       if(z0.le.4..and.iz.eq.0)then
               ti=1.-alog(4.)/sqrt(0.73)
               iz=1
               print *,ti,'time at z=3'
       end if
      if(x(i).lt. 0.144+d   .and.lflag.eq.0)then
              lflag=1
      j=i
      k(m)=i
      print *,'DE',i,xm(i,m)
      write(11,*)'DE',i,xm(i,m),x(i),y(i)           
      end if
      if(x(i).lt. 0.144+d)then
              if(mflag.eq.1) go to 83 
               dltde=dla/10.!increase time resolution 
              dlmdltde=dlmdlt
              xm0=xm(i,m)
              tt=ti
  83          dla=dltde
              dt=-dla/sqrt(0.73)
              tt=tt-dla
              dlt=-dt/(tt-ti)
              dm=dlmdltde*dlt+3.*(xm0-xm(i,m))
              if(dm.lt.0.)then
              print *,dlmdltde,dlt,'mass increase!'
              call pgend
              stop
      end if
      xm(i,m)=xm(i-1,m)-dm   
      y(i)=y(i-1)+dm
      if(iz.eq.1.and.dltde.ne.0.)x(i)=x(i-1)+dltde!incr time resolution in DE
      mflag=1
      end if
      if(x(i).lt.xt8.  and.it8.eq.0)then
      write(11,*)'reached t8',i,start(m),xm(i,m)!,'log mass lost'
      print *, start0,xm(i,m),'mass lost by t8',start0  -xm(i,m)
      if(m.eq.1)write(11, *)'mass lost by t8',start0-xm(i,m)
      if(m.gt.1)write(11, *)'mass lost by t8',xm(1,m)-xm(i,m)
      izflag=1
      it8=1
      write(11,*)pbh(m),'mass loss rate',x(i  ),(xm(i  ,m)-xm(i-2,m))/
     1(x(i  )-x(i-2))
      if(xm(i,m).ne.0.and.ian.eq.0.and.itg.eq.0)write(11,*)i,xm(i,m)
      if(ian.eq.0.and.xm(i,m)-start(m).lt.9.)then
      ian=1
      end if
      end if!uncertain about this
      write(11,100)i,'xm',xm(i,m),'x:',x(i),y(i),jflag,kflag,lflag,!mflag
     1 im(m),iz,int(vgt),dlmdlt,dlt
  100 format(i5,2x,a2,f8.3,2x,a2,f8.3,f8.3,6i5,1x,e10.3,f8.3)
      if(x(i).lt.d.and.mmflag.eq.0)then
      write(11,*)'mass lost',tart(m)-xm(i,m),'redshift 0'
      print *,'mass lost',tart(m)-xm(i,m),'redshift zero new
     1mass',xm(i,m)
      xl=(xm(i,m)-xm(i-1,m))-(x(i)-x(i-1))
     1+2.*.477+20.+xm(i,m)-alog10(13.7*3.14/2.)-16.
      print *,'Luminosity now',xl,'log Lsun',xm(i,m)+33.3+xl,'log M/L'
      mflag=1
      mmflag=1
      end if
      if(xm(i,m).lt.-39.)then
      print *,'m < -39',i,xm(i-1,m)
      go to 9
      end if
      call pgsci(1)
c       if(im(m).gt.0.and.icmbflag.eq.0)then
c       xc(m)=x(i)
c       yc(m)=y(i)
c       icmbflag=1
c       end if
      if(x(i).lt.-6.)go to 9!don't track off into the future
      end do!end of i loop
      if(i.gt.j)x(i)=x1-float(i-j+1)
 9    print *,'exit track loop',m,l(m),n,j,i
      i=min0(i,1000)
      kend(m)=i
      l(m)=min0(i,l(m))
      n=min0(i,n)
c     write(23,*)xl,y(i),pbh(m)
      call pgsch(1.2)
      call pgsci(m+1)
      pos=start(m)!7.-.5*tart(m)
      if(tart(m).le.-20.)then!these masses evaporate before end of RE
             n=i
             if(m.eq.1)call pgtext(xart(m)-3.8,start(m)+2.0,pbh(1))
             if(m.eq.2)call pgtext(3.7,18.,pbh(2))
      call pgline(n,x,y)
              else
              call pgtext(1.5+xart(m),pos-2.0,pbh(m))
      end if
  55  call pgsls(1)
      call pgline(n,x,y)!dashed line for MD era
      call pgsls(2)
      l(m)=n+33
      if(tart(m).gt.-19.) then
      call pgline(l(m)-n,x(n+1),y(n+1))
      else
      if(tart(m).eq.-19.)call pgline(10,x(n+1),y(n+1))
      end if
      call pgslw(3)
      call pgsls(4)!dotted line for DE era
      call pgline(kend(m)-l(m),x(l(m)+1),y(l(m)+1))
      end do!end of m loop
      call pgend
      end
      subroutine alf(xm,alpha,i,m)!from Mosbech & Picker
              dimension xm(1000,10)
              if(xm(i,m).lt.-15.3)then
              alpha=-0.3015+0.3113*10.**((xm(i,m)+33.3)*(-0.0008))
              alpha=alpha/0.0002
              else
              alpha=2.011/2.
              end if
              return
              end
      subroutine labels0(b,c,d,s,t,xt8)
      dimension c(2),s(12),t(12)
      c(1)=b-5.+alog10(2.7)
      d=c(1)
      call pgtext(10.+1.+b,25.9,'10\u15\dK')
      call pgtext(5.+1.+b,25.9,'10\u10\dK')
      call pgtext(b-.5,25.9,'CMB')
      call pgtext(-6.4,10.-.5+b,'10\u15\dK')
      call pgtext(-6.4,5.-5.5+b,'10\u5\dK')
      call pgtext(-6.4,b-10.5,'10\u-5\dK')
      call pgtext(d+1.,25.9,'2.7K')
      s(1)=alog10(1301.)-3.633!CMB
      s(2)=s(1)
      t(1)=-30
      t(2)=11.5
      call pgsls(4)
      call pgline(2,s,t)!CMB line
      s(1)=6.+alog10(220.)
      s(2)=s(1)
      t(1)=30
      t(2)=-10
      call pgline(2,s,t)!220 line QCD transition
c mark the  temperature axes      
      print *,'2.7K',c(1)
      s(1)=c(1)
      s(2)=s(1)
      t(2)=-30
      call pgsls(4)
      call pgline(2,s,t)!now line
      call pglabel('log kT eV universe','log kT eV PBH',' ')
      s(1)=13
      t(1)=s(1)
      s(2)=-5
      t(2)=s(2)
      call pgsci(3)
      call pgsls(1)
       call pgslw(4)
c     call pgline(2,s,t)!birth line
      s(1)=alog10(26.)-3.633
      s(2)=s(1)
      xt8=s(1)
      t(1)=-30
      t(2)=11.5
      call pgsls(2)
      call pgsch(1.2)
      call pgtext(s(1)-.2,10.5,'t8')
      call pgline(2,s,t)!t8
      call pgslw(1)
      call pgsls(4)
      s(1)=alog10(11.)-3.633
      s(2)=s(1)
      t(1)=20
      t(2)=12
c     call pgline(2,s,t)!JWST look back line
      call pgsch(1.)
      s(1)=-.144-3.633
      s(2)=s(1)
      t(1)=-10
      t(2)=20
c     call pgsls(4)
      call pgsci(1)
      call pgsls(1)
      call pgsah(2,45.,0.3)
      do i=1,5
      call pgsch(1.3)
      s(i)=.15+6.+(i-1)*3.-alog10(1.31)/2.
      t(i)=-14
      if(i.eq.1)call pgtext(s(i),t(i)+.5,'1s')
      if(i.eq.2)call pgtext(s(i),t(i)+.5,'1\gms')
      if(i.eq.3)call pgtext(s(i),t(i)+.5,'1ps')
      if(i.eq.4)call pgtext(s(i),t(i)+.5,'1as')
      if(i.eq.5)call pgtext(s(i),t(i)+.5,'1ys')
      call pgsch(1.1)
      if(i.eq.1)call pgtext(1.-3.633,22.5,'1')
      if(i.eq.2)call pgtext(7.8-3.633,22.5,'log 1+z =')
      if(i.eq.2)call pgtext(-.4,22.5,'3')
      if(i.eq.3)call pgtext(2.2-3.633,22.5,'2')
      w=i
      w=w-3.633
      call pgsch(1.1)
      if(i.le.3)call pgarro(w,21.,w,19.)
      end do
      call pgpoint(5,s,t,2)
      s(1)=20
      s(2)=3
      call pgline(2,s,t)!time line
      return
      end
      subroutine setzeros(x,y,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)
      dimension x(1000),y(1000)!,xm(1000,10),iflags(15)
      do i=1,1000!clear plotting arrays
      x(i)=0
      y(i)=0
      end do
      call pgsls(1)
      i1=0
      i2=0
      i3=0
      i4=0
      i5=0
      i6=0
      i7=0
      i8=0
      i9=0
      i10=0
      return
      end
