  integer,parameter:: mx=201,my=101,ma=max(mx,my),mi=1

  real, dimension(mx,my)::tet,tets,tet1,tet2,vi,vis,vi1,vi2,ft,ft1,ft2,fts,ux,uy
  real, dimension(mi:ma)::a,b,c,d,e
  character(100):: sim ,name

  namelist/Kanal_plosk/ Re,Pr, dt, dl_trub, dx, dy, sumft, sumvi, sumtet

	!param----------------
  dt=0.002
  Gr=0.
  Re=20.
  Pr=1.
  y1=0.4
  y2=0.6
  x1=0.5

  dl_trub=4.

  dy=1./(my-1)
  dy2=dy*dy

  dx=dl_trub/(mx-1)
  dx2=dx*dx
  j1=nint(y1/dy)+1
  j2=nint(y2/dy)+1
  i1=nint(x1/dx)+1
  x1=(i1-1)*dx
  y1=(j1-1)*dy
  y2=(j2-1)*dy
                     
     !-------------------
     
  do j=2,j1-1
    tet(1,j)=1.
  enddo

  do j=j2+1,my-1
    tet(1,j)=0.
  enddo

  tets=tet

   !----------------------  
  do j=2,j1-1
    ux(1,j)=1.
  enddo
  do j=j2+1,my-1
    ux(1,j)=1.
  enddo
  !--------------------  
  do j=1,j1
    y=(j-1)*dy
    ft(1,j)=y
  enddo
  do i=2,i1
    ft(i,j1)=ft(1,j1) !!!������������ ������
  enddo
    
  do j=j1,j2
    ft(i1,j)=ft(1,j1)
  enddo
  do i=1,i1
    ft(i,j2)=ft(1,j1)
  enddo
  
  do j=j2+1,my
    y=(j-1)*dy
    ft(1,j)=ft(1,j1)+y-y2 !
  enddo
  do i=2,mx
    ft(i,my)=ft(1,my)
    ft(i,1)=0.
  enddo
  fts=ft
    
!-----------------------
  
    
  write(*,"(' kolichestvo itk='/,)")
  read(*,*) itk
  write(*,"(' shag vivoda=',)")
  read(*,*) ish
  isht=ish
  it=0
            
  !���� �� ������� !--------------------------   
            
  do while (it<itk)
    it=it+1        

    ! �������� � 1         
    do j=2, j1-1

      do i=2,mx-1
        a(i)=1/dx2
        c(i)=1/dx2
        b(i)=1/dt+2/dx2
        d(i)=fts(i,j)/dt-vi(i,j)
      enddo !i

      a(1)=0.;  b(1)=1.;   c(1)=0.;   d(1)=fts(1,j)
      a(mx)=1.; b(mx)=1.;  c(mx)=0.;  d(mx)=0.

      call Tom(1,mx,a,b,c,d,e,mi,ma)
      do i=1,mx
        ft1(i,j)=e(i)
      enddo !i
    enddo
       
       
    ! �������� � 2         
    do j=j1, j2

      do i=i1+1,mx-1
        a(i)=1/dx2
        c(i)=1/dx2
        b(i)=1/dt+2/dx2
        d(i)=fts(i,j)/dt-vi(i,j)
      enddo !i

      a(i1)=0.;  b(i1)=1.;   c(i1)=0.;   d(i1)=fts(i1,j)
      a(mx)=1.; b(mx)=1.;  c(mx)=0.;  d(mx)=0.

      call Tom(i1,mx,a,b,c,d,e,mi,ma)
      do i=i1,mx
        ft1(i,j)=e(i)
      enddo !i
    enddo
    ! �������� � 3         
    do j=j2+1, my-1
      do i=2,mx-1
        a(i)=1/dx2
        c(i)=1/dx2
        b(i)=1/dt+2/dx2
        d(i)=fts(i,j)/dt-vi(i,j)
      enddo !i

      a(1)=0.;  b(1)=1.;   c(1)=0.;   d(1)=fts(1,j)
      a(mx)=1.; b(mx)=1.;  c(mx)=0.;  d(mx)=0.

      call Tom(1,mx,a,b,c,d,e,mi,ma)
      do i=1,mx
      ft1(i,j)=e(i)
      enddo !i
    enddo   
    ! �������� � 1 
    do i=2, i1

      do j=2,j1-1                                                                     
        a(j)=1/dy2
        c(j)=1/dy2
        b(j)=1/dt+2/dy2
        d(j)=ft1(i,j)/dt
      enddo !i

      a(1)=0.;  b(1)=1.;   c(1)=0.;   d(1)=0.
      a(j1)=0.; b(j1)=1.;  c(j1)=0.;  d(j1)=fts(i,j1)

      call Tom(1,j1,a,b,c,d,e,mi,ma)
      do j=1,j1
        ft(i,j)=e(j)
      enddo !i
    enddo 
    ! �������� � 2
    do i=2, i1
      do j=j2+1,my-1                                                                    
        a(j)=1/dy2
        c(j)=1/dy2
        b(j)=1/dt+2/dy2
        d(j)=ft1(i,j)/dt
      enddo !i

      a(j2)=0.;  b(j2)=1.;   c(j2)=0.;   d(j2)=fts(i,j2)
      a(my)=0.; b(my)=1.;  c(my)=0.;  d(my)=fts(i,my)

      call Tom(j2,my,a,b,c,d,e,mi,ma)
      do j=j2,my
        ft(i,j)=e(j)
      enddo !i
    enddo 
    ! �������� � 3
    do i=i1+1, mx-1
      do j=2,my-1
        a(j)=1/dy2
        c(j)=1/dy2
        b(j)=1/dt+2/dy2
        d(j)=ft1(i,j)/dt
      enddo !i

      a(1)=0.;  b(1)=1.;   c(1)=0.;   d(1)=fts(i,1)
      a(my)=0.; b(my)=1.;  c(my)=0.;  d(my)=fts(i,my)

      call Tom(1,my,a,b,c,d,e,mi,ma)
      do j=1,my
        ft(i,j)=e(j)
      enddo !i
    enddo 
 
    ! 
    !!������� �� ���.     
    do j=2, my-1
      ft(mx,j)=ft(mx-1,j)
    enddo
    ft2=ft-fts   
  
    !-------------------------------------------  
    do j=2, my-1
      do i=2,mx-1
        if(.not.(i<=i1.and.j>=j1.and.j<=j2)) then
          ux(i,j)=(ft(i,j+1)-ft(i,j-1))/2/dy     
          uy(i,j)=-(ft(i+1,j)-ft(i-1,j))/2/dx
        endif
      enddo
    enddo

    !+++++++++++++++++++++++++++++++++++++++++++++++++++
    do j=2,my-1
        ux(mx,j)=ux(mx-1,j)
        uy(mx,j)=uy(mx-1,j)
    enddo
   
 !-----------------------------------   
    do i=1, i1
      ux(i,j1)=0.
      uy(i,j1)=0.
      ux(i,j2)=0.
      uy(i,j2)=0.
    enddo
    do j=j1, j2
      ux(i1,j)=0.
      uy(i1,j)=0.
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!
    do j=j2+1, my-1
      uy(1,j)=uy(2,j)
    enddo
    do j=2, j1-1
      uy(1,j)=uy(2,j)
    enddo
       
    !  VI
    !!! x1
    do  j=2, j1-1
      do i=2,mx-1
        aux=abs(ux(i,j))
        a(i)=(aux+ux(i,j))/2/dx+1/Re/dx2
        c(i)=(aux-ux(i,j))/2/dx+1./Re/dx2
        b(i)=1/dt+2/Re/dx2+aux/dx            
        !  d(i)=vis(i,j)/dt+Gr/Re/Re*(tet(i,j+1)-tet(i,j-1))/2/dy
        d(i)=vis(i,j)/dt-Gr/Re/Re*(tet(i+1,j)-tet(i-1,j))/2/dx
      enddo !i

      a(1)=0.;  b(1)=1.;   c(1)=0.;   d(1)=0.
      a(mx)=1.; b(mx)=1.;  c(mx)=0.;  d(mx)=0.

      call Tom(1,mx,a,b,c,d,e,mi,ma)

      do i=1,mx
        vi1(i,j)=e(i)
      enddo !i
    enddo
    !---------x2
    do j=j1, j2
      do i=i1+1,mx-1
        aux=abs(ux(i,j))
        a(i)=(aux+ux(i,j))/2/dx+1/Re/dx2
        c(i)=(aux-ux(i,j))/2/dx+1./Re/dx2
        b(i)=1/dt+2/Re/dx2+aux/dx            
        !  d(i)=vis(i,j)/dt+Gr/Re/Re*(tet(i,j+1)-tet(i,j-1))/2/dy
        d(i)=vis(i,j)/dt-Gr/Re/Re*(tet(i+1,j)-tet(i-1,j))/2/dx
      enddo !i

      a(i1)=0.;  b(i1)=1.;   c(i1)=0.;   d(i1)=2*(ft(i1+1,j)-ft(i1,j))/dx2                        
      a(mx)=1.; b(mx)=1.;  c(mx)=0.;  d(mx)=0.

      call Tom(i1,mx,a,b,c,d,e,mi,ma)                                                              

      do i=i1,mx
        vi1(i,j)=e(i)
      enddo !i
    enddo
    !!!!!!!!!x3
    do j=j2+1, my-1
      do i=2,mx-1
        aux=abs(ux(i,j))
        a(i)=(aux+ux(i,j))/2/dx+1/Re/dx2
        c(i)=(aux-ux(i,j))/2/dx+1./Re/dx2
        b(i)=1/dt+2/Re/dx2+aux/dx            
        !  d(i)=vis(i,j)/dt+Gr/Re/Re*(tet(i,j+1)-tet(i,j-1))/2/dy
        d(i)=vis(i,j)/dt-Gr/Re/Re*(tet(i+1,j)-tet(i-1,j))/2/dx
      enddo !i

      a(1)=0.;  b(1)=1.;   c(1)=0.;   d(1)=0.
      a(mx)=1.; b(mx)=1.;  c(mx)=0.;  d(mx)=0.

      call Tom(1,mx,a,b,c,d,e,mi,ma)

      do i=1,mx
        vi1(i,j)=e(i)
      enddo !i
    enddo
    !!!!!!!!!!!y1
    do i=2, i1
      do j=2,j1-1 
        auy=abs(uy(i,j))
        a(j)=(auy+uy(i,j))/2/dy+1/Re/dy2
        c(j)=(auy-uy(i,j))/2/dy+1/dy2/Re
        b(j)=2./dy2/re+1./dt+auy/dy
        d(j)=vi1(i,j)/dt
      enddo !i

      a(1)=0;  b(1)=1.;   c(1)=0;   d(1)=2*(ft(i,2)-ft(i,1))/dy2                         
      a(j1)=0.; b(j1)=1.;  c(j1)=0.;  d(j1)=2*(ft(i,j1-1)-ft(i,j1))/dy2                 

      call Tom(1,j1,a,b,c,d,e,mi,ma)

      do j=1,j1
        vi(i,j)=e(j)
      enddo !i
      !!!!!!!!!!y2
    enddo
    do i=2, i1
      do j=j2+1,my-1
        auy=abs(uy(i,j))
        a(j)=(auy+uy(i,j))/2/dy+1/Re/dy2
        c(j)=(auy-uy(i,j))/2/dy+1/dy2/Re
        b(j)=2./dy2/re+1./dt+auy/dy
        d(j)=vi1(i,j)/dt
      enddo !i

      a(j2)=0;  b(j2)=1.;   c(j2)=0;   d(j2)=2*(ft(i,j2+1)-ft(i,j2))/dy2                               
      a(my)=0.; b(my)=1.;  c(my)=0.;  d(my)=2*(ft(i,my-1)-ft(i,my))/dy2

      call Tom(j2,my,a,b,c,d,e,mi,ma)                                                             

      do j=j2,my
        vi(i,j)=e(j)
      enddo !i
    enddo
    !!!!!!!!y3
    do i=i1+1, mx-1
      do j=2,my-1
        auy=abs(uy(i,j))
        a(j)=(auy+uy(i,j))/2/dy+1/Re/dy2
        c(j)=(auy-uy(i,j))/2/dy+1/dy2/Re
        b(j)=2./dy2/re+1./dt+auy/dy
        d(j)=vi1(i,j)/dt
      enddo !i

      a(1)=0;  b(1)=1.;   c(1)=0;   d(1)=2*(ft(i,2)-ft(i,1))/dy2
      a(my)=0.; b(my)=1.;  c(my)=0.;  d(my)=2*(ft(i,my-1)-ft(i,my))/dy2

      call Tom(1,my,a,b,c,d,e,mi,ma)

      do j=1,my
        vi(i,j)=e(j)
      enddo !i
    enddo

    !��� �� ���.
    do j=2,my-1
      vi(mx,j)=vi(mx-1,j)
    enddo
    do j=j1, j2
      vi(i1,j)=2*(ft(i1+1,j)-ft(i1,j))/dx2
    enddo
    vi2=vi-vis  

    !------------�����������
               
    do ll=1,3

      !!!!!x1
      do  j=2, j1-1
        do i=2,mx-1
          aux=abs(ux(i,j))
          a(i)=(aux+ux(i,j))/2/dx+1/Re/Pr/dx2
          c(i)=(aux-ux(i,j))/2/dx+1./Re/Pr/dx2
          b(i)=1/dt+2/Re/Pr/dx2+aux/dx
          d(i)=tets(i,j)/dt
        enddo !i

        a(1)=0.;  b(1)=1.;   c(1)=0.;   d(1)=1.                
        a(mx)=1.; b(mx)=1.;  c(mx)=0.;  d(mx)=0. 

        call Tom(1,mx,a,b,c,d,e,mi,ma)

        do i=1,mx
          tet1(i,j)=e(i)
        enddo !i
      enddo
      !------------x2
      do j=j1, j2
        do i=i1+1,mx-1
          aux=abs(ux(i,j))
          a(i)=(aux+ux(i,j))/2/dx+1/Re/Pr/dx2
          c(i)=(aux-ux(i,j))/2/dx+1./Re/Pr/dx2
          b(i)=1/dt+2/Re/Pr/dx2+aux/dx
          d(i)=tets(i,j)/dt
        enddo !i

        a(i1)=0.;  b(i1)=1.;   c(i1)=1.;   d(i1)=0.                 
        a(mx)=1.; b(mx)=1.;  c(mx)=0.;  d(mx)=0. 

        call Tom(i1,mx,a,b,c,d,e,mi,ma)

        do i=i1,mx
          tet1(i,j)=e(i)
        enddo !i
      enddo
      !!!x3
      do j=j2+1, my-1
        do i=2,mx-1
          aux=abs(ux(i,j))
          a(i)=(aux+ux(i,j))/2/dx+1/Re/Pr/dx2
          c(i)=(aux-ux(i,j))/2/dx+1./Re/Pr/dx2
          b(i)=1/dt+2/Re/Pr/dx2+aux/dx
          d(i)=tets(i,j)/dt
        enddo !i

        a(1)=0.;  b(1)=1.;   c(1)=0.;   d(1)=0. 
        a(mx)=1.; b(mx)=1.;  c(mx)=0.;  d(mx)=0. 

        call Tom(1,mx,a,b,c,d,e,mi,ma)

        do i=1,mx
          tet1(i,j)=e(i)
        enddo !i
      enddo
      !!!y1
      do i=2, i1
        do j=2,j1-1 
          auy=abs(uy(i,j))
          a(j)=(auy+uy(i,j))/2/dy+1/Re/Pr/dy2
          c(j)=(auy-uy(i,j))/2/dy+1/dy2/Re/Pr
          b(j)=2./dy2/re/Pr+1./dt+auy/dy
          d(j)=tet1(i,j)/dt
        enddo !i

        a(1)=0.;  b(1)=1.;   c(1)=1.;   d(1)=0.                      
        a(j1)=1.; b(j1)=1.;  c(j1)=0.;  d(j1)=0. 

        call Tom(1,j1,a,b,c,d,e,mi,ma)

        do j=1,j1
          tet(i,j)=e(j)
        enddo !i
      enddo
      !!y2
      do i=2, i1
        do j=j2+1,my-1
          auy=abs(uy(i,j))
          a(j)=(auy+uy(i,j))/2/dy+1/Re/Pr/dy2
          c(j)=(auy-uy(i,j))/2/dy+1/dy2/Re/Pr
          b(j)=2./dy2/re/Pr+1./dt+auy/dy
          d(j)=tet1(i,j)/dt
        enddo !i

        a(j2)=0.;  b(j2)=1.;   c(j2)=1.;   d(j2)=0. 
        a(my)=1.; b(my)=1.;  c(my)=0.;  d(my)=0. 

        call Tom(j2,my,a,b,c,d,e,mi,ma)

        do j=j2,my
          tet(i,j)=e(j)   
        enddo !i
      enddo


      !-----------y3  
      do i=i1+1, mx-1
        do j=2,my-1
          auy=abs(uy(i,j))
          a(j)=(auy+uy(i,j))/2/dy+1/Re/Pr/dy2
          c(j)=(auy-uy(i,j))/2/dy+1/dy2/Re/Pr
          b(j)=2./dy2/re/Pr+1./dt+auy/dy
          d(j)=tet1(i,j)/dt
        enddo !i

        a(1)=0.;  b(1)=1.;   c(1)=1.;   d(1)=0.                                        
        a(my)=1.; b(my)=1.;  c(my)=0.;  d(my)=0.                       

        call Tom(1,my,a,b,c,d,e,mi,ma)

        do j=1,my
          tet(i,j)=e(j)
        enddo !i
      enddo
      !!
      do j=2,my-1
        tet(mx,j)=tet(mx-1,j)
      enddo

      do i=2,mx-1
        tet(i,1)=tet(i,2)
      enddo                                                
      do i=2,mx-1
        tet(i,my)=tet(i,my-1)
      enddo
      do i=2,i1
        tet(i,j1)=tet(i,j1-1)
        tet(i,j2)=tet(i,j2+1)
      enddo
      do j=j1,j2
        tet(i1,j)=tet(i1+1,j)
      enddo

      tet2=tet-tets
      tets=tet
    enddo


    vis=vi
    fts=ft
    tets=tet
 
! ������ ������
 
    if(mod(it,isht)==0) then
      isht=isht+ish

      sumft=0.
      sumvi=0.
      sumtet=0.

      do i=2, mx-1
        do j=2, my-1
          sumft= sumft+abs(ft2(i,j))/dt
          sumvi= sumvi+abs(vi2(i,j))/dt
          sumtet= sumtet+abs(tet2(i,j))/dt
        enddo
      enddo

      write(*,'(I7,3es14.4)') it, sumft, sumvi, sumtet
    endif 
 
 
    if(it<itk) cycle

    print*
    write(*,"('vvod dopol iter=',\)")
    read(*,*) itdop
    itk=itk+itdop

    if(itdop>0)then
      isht=isht-ish
      write(*,"('vvod shaga cxod=',\)")
      read(*,*) ish
      isht=isht+ish
    endif

    if(itdop==0) exit
  enddo !while
    
    
  print*,'Name varianta'

  read(*,*)name


  ! sim=Adjustl(sim)
  open(23,file=trim(name)//' izo_tem_ft_ '//' .dat')
  do i=1,mx
    do j=1,my
      x=(i-1)*dx
      y=(j-1)*dy
      write(23,'(7es14.4)') x,y, tet (i,j), ft (i,j)
    enddo
  enddo
  !!!!
  open(52,file=trim(name)//' izo_ux_uy_vi '//' .dat')

  do i=1,mx
    do j=1,my-1
      x=(i-1)*dx
      y=(j-1)*dy
      write(52,'(7es14.4)') x,y, ux(i,j), uy(i,j) ,vi(i,j)
    enddo
  enddo

  !!!
  open(24,file=trim(name)//' ux nizhnee sechenie'//' .dat')
  do j=1,j1
    y=(j-1)*dy
    write(24,'(8es14.4)') y, ux(i1/2,j), ux (i1-1,j)
  enddo

  open(204,file=trim(name)//' ux verhnee sechenie'//' .dat')
  do j=j2,my
    y=(j-1)*dy
    write(204,'(8es14.4)') y, ux(i1/2,j), ux (i1-1,j)
  enddo

  open(214,file=trim(name)//' ux obshaya chast' //' .dat')
  do j=1,my
    y=(j-1)*dy
    write(214,'(8es14.4)') y, ux((i1+mx)/2,j),ux(2*(i1+mx)/3,j), ux(mx,j) 
  enddo

  open(26,file=trim(name)//' tet nizhnee sechenie'//' .dat')
  do j=1,j1
    y=(j-1)*dy
    write(26,'(8es14.4)') y, tet(i1/2,j), tet (i1-1,j)
  enddo

  open(264,file=trim(name)//' tet verhnee sechenie'//' .dat')
  do j=j2,my
    y=(j-1)*dy
    write(264,'(8es14.4)') y,tet(i1/2,j), tet (i1-1,j)
  enddo

  open(216,file=trim(name)//' tet obshaya chast' //' .dat')
  do j=1,my
    y=(j-1)*dy
    write(216,'(8es14.4)') y, tet((i1+mx)/2,j),tet(2*(i1+mx)/3,j),tet(mx,j) 
  enddo

  open(27,file=trim(name)//' uy nizhnee sechenie'//' .dat')
  do j=1,j1
    y=(j-1)*dy
    write(27,'(8es14.4)') y, uy(i1/2,j), uy (i1-1,j)
  enddo

  open(274,file=trim(name)//' uy verhnee sechenie'//' .dat')
  do j=j2,my
    y=(j-1)*dy
    write(274,'(8es14.4)') y, uy(i1/2,j), uy (i1-1,j)
  enddo

  open(217,file=trim(name)//' uy obshaya chast' //' .dat')
  do j=1,my
    y=(j-1)*dy
    write(217,'(8es14.4)') y, uy((i1+mx)/2,j),uy(2*(i1+mx)/3,j), uy(mx,j) 
  enddo 

  open (77,file=trim(name)//'_Namelist.dat')
  write(77,Kanal_plosk) 

end Program

! Прогонка
! A -B C D
! mi - minimum
! ma - maximum
SUBROUTINE Tom(i0, iN, a, b, c, d, e, mi, ma)
  integer i0,iN,mi,ma,i
  real,dimension (mi:ma):: a,b,c,d,e,alf,bet

  alf(i0) = c(i0)/b(i0);
  bet(i0) = d(i0)/b(i0)
  do i=i0+1, iN-1
    alf(i)= c(i) / ( b(i)-a(i)*alf(i-1) )
  enddo

  do i=i0+1, iN
    bet(i)= ( d(i)+a(i)*bet(i-1) ) / ( b(i)-a(i)*alf(i-1) )
  enddo
  
  e(iN) = bet(iN)
  do i=iN-1,i0,-1
    e(i)= alf(i) * e(i+1) + bet(i)
  enddo

  end
   

