program vihrToka

    character(100):: name


    !        y   ^
    !            '
    !            '
    !       m    --------------
    !   Ux2->       '
    ! m * y2/y   ---'
    !               |
    ! m * y1/y+1 ---'    
    !   Ux1->       '
    !            --------------  --> x
    !           0              n
    ! Соотношений 1:1(ширина:высота)

    ! Объявление констант
    ! n - кол-во точек по x
    ! m - кол-во точек по y
    !               m --------------
    !
    !    topWallPoint ---'
    !                    |
    ! bottomWallPoint ---'    
    !                    '
    !                 ---'----------
    !                    ^-rightWallPoint
    integer, parameter:: n=12, m=12
    real, parameter:: y=3, x=3, y1=1, y2=2, x1=1, Re=10
    integer, parameter::bottomWallPoint= nint(m * y1/y) + 1 ! j1
    integer, parameter::topWallPoint= nint(m * y2/y) ! j2
    integer, parameter::rightWallPoint = nint(n * x1/x) ! i2
     
    real, parameter::dt=0.01

    ! Объявление массивов
    ! tok - функция тока на n
    ! tokTemp - промежуточный слой 
    ! tokn1 -  функция тока на n+1 слое
    real, dimension(n, m):: tok=0, tokTemp=0, tokn1=0

    ! Объявление массивов
    ! vihr - вихря на n
    ! vihrTemp - вихря промежуточный слой 
    ! vihrn1 -  вихря на n+1 слое
    real, dimension(n, m):: vihr=0, vihrTemp=0, vihrn1=0
    ! real, dimension(n):: a=0, b=0, c=0, d=0, e=0
    real tokConvergence, vihrConvergence, Ux1, Ux2, current_x, current_y
    ! integer lowerBoundary, upperBoundary, leftBoundary, rightBoundary
    tokConvergence = 0; vihrConvergence = 0;
    Ux1=1
    Ux2=1
    dy = y/(m-1)
    dx = x/(n-1)

    ! Начальные/граничные условия условия левая стенка
    call setBoundaryTokValue(n, m, tok, bottomWallPoint, topWallPoint, rightWallPoint, dy, Ux1, Ux2);
    call setBoundaryVihrValue(n, m, vihr, tok, bottomWallPoint, topWallPoint, rightWallPoint, dy, dx);

    tokConvergence = calcConvergence(n, m, tok, tokn1, dt) ! получаем начальную сходимость, чтобы войти в цикл
    vihrConvergence = calcConvergence(n, m, vihr, vihrn1, dt) ! получаем начальную сходимость, чтобы войти в цикл
    ! call PrintArray(n, m, tok) ! это для отладки
    call PrintArray(n, m, vihrn1)

    ! Вывод на консоль
    print*,'Enter name of the file'
    ! Считать значение с консоли
    read(*,*) name

    ! start while
    do while (tokConvergence > 0.001 .AND. vihrConvergence > 0.001)

        ! по X
        ! c j = 2 до j = bottomWallPoint - 1
        !   i = 1 до i = n
        call calcTokX(2, bottomWallPoint-1, 1, n, n, m, tokTemp, tok, vihr, dx, dt)

        ! c j = bottomWallPoint    до j = topWallPoint
        !   i = rightWallPoint     до i = n
        call calcTokX(bottomWallPoint, topWallPoint, rightWallPoint, n, n, m, tokTemp, tok, vihr, dx, dt)

        ! c j = topWallPoint+1 до j = m-1
        !   i = 1              до i = n
        call calcTokX(topWallPoint+1, m-1, 1, n, n, m, tokTemp, tok, vihr, dx, dt)

        ! по Y
        ! c i = 2 до i = rightWallPoint
        ! c j = 1 до j = bottomWallPoint 
        call calcTokY(2, rightWallPoint, 1, bottomWallPoint, n, m, tokn1, tokTemp, tok, dy, dt)

        ! c i = rightWallPoint + 1 до i = n-1
        ! c j = 1                  до j = m
        call calcTokY(rightWallPoint + 1, n-1, 1, m, n, m, tokn1, tokTemp, tok, dy, dt)

        ! c i = 2            до i = rightWallPoint
        ! c j = topWallPoint до j = m
        call calcTokY(2, rightWallPoint, topWallPoint, m, n, m, tokn1, tokTemp, tok, dy, dt)

        ! Начальные/граничные условия условия левая стенка
        call setBoundaryTokValue(n, m, tokn1, bottomWallPoint, topWallPoint, rightWallPoint, dy, Ux1, Ux2);


        ! по X
        ! c j = 2 до j = bottomWallPoint - 1
        !   i = 1 до i = n
        call calcVihrX(2, bottomWallPoint-1, 1, n, n, m, vihrTemp, tokn1, vihr, dx, dy, dt, Re, .false.)

        ! c j = bottomWallPoint    до j = topWallPoint
        !   i = rightWallPoint     до i = n
        call calcVihrX(bottomWallPoint, topWallPoint, rightWallPoint, n, n, m, vihrTemp, tokn1, vihr, dx, dy, dt, Re, .true.)

        ! c j = topWallPoint+1 до j = m-1
        !   i = 1              до i = n
        call calcVihrX(topWallPoint+1, m-1, 1, n, n, m, vihrTemp, tokn1, vihr, dx, dy, dt, Re, .false.)

        ! по Y
        ! c i = 2 до i = rightWallPoint
        ! c j = 1 до j = bottomWallPoint
        call calcVihrY(2, rightWallPoint, 1, bottomWallPoint, n, m, vihrn1, vihrTemp, tokn1, dx, dy, dt, Re)

        ! c i = rightWallPoint + 1 до i = n-1
        ! c j = 1                  до j = m
        call calcVihrY(rightWallPoint + 1, n-1, 1, m, n, m, vihrn1, vihrTemp, tokn1, dx, dy, dt, Re)

        ! c i = 2            до i = rightWallPoint
        ! c j = topWallPoint до j = m
        call calcVihrY(2, rightWallPoint, topWallPoint, m, n, m, vihrn1, vihrTemp, tokn1, dx, dy, dt, Re)

        ! Начальные/граничные условия условия левая стенка
        call setBoundaryVihrValue(n, m, vihrn1, tokn1, bottomWallPoint, topWallPoint, rightWallPoint, dy, dx);

        vihrConvergence = calcConvergence(n, m, vihr, vihrn1, dt)
        tokConvergence = calcConvergence(n, m, tok, tokn1, dt)

        print *, tokConvergence;
        print *, vihrConvergence;
        print *, "-------------";

        vihr = vihrn1;
        tok = tokn1;
        
    enddo ! end while

    ! call PrintArray(n, m, tokn1)
    call PrintArray(n, m, vihrn1)

    ! Запись в файл
    open(23,file=trim(name)//'.dat')
    do i=1, n
        do j=1, m
            current_x = (i-1)*dx
            current_y = (j-1)*dy
            write(23,'(3f9.5)') current_x, current_y, tok(i, j)
        enddo
    enddo

end program vihrToka

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

! это для отладки
subroutine PrintArray(n, m, arr)
    integer n, m
    real,dimension (n, m):: arr

    do j=m, 1, -1
        do i=1, n
            write (*,20, advance='no') arr(i,j)
            20 format(3f9.5)
        enddo
        print *, ''
    enddo

    end

function calcConvergence(n, m, arr, arrn1, dt)
    integer n, m
    real,dimension (n, m):: arr, arrn1
    real calcConvergence, dt
    calcConvergence = 0

    do j=2, m-1
        do i=2, n-1
            calcConvergence = calcConvergence + abs(arrn1(i, j) - arr(i, j)) / dt
        enddo
    enddo

    end function calcConvergence

subroutine setBoundaryTokValue(n, m, tok, bottomWallPoint, topWallPoint, rightWallPoint, dy, Ux1, Ux2)
    integer n, m, bottomWallPoint, topWallPoint, rightWallPoint
    real,dimension (n, m):: tok
    real dy, Ux1, Ux2

    tok(1:n, 1) = 0 ! нижняя стенка
 
    do j = 1, bottomWallPoint
        tok(1, j) = (j - 1) * dy * Ux1 ! нижний вход граничные
    enddo

    tok(1:rightWallPoint, bottomWallPoint) = tok(1, bottomWallPoint); ! нижняя внутренняя стенка
    tok(rightWallPoint, bottomWallPoint:topWallPoint) = tok(1, bottomWallPoint) ! средний проход граничные
    tok(1:rightWallPoint, topWallPoint) = tok(1, bottomWallPoint); ! верхняя внутренняя стенка

    do j = topWallPoint + 1, m
        tok(1, j) = tok(1, bottomWallPoint) + (j - topWallPoint) * dy * Ux2 ! верхний вход граничные
    enddo

    tok(2:n, m) = tok(1, m) ! верхняя стенка
    
    tok(n, 1:m) = tok(n-1, 1:m); ! правый выход

    end

subroutine setBoundaryVihrValue(n, m, vihr, tok, bottomWallPoint, topWallPoint, rightWallPoint, dy, dx)
    integer n, m, bottomWallPoint, topWallPoint, rightWallPoint
    real,dimension (n, m):: vihr, tok
    real dy, dx

    vihr(1, 2:bottomWallPoint - 1) = 0; ! нижний вход граничные

    do j = bottomWallPoint, topWallPoint
        vihr(rightWallPoint, j) = 2 * ( tok(rightWallPoint + 1, j) - tok(rightWallPoint, j) ) / dx**2; ! средний проход граничные
    enddo

    vihr(1, topWallPoint + 1:m - 1) = 0; ! верхний вход граничные

    do i = 1, n
        vihr(i, 1) = 2 * ( tok(i, 2) - tok(i, 1) ) / dy**2; ! нижняя стенка
        vihr(i, m) = 2 * ( tok(i, m - 1) - tok(i, m) ) / dy**2; ! верхняя стенка
    enddo

    do i = 1, rightWallPoint
        vihr(i, topWallPoint)    = 2 * ( tok(i, topWallPoint - 1)    - tok(i, topWallPoint)    ) / dy**2; ! верхняя внутренняя стенка
        vihr(i, bottomWallPoint) = 2 * ( tok(i, bottomWallPoint - 1) - tok(i, bottomWallPoint) ) / dy**2; ! нижняя внутренняя стенка
    enddo

    vihr(n, 1:m) = vihr(n-1, 1:m); ! правый выход

    end

subroutine calcTokX(jStart, jEnd, leftBoundary, rightBoundary, n, m, tokTemp, tok, vihr, dx, dt)
    integer n, m, jStart, jEnd, leftBoundary, rightBoundary
    real, dimension(n):: a, b, c, d, e
    real, dimension(n, m):: tokTemp, tok, vihr

    do j=jStart, jEnd

        do i=leftBoundary+1, rightBoundary-1
            a(i) = 1./dx**2
            c(i) = 1./dx**2
            b(i) = 1./dt + 2./dx**2
            d(i) = tok(i, j)/dt + vihr(i, j)
        enddo
        
        a(leftBoundary)  = 0; b(leftBoundary)   = 1; c(leftBoundary)  = 0; d(leftBoundary)  = tok(leftBoundary, j); 
        a(rightBoundary) = 1; b(rightBoundary)  = 1; c(rightBoundary) = 0; d(rightBoundary) = 0;

        ! Вызов прогонки
        call Tom(leftBoundary, rightBoundary, a, b, c, d, e, 1, max(n,m))
        ! Присваивам в промежуточный слой
        do i=leftBoundary, rightBoundary
            tokTemp(i, j) = e(i)
        enddo
    enddo

    end

subroutine calcTokY(iStart, iEnd, lowerBoundary, upperBoundary, n, m, tokn1, tokTemp, tok, dy, dt)
    integer n, m, iStart, iEnd, lowerBoundary, upperBoundary
    real, dimension(n):: a, b, c, d, e
    real, dimension(n, m):: tokn1, tokTemp, tok

    do i=iStart, iEnd

        do j= lowerBoundary + 1, upperBoundary - 1 
            a(j) = 1./dy**2
            c(j) = 1./dy**2
            b(j) = 1./dt + 2./dy**2
            d(j) = tokTemp(i, j)/dt
        enddo
        
        a(lowerBoundary) = 0; b(lowerBoundary) = 1; c(lowerBoundary) = 0; d(lowerBoundary) = tok(i, lowerBoundary);
        a(upperBoundary) = 0; b(upperBoundary) = 1; c(upperBoundary) = 0; d(upperBoundary) = tok(i, upperBoundary);

        ! Вызов прогонки
        call Tom(lowerBoundary, upperBoundary, a, b, c, d, e, 1, max(n,m))

        ! Присваивам в n+1 слой
        do j= lowerBoundary, upperBoundary 
            tokn1(i, j) = e(j)
        enddo
    enddo

    end

subroutine calcVihrX(jStart, jEnd, leftBoundary, rightBoundary, n, m, vihrTemp, tokn1, vihr, dx, dy, dt, Re, isMiddle)
    integer n, m, jStart, jEnd, leftBoundary, rightBoundary
    real, dimension(n):: a, b, c, d, e
    real, dimension(n, m):: vihrTemp, tokn1, vihr
    logical:: isMiddle
    real Re

    real:: Ux = 0

    do j=jStart, jEnd

        do i=leftBoundary+1, rightBoundary-1
            Ux = ( tokn1(i, j + 1) - tokn1(i, j - 1) ) / dy
            a(i) = -Ux/dx - 1./Re/dx**2
            c(i) = Ux/dx - 1./Re/dx**2
            b(i) = 2./Re/dx**2 + 1./dt
            d(i) = vihr(i, j)/dt
        enddo

        a(leftBoundary) = 0; b(leftBoundary) = 1; c(leftBoundary) = 0;
        a(rightBoundary) = 1; b(rightBoundary)  = 1; c(rightBoundary) = 0;

        if (isMiddle) then
            d(leftBoundary) = 2 * (tokn1(leftBoundary + 1, j) - tokn1(leftBoundary, j)) / dx**2;
        else
            d(leftBoundary) = 0;
        endif

         d(rightBoundary) = 0;

        ! Вызов прогонки
        call Tom(leftBoundary, rightBoundary, a, b, c, d, e, 1, max(n,m))

        ! Присваивам в промежуточный слой
        do i=leftBoundary, rightBoundary
            vihrTemp(i, j) = e(i)
        enddo
    enddo

    end

subroutine calcVihrY(iStart, iEnd, lowerBoundary, upperBoundary, n, m, vihrn1, vihrTemp, tokn1, dx, dy, dt, Re)
    integer n, m, iStart, iEnd, lowerBoundary, upperBoundary
    real, dimension(n, m):: vihrn1, vihrTemp, tokn1
    real, dimension(n):: a, b, c, d, e
    real Re

    real:: Uy = 0

    do i=iStart, iEnd

        do j= lowerBoundary + 1, upperBoundary - 1
            Uy = ( tokn1(i + 1, j) - tokn1(i - 1, j) ) / dx
            a(i) = -Uy/dy - 1./Re/dy**2
            c(i) = Uy/dy - 1./Re/dy**2
            b(i) = 2./Re/dy**2 + 1./dt
            d(j) = vihrTemp(i, j)/dt
        enddo

        a(lowerBoundary) = 0; b(lowerBoundary) = 1; c(lowerBoundary) = 0;
        a(upperBoundary) = 0; b(upperBoundary) = 1; c(upperBoundary) = 0;

        d(lowerBoundary) = 2 * ( tokn1(i, lowerBoundary + 1) - tokn1(i, lowerBoundary) ) / dy**2;
        d(upperBoundary) = 2 * ( tokn1(i, upperBoundary - 1) - tokn1(i, upperBoundary) ) / dy**2;

        ! Вызов прогонки
        call Tom(lowerBoundary, upperBoundary, a, b, c, d, e, 1, max(n,m))

        ! Присваивам в n+1 слой
        do j= lowerBoundary, upperBoundary
            vihrn1(i, j) = e(j)
        enddo
    enddo

    end
