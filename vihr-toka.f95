program vihrToka

    character(100):: name


    !        y  ^
    !           '
    !           '
    !       m   --------------
    !   Ux2->      '
    ! m * 2/3   ---'
    !              |
    ! m * 1/3+1 ---'    
    !   Ux1->      '
    !           --------------  --> x
    !          0              n
    ! Соотношений 4:3 (ширина:высота)

    ! Объявление констант
    ! n - кол-во точе по x
    ! m - кол-во точе по y
    !               m --------------
    !
    !    topWallPoint ---'
    !                    |
    ! bottomWallPoint ---'    
    !                    '
    !                 ---'----------
    !                    ^-rightWallPoint
    integer, parameter:: n=12, m=12
    integer, parameter::bottomWallPoint= nint(m/3.) + 1 ! j1
    integer, parameter::topWallPoint= nint(m * 2./3) ! j2
    integer, parameter::rightWallPoint = nint(n/3.) ! i2
     
    real, parameter::dt=0.01

    ! Объявление массивов
    ! tok - функция тока на n
    ! tokTemp - промежуточный слой 
    ! tokn1 -  функция тока на n+1 слоей
    real, dimension(n, m):: tok=0, tokTemp=0, tokn1=0

    ! Объявление массивов
    ! vihr - вихря на n
    ! vihrTemp - вихря промежуточный слой 
    ! vihrn1 -  вихря на n+1 слоей
    real, dimension(n, m):: vihr=0, vihrTemp=0, vihrn1=0
    real, dimension(n):: a=0, b=0, c=0, d=0, e=0
    real x,y, tokConvergence, Ux1, Ux2, current_x, current_y
    integer lowerBoundary, upperBoundary, leftBoundary, rightBoundary
    tokConvergence = 0
    x=3
    y=3
    Ux1=1
    Ux2=1
    dy = y/(m-1)
    dx = x/(n-1)

    ! Начальные/граничные условия условия левая стенка
    call setBoundaryTokValue(n, m, tok, bottomWallPoint, topWallPoint, rightWallPoint, dy, Ux1, Ux2, y);

    tokConvergence = calcConvergence(n, m, tok, tokn1) ! получаем начальную сходимость, чтобы войти в цикл
    call PrintArray(n, m, tok) ! это для отладки

    ! Вывод на консоль
    print*,'Enter name of the file'
    ! Считать значение с консоли
    read(*,*) name

    ! start while
    do while (tokConvergence > 0.001)

        ! по X
        ! c j = 2 до j = bottomWallPoint - 1
        !   i = 1 до i = n
        do j=2, bottomWallPoint-1
            leftBoundary = 1
            rightBoundary = n

            do i=leftBoundary+1, rightBoundary-1
                a(i) = 1./dx**2
                c(i) = 1./dx**2
                b(i) = 1./dt + 2./dx**2
                d(i) = tok(i, j)/dt + vihr(i, j)
            enddo
            
            a(leftBoundary)  = 0; b(leftBoundary)   = 1; c(leftBoundary)  = 0; d(leftBoundary)  = tok(leftBoundary, j); 
            a(rightBoundary) = 1; b(rightBoundary)  = 1; c(rightBoundary) = 0; d(rightBoundary) = 0;

            ! Вызов прогонки
            call Tom(leftBoundary, rightBoundary, a, b, c, d, e, leftBoundary, rightBoundary)

            ! Присваивам в промежуточный слой
            do i=leftBoundary, rightBoundary
                tokTemp(i, j) = e(i)
            enddo
        enddo

        ! c j = bottomWallPoint    до j = topWallPoint
        !   i = rightWallPoint     до i = n
        do j=bottomWallPoint, topWallPoint
            leftBoundary  = rightWallPoint
            rightBoundary = n

            do i= leftBoundary + 1, rightBoundary - 1
                a(i) = 1./dx**2
                c(i) = 1./dx**2
                b(i) = 1./dt + 2./dx**2
                d(i) = tok(i, j)/dt + vihr(i, j)
            enddo
            
            a(leftBoundary)  = 0; b(leftBoundary)  = 1; c(leftBoundary)  = 0; d(leftBoundary)  = tok(leftBoundary, j); 
            a(rightBoundary) = 1; b(rightBoundary) = 1; c(rightBoundary) = 0; d(rightBoundary) = 0;

            ! Вызов прогонки
            call Tom(leftBoundary, rightBoundary, a, b, c, d, e, leftBoundary, rightBoundary)

            ! Присваивам в промежуточный слой
            do i=leftBoundary, rightBoundary
                tokTemp(i, j) = e(i)
            enddo
        enddo

        ! c j = topWallPoint+1 до j = m-1
        !   i = 1              до i = n
        do j=topWallPoint+1, m-1
            leftBoundary = 1
            rightBoundary = n

            do i= leftBoundary + 1, rightBoundary - 1
                a(i) = 1./dx**2
                c(i) = 1./dx**2
                b(i) = 1./dt + 2./dx**2
                d(i) = tok(i, j)/dt + vihr(i, j)
            enddo

            a(leftBoundary)  = 0; b(leftBoundary)  = 1; c(leftBoundary)  = 0; d(leftBoundary)  = tok(leftBoundary, j); 
            a(rightBoundary) = 1; b(rightBoundary) = 1; c(rightBoundary) = 0; d(rightBoundary) = 0;

            ! Вызов прогонки
            call Tom(leftBoundary, rightBoundary, a, b, c, d, e, leftBoundary, rightBoundary)

            ! Присваивам в промежуточный слой
            do i=leftBoundary, rightBoundary
                tokTemp(i, j) = e(i)
            enddo
        enddo

        ! по Y
        ! c i = 2 до i = rightWallPoint
        ! c j = 1 до j = bottomWallPoint 
        do i=2, rightWallPoint
            lowerBoundary = 1
            upperBoundary = bottomWallPoint;

            do j= lowerBoundary + 1, upperBoundary - 1 
                a(j) = 1./dy**2
                c(j) = 1./dy**2
                b(j) = 1./dt + 2./dy**2
                d(j) = tokTemp(i, j)/dt
            enddo
            
            a(lowerBoundary) = 0; b(lowerBoundary) = 1; c(lowerBoundary) = 0; d(lowerBoundary) = tok(i, lowerBoundary);
            a(upperBoundary) = 0; b(upperBoundary) = 1; c(upperBoundary) = 0; d(upperBoundary) = tok(i, upperBoundary);

            ! Вызов прогонки
            call Tom(lowerBoundary, upperBoundary, a, b, c, d, e, lowerBoundary, upperBoundary)

            ! Присваивам в промежуточный слой
            do j= lowerBoundary, upperBoundary 
                tokn1(i, j) = e(i)
            enddo
        enddo

        ! c i = rightWallPoint + 1 до i = n-1
        ! c j = 1                  до j = m
        do i= rightWallPoint + 1, n-1
            lowerBoundary = 1
            upperBoundary = m

            do j= lowerBoundary + 1, upperBoundary - 1 
                a(j) = 1./dy**2
                c(j) = 1./dy**2
                b(j) = 1./dt + 2./dy**2
                d(j) = tokTemp(i, j)/dt
            enddo
            
            a(lowerBoundary) = 0; b(lowerBoundary) = 1; c(lowerBoundary) = 0; d(lowerBoundary) = tok(i, lowerBoundary);
            a(upperBoundary) = 0; b(upperBoundary) = 1; c(upperBoundary) = 0; d(upperBoundary) = tok(i, upperBoundary);

            ! Вызов прогонки
            call Tom(lowerBoundary, upperBoundary, a, b, c, d, e, lowerBoundary, upperBoundary)

            ! Присваивам в промежуточный слой
            do j= lowerBoundary, upperBoundary 
                tokn1(i, j) = e(i)
            enddo
        enddo

        ! c i = 2            до i = rightWallPoint
        ! c j = topWallPoint до j = m
        do i=2, rightWallPoint
            lowerBoundary = topWallPoint
            upperBoundary = m

            do j= lowerBoundary + 1, upperBoundary - 1 
                a(j) = 1./dy**2
                c(j) = 1./dy**2
                b(j) = 1./dt + 2./dy**2
                d(j) = tokTemp(i, j)/dt
            enddo
            
            a(lowerBoundary) = 0; b(lowerBoundary) = 1; c(lowerBoundary) = 0; d(lowerBoundary) = tok(i, lowerBoundary);
            a(upperBoundary) = 0; b(upperBoundary) = 1; c(upperBoundary) = 0; d(upperBoundary) = tok(i, upperBoundary);

            ! Вызов прогонки
            call Tom(lowerBoundary, upperBoundary, a, b, c, d, e, lowerBoundary, upperBoundary)

            ! Присваивам в промежуточный слой
            do j= lowerBoundary, upperBoundary 
                tokn1(i, j) = e(i)
            enddo
        enddo

        ! Начальные/граничные условия условия левая стенка
        call setBoundaryTokValue(n, m, tokn1, bottomWallPoint, topWallPoint, rightWallPoint, dy, Ux1, Ux2, y);

        tokConvergence = calcConvergence(n, m, tok, tokn1)
        print *, tokConvergence;

        vihr = vihrn1;
        tok = tokn1;
        
    enddo ! end while

    call PrintArray(n, m, tokn1)

    ! Запись в файл
    open(23,file=trim(name)//'.dat')
    do i=1, n
        do j=1, m
            current_x = (i-1)*dx
            current_y = (j-1)*dy
            write(23,'(3f9.5)') current_x, current_y, tok(i, j)
            ! 23 format(3f9.5)
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

function calcConvergence(n, m, arr, arrn1)
    integer n, m
    real,dimension (n, m):: arr, arrn1
    real calcConvergence
    calcConvergence = 0

    do j=2, m-1
        do i=2, n-1
            calcConvergence = calcConvergence + abs(arrn1(i, j) - arr(i, j))
        enddo
    enddo

    end function calcConvergence

subroutine setBoundaryTokValue(n, m, tok, bottomWallPoint, topWallPoint, rightWallPoint, dy, Ux1, Ux2, y)
    integer n, m, bottomWallPoint, topWallPoint, rightWallPoint
    real,dimension (n, m):: tok
    real dy, Ux1, Ux2

    do j=1, m

        if (j > 1 .AND. j < bottomWallPoint) then
             ! Нижний вход Ux1
            tok(1, j) = (j - 1) * dy * Ux1

        else if (j >= bottomWallPoint .AND. j <= topWallPoint) then
             ! Установка значений внутренних стенок и ее внутренней области
            tok(1:rightWallPoint, j) = (j - 1) * dy

        else if (j > topWallPoint .AND. j < m) then 
             ! Верхний вход Ux2
            tok(1, j) = (j - 1) * dy * Ux2

        else
             ! Левая граница верхняя и нижняя стенка (точка)
            tok(1, j) = (j - 1) * dy
        endif

        tok(n, j) = tok(n-1, j) ! правая граница
    enddo
    tok(1:n, m) = y ! верхняя стенка
    tok(1:n, 1) = 0 ! нижняя стенка

    end