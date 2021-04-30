program vihrToka
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
    integer, parameter:: n=69, m=36
    real, parameter:: x1=1, x=3
    real, parameter:: y1=.4, y2=.6, y=1
    integer, parameter:: bottomWallPoint = nint(m * y1/y) + 1 ! j1
    integer, parameter:: topWallPoint    = nint(m * y2/y) ! j2
    integer, parameter:: rightWallPoint  = nint(n * x1/x) ! i2
     
    real, parameter:: dt=0.001
    real, parameter:: Gr=0.
    real, parameter:: Re=20.
    real, parameter:: Pr=.1

    ! Объявление Тока
    ! tok - функция тока на n
    ! tokTemp - промежуточный слой 
    ! tokn1 - функция тока на n+1 слое
    real, dimension(n, m):: tok=0, tokTemp=0, tokn1=0

    ! Объявление Вихря
    ! vihr - вихря на n
    ! vihrTemp - вихря промежуточный слой 
    ! vihrn1 - вихря на n+1 слое
    real, dimension(n, m):: vihr=0, vihrTemp=0, vihrn1=0

    ! Объявление Температуры
    ! theta - на n
    ! thetaTemp - промежуточный слой 
    ! thetaN1 - n+1 слое
    ! thetaConvergence - разница между thetaN1 и theta
    real, dimension(n, m):: theta=0, thetaTemp=0, thetaN1=0, thetaConvergence=0

    ! Объявление Поля скорости
    real, dimension(n, m):: Ux = 0, Uy = 0

    real, dimension(n):: a=0, b=0, c=0, d=0, e=0
    real Ux1, Ux2, UxAbs, UyAbs
    integer lowerBoundary, upperBoundary, leftBoundary, rightBoundary
    integer additionalIterationNum, curIterationNum, printStepNumber, curPrintIterationNumber

    character(100):: filePrefix
    write (filePrefix,"('Re', i0, 'Pr', i0, 'Gr', i0)") int(Re), int(Pr), int(Gr)
    print*,"file name is ", trim(filePrefix)

    additionalIterationNum = 0
    Ux1=1
    Ux2=1
    dy = y/(m-1)
    dx = x/(n-1)

    ! Начальные/граничные условия
    call setBoundaryTokValue(n, m, tok, bottomWallPoint, topWallPoint, rightWallPoint, dy, Ux1, Ux2);

    ! Считывание итераций с консоли
    write(*,"('Number of iteration:')")
    read(*,*) iterationNum
    write(*,"('Step print:')")
    read(*,*) printStepNumber
    curPrintIterationNumber = printStepNumber
    curIterationNum = 0

    ! start while
    do while (curIterationNum < iterationNum)
        curIterationNum = curIterationNum + 1

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
        ! Начальные/граничные условия условия
        call setBoundaryTokValue(n, m, tokn1, bottomWallPoint, topWallPoint, rightWallPoint, dy, Ux1, Ux2);
        
        do j=2, m-1
            do i=2, n-1
                ! Кроме внутреннего блока из стенок
                ! Расчитываем скорость
                if(.not.(i <= rightWallPoint .and. j >= bottomWallPoint .and. j <= topWallPoint)) then
                    Ux(i,j)=(tokN1(i,j+1)-tokN1(i,j-1))/2/dy     
                    Uy(i,j)=-(tokN1(i+1,j)-tokN1(i-1,j))/2/dx
                endif
            enddo
        enddo

        ux(n, 2:m-1)=ux(n-1, 2:m-1) ! правый выход
        uy(n, 2:m-1)=uy(n-1, 2:m-1) ! правый выход
        ux(1:rightBoundary, bottomWallPoint) = 0 ! Ниэняя внутренняя стенка 
        uy(1:rightBoundary, bottomWallPoint) = 0 ! Ниэняя внутренняя стенка
        ux(1:rightBoundary, topWallPoint)    = 0 ! Верхняя внутренняя стенка
        uy(1:rightBoundary, topWallPoint)    = 0 ! Верхняя внутренняя стенка
        ux(rightBoundary, bottomWallPoint:topWallPoint) = 0 ! Внутренняя вертикальная стенка
        uy(rightBoundary, bottomWallPoint:topWallPoint) = 0 ! Внутренняя вертикальная стенка
        uy(1, 2:bottomWallPoint - 1) = uy(2, 2:bottomWallPoint - 1) ! Вход нижний
        uy(1, topWallPoint+1:m-1)    = uy(2, topWallPoint+1:m-1) ! Вход верхний

        ! по X
        ! c j = 2 до j = bottomWallPoint - 1
        !   i = 1 до i = n
        do j=2, bottomWallPoint-1
            leftBoundary = 1
            rightBoundary = n

            do i=leftBoundary+1, rightBoundary-1
                UxAbs = abs(Ux(i, j))
                a(i) = (UxAbs + Ux(i, j))/dx/2 + 1./Re/dx**2
                c(i) = (UxAbs - Ux(i, j))/dx/2 + 1./Re/dx**2
                b(i) = UxAbs/dx + 2./Re/dx**2 + 1./dt
                d(i) = vihr(i, j)/dt - Gr/Re/Re*(thetaN1(i+1,j)-thetaN1(i-1,j))/2/dx
            enddo

            a(leftBoundary) = 0;  b(leftBoundary)   = 1; c(leftBoundary)  = 0;  d(leftBoundary)  = 0;
            a(rightBoundary) = 1; b(rightBoundary)  = 1; c(rightBoundary) = 0;  d(rightBoundary) = 0;

            ! Вызов прогонки
            call Tom(leftBoundary, rightBoundary, a, b, c, d, e, 1, max(n,m))

            ! Присваивам в промежуточный слой
            do i=leftBoundary, rightBoundary
                vihrTemp(i, j) = e(i)
            enddo
        enddo
        ! c j = bottomWallPoint    до j = topWallPoint
        !   i = rightWallPoint     до i = n
        do j=bottomWallPoint, topWallPoint
            leftBoundary = rightWallPoint
            rightBoundary = n

            do i=leftBoundary+1, rightBoundary-1
                UxAbs = abs(Ux(i, j))
                a(i) = (UxAbs + Ux(i, j))/dx/2 + 1./Re/dx**2
                c(i) = (UxAbs - Ux(i, j))/dx/2 + 1./Re/dx**2
                b(i) = UxAbs/dx + 2./Re/dx**2 + 1./dt
                d(i) = vihr(i, j)/dt - Gr/Re/Re*(thetaN1(i+1,j)-thetaN1(i-1,j))/2/dx
            enddo

            a(leftBoundary) = 0;  b(leftBoundary)   = 1; c(leftBoundary)  = 0;  
            a(rightBoundary) = 1; b(rightBoundary)  = 1; c(rightBoundary) = 0;

            d(leftBoundary) = 2 * (tok(leftBoundary + 1, j) - tok(leftBoundary, j)) / dx**2;
            d(rightBoundary) = 0;
            ! Вызов прогонки
            call Tom(leftBoundary, rightBoundary, a, b, c, d, e, 1, max(n,m))

            ! Присваивам в промежуточный слой
            do i=leftBoundary, rightBoundary
                vihrTemp(i, j) = e(i)
            enddo
        enddo
        ! c j = topWallPoint+1 до j = m-1
        !   i = 1              до i = n
        do j=topWallPoint+1, m-1
            leftBoundary = 1
            rightBoundary = n

            do i=leftBoundary+1, rightBoundary-1
                UxAbs = abs(Ux(i, j))
                a(i) = (UxAbs + Ux(i, j))/dx/2 + 1./Re/dx**2
                c(i) = (UxAbs - Ux(i, j))/dx/2 + 1./Re/dx**2
                b(i) = UxAbs/dx + 2./Re/dx**2 + 1./dt
                d(i) = vihr(i, j)/dt - Gr/Re/Re*(thetaN1(i+1,j)-thetaN1(i-1,j))/2/dx
            enddo

            a(leftBoundary) = 0;  b(leftBoundary)   = 1; c(leftBoundary)  = 0;  d(leftBoundary)  = 0;
            a(rightBoundary) = 1; b(rightBoundary)  = 1; c(rightBoundary) = 0;  d(rightBoundary) = 0;

            ! Вызов прогонки
            call Tom(leftBoundary, rightBoundary, a, b, c, d, e, 1, max(n,m))

            ! Присваивам в промежуточный слой
            do i=leftBoundary, rightBoundary
                vihrTemp(i, j) = e(i)
            enddo
        enddo
        ! по Y
        ! c i = 2 до i = rightWallPoint
        ! c j = 1 до j = bottomWallPoint 
        do i=2, rightWallPoint
            lowerBoundary = 1
            upperBoundary = bottomWallPoint;

            do j= lowerBoundary + 1, upperBoundary - 1
                UyAbs = abs(Uy(i, j))
                a(j) = (UyAbs + Uy(i, j))/dy/2 + 1./Re/dy**2
                c(j) = (UyAbs - Uy(i, j))/dy/2 + 1./Re/dy**2
                b(j) = UyAbs/dy + 2./Re/dy**2 + 1./dt
                d(j) = vihrTemp(i, j)/dt
            enddo

            a(lowerBoundary) = 0; b(lowerBoundary) = 1; c(lowerBoundary) = 0;
            a(upperBoundary) = 0; b(upperBoundary) = 1; c(upperBoundary) = 0;

            d(lowerBoundary) = 2 * ( tok(i, lowerBoundary + 1) - tok(i, lowerBoundary) ) / dy**2;
            d(upperBoundary) = 2 * ( tok(i, upperBoundary - 1) - tok(i, upperBoundary) ) / dy**2;

            ! Вызов прогонки
            call Tom(lowerBoundary, upperBoundary, a, b, c, d, e, 1, max(n,m))

            ! Присваивам в промежуточный слой
            do j= lowerBoundary, upperBoundary 
                vihrn1(i, j) = e(j)
            enddo
        enddo
        ! c i = rightWallPoint + 1 до i = n-1
        ! c j = 1                  до j = m
        do i=rightWallPoint + 1, n-1
            lowerBoundary = 1
            upperBoundary = m;

            do j= lowerBoundary + 1, upperBoundary - 1 
                UyAbs = abs(Uy(i, j))
                a(j) = (UyAbs + Uy(i, j))/dy/2 + 1./Re/dy**2
                c(j) = (UyAbs - Uy(i, j))/dy/2 + 1./Re/dy**2
                b(j) = UyAbs/dy + 2./Re/dy**2 + 1./dt
                d(j) = vihrTemp(i, j)/dt
            enddo

            a(lowerBoundary) = 0; b(lowerBoundary) = 1; c(lowerBoundary) = 0;
            a(upperBoundary) = 0; b(upperBoundary) = 1; c(upperBoundary) = 0;

            d(lowerBoundary) = 2 * ( tok(i, lowerBoundary + 1) - tok(i, lowerBoundary) ) / dy**2;
            d(upperBoundary) = 2 * ( tok(i, upperBoundary - 1) - tok(i, upperBoundary) ) / dy**2;

            ! Вызов прогонки
            call Tom(lowerBoundary, upperBoundary, a, b, c, d, e, 1, max(n,m))

            ! Присваивам в промежуточный слой
            do j= lowerBoundary, upperBoundary 
                vihrn1(i, j) = e(j)
            enddo
        enddo
        ! c i = 2            до i = rightWallPoint
        ! c j = topWallPoint до j = m
        do i=2, rightWallPoint
            lowerBoundary = topWallPoint
            upperBoundary = m;

            do j= lowerBoundary + 1, upperBoundary - 1 
                UyAbs = abs(Uy(i, j))
                a(j) = (UyAbs + Uy(i, j))/dy/2 + 1./Re/dy**2
                c(j) = (UyAbs - Uy(i, j))/dy/2 + 1./Re/dy**2
                b(j) = UyAbs/dy + 2./Re/dy**2 + 1./dt
                d(j) = vihrTemp(i, j)/dt
            enddo

            a(lowerBoundary) = 0; b(lowerBoundary) = 1; c(lowerBoundary) = 0;
            a(upperBoundary) = 0; b(upperBoundary) = 1; c(upperBoundary) = 0;

            d(lowerBoundary) = 2 * ( tok(i, lowerBoundary + 1) - tok(i, lowerBoundary) ) / dy**2;
            d(upperBoundary) = 2 * ( tok(i, upperBoundary - 1) - tok(i, upperBoundary) ) / dy**2;

            ! Вызов прогонки
            call Tom(lowerBoundary, upperBoundary, a, b, c, d, e, 1, max(n,m))

            ! Присваивам в промежуточный слой
            do j= lowerBoundary, upperBoundary 
                vihrn1(i, j) = e(j)
            enddo
        enddo

        ! Начальные/граничные условия условия левая стенка
        call setBoundaryVihrValue(n, m, vihrn1, tok, bottomWallPoint, topWallPoint, rightWallPoint, dy, dx);

        ! Температура -------------
        do k = 1,3
            ! по X
            ! c j = 2 до j = bottomWallPoint - 1
            !   i = 1 до i = n
            do j=2, bottomWallPoint-1
                leftBoundary = 1
                rightBoundary = n

                do i=leftBoundary+1, rightBoundary-1
                    UxAbs = abs(Ux(i, j))
                    a(i) = (UxAbs + Ux(i, j))/dx/2 + 1/Re/Pr/dx**2
                    c(i) = (UxAbs - Ux(i, j))/dx/2 + 1./Re/Pr/dx**2
                    b(i) = UxAbs/dx + 2/Re/Pr/dx**2 + 1/dt
                    d(i) = theta(i,j)/dt
                enddo

                a(leftBoundary)  = 0; b(leftBoundary)  =1; c(leftBoundary) =  0; d(leftBoundary)  = 1;                
                a(rightBoundary) = 1; b(rightBoundary) =1; c(rightBoundary) = 0; d(rightBoundary) = 0;

                call Tom(leftBoundary, rightBoundary, a, b, c, d, e, 1, max(n,m))

                ! Присваивам в промежуточный слой
                do i=leftBoundary, rightBoundary
                    thetaTemp(i, j) = e(i)
                enddo
            enddo
            ! c j = bottomWallPoint    до j = topWallPoint
            !   i = rightWallPoint     до i = n
            do j=bottomWallPoint, topWallPoint
                leftBoundary = rightWallPoint
                rightBoundary = n

                do i=leftBoundary+1, rightBoundary-1
                    UxAbs = abs(Ux(i, j))
                    a(i) = (UxAbs + Ux(i, j))/dx/2 + 1/Re/Pr/dx**2
                    c(i) = (UxAbs - Ux(i, j))/dx/2 + 1./Re/Pr/dx**2
                    b(i) = UxAbs/dx + 2/Re/Pr/dx**2 + 1/dt
                    d(i) = theta(i,j)/dt
                enddo

                a(leftBoundary)  = 0; b(leftBoundary)  =1; c(leftBoundary) =  1; d(leftBoundary)  = 1;                
                a(rightBoundary) = 1; b(rightBoundary) =1; c(rightBoundary) = 0; d(rightBoundary) = 0;

                call Tom(leftBoundary, rightBoundary, a, b, c, d, e, 1, max(n,m))

                ! Присваивам в промежуточный слой
                do i=leftBoundary, rightBoundary
                    thetaTemp(i, j) = e(i)
                enddo
            enddo
            ! c j = topWallPoint+1 до j = m-1
            !   i = 1              до i = n
            do j=topWallPoint+1, m-1
                leftBoundary = 1
                rightBoundary = n

                do i=leftBoundary+1, rightBoundary-1
                    UxAbs = abs(Ux(i, j))
                    a(i) = (UxAbs + Ux(i, j))/dx/2 + 1/Re/Pr/dx**2
                    c(i) = (UxAbs - Ux(i, j))/dx/2 + 1./Re/Pr/dx**2
                    b(i) = UxAbs/dx + 2/Re/Pr/dx**2 + 1/dt
                    d(i) = theta(i,j)/dt
                enddo

                a(leftBoundary)  = 0; b(leftBoundary)  =1; c(leftBoundary) =  0; d(leftBoundary)  = 0;                
                a(rightBoundary) = 1; b(rightBoundary) =1; c(rightBoundary) = 0; d(rightBoundary) = 0;

                call Tom(leftBoundary, rightBoundary, a, b, c, d, e, 1, max(n,m))

                ! Присваивам в промежуточный слой
                do i=leftBoundary, rightBoundary
                    thetaTemp(i, j) = e(i)
                enddo
            enddo
            ! по Y
            ! c i = 2 до i = rightWallPoint
            ! c j = 1 до j = bottomWallPoint 
            do i=2, rightWallPoint
                lowerBoundary = 1
                upperBoundary = bottomWallPoint;

                do j= lowerBoundary + 1, upperBoundary - 1
                    UyAbs = abs(Uy(i, j))
                    a(j)= (UyAbs + Uy(i, j))/dy/2 + 1/Re/Pr/dy**2
                    c(j)= (UyAbs - Uy(i, j))/dy/2 + 1/Re/Pr/dy**2
                    b(j)= UyAbs/dy + 2./Re/Pr/dy**2 + 1./dt
                    d(j)= thetaTemp(i,j)/dt
                enddo

                a(lowerBoundary) = 0; b(lowerBoundary) = 1; c(lowerBoundary) = 1; d(lowerBoundary) = 0;                     
                a(upperBoundary) = 1; b(upperBoundary) = 1; c(upperBoundary) = 0; d(upperBoundary) = 0; 

                call Tom(lowerBoundary, upperBoundary, a, b, c, d, e, 1, max(n,m))

                do j= lowerBoundary, upperBoundary
                    thetaN1(i,j) = e(j)
                enddo
            enddo
            ! c i = rightWallPoint + 1 до i = n-1
            ! c j = 1                  до j = m
            do i=rightWallPoint + 1, n-1
                lowerBoundary = 1
                upperBoundary = m;

                do j= lowerBoundary + 1, upperBoundary - 1 
                    UyAbs = abs(Uy(i, j))
                    a(j)= (UyAbs + Uy(i, j))/dy/2 + 1/Re/Pr/dy**2
                    c(j)= (UyAbs - Uy(i, j))/dy/2 + 1/Re/Pr/dy**2
                    b(j)= UyAbs/dy + 2./Re/Pr/dy**2 + 1./dt
                    d(j)= thetaTemp(i,j)/dt
                enddo

                a(lowerBoundary) = 0; b(lowerBoundary) = 1; c(lowerBoundary) = 1; d(lowerBoundary) = 0; 
                a(upperBoundary) = 1; b(upperBoundary) = 1; c(upperBoundary) = 0; d(upperBoundary) = 0; 

                call Tom(lowerBoundary, upperBoundary, a, b, c, d, e, 1, max(n,m))

                do j= lowerBoundary, upperBoundary
                    thetaN1(i,j) = e(j)
                enddo
            enddo
            ! c i = 2            до i = rightWallPoint
            ! c j = topWallPoint до j = m
            do i=2, rightWallPoint
                lowerBoundary = topWallPoint
                upperBoundary = m;

                do j= lowerBoundary + 1, upperBoundary - 1
                    UyAbs = abs(Uy(i, j))
                    a(j)= (UyAbs + Uy(i, j))/dy/2 + 1/Re/Pr/dy**2
                    c(j)= (UyAbs - Uy(i, j))/dy/2 + 1/Re/Pr/dy**2
                    b(j)= UyAbs/dy + 2./Re/Pr/dy**2 + 1./dt
                    d(j)= thetaTemp(i,j)/dt
                enddo

                a(lowerBoundary) = 0; b(lowerBoundary) = 1; c(lowerBoundary) = 1; d(lowerBoundary) = 0; 
                a(upperBoundary) = 1; b(upperBoundary) = 1; c(upperBoundary) = 0; d(upperBoundary) = 0; 

                call Tom(lowerBoundary, upperBoundary, a, b, c, d, e, 1, max(n,m))

                do j= lowerBoundary, upperBoundary
                    thetaN1(i,j) = e(j)
                enddo
            enddo

            ! Граничные
            thetaN1(n, 2:m-1) = thetaN1(n-1, 2:m-1) ! правый выход
            thetaN1(2:n-1, 1) = thetaN1(2:n-1, 2) ! нижняя стенка
            thetaN1(2:n-1, m) = thetaN1(2:n-1, m-1) ! Верхняя стенка
            thetaN1(2:rightBoundary, bottomWallPoint) = thetaN1(2:rightBoundary, bottomWallPoint - 1) ! нижняя внутренняя стенка
            thetaN1(2:rightBoundary, topWallPoint)    = thetaN1(2:rightBoundary, topWallPoint + 1) ! верхняя внутренняя стенка

            do j=bottomWallPoint, topWallPoint
                thetaN1(rightWallPoint, j) = thetaN1(rightWallPoint+1, j) ! внутренняя вертикальная стенка
            enddo

            thetaConvergence=thetaN1-theta
            theta=thetaN1
        enddo
        ! Температура -------------


        ! Выводим результаты на определенной введенной итерации
        if (mod(curIterationNum,curPrintIterationNumber) == 0) then
            curPrintIterationNumber = curPrintIterationNumber + printStepNumber
            call calcConvergence(n, m, dt, curIterationNum, tok, tokn1, vihr, vihrn1, thetaConvergence);
        endif

        ! Переприсваиваем слои
        vihr = vihrn1;
        tok = tokn1;

        ! Условие продолжения
        if (curIterationNum < iterationNum) cycle

        print*
        write(*,"('Additional number iteration or type 0 for end:')")
        read(*,*) additionalIterationNum
        iterationNum = iterationNum + additionalIterationNum

        if (additionalIterationNum > 0) then
            curPrintIterationNumber = curPrintIterationNumber - printStepNumber
            write(*,"('New print step number:')")
            read(*,*) printStepNumber
            curPrintIterationNumber = curPrintIterationNumber + printStepNumber
        else
           exit ! Условие выхода из цикла
        endif
    enddo ! end while
    
    open(23, file=trim(filePrefix)//'_iso_theta_tok_ux.dat')
    do i=1, n
        do j=1, m
            tempX = (i-1)*dx
            tempY = (j-1)*dy
            write(23, '(7es14.4)') tempX, tempY, thetaN1(i,j), tokN1(i,j), ux(i,j)
        enddo
    enddo

    open(214, file=trim(filePrefix)//'_theta_ux.dat')
    do j= 1, m
        tempY = (j-1)*dy
        write(214,'(8es14.4)') tempY, thetaN1(n - 1, j), ux(n - 1, j) 
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

endsubroutine

subroutine calcConvergence(n, m, dt, curIterationNum, tok, tokn1, vihr, vihrn1, thetaConvergence)
    integer n, m, curIterationNum
    real, dimension(n, m):: tok, tokn1, vihr, vihrn1, thetaConvergence
    real, dimension(n, m):: tokConvergence, vihrConvergence
    real sumTokConvergence, sumVihrConvergence, sumThetaConvergence
    tokConvergence = tokn1 - tok
    vihrConvergence = vihrn1 - vihr

    sumTokConvergence = 0.
    sumVihrConvergence = 0.
    sumThetaConvergence=0.

    do i=2, n-1
        do j=2, m-1
        sumTokConvergence = sumTokConvergence + abs(tokConvergence(i,j))/dt
        sumVihrConvergence = sumVihrConvergence + abs(vihrConvergence(i,j))/dt
        sumThetaConvergence = sumThetaConvergence + abs(thetaConvergence(i,j))/dt
        enddo
    enddo 

    write(*,'(I7,3es14.4)') curIterationNum, sumTokConvergence, sumVihrConvergence, sumThetaConvergence
endsubroutine

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

endsubroutine

subroutine setBoundaryVihrValue(n, m, vihr, tok, bottomWallPoint, topWallPoint, rightWallPoint, dy, dx)
    integer n, m, bottomWallPoint, topWallPoint, rightWallPoint
    real,dimension (n, m):: vihr, tok
    real dy, dx

    vihr(1, 2:bottomWallPoint - 1) = 0; ! нижний вход граничные

    do j = bottomWallPoint, topWallPoint
        vihr(rightWallPoint, j) = 2 * ( tok(rightWallPoint + 1, j) - tok(rightWallPoint, j) ) / dx**2; ! средний проход граничные
    enddo

    vihr(1, topWallPoint + 1:m - 1) = 0; ! верхний вход граничные

    do i = 2, n
        vihr(i, 1) = 2 * ( tok(i, 2) - tok(i, 1) ) / dy**2; ! нижняя стенка
        vihr(i, m) = 2 * ( tok(i, m - 1) - tok(i, m) ) / dy**2; ! верхняя стенка
    enddo

    do i = 2, rightWallPoint
        vihr(i, topWallPoint)    = 2 * ( tok(i, topWallPoint + 1)    - tok(i, topWallPoint)    ) / dy**2; ! верхняя внутренняя стенка
        vihr(i, bottomWallPoint) = 2 * ( tok(i, bottomWallPoint - 1) - tok(i, bottomWallPoint) ) / dy**2; ! нижняя внутренняя стенка
    enddo

    vihr(n, 1:m) = vihr(n-1, 1:m); ! правый выход

endsubroutine

subroutine calcTokX(jStart, jEnd, leftBoundary, rightBoundary, n, m, tokTemp, tok, vihr, dx, dt)
    integer n, m, jStart, jEnd, leftBoundary, rightBoundary
    real, dimension(n):: a, b, c, d, e
    real, dimension(n, m):: tokTemp, tok, vihr

    do j=jStart, jEnd

        do i=leftBoundary+1, rightBoundary-1
            a(i) = 1./dx**2
            c(i) = 1./dx**2
            b(i) = 1./dt + 2./dx**2
            d(i) = tok(i, j)/dt - vihr(i, j)
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

endsubroutine

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
 
endsubroutine

elemental subroutine my_incr( var, incr )
  implicit none
  integer,intent(inout) :: var
  integer,intent(in)    :: incr

  var = var + incr
end subroutine
