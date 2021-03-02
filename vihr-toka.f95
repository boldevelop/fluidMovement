program vihrToka
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
    integer, parameter:: n=12, m=9
    integer, parameter::topWallPoint= nint(m * 2./3)
    integer, parameter::bottomWallPoint= nint(m/3.) + 1
    integer, parameter::rightWallPoint = nint(n/3.) + 1
     
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
    real x,y, tokConvergence, Ux1, Ux2
    integer lowerBoundary, upperBoundary, leftBoundary, rightBoundary
    tokConvergence = 0
    x=4
    y=3
    Ux1=0.2
    Ux2=0.6

    dy = y/(m-1)
    dx = x/(n-1)

    ! Начальные условия левая стенка
    do j=1, m

        if (j > 1 .AND. j <= m * 1/3) then
            tok(1, j) = (j - 1) * dy * Ux1  ! Нижний вход Ux1
        else if (j > m * 1/3 .AND. j <= topWallPoint) then 
            tok(1:(n/4 + 1), j) = (j - 1) * dy ! Установка внутренней стенки
        else if (j > topWallPoint .AND. j < m) then 
            tok(1, j) = (j - 1) * dy * Ux2 ! Верхний вход Ux2
        else
            tok(1, j) = (j - 1) * dy ! Стенки (границы)
        endif

    enddo
    tok(1:n, m) = y

    tokConvergence = calcConvergence(n, m, tok, tokn1)
    call PrintArray(n, m, tok)
    do while (tokConvergence .GT. 0.001)

        ! по X
        ! c j = 2 до j = m*1/3 | i = 2 до i = n - 1
        do j=2, m*1/3
            leftBoundary = 1
            rightBoundary = n

            do i=leftBoundary+1, rightBoundary-1
                a(i) = 1./dx**2
                c(i) = 1./dx**2
                b(i) = 1./dt + 2./dx**2
                d(i) = -tok(i, j)/dt + vihr(i, j)
            enddo
            
            a(leftBoundary) = 0;  b(leftBoundary)  = 1; c(leftBoundary)  = 0; d(leftBoundary)  = (j - 1) * dy * Ux1; 
            a(rightBoundary) = -1; b(rightBoundary) = 1; c(rightBoundary) = 0; d(rightBoundary) = 0;

            ! Вызов прогонки
            call Tom(leftBoundary, rightBoundary, a, b, c, d, e, leftBoundary, rightBoundary)

            ! Присваивам в промежуточный слой
            do i=leftBoundary, rightBoundary
                tokTemp(i, j) = e(i)
            enddo
        enddo

        ! c j = m*1/3+1 до j = m*2/3 | i = n/4 + 2 до i = n - 1
        do j=m*1/3+1, topWallPoint
            leftBoundary = n/4 + 1
            rightBoundary = n

            do i= leftBoundary + 1, rightBoundary - 1
                a(i) = 1./dx**2
                c(i) = 1./dx**2
                b(i) = 1./dt + 2./dx**2
                d(i) = -tok(i, j)/dt + vihr(i, j)
            enddo
            
            a(leftBoundary) = 0;  b(leftBoundary)  = 1; c(leftBoundary)  = 0; d(leftBoundary)  = 0; 
            a(rightBoundary) = -1; b(rightBoundary) = 1; c(rightBoundary) = 0; d(rightBoundary) = 0;

            ! Вызов прогонки
            call Tom(leftBoundary, rightBoundary, a, b, c, d, e, leftBoundary, rightBoundary)

            ! Присваивам в промежуточный слой
            do i=leftBoundary, rightBoundary
                tokTemp(i, j) = e(i)
            enddo
        enddo

        ! c j = m*2/3+1 до j = m-1 | i = 2 до i = n - 1
        do j=topWallPoint+1, m-1
            leftBoundary = 1
            rightBoundary = n

            do i= leftBoundary + 1, rightBoundary - 1
                a(i) = 1./dx**2
                c(i) = 1./dx**2
                b(i) = 1./dt + 2./dx**2
                d(i) = -tok(i, j)/dt + vihr(i, j)
            enddo

            a(leftBoundary) = 0;  b(leftBoundary)  = 1; c(leftBoundary)  = 0; d(leftBoundary)  = (j - 1) * dy * Ux2; 
            a(rightBoundary) = -1; b(rightBoundary) = 1; c(rightBoundary) = 0; d(rightBoundary) = 0;

            ! Вызов прогонки
            call Tom(leftBoundary, rightBoundary, a, b, c, d, e, leftBoundary, rightBoundary)

            ! Присваивам в промежуточный слой
            do i=leftBoundary, rightBoundary
                tokTemp(i, j) = e(i)
            enddo
        enddo

        ! по Y
        ! c i = 1 до i = n/4 + 1 | c j = 2 до j = m*1/3 
        do i=2, n/4 + 1
            lowerBoundary = 1
            upperBoundary = m*1/3 + 1;

            do j= lowerBoundary + 1, upperBoundary - 1 
                a(j) = 1./dy**2
                c(j) = 1./dy**2
                b(j) = 2./dy**2
                d(j) = tokTemp(i, j)
            enddo
            
            a(lowerBoundary) = 0; b(lowerBoundary) = 1; c(lowerBoundary) = 0; d(lowerBoundary) = 0;
            a(upperBoundary) = 0; b(upperBoundary) = 1; c(upperBoundary) = 0; d(upperBoundary) = 0;

            ! Вызов прогонки
            call Tom(lowerBoundary, upperBoundary, a, b, c, d, e, lowerBoundary, upperBoundary)

            ! Присваивам в промежуточный слой
            do j= lowerBoundary, upperBoundary 
                tokn1(i, j) = e(i)
            enddo
        enddo
        ! c i = n/4 + 2 до i = n-1 | c j = 2 до j = m-1
        do i= n/4 + 2, n-1
            lowerBoundary = 1
            upperBoundary = m

            do j= lowerBoundary + 1, upperBoundary - 1 
                a(j) = 1./dy**2
                c(j) = 1./dy**2
                b(j) = 2./dy**2
                d(j) = tokTemp(i, j)
            enddo
            
            a(lowerBoundary) = 0; b(lowerBoundary) = 1;  c(lowerBoundary) = 0; d(lowerBoundary) = 0;
            a(upperBoundary) = 0; b(upperBoundary) = 1; c(upperBoundary) = 0; d(upperBoundary) = 0;

            ! Вызов прогонки
            call Tom(lowerBoundary, upperBoundary, a, b, c, d, e, lowerBoundary, upperBoundary)

            ! Присваивам в промежуточный слой
            do j= lowerBoundary, upperBoundary 
                tokn1(i, j) = e(i)
            enddo
        enddo

        ! c i = 2 до i = n/4 + 1 | c j = m*2/3+1 до j = m-1
        do i= 2, n/4 + 1
            lowerBoundary = topWallPoint
            upperBoundary = m

            do j= lowerBoundary + 1, upperBoundary - 1 
                a(j) = 1./dy**2
                c(j) = 1./dy**2
                b(j) = 2./dy**2
                d(j) = tokTemp(i, j)
            enddo
            
            a(lowerBoundary) = 0; b(lowerBoundary) = 1;  c(lowerBoundary) = 0; d(lowerBoundary) = 0;
            a(upperBoundary) = 0; b(upperBoundary) = 1; c(upperBoundary) = 0; d(upperBoundary) = 0;

            ! Вызов прогонки
            call Tom(lowerBoundary, upperBoundary, a, b, c, d, e, lowerBoundary, upperBoundary)

            ! Присваивам в промежуточный слой
            do j= lowerBoundary, upperBoundary 
                tokn1(i, j) = e(i)
            enddo
        enddo

        ! Начальные условия левая стенка
        do j=1, m

            if (j > 1 .AND. j <= m * 1/3) then
                tokn1(1, j) = (j - 1) * dy * Ux1 ! Нижний вход Ux1
            else if (j > m * 1/3 .AND. j <= topWallPoint) then
                tokn1(1:(n/4 + 1), j) = (j - 1) * dy ! Установка внутренней стенки
            else if (j > topWallPoint .AND. j < m) then 
                tokn1(1, j) = (j - 1) * dy * Ux2 ! Верхний вход Ux2
            else
                tokn1(1, j) = (j - 1) * dy ! Стенки (границы)
            endif
    
            tokn1(n, j) = tokn1(n-1, j)
        enddo
        tokn1(1:n, m) = y

        tokConvergence = calcConvergence(n, m, tok, tokn1)
        print *, tokConvergence;

        vihr = vihrn1;
        
    enddo ! end while

    call PrintArray(n, m, tokn1)
    

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
