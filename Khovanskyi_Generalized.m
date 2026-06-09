function Khovanskyi_Generalized()
    % Реалізація узагальненого методу Хованського для матричних рівнянь
    clear; clc;

    function [X, n_iter, residuals] = khovanskii_solver(A2, A1, A0, X0, l, k, eps, max_iter, step_size)
        % Обчислюємо X_n = (L*A2*X_{n-1} + L*A1 + K)^-1 * (K*X_{n-1} - L*A0)

        m = size(A2, 1);
        E = eye(m);
        L = l * E;
        K = k * E;

        X = X0;
        residuals = [];

        fprintf('\n Початок ітераційного процесу \n');

        for n = 1:max_iter
            X_old = X;

            % Перший множник у формулі
            Denominator = L * A2 * X_old + L * A1 + K;

            % Другий множник у формулі
            Numerator = K * X_old - L * A0;

            % Обчислення формули X_n = (L*A2*X_{n-1} + L*A1 + K)^-1 * (K*X_{n-1} - L*A0)
            X = Denominator \ Numerator;

            % Контроль приросту ||X_n - X_{n-1}||
            res = norm(X - X_old);
            residuals(n) = res;


            if mod(n, step_size) == 0
                fprintf('\n>>> Ітерація №%d\n', n);
                fprintf('Поточна норма приросту: %.6f\n', res);
                disp('Матриця X_n:');
                disp(X);
                fprintf('------------------------------------\n');
            end

            if res < eps
                n_iter = n;
                return;
            end
        end
        n_iter = max_iter;
        warning('Метод не збігся за максимальну кількість ітерацій');
    end


    function [X, n_iter, residuals] = khovanskii_solver_cubic(A3, A2, A1, A0, X0, l, k, eps, max_iter, step_size)
        m = size(A3, 1);
        E = eye(m);
        L = l * E;
        K = k * E;

        X = X0;
        Y0 = inv(X0); % Нам потрібно лиш Y0
        residuals = [];

        fprintf('\n Початок ітераційного процесу (Кубічне рівняння) \n');

        for iter = 1:max_iter
            X_old = X;
            Y0_old = Y0;

            % Загальний знаменник (L*An*X_{n-1} + K)
            Denominator = L * A3 * X_old + K;

            % Обчислення Y_0^{(n)}
            Num_Y0 = L * A3 + K * Y0_old;

            Y0 = Denominator \ Num_Y0;

            % Обчислення X^{(n)} (Для n=3)
            Num_X = (K - L * A2) * X_old - L * A1 - L * A0 * Y0;
            X = Denominator \ Num_X;

            res = norm(X - X_old);
            residuals(iter) = res;

            if mod(iter, step_size) == 0
                fprintf('\n>>> Ітерація №%d\n', iter);
                fprintf('Поточна норма приросту: %.6f\n', res);
                disp('Матриця X_n:');
                disp(X);
                fprintf('------------------------------------\n');
            end

            if res < eps
                n_iter = iter;
                return;
            end
        end
        n_iter = max_iter;
        warning('Метод не збігся за максимальну кількість ітерацій');
    end



    fprintf('УНІВЕРСАЛЬНИЙ МЕТОД ХОВАНСЬКОГО\n');
    fprintf('0. Приклад з курсової 1 (2x2)\n');
    fprintf('1. Приклад з курсової 2 (3x3)\n');
    fprintf('2. Приклад 1 (3x3)\n');
    fprintf('3. Приклад 2 (4x4)\n');
    fprintf('4. Власний приклад \n');
    fprintf('5. Приклад з кубічним многочленом (3x3)\n');
    choice = input('Оберіть опцію: ');

    if choice == 4
        fprintf('\n НАЛАШТУВАННЯ ВЛАСНОГО ПРИКЛАДУ\n');
        m = input('1. Розмірність матриць (m): ');

        fprintf('\nПІДКАЗКА: \n');
        fprintf('Приклад для 2x2: [1 2; 3 4]\n\n');

        A2 = input('2. Матриця A2: ');
        A1 = input('3. Матриця A1: ');
        A0 = input('4. Матриця A0: ');

        X0 = eye(m) ;

        l = input('6. Параметр l (скаляр): ');
        k = input('7. Параметр k (скаляр): ');
        eps = input('8. Точність eps (напр., 0.0001): ');
        max_iter = input('9. Максимальна кількість ітерацій (напр., 500): ');
        step_size = input('10. Частота виводу (через скільки ітерацій показувати X): ');

        name = 'Власний приклад';

    elseif choice == 5

        step_size = input('Через скільки ітерацій виводити матрицю: ');
        eps = 0.001;
        max_iter = 200;
        l = 1;
        k = 1;


        A3 = eye(3);
        A2 = -[39 3 3; 3 42 3; 3 3 45];
        A1 = [501 81 84; 81 582 87; 84 87 669];
        A0 = -[2109 526 565; 529 2647 607; 571 610 3269];

        X0 = [15.9 1 1; 1 16.9 1; 1 1 17.9];
        name = 'Приклад з кубічним характеристичним многочленом (3x3)';

        [X_final, n, res_history] = khovanskii_solver_cubic(A3, A2, A1, A0, X0, l, k, eps, max_iter, step_size);

        fprintf('\nОСТАТОЧНИЙ РЕЗУЛЬТАТ: %s \n', name);
        fprintf('Успішно завершено на ітерації: %d\n', n);
        fprintf('Фінальна норма приросту: %g\n', res_history(end));


        Residual_Matrix = A3 * (X_final^3) + A2 * (X_final^2) + A1 * X_final + A0;
        fprintf('Норма матричного рівняння: %g\n', norm(Residual_Matrix));

        disp('Фінальна матриця X:');
        disp(X_final);

        return;

    else
        step_size = input('Через скільки ітерацій виводити матрицю: ');
        eps = 0.001; max_iter = 1000; l = 1; k = 1;

        switch choice
            case 0
                A2 = eye(2);
                A1 = [-12 -2; -2 -14];
                A0 = [34 11; 11 47];
                X0 = [1 -3; 8 3];
                l = 2;
                k = 25;
                eps = 0.0001;
                name = 'Приклад з характеристичного многочлена (2x2)';
            case 1
                A2 = eye(3);
                A1 = [-5 -1 -1; -1 -4 -1; -1 -1 -3];
                A0 = [4 1 1; 1 3 1; 1 1 2];
                X0 = [3 2 1; 1 2 1; 2 1 2];
                l = 2;
                k = 10;
                eps = 0.0001;
                name = 'Приклад з характеристичного многочлена (3x3)';
            case 2
                A2 = eye(3);
                A1 = [1 2 3; 2 3 4; 3 4 5];
                A0 = [-13 -13 -14; -16 -18 -18; -20 -21 -23];
                X0 = eye(3); name = 'Приклад 1';
            case 3
                A2 = eye(4);
                A1 = [-1 0 2 1; 0 1 0 2; 0 0 4 1; 0 0 0 -5];
                A0 = [-8 -8 -10 -9; -9 -11 -9 -11; -11 -11 -16 -12; -1 -1 -1 3];
                X0 = eye(4); name = 'Приклад 2';
        end
    end


    [X_final, n, res_history] = khovanskii_solver(A2, A1, A0, X0, l, k, eps, max_iter, step_size);


    fprintf('\nОСТАТОЧНИЙ РЕЗУЛЬТАТ: %s \n', name);
    fprintf('Успішно завершено на ітерації: %d\n', n);
    fprintf('Фінальна норма приросту: %g\n', res_history(end));


    Residual_Matrix = A2 * (X_final^2) + A1 * X_final + A0;
    fprintf('Норма матричного рівняння: %g\n', norm(Residual_Matrix));

    disp('Фінальна матриця X:');
    disp(X_final);

end
