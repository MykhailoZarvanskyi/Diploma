function generalized_scalar_product_method()
  clc;
  % Головна функція для узагальненого методу скалярних добутків



  fprintf('\nУЗАГАЛЬНЕНИЙ МЕТОД СКАЛЯРНИХ ДОБУТКІВ \n\n');
  fprintf('Програма дозволяє знаходити власні значення та вектори блочних матриць.\n');
  fprintf('Матриця A розбита на n×n блоків, кожен блок має розмір p×p.\n\n');

  n = input('Введіть кількість блоків n: ');
  p = input('Введіть розмірність кожного блоку p: ');


  A = cell(n, n);


  fprintf('\nВиберіть спосіб введення матриці A:\n');
  fprintf('1 - Ввести вручну (покоординатно для кожного блоку)\n');
  fprintf('2 - Створити випадкову матрицю\n');
  fprintf('3 - Використати тестовий приклад\n');
  input_method = input('Ваш вибір (1/2/3): ');

  if input_method == 1

    fprintf('\nВВЕДЕННЯ БЛОЧНОЇ МАТРИЦІ A \n');
    fprintf('Потрібно ввести %d блоків, кожен розміром %d×%d.\n', n*n, p, p);
    fprintf('Введення виконується по рядках: спочатку елементи першого рядка блоку, потім другого і т.д.\n\n');

    for i = 1:n
      for j = 1:n
        fprintf('Блок A{%d,%d} розміром %d×%d:\n', i, j, p, p);
        A{i,j} = zeros(p, p);


        fprintf('Виберіть спосіб введення цього блоку:\n');
        fprintf('1 - Ввести покоординатно\n');
        fprintf('2 - Нульова матриця\n');
        fprintf('3 - Одинична матриця\n');
        fprintf('4 - Діагональна матриця зі сталим значенням\n');
        block_input = input('Спосіб (1/2/3/4): ');

        if block_input == 1

          for r = 1:p
            fprintf('Рядок %d (введіть %d чисел, розділених пробілами): ', r, p);
            row_input = input('', 's');

            row_values = str2num(row_input);

            if length(row_values) ~= p
               fprintf('Помилка! Введено %d елементів замість %d. Спробуйте ще раз.\n', length(row_values), p);
               r = r - 1;
               continue;
            end

        A{i,j}(r, :) = row_values;
       end

        elseif block_input == 2

          A{i,j} = zeros(p, p);
        elseif block_input == 3

          A{i,j} = eye(p, p);
        elseif block_input == 4

          diagonal_value = input('Введіть значення для діагоналі: ');
          A{i,j} = diagonal_value * eye(p, p);
        else
          error('Невідомий спосіб введення блоку.');
        end
      end
    end
  elseif input_method == 2

    for i = 1:n
      for j = 1:n
        A{i,j} = rand(p, p);
      end
    end
    fprintf('Створено випадкову матрицю.\n');
  elseif input_method == 3

    fprintf('\nДоступні тестові приклади:\n');

    if n == 3 && p == 2
      fprintf('1 - Приклад з general_scalar_product_method1.m (3×3 блоків по 2×2)\n');
      test_example = 1;
    elseif n == 2 && p == 3
      fprintf('2 - Приклад з general_scalar_product_method2.m (2×2 блоків по 3×3)\n');
      test_example = 2;
    else
      fprintf('3 - Проста матриця з одиничною діагоналлю (%d×%d блоків по %d×%d)\n', n, n, p, p);
      test_example = 3;
    end

    if n == 3 && p == 2 || n == 2 && p == 3
      test_example = input('Виберіть тестовий приклад: ');
    end

    if test_example == 1 && n == 3 && p == 2
      % Перший тестовий приклад з general_scalar_product_method1.m
      A{1,1} = [ 4  0; 0 6 ];
      A{1,2} = [ 1 0; 0 -11 ];
      A{1,3} = [ -4  0; 0 6 ];
      A{2,1} = [ 1  0; 0 1 ];
      A{2,2} = [ 0  0; 0 0 ];
      A{2,3} = [ 0  0; 0 0 ];
      A{3,1} = [ 0  0; 0 0 ];
      A{3,2} = [ 1 0; 0 1 ];
      A{3,3} = [ 0  0; 0 0 ];
      fprintf('Використано тестовий приклад 1 (general_scalar_product_method1.m).\n');
    elseif test_example == 2 && n == 2 && p == 3
      % Другий тестовий приклад з general_scalar_product_method2.m
      A6 = [
         5   1   1  -4  -1  -1;
         1   4   1  -1  -3  -1;
         1   1   3  -1  -1  -2;
         1   0   0   0   0   0;
         0   1   0   0   0   0;
         0   0   1   0   0   0
      ];
      A{1,1} = A6(1:3, 1:3);
      A{1,2} = A6(1:3, 4:6);
      A{2,1} = A6(4:6, 1:3);
      A{2,2} = A6(4:6, 4:6);
      fprintf('Використано тестовий приклад 2 (general_scalar_product_method2.m).\n');
    else

      for i = 1:n
        for j = 1:n
          if i == j
            A{i,j} = eye(p);
          else
            A{i,j} = zeros(p, p);
          end
        end
      end
      fprintf('Використано простий тестовий приклад з одиничною діагоналлю.\n');
    end
  else
    error('Невідомий метод введення.');
  end


  display_block_matrix(A, 'Введена блочна матриця A');


  X0 = cell(n, 1);

  fprintf('\nВВЕДЕННЯ ПОЧАТКОВОГО НАБЛИЖЕННЯ X0 \n');
  fprintf('Початкове наближення X0 - блочний вектор розміром %d×1, кожен блок %d×%d.\n\n', n, p, p);
  fprintf('Виберіть спосіб завдання початкового наближення X0:\n');
  fprintf('1 - Ввести вручну\n');
  fprintf('2 - Використати стандартне наближення (одиничні/нульові блоки)\n');
  fprintf('3 - Випадкові значення\n');

  if (n == 3 && p == 2) || (n == 2 && p == 3)
    fprintf('4 - Використати початкове наближення з тестового прикладу\n');
  end

  x0_method = input('Ваш вибір: ');

  if x0_method == 1

    fprintf('\nПотрібно ввести %d блоків, кожен розміром %d×%d.\n', n, p, p);
    for i = 1:n
      fprintf('Блок X0{%d} розміром %d×%d:\n', i, p, p);
      fprintf('Виберіть спосіб введення цього блоку:\n');
      fprintf('1 - Ввести покоординатно\n');
      fprintf('2 - Нульова матриця\n');
      fprintf('3 - Одинична матриця\n');
      fprintf('4 - Діагональна матриця зі сталим значенням\n');

      block_input = input('Спосіб (1/2/3/4): ');

      if block_input == 1

        X0{i} = zeros(p, p);
        for r = 1:p
          fprintf('Рядок %d (введіть %d чисел, розділених пробілами): ', r, p);
          row_input = input('', 's');


          row_values = str2num(row_input);


          if length(row_values) ~= p
            fprintf('Помилка! Введено %d елементів замість %d. Спробуйте ще раз.\n', length(row_values), p);
            r = r - 1;
            continue;
          end

          X0{i}(r, :) = row_values;
        end

      elseif block_input == 2

        X0{i} = zeros(p, p);
      elseif block_input == 3

        X0{i} = eye(p, p);
      elseif block_input == 4

        diagonal_value = input('Введіть значення для діагоналі: ');
        X0{i} = diagonal_value * eye(p, p);
      else
        error('Невідомий спосіб введення блоку.');
      end
    end
  elseif x0_method == 2

    for i = 1:n
      if i == 1
        X0{i} = eye(p);
      else
        X0{i} = zeros(p, p);
      end
    end
    fprintf('Використано одиничний блок для X0{1} та нульові для інших.\n');
  elseif x0_method == 3

    for i = 1:n
      X0{i} = rand(p, p);
    end
    fprintf('Створено випадкові блоки для X0.\n');
  elseif x0_method == 4 && n == 3 && p == 2

    X0{1} = [ 1  0; 0  3 ];
    X0{2} = [ 9  0; 0  1 ];
    X0{3} = [ 10  0; 0  8 ];
    fprintf('Використано початкове наближення з general_scalar_product_method1.m\n');
  elseif x0_method == 4 && n == 2 && p == 3

    X0{1} = 2*eye(p);
    X0{2} = zeros(p, p);
    fprintf('Використано початкове наближення з general_scalar_product_method2.m\n');
  else
    error('Невідомий метод введення X0.');
  end


  display_block_vector(X0, 'Початкове наближення X0');


  fprintf('\nПАРАМЕТРИ МЕТОДУ \n');
  tol = input('Введіть точність (напр., 1.0e-6): ');
  max_iter = input('Введіть максимальну кількість ітерацій (напр., 100): ');


  [Lambda, X, iter, convergence_history] = block_scalar_product_method(A, X0, tol, max_iter);


  fprintf('\nРЕЗУЛЬТАТИ УЗАГАЛЬНЕНОГО МЕТОДУ СКАЛЯРНИХ ДОБУТКІВ \n');
  fprintf('Кількість ітерацій: %d\n', iter);
  fprintf('Власне значення Lambda (%d×%d) =\n', p, p);
  disp(Lambda);

  fprintf('Власний блок-вектор X (%d×1), кожен блок %d×%d:\n', n, p, p);
  display_block_vector(X, 'X');


  AX = block_matvec(A, X);


  LambdaX = cell(n, 1);
  for i = 1:n
    LambdaX{i} = X{i} * Lambda;
  end


  diff = 0;
  for i = 1:n
    diff = diff + norm(AX{i} - LambdaX{i}, 'fro')^2;
  end
  diff = sqrt(diff);

  fprintf('Перевірка: ||A*X - X*Lambda|| = %g\n', diff);


  if ~isempty(convergence_history)
    figure;
    semilogy(1:length(convergence_history), convergence_history, 'b-o');
    grid on;
    xlabel('Ітерація');
    ylabel('||Lambda^{(k)} - Lambda^{(k-1)}||_F');
    title('Збіжність методу скалярних добутків');
  end


  save_results = input('\nЗберегти результати у файл? (1 - Так, 0 - Ні): ');
  if save_results == 1
    save_results_to_file(A, X0, X, Lambda, iter, diff, convergence_history, n, p, tol);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Lambda_k, X_k, iter_count, convergence_history] = block_scalar_product_method(A, X0, tol, max_iter)
  % Реалізація узагальненого методу скалярних добутків
  % для блочної матриці A (n×n), кожен блок p×p
  %
  % A   - (n×n) cell-масив із блоками p×p
  % X0  - (n×1) cell-масив, початковий блочний вектор
  % tol - точність
  % max_iter - макс. кількість ітерацій

  % Повертає:
  %   Lambda_k - власне значення (p×p матриця)
  %   X_k - власний блок-вектор (n×1 cell-масив із p×p блоками)
  %   iter_count - кількість виконаних ітерацій
  %   convergence_history - масив різниць норм Lambda_k між ітераціями

  n = size(A, 1);
  p = size(X0{1}, 1);

  % Перевірка, щоб X0 не був нульовим
  if block_fro_norm(X0) < 1.0e-15
    error('Початковий вектор X0 занадто малий або нульовий.');
  end

  % Початкове Lambda^(0) - нульова p×p матриця
  Lambda_prev = zeros(p, p);

  % Крок 1: нормування X0
  S0 = block_scalar(X0, X0);   % (X0, X0) -> p×p

  % Перевірка на близьку до сингулярної S0
  if rcond(S0) < 1.0e-12
    shift = 1.0e-12;
    fprintf('Увага: S0 майже сингулярна, додаємо shift=%g до діагоналі\n', shift);
    S0 = S0 + shift*eye(p);
  end

  M0 = sqrtm(S0);              % корінь від p×p матриці

  % Перевірка на близьку до сингулярної M0
  if rcond(M0) < 1.0e-12
    shift = 1.0e-12;
    fprintf('Увага: M0 майже сингулярна, додаємо shift=%g до діагоналі\n', shift);
    M0 = sqrtm(S0 + shift*eye(p));
  end

  M0_inv = inv(M0);

  for i = 1:n
    X0{i} = M0_inv * X0{i};
  end

  X_k = X0;
  iter_count = 0;
  convergence_history = [];

  fprintf('\nПочаток ітераційного процесу:\n');
  fprintf('-------------------------------\n');

  for k = 1:max_iter
    iter_count = k;

    % 2) Y^(k) = A * X^(k-1)
    Y_k = block_matvec(A, X_k);

    % 3) S^(k) = (Y^(k), Y^(k)),   T^(k) = (Y^(k), X^(k-1))
    S_k = block_scalar(Y_k, Y_k);   % p×p
    T_k = block_scalar(Y_k, X_k);   % p×p

    % 3.1) M^(k) = sqrtm(S_k). Перевірка на виродженість S_k
    if rcond(S_k) < 1.0e-12
      shift = 1.0e-12;
      fprintf('Увага: S^(%d) майже сингулярна, додаємо shift=%g\n', k, shift);
      S_k = S_k + shift*eye(p);
    end

    M_k = sqrtm(S_k);

    % Перевірка на виродженість M_k
    if rcond(M_k) < 1.0e-12
      shift = 1.0e-12;
      fprintf('Увага: M^(%d) майже сингулярна, додаємо shift=%g\n', k, shift);
      M_k = sqrtm(S_k + shift*eye(p));
    end

    M_k_inv = inv(M_k);

    % 3.2) X^(k) = M_k_inv * Y^(k)
    X_new = cell(n, 1);
    for i = 1:n
      X_new{i} = M_k_inv * Y_k{i};
    end

    % 3.3) Lambda^(k) = S_k * inv(T_k). Перевірка на виродженість T_k перед інверсією
    if rcond(T_k) < 1.0e-12
      shift = 1.0e-12;
      fprintf('Увага: T^(%d) майже сингулярна, додаємо shift=%g\n', k, shift);
      T_k = T_k + shift*eye(p);
    end

    Lambda_k = S_k * inv(T_k);

    % 4) Перевірка збіжності за Lambda
    diff_mat = Lambda_k - Lambda_prev;
    diff_val = norm(diff_mat, 'fro');  % Фробеніусова норма p×p


    convergence_history(k) = diff_val;

    fprintf('Ітерація %d: ||Lambda^(%d) - Lambda^(%d)|| = %g\n', k, k, k-1, diff_val);

    if diff_val < tol
      fprintf('Збіжність досягнута на ітерації %d\n', k);
      X_k = X_new;
      break;
    end

    Lambda_prev = Lambda_k;
    X_k = X_new;
  end

  if iter_count == max_iter && diff_val >= tol
    fprintf('Увага: Досягнуто максимальну кількість ітерацій без збіжності.\n');
    fprintf('Остання різниця: %g (> %g)\n', diff_val, tol);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = block_fro_norm(X)
  % Фробеніусова норма блочного вектора X (n×1), де X{i} - матриця p×p
  n = numel(X);
  s = 0;
  for i = 1:n
    s = s + sum(sum(X{i}.^2));
  end
  val = sqrt(s);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = block_scalar(X, Y)
  % "Скалярний добуток" (X, Y) = сума X_i * Y_i для i=1..n. Результат - p×p
  n = numel(X);
  p = size(X{1}, 1);
  S = zeros(p, p);
  for i = 1:n
    S = S + X{i} * Y{i};
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = block_matvec(A, X)
  % Y = A * X для блочної матриці A (n×n) та вектора X (n×1). Кожен блок і кожен елемент вектора - p×p
  n = size(A, 1);
  p = size(X{1}, 1);
  Y = cell(n, 1);
  for i = 1:n
    tmp = zeros(p, p);
    for j = 1:n
      tmp = tmp + A{i,j} * X{j};
    end
    Y{i} = tmp;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_block_matrix(A, name)
  % Виведення блочної матриці з форматуванням
  n = size(A, 1);

  fprintf('\n%s (%d×%d блоків):\n', name, n, n);
  fprintf('-------------------------------\n');

  for i = 1:n
    for j = 1:n
      fprintf('Блок A{%d,%d}:\n', i, j);
      disp(A{i,j});
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_block_vector(X, name)
  % Виведення блочного вектора з форматуванням
  n = numel(X);

  fprintf('\n%s (%d×1 блоків):\n', name, n);
  fprintf('-------------------------------\n');

  for i = 1:n
    fprintf('Блок %s{%d}:\n', name, i);
    disp(X{i});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_results_to_file(A, X0, X, Lambda, iter, diff, convergence_history, n, p, tol)
  % Зберігає результати розрахунку у текстовий файл


  timestamp = datestr(now, 'yyyymmdd_HHMMSS');
  filename = ['results_', timestamp, '.txt'];

  fid = fopen(filename, 'w');

  fprintf(fid, '====== РЕЗУЛЬТАТИ РОЗРАХУНКУ МЕТОДОМ СКАЛЯРНИХ ДОБУТКІВ ======\n');
  fprintf(fid, 'Дата і час: %s\n\n', datestr(now));

  fprintf(fid, 'Розмірність: %d×%d блоків, кожен блок %d×%d\n\n', n, n, p, p);

  fprintf(fid, 'Параметри методу:\n');
  fprintf(fid, '- Точність: %g\n', tol);
  fprintf(fid, '- Кількість ітерацій: %d\n\n', iter);

  fprintf(fid, 'Вхідна матриця A:\n');
  for i = 1:n
    for j = 1:n
      fprintf(fid, 'A{%d,%d} =\n', i, j);
      for r = 1:p
        for c = 1:p
          fprintf(fid, '%12.6g ', A{i,j}(r,c));
        end
        fprintf(fid, '\n');
      end
      fprintf(fid, '\n');
    end
  end

  fprintf(fid, 'Початкове наближення X0:\n');
  for i = 1:n
    fprintf(fid, 'X0{%d} =\n', i);
    for r = 1:p
      for c = 1:p
        fprintf(fid, '%12.6g ', X0{i}(r,c));
      end
      fprintf(fid, '\n');
    end
    fprintf(fid, '\n');
  end

  fprintf(fid, 'Результати:\n');
  fprintf(fid, 'Власне значення Lambda =\n');
  for r = 1:p
    for c = 1:p
      fprintf(fid, '%12.6g ', Lambda(r,c));
    end
    fprintf(fid, '\n');
  end
  fprintf(fid, '\n');

  fprintf(fid, 'Власний вектор X:\n');
  for i = 1:n
    fprintf(fid, 'X{%d} =\n', i);
    for r = 1:p
      for c = 1:p
        fprintf(fid, '%12.6g ', X{i}(r,c));
      end
      fprintf(fid, '\n');
    end
    fprintf(fid, '\n');
  end

  fprintf(fid, 'Точність: ||A*X - X*Lambda|| = %g\n\n', diff);

  fprintf(fid, 'Історія збіжності:\n');
  fprintf(fid, 'Ітерація\tРізниця\n');
  for k = 1:length(convergence_history)
    fprintf(fid, '%d\t%g\n', k, convergence_history(k));
  end

  fclose(fid);
  fprintf('Результати збережено у файл %s\n', filename);
end
