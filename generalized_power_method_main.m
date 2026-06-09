function generalized_power_method_main()

  % Узагальнений степеневий метод для блокових матриць



  n = input('Введіть кількість блоків n (розмір матриці n×n блоків): ');
  p = input('Введіть розмір кожного блоку p (блоки p×p): ');


  total_size = n * p;


  disp('Введення блокової матриці:');
  A = cell(n, n);
  for i = 1:n
    for j = 1:n
      fprintf('Введення блоку A{%d,%d} (розмір %d×%d):\n', i, j, p, p);
      A{i,j} = input_matrix(p, p);
    end
  end


  disp('Введення початкового наближення X(0):');
  X0 = cell(n, 1);
  for i = 1:n
    fprintf('Введення блоку X0{%d} (розмір %d×%d):\n', i, p, p);
    user_choice = input('Ввести вручну (1) чи використати одиничну/нульову матрицю (2)? ');

    if user_choice == 1
      X0{i} = input_matrix(p, p);
    else
      block_type = input('Введіть тип блоку (1 - одинична матриця, 0 - нульова матриця): ');
      if block_type == 1
        X0{i} = eye(p);
      else
        X0{i} = zeros(p, p);
      end
    end
  end


  Lambda0 = cell(n, 1);
  for i = 1:n
    Lambda0{i} = zeros(p, p);
  end


  tol = input('Введіть точність обчислень (наприклад, 1e-6): ');
  max_iter = input('Введіть максимальну кількість ітерацій: ');


  [Lambda, X, iter, convergence_history] = generalized_power_method_core(A, X0, Lambda0, tol, max_iter);


  fprintf('\nРЕЗУЛЬТАТИ УЗАГАЛЬНЕНОГО СТЕПЕНЕВОГО МЕТОДУ\n');
  fprintf('Кількість виконаних ітерацій: %d\n', iter);
  fprintf('Точність на останній ітерації: %e\n', convergence_history(end));


  fprintf('\nЗнайдене власне значення (Lambda):\n');
  disp(Lambda);


  fprintf('Нормований власний блок-вектор X:\n');
  for i = 1:n
    fprintf('X{%d} =\n', i);
    disp(X{i});
  end


  fprintf('\nПеревірка точності обчислень:\n');




  figure;
  semilogy(1:iter, convergence_history, 'b-o');
  grid on;
  xlabel('Ітерація');
  ylabel('Похибка');
  title('Графік збіжності методу');
end

function [Lambda, X, iter_count, convergence_history] = generalized_power_method_core(A, X0, Lambda0, tol, max_iter)

  % Реалізація узагальненого степеневого методу для блокової матриці
  % Вхідні параметри:
  %   A - клітинна матриця n×n блоків розміром p×p
  %   X0 - початкове наближення блок-вектора (n блоків розміром p×p)
  %   Lambda0 - початкове наближення для власного значення (n блоків)
  %   tol - точність обчислень
  %   max_iter - максимальна кількість ітерацій

  % Вихідні параметри:
  %   Lambda - наближене власне значення
  %   X - наближений власний блок-вектор
  %   iter_count - кількість виконаних ітерацій
  %   convergence_history - історія збіжності


  n = size(A, 1);
  convergence_history = [];

  % Початкова нормалізація X(0)
  normX0 = block_fro_norm(X0);
  if normX0 < 1e-15
    error('Початковий вектор X(0) є нульовим або майже нульовим.');
  end

  for i = 1:n
    X0{i} = X0{i} / normX0;
  end

  Lold = Lambda0;
  iter_count = 0;

  for k = 1:max_iter
    iter_count = k;

    % КРОК 2: Y^(k) = A * X^(k-1)
    Y = block_matvec(A, X0);

    % КРОК 3: Нормування Y^(k) для отримання X^(k)
    normY = block_fro_norm(Y);
    if normY < 1e-15
      warning('Y^(k) став майже нульовим. Метод не може продовжуватись.');
      break;
    end

    Xk = cell(n, 1);
    for i = 1:n
      Xk{i} = Y{i} / normY;
    end

    % КРОК 4: Lambda_i^(k) = Y_i^(k) * [X_i^(k-1)]^-1
    Lnew = cell(n, 1);
    for i = 1:n
      % Перевірка на виродженість X_i^(k-1)
      d = det(X0{i});
      if abs(d) < 1e-15
        % Якщо матриця вироджена, додаємо малий зсув до діагоналі
        warning('Виродженість: det(X0{%d}) ~ 0. Додаємо зсув...', i);
        X0{i} = X0{i} + 1e-12 * eye(size(X0{i}));
      end

      % Обчислення Lambda_i^(k), яке еквівалентно Y{i} * inverse(X0{i})
      Lnew{i} = Y{i} / X0{i};
    end

    % КРОК 5: Перевірка збіжності
    diff_sum = 0;
    for i = 1:n
      diff_block = Lnew{i} - Lold{i};
      diff_sum = diff_sum + sumsq(diff_block(:));
    end
    diff_val = sqrt(diff_sum);


    convergence_history(k) = diff_val;


    fprintf('Ітерація %d: похибка = %e\n', k, diff_val);

    if diff_val < tol

      X0 = Xk;
      Lold = Lnew;
      break;
    end


    X0 = Xk;
    Lold = Lnew;
  end

  % Обчислення фінального Lambda як середнього з Lambda_i^(k)
  Lsum = zeros(size(Lold{1}));
  for i = 1:n
    Lsum = Lsum + Lold{i};
  end
  Lambda = Lsum / n;
  X = X0;
end

function Y = block_matvec(A, X)
  % Множення блочної матриці A на блочний вектор X

  n = size(A, 1);
  Y = cell(n, 1);

  for i = 1:n
    S = zeros(size(X{1}));
    for j = 1:n
      S = S + A{i,j} * X{j};
    end
    Y{i} = S;
  end
end

function val = block_fro_norm(X)
  % Обчислення блочної норми Фробеніуса для вектора X

  n = numel(X);
  s = 0;

  for i = 1:n
    s = s + sumsq(X{i}(:));
  end

  val = sqrt(s);
end

function M = input_matrix(rows, cols)
  % Функція для зручного введення матриці

  fprintf('Введіть матрицю %d×%d (розділяйте рядки символом ;, а елементи пробілами):\n', rows, cols);
  fprintf('Наприклад: [1 2 3; 4 5 6; 7 8 9] для матриці 3×3\n');

  valid_input = false;
  while ~valid_input
    try
      input_str = input('Матриця: ', 's');
      if strncmp(input_str, '[', 1) && strcmp(input_str(end), ']')

        M = eval(input_str);
      else

        M = eval(['[' input_str ']']);
      end

      [r, c] = size(M);
      if r == rows && c == cols
        valid_input = true;
      else
        fprintf('Помилка: розмір введеної матриці (%d×%d) не відповідає очікуваному (%d×%d)\n', r, c, rows, cols);
      end
    catch
      fprintf('Помилка: неправильний формат введення. Спробуйте знову.\n');
    end
  end
end

