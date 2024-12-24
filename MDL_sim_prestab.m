function [u, x, y] = MDL_sim_prestab(sys, u_init, y_init, K, noise_max, MULT_NOISE, N)
    % 参数
    n = size(y_init, 2);
    m = size(u_init, 1);
    p = size(y_init, 1);

    % 构建马尔可夫参数
    markov = zeros(n*p, n*m);
    for i = 1:n
        for j = 1:n
            if i > j
                markov((i-1)*p+1:i*p, (j-1)*m+1:j*m) = sys.C * sys.A^(i-j-1) * sys.B;
            elseif i == j
                markov((i-1)*p+1:i*p, (j-1)*m+1:j*m) = sys.D;
            end
        end
    end

    % 初始状态
    Obs = observability_matrix(sys.A, sys.C);
    x0 = Obs \ (y_init(:) - markov * u_init(:));

    % 检查x0的形状
    if length(x0) ~= n
        error('Initial state x0 must be of length %d', n);
    end

    % 基于模型的仿真
    x = zeros(n, N);
    x(:, 1) = x0;
    u = 2 * rand(m, N) - 1;  % 为什么乘二
    u(:, 1:n) = u_init;
    y = zeros(p, N);
    y(:, 1:n) = y_init;
    eps = rand(p, N);
    % 仿真循环
    for i = 1:N-1
        if i >= n
            u(:, i) = u(:, i) + K * y(:, i);
        end
        x(:, i+1) = sys.A * x(:, i) + sys.B * u(:, i);
        if MULT_NOISE
            y(:, i+1) = (sys.C * x(:, i+1) + sys.D * u(:, i+1)) * (1 + noise_max * (-1 + 2 * eps(:, i+1)));
        else
            y(:, i+1) = sys.C * x(:, i+1) + sys.D * u(:, i+1) + noise_max * (-1 + 2 * eps(:, i+1));
        end
    end
    
    % 为输出平铺u和y
    u = u(:);
    y = y(:);
    end

function Obs = observability_matrix(A, C)
    % 计算系统A和C的可观测性矩阵。
    n = size(A, 1);
    Obs = C;
    for i = 1:n-1
        Obs = [Obs; C * A^i];
    end
end

