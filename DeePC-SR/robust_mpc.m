close all
clear
clc

rng(0);

%% System Setup
data = load('tank_sys2.mat');
A_d = data.A_d;
B_d = data.B_d;
C = data.C;
D = data.D;
T_s = data.T_s; % 假设为标量
sys = ss(A_d, B_d, C, D, T_s);
%% Params and Buffers
% 设置选项
robust = true;
TEC = true;
tol_opt = 1e-4;
opt_settings = struct('OptimalityTolerance', 1e-9, 'MaxIterations', 20000, 'ConstraintTolerance', 1e-9);

% 创建数据
n = 4; % 真实系统阶数
nu = 4; % 估计系统阶数
M = 10; % 连续应用最优输入的次数（多步）
noise_max = 0.002; % 测量噪声对y的影响
N = 400; % 用于预测的数据长度
L_true = 30; % 实际预测范围（不包括初始条件）
L = L_true + nu; % 完整预测范围（包括初始条件）
T = 600; % “闭环范围”（模拟长度）
m = 2; p = 2; % 输入和输出的数量

% 初始I/O轨迹用于模拟以获得一些数据
ui = rand(m, n); % 初始输入
xi0 = rand(n, 1);
xi = zeros(n, n);
xi(:, 1) = xi0;

for i = 1:n-1
    xi(:, i+1) = sys.A * xi(:, i) + sys.B * ui(:, i);
end

yi = sys.C * xi + sys.D * ui;
K = eye(m); % 预稳定控制器

% 终端约束
u_term = ones(m, 1);
y_term = sys.C * inv(eye(n) - sys.A) * sys.B * u_term;

% 设置MPC
% 成本矩阵
R = 1e-4 * eye(m); % 输入加权
Q = 3 * eye(p); % 输出加权
S = zeros(m, p);

% 构建Pi矩阵 
Pi = blkdiag(kron(eye(L),R), kron(eye(L),Q)); %注意检查这里

% 成本QP
if robust
    cost_sigma = 1e3;
    cost_alpha = 1e-1;
    H = 2 * [cost_alpha * eye(N-L+1), zeros(N-L+1, (m+p)*L), zeros(N-L+1, p*L); ...
            zeros((m+p)*L, N-L+1), Pi, zeros((m+p)*L, p*L); ...
            zeros(p*L, N-L+1 + (m+p)*L), cost_sigma * eye(p*L)];
    f = [zeros(N-L+1, 1); ...
         -2 * kron(eye(L),R) * repmat(u_term, L, 1); ...
         -2 * kron(eye(L),Q) * repmat(y_term, L, 1); ...
         zeros(p*L, 1)];
else
    cost_alpha = 0;
    H = 2 * [cost_alpha * eye(N-L+1), zeros(N-L+1, (m+p)*L); ...
            zeros((m+p)*L, N-L+1), Pi];
    f = [zeros(N-L+1, 1); ...
         -2 * kron(eye(L),R) * repmat(u_term, L, 1); ...
         -2 * kron(eye(L),Q) * repmat(y_term, L, 1)];
end

% 不等式约束
u_max = inf * ones(m, 1);
u_min = -inf * ones(m, 1);
y_max = inf * ones(p, 1);
y_min = -inf * ones(p, 1);

if robust
    sigma_max = inf * ones(p, 1);
    sigma_min = -sigma_max;
    ub = [ones(N-L+1, 1)*inf; ...
         repmat(u_max, L, 1); ...
         repmat(y_max, L, 1); ...
         repmat(sigma_max, L, 1)];
    lb = [-ones(N-L+1, 1)*inf; ...
         repmat(u_min, L, 1); ...
         repmat(y_min, L, 1); ...
         repmat(sigma_min, L, 1)];
else
    ub = [ones(N-L+1, 1)*1; ...
         repmat(u_max, L, 1); ...
         repmat(y_max, L, 1)];
    lb = [-ones(N-L+1, 1)*1; ...
         repmat(u_min, L, 1); ...
         repmat(y_min, L, 1)];
end

[u, x, y] = MDL_sim_prestab(sys, ui, yi, K, noise_max, 0, N);
% 构建Hankel矩阵
Hu = hankel(u(:), L, N-L+1, m); % 为什么用MDL_sim_prestab的输入
Hy = hankel(y(:), L, N-L+1, p);

% 构建等式约束
if robust
    if TEC % 带终端约束
        B = [Hu, -eye(m*L), zeros(p*L, p*L), zeros(p*L, p*L); ...
             Hy, zeros(p*L, m*L), -eye(p*L), -eye(p*L); ...
             zeros(m*nu, N-L+1), [eye(m*nu), zeros(m*nu, m*(L-nu))], zeros(m*nu, p*L), zeros(m*nu, p*L); ...
             zeros(p*nu, N-L+1), zeros(p*nu, m*L), [eye(p*nu), zeros(p*nu, p*(L-nu))], zeros(p*nu, p*L);...
             zeros(m*nu, N-L+1), [zeros(m*nu, m*(L-nu)), eye(m*nu)], zeros(m*nu, p*L), zeros(m*nu, p*L);...
             zeros(p*nu, N-L+1), zeros(p*nu, m*L), [zeros(p*nu, p*(L-nu)), eye(p*nu)], zeros(p*nu, p*L)];
    else % 不带终端约束
        B = [Hu, -eye(m*L), zeros(p*L), zeros(p*L); ...
             Hy, zeros(p*L, m*L), -eye(p*L), -eye(p*L); ...
             zeros(m*nu, N-L+1), [eye(m*nu), zeros(m*nu, m*(L-nu))], zeros(m*nu, p*L), zeros(m*nu, p*L); ...
             zeros(p*nu, N-L+1), zeros(p*nu, m*L), [eye(p*nu), zeros(p*nu, p*(L-nu))], zeros(p*nu, p*L)];
    end
else
    if TEC % 带终端约束
        B = [Hu, -eye(m*L), zeros(m*L, p*L); ...
             Hy, zeros(p*L, m*L), -eye(p*L); ...
             zeros(m*nu, N-L+1), eye(m*nu), zeros(m*nu, m*(L-nu)), zeros(m*nu, p*L); ...
             zeros(p*nu, N-L+1), zeros(p*nu, m*L), eye(p*nu), zeros(p*nu, p*(L-nu))];
    else % 不带终端约束
        B = [Hu, -eye(mL), zeros(mL, pL); ...
        Hy, zeros(pL, mL), -eye(pL); ...
        zeros(mnu, N-L+1), eye(mnu), zeros(mnu, m(L-nu)); ...
        zeros(pnu, N-L+1), zeros(pnu, mL), eye(pnu), zeros(pnu, p(L-nu))];
    end
end

% 初始I/O轨迹
u_init = 0.8 * ones(m, nu); % 初始输入
x0 = [0.4, 0.4, 0.5, 0.5]; % 初始状态
x_init = zeros(n, nu);
x_init(:, 1) = x0;

for i = 1:nu-1
    x_init(:, i+1) = sys.A * x_init(:, i) + sys.B * u_init(:, i);
end

y_init = sys.C * x_init + sys.D * u_init;

% 闭环存储变量
u_cl = zeros(m, T);
u_cl(:, 1:nu) = u_init; % 设置初始输入

y_cl = zeros(p, T);
y_cl(:, 1:nu) = y_init; % 设置初始输出

y_cl_noise = y_cl; % 用于噪声输出的y_cl副本

x_cl = zeros(n, T);
x_cl(:, 1) = x0; % 设置初始状态

% 模拟前nu步
for j = 1:nu
    x_cl(:, j+1) = sys.A * x_cl(:, j) + sys.B * u_cl(:, j);
end
% 开环存储变量
u_ol = zeros(m * L, T);
y_ol = zeros(p * L, T);
sigma_ol = zeros(p * L, T);
alpha_ol = zeros(N - L + 1, T);
u_init_store = zeros(m * nu, T);
y_init_store = zeros(p * nu, T);

% 候选解存储变量
% u_cand = u_ol;
% y_cand = y_ol;
% alpha_cand = alpha_ol;
% fval_cand = zeros(1, T);

if robust
    sol_store = zeros((m + 2 * p) * L + N - L + 1, T);
else
    sol_store = zeros((m + p) * L + N - L + 1, T);
end

sol_cand = sol_store;
fval = zeros(T);

% MPC循环
for j = nu + 1:M:T
    disp(j);
    % 更新等式约束
    if TEC
        c = [zeros((m + p) * L, 1); ...
            u_init(:); ...
            y_init(:); ...
            repmat(u_term, nu, 1); ...
            repmat(y_term, nu, 1)];
    else
        c = [zeros((m + p) * L, 1); ...
            u_init(:); ...
            y_init(:)];
    end
    % 定义并解决QP问题
    x = optimvar('x', length(H));
    objective = 0.5 * x' * H * x + f' * x;
    constraints = [B * x == c];
    prob = optimproblem('Objective', objective, 'Constraints', constraints);
    opts = optimoptions('quadprog', 'Algorithm', 'interior-point-convex');
    % options = optimoptions('quadprog','Display','off','OptimalityTolerance', 1e-9, 'MaxIterations', 50000);
    [x_result, fval_result, exitflag, output] = quadprog(H, f, [], [], B, c, [], [], [], opts);

    sol = x_result;
    fval(j) = fval_result + dot(repmat(y_term, L, 1).', kron(eye(L), Q) * repmat(y_term, L, 1)) + dot(repmat(u_term, L, 1).', kron(eye(L), R) * repmat(u_term, L, 1));
    sol_store(:, j) = sol;
    alpha_ol(:, j) = sol(1:N-L+1);
    u_ol(:, j) = sol(N-L+2:N-L+1+m*L);
    y_ol(:, j) = sol(N-L+2+m*L:N-L+1+(m+p)*L);
    if robust
        sigma_ol(:, j) = sol(N-L+2+(m+p)*L:N-L+1+(m+2*p)*L);
    end

    u_init_store(:, j-nu) = u_init(:);
    y_init_store(:, j-nu) = y_init(:);

    % 模拟闭环
    for k = j:min(j + M, T - 1)
        u_cl(:, k) = u_ol(m*n + (k - j)*m + 1:m*n + m + (k - j)*m, j);
        if k + 1 < T
            x_cl(:, k + 1) = sys.A * x_cl(:, k) + sys.B * u_cl(:, k);
        end
        y_cl(:, k) = sys.C * x_cl(:, k) + sys.D * u_cl(:, k);
        y_cl_noise(:, k) = y_cl(:, k) + noise_max * (-1 + 2 * rand(p, 1));

        % 设置下一次迭代的新初始条件
        u_init = [u_cl(:, k-n+1:k)];
        y_init = [y_cl_noise(:, k-n+1:k)];
        disp(u_init);
    end
end

figure;

% u_1子图
subplot(2, 2, 1);
plot(1:T, u_cl(1, :), 'DisplayName', 'u_1');
hold on;
plot(1:T, repmat(u_term(1), T, 1), 'DisplayName', 'u_{1,eq}');
legend;

% u_2子图
subplot(2, 2, 2);
plot(1:T, u_cl(2, :), 'DisplayName', 'u_2');
hold on;
plot(1:T, repmat(u_term(2), T, 1), 'DisplayName', 'u_{2,eq}');
legend;

% y_1子图
subplot(2, 2, 3);
plot(1:T, y_cl(1, :), 'DisplayName', 'y_1');
hold on;
plot(1:T, repmat(y_term(1), T, 1), 'DisplayName', 'y_{1,eq}');
legend;

% y_2子图
subplot(2, 2, 4);
plot(1:T, y_cl(2, :), 'DisplayName', 'y_2');
hold on;
plot(1:T, repmat(y_term(2), T, 1), 'DisplayName', 'y_{2,eq}');
legend;

% 显示图表
hold off;