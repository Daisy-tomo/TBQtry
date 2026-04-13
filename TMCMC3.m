clear; clc; close all;
rng(206, 'twister');

%% 1. 参数定义 / Определение параметров
param_labels       = {'eta_к', 'eta_т', 'eta_в', 'eta_тв', 'eta_с2', 'lambda'};
param_labels_latex = {'eta_к', 'eta_т', 'eta_в', 'eta_тв', 'eta_с2', '$lambda$'};
n_params = 6;   % 待辨识参数个数 / Количество идентифицируемых параметров

% 均匀先验的上下界 / Нижняя и верхняя границы равномерного априорного распределения 
lb = [0.80,  0.80,  0.75,  0.85,  0.83, 0.9];   % 下界 / Нижняя граница
ub = [0.90,  0.99,  0.95,  0.97,  0.99, 1.1];   % 上界 / Верхняя граница

% “真值”（用于生成虚拟观测数据） / Истинные значения (для генерации виртуальных наблюдений)
theta_true = [0.86, 0.92, 0.86, 0.91, 0.93, 1.030];

% 固定参数向量 / Вектор фиксированных параметров
% eta_m, eta_c1, sigma_cc, sigma_kan, sigma_kask, sigma_ks / eta_m, eta_c1, sigma_cc, sigma_kan, sigma_kask, sigma_ks
theta_fixed = [0.995, 0.95, 1.000, 0.98, 0.99, 0.960];

% 工况参数（地面静态条件） / Параметры режима (наземное статическое состояние)
cond.T_H      = 288;       % 大气总温 (K) / Полная температура атмосферы (K)
cond.M_flight = 0.0;       % 飞行马赫数 / Число Маха полета
cond.m        = 10.0;      % 涵道比 / Степень двухконтурности
cond.pi_k     = 33.0;      % 压气机增压比 / Степень повышения давления компрессора
cond.T_g      = 1700.0;    % 涡轮前总温 (K) / Полная температура перед турбиной (K)

% 检查真值是否在先验范围内 / Проверка, попадает ли истинное значение в априорный диапазон
assert(all(theta_true >= lb) && all(theta_true <= ub), ...
    'Ошибка: theta_true выходит за пределы априорного диапазона [lb, ub]');

fprintf('Параметры режима: T_H=%g K, M=%g, m=%g, pi_k=%g, T_g=%g K\n\n', ...
    cond.T_H, cond.M_flight, cond.m, cond.pi_k, cond.T_g);

%% 2. 前向模型验证 / Проверка прямой модели
fprintf('Проверка прямой модели\n');
[y_true, ~] = engine_forward(theta_true, theta_fixed, cond);
if ~all(isfinite(y_true))
    error('Прямая модель возвращает нечисловой результат в точке theta_true; проверьте настройки параметров');
end
fprintf('  R_ud = %.4f  [Н·с/кг]\n',   y_true(1));
fprintf('  C_ud = %.6f  [кг/(Н·ч)]\n\n', y_true(2));

%% 3. 生成虚拟观测数据 / Генерация виртуальных наблюдений
% 噪声水平：相对标准差 0.1% / Уровень шума: относительное стандартное отклонение 0.1% (как в исходном коде)
data = generate_virtual_data(theta_true, theta_fixed, cond, 0.001, 0.001);

fprintf('Виртуальные наблюдения \n');
fprintf('  R: истинное=%.4f, наблюдение=%.4f, sigma_R=%.4f\n', ...
    data.R_true, data.R_obs, data.sigma_R);
fprintf('  C: истинное=%.6f, наблюдение=%.6f, sigma_C=%.6f\n\n', ...
    data.C_true, data.C_obs, data.sigma_C);

%% 4. TMCMC 算法配置 / Настройка алгоритма TMCMC
opts.N          = 3000;    % 粒子数 / Число частиц
opts.COV_target = 1.0;      % 目标权重变异系数 / Целевой коэффициент вариации весов
opts.N_MH       = 8;        % 每粒子每阶段 Metropolis-Hastings 步数 / Число шагов Metropolis-Hastings на частицу на каждом этапе
opts.scale      = 0.80;     % 提议协方差缩放因子 / Коэффициент масштабирования ковариации предложения
opts.max_stages = 60;       % 最大阶段数上限 / Максимально допустимое число этапов

fprintf('===== Параметры TMCMC =====\n');
fprintf('  Число частиц N=%d, COV_target=%.1f, N_MH=%d, scale=%.2f\n\n', ...
    opts.N, opts.COV_target, opts.N_MH, opts.scale);

%% 5. 运行 TMCMC / Запуск TMCMC
fprintf('===== Запуск TMCMC =====\n');
results = run_tmcmc(data, theta_fixed, cond, lb, ub, opts);
fprintf('\nTMCMC завершен, всего %d переходных этапов\n\n', results.n_stages);

%% 6. 后验统计结果 / Апостериорная статистика
fprintf('===== Апостериорная статистика =====\n');
fprintf('%-12s %8s %10s %10s %10s %10s\n', ...
    'Параметр', 'Истина', 'Среднее', 'MAP', 'CI95_низ', 'CI95_верх');
fprintf('%s\n', repmat('-', 1, 66));
for i = 1:n_params
    fprintf('%-12s %8.4f %10.4f %10.4f %10.4f %10.4f\n', ...
        param_labels{i}, theta_true(i), results.theta_mean(i), ...
        results.theta_map(i), results.theta_ci95(i,1), results.theta_ci95(i,2));
end
fprintf('\n');

%% 7. 后验预测对比 / Сравнение апостериорных прогнозов
[y_mn, ~] = engine_forward(results.theta_mean, theta_fixed, cond);
[y_mp, ~] = engine_forward(results.theta_map,  theta_fixed, cond);

fprintf('===== Сравнение апостериорных прогнозов =====\n');
fprintf('  %-24s %14s %16s\n', '', 'R_ud [Н·с/кг]', 'C_ud [кг/(Н·ч)]');
fprintf('  %-24s %14.4f %16.6f\n', 'Истинное значение', data.R_true, data.C_true);
fprintf('  %-24s %14.4f %16.6f\n', 'Наблюдение',        data.R_obs,  data.C_obs);
fprintf('  %-24s %14.4f %16.6f\n', 'Прогноз по среднему', y_mn(1), y_mn(2));
fprintf('  %-24s %14.4f %16.6f\n', 'MAP-прогноз',        y_mp(1), y_mp(2));
fprintf('\n');
fprintf('  Относительная ошибка апостериорного среднего: R=%.3f%%, C=%.3f%%\n', ...
    abs(y_mn(1) - data.R_true) / data.R_true * 100, ...
    abs(y_mn(2) - data.C_true) / data.C_true * 100);
fprintf('  Относительная ошибка MAP: R=%.3f%%, C=%.3f%%\n\n', ...
    abs(y_mp(1) - data.R_true) / data.R_true * 100, ...
    abs(y_mp(2) - data.C_true) / data.C_true * 100);

%% 8. 绘图 / Построение графиков
plot_results(results, theta_true, lb, ub, param_labels, param_labels_latex);


%% 分段绝热指数 / Кусочно-заданный показатель адиабаты
function kT = piecewise_kT(T_g)
    if T_g > 1600
        kT = 1.25;
    elseif T_g > 1400 && T_g <= 1600
        kT = 1.30;
    else
        kT = 1.33;
    end
end

%% 分段气体常数 / Кусочно-заданная газовая постоянная
function RT = piecewise_RT(T_g)
    if T_g > 1600
        RT = 288.6;
    elseif T_g > 1400 && T_g <= 1600
        RT = 288.0;
    else
        RT = 287.6;
    end
end

%% 冷却引气系数 / Коэффициент отбора воздуха на охлаждение
function d = delta_cooling(T_g)
    d = max(0.0, min(0.15, 0.02 + (T_g - 1200) / 100 * 0.02));
end

%% 前向模型 / Прямая модель
function [y, aux] = engine_forward(theta_s, theta_f, cond)

    % 待辨识参数 / Идентифицируемые параметры
    eta_k  = theta_s(1);    % 压气机效率 / КПД компрессора
    eta_t  = theta_s(2);    % 涡轮效率 / КПД турбины
    eta_v  = theta_s(3);    % 进气道效率修正 / КПД вентилятора
    eta_tv = theta_s(4);    % 推力效率修正 / КПД турбины вентилятора
    eta_c2 = theta_s(5);    % 外涵喷管速度系数/КПД сопла второго контура
    lambda = theta_s(6);    % 过量空气系数 / Коэффициент избытка воздуха

    % 固定参数 / Фиксированные параметры
    eta_m      = theta_f(1);    % 机械效率 / механический  КПД  турбокомпрессора 
    eta_c1     = theta_f(2);    % 内涵喷管速度系数 / КПД сопла первого контура
    sigma_cc   = theta_f(3);    % 燃烧室总压恢复/Коэф. восст. давления в скачках уплотнения сверхзвукового входного устройства
    sigma_kan  = theta_f(4);    % 进气道总压恢复 / Коэф. восст. давления в канале входного устройства
    sigma_kask = theta_f(5);    % 排气管道总压恢复 / Коэф. восст. давления между каскадами компрессора
    sigma_ks   = theta_f(6);    % 喷管总压恢复 / Коэф. восст. давления в камере сгорания

    % 工况 / Рабочий режим
    T_H      = cond.T_H;        % 大气总温 / Полная температура атмосферы
    M_flight = cond.M_flight;   % 飞行马赫数 / Число Маха полета
    m        = cond.m;          % 涵道比 / Степень двухконтурности
    pi_k     = cond.pi_k;       % 压气机增压比 / Степень повышения давления компрессора
    T_g      = cond.T_g;        % 涡轮前总温 / Полная температура перед турбиной

    % 初始化输出 / Инициализация выходных величин
    y   = [NaN, NaN];
    aux = struct();

    try
        % 空气物性常数 / Константы свойств воздуха
        k_air = 1.4;        % 空气绝热指数 / Показатель адиабаты воздуха
        R_air = 287.3;      % 空气气体常数 / Газовая постоянная воздуха 

        % (1) 声速与飞行速度 / Скорость звука и скорость полета
        a = sqrt(k_air * R_air * T_H);
        if ~isfinite(a) || a <= 0, return; end
        V_flight = a * M_flight;

        % (2) 燃气物性 / Свойства газа
        kT = piecewise_kT(T_g);     % 燃气绝热指数 / Показатель адиабаты газа
        RT = piecewise_RT(T_g);     % 燃气气体常数 / Газовая постоянная газа
        d  = delta_cooling(T_g);    % 冷却引气系数 / Коэффициент отбора воздуха на охлаждение

        % (3) степень повышения давления во входном устройства,температура на входе компрессора
        inner = 1 + V_flight^2 / (2 * (k_air / (k_air - 1)) * R_air * T_H);
        if inner <= 0, return; end
        tau_v = inner^(k_air / (k_air - 1));   % 进气道总压比 / степень повышения давления во входном устройства
        T_B   = T_H * (inner^k_air);           % 压气机入口总温 / температура на входе компрессора
        if ~isfinite(T_B) || T_B <= 0, return; end

        % (4) 压气机出口总温 / температура воздуха за компрессором
        pi_k_ratio = pi_k^((k_air - 1) / k_air);
        if ~isfinite(pi_k_ratio) || pi_k_ratio < 1, return; end
        T_k = T_B * (1 + (pi_k_ratio - 1) / eta_k);
        if ~isfinite(T_k) || T_k <= 0, return; end
        g_T = 3e-5 * T_g - 2.69e-5 * T_k - 0.003;
        if ~isfinite(g_T) || g_T <= 0, return; end

        % (5) 燃烧室前总压恢复系数 / коэффициент восстановления давления во входном устройстве
        sigma_bx = sigma_cc * sigma_kan;

        % (6) 单位自由能 удельная свободная энергия
        exp_T = (kT - 1) / kT;
        expansion_pr_denom = tau_v * sigma_bx * pi_k * sigma_kask * sigma_ks;
        if expansion_pr_denom <= 0, return; end
        expansion_term = (1.0 / expansion_pr_denom)^exp_T;
        if ~isfinite(expansion_term), return; end
        term1 = (kT / (kT - 1)) * RT * T_g * (1 - expansion_term);
        compress_work = (k_air / (k_air - 1)) * R_air * T_B * ...
                        (pi_k^((k_air - 1) / k_air) - 1);
        denom2 = (1 + g_T) * eta_k * eta_t * eta_m * (1 - d);
        if abs(denom2) < 1e-10, return; end
        term2 = compress_work / denom2;

        L_sv = lambda * (term1 - term2);
        if ~isfinite(L_sv) || L_sv <= 0, return; end

        % (7) 最优自由能分配系数 Оптимальное значение коэффициента распределения 
        % свободной энергии между контурами на заданном режиме полёта
        % с раздельными соплами
        V2_term = m * V_flight^2;
        num_xpc = 1 + V2_term / (2 * L_sv * eta_tv * eta_v * eta_c2);
        den_xpc = 1 + (m * eta_tv * eta_v * eta_c2) / (eta_c1 * lambda);
        if abs(den_xpc) < 1e-10, return; end
        x_pc = num_xpc / den_xpc;
        if ~isfinite(x_pc) || x_pc <= 0, return; end

        % (8) 内外涵喷流速度与比推力 / Скорости струй внутреннего и внешнего контуров и удельная тяга
        inner_sq1 = 2 * eta_c1 * lambda * x_pc * L_sv;
        if inner_sq1 < 0, return; end
        V_j1 = (1 + g_T) * sqrt(inner_sq1) - V_flight;    % 内涵喷流 / Струя внутреннего контура

        inner_sq2 = 2 * (1 - x_pc) / m * L_sv * eta_tv * eta_v * eta_c2 ...
                    + V_flight^2;
        if inner_sq2 < 0, return; end
        V_j2 = sqrt(inner_sq2) - V_flight;                 % 外涵喷流 / Струя внешнего контура

        R_ud = (1 / (1 + m)) * V_j1 + (m / (1 + m)) * V_j2;  % 比推力 / Удельная тяга
        if ~isfinite(R_ud) || R_ud <= 0, return; end

        % (9) 比油耗 / Удельный расход топлива
        denom_C = R_ud * (1 + m);
        if abs(denom_C) < 1e-10, return; end
        C_ud = 3600 * g_T / denom_C;
        if ~isfinite(C_ud) || C_ud <= 0, return; end

        y = [R_ud, C_ud];

        aux.T_B         = T_B;
        aux.T_k         = T_k;
        aux.tau_v       = tau_v;
        aux.g_T         = g_T;
        aux.sigma_bx    = sigma_bx;
        aux.L_sv        = L_sv;
        aux.x_pc        = x_pc;
        aux.kT          = kT;
        aux.RT          = RT;
        aux.delta       = d;

    catch ME
        warning('Исключение в engine_forward: %s', ME.message);
    end
end

%% ---------- 生成虚拟观测数据 / Генерация виртуальных наблюдений ----------
function data = generate_virtual_data(theta_true, theta_f, cond, nR, nC)
    [yt, ~] = engine_forward(theta_true, theta_f, cond);
    if ~all(isfinite(yt))
        error('generate_virtual_data: прямая модель неработоспособна в точке истинных значений');
    end

    % 绝对噪声标准差 = 相对噪声 × 真值 / Абсолютное стандартное отклонение шума = относительный шум x истинное значение
    sR = nR * abs(yt(1));
    sC = nC * abs(yt(2));

    data.R_true  = yt(1);
    data.C_true  = yt(2);
    data.R_obs   = yt(1) + sR * randn();   % 加高斯噪声 / Добавление гауссова шума
    data.C_obs   = yt(2) + sC * randn();
    data.sigma_R = sR;
    data.sigma_C = sC;
end

%% 对数似然函数 / Функция логарифма правдоподобия
% 假设两个观测量独立高斯噪声： / Предполагается, что два наблюдаемых значения имеют независимый гауссов шум:
% ln L = -0.5 * [ ((R_obs - R_pred)/sigma_R)^2 + ((C_obs - C_pred)/sigma_C)^2 ] / ln L = -0.5 * [ ((R_obs - R_pred)/sigma_R)^2 + ((C_obs - C_pred)/sigma_C)^2 ]
function ll = log_likelihood(theta, theta_f, data, cond)
    ll = -Inf;
    [yp, ~] = engine_forward(theta, theta_f, cond);
    if ~all(isfinite(yp)), return; end

    rR = (data.R_obs - yp(1)) / data.sigma_R;
    rC = (data.C_obs - yp(2)) / data.sigma_C;
    ll = -0.5 * (rR^2 + rC^2);

    if ~isfinite(ll), ll = -Inf; end
end

%% 自适应选择下一个 beta / Адаптивный выбор следующего beta
% 通过二分法寻找使权重变异系数 <= cov_tgt 的最大 beta_next / Бинарным поиском находится максимальное beta_next, при котором коэффициент вариации весов <= cov_tgt
function beta_next = find_next_beta(log_likes, beta_prev, cov_tgt)
    valid = isfinite(log_likes);
    if sum(valid) < 2
        beta_next = 1.0;
        return;
    end
    ll = log_likes(valid);

    % 如果直接跳到 beta=1 也满足 COV 约束，则直接到达 / Если переход сразу к beta=1 удовлетворяет ограничению COV, выполняется прямой переход
    if compute_cov(ll, 1.0 - beta_prev) <= cov_tgt
        beta_next = 1.0;
        return;
    end

    % 二分法搜索 / Бинарный поиск
    lo = beta_prev;
    hi = 1.0;
    for iter = 1:60
        mid = 0.5 * (lo + hi);
        if compute_cov(ll, mid - beta_prev) <= cov_tgt
            lo = mid;       % 还能继续增大 / Можно еще увеличить
        else
            hi = mid;       % 需要缩小 / Нужно уменьшить
        end
        if hi - lo < 1e-7, break; end
    end
    beta_next = min(max(lo, beta_prev + 1e-6), 1.0);
end

%% 权重变异系数计算 / Вычисление коэффициента вариации весов
function cv = compute_cov(ll_vec, db)
    if db <= 0, cv = 0; return; end
    lw = db * ll_vec;           % 对数权重 / Логарифмы весов
    lw = lw - max(lw);          % 减去最大值防止溢出 / Вычитание максимума для предотвращения переполнения
    w  = exp(lw);
    mu = mean(w);
    if mu <= 0, cv = Inf; return; end
    cv = std(w) / mu;
end

%% TMCMC 主算法 / Основной алгоритм TMCMC
function results = run_tmcmc(data, theta_f, cond, lb, ub, opts)

    N       = opts.N;               % 粒子数 / Число частиц
    n_p     = length(lb);           % 参数维度 / Размерность параметров
    cov_tgt = opts.COV_target;      % 权重 COV 目标 / Целевое значение COV весов
    nMH     = opts.N_MH;            % 每阶段 MCMC 步数 / Число шагов MCMC на этапе
    sc      = opts.scale;           % 提议协方差缩放因子 / Масштаб ковариации предложения

    % ---- 从均匀先验采样初始粒子 / Выборка начальных частиц из равномерного априора ----
    particles = bsxfun(@plus, lb, bsxfun(@times, ub - lb, rand(N, n_p)));

    % ---- 计算所有粒子的初始对数似然 / Вычисление начального логарифма правдоподобия для всех частиц ----
    log_L = zeros(N, 1);
    fprintf('  Инициализация правдоподобий частиц...');
    for i = 1:N
        log_L(i) = log_likelihood(particles(i,:), theta_f, data, cond);
    end
    fprintf(' завершено (допустимых частиц %d/%d)\n', sum(isfinite(log_L)), N);

    % ---- 初始化状态变量 / Инициализация переменных состояния ----
    beta_cur = 0;           % 当前温度参数 / Текущий температурный параметр
    betas    = 0;           % 记录所有 beta 值 / Запись всех значений beta
    log_evid = 0;           % 累积对数模型证据 / Накопленный логарифм свидетельства модели
    stage    = 0;           % 阶段计数 / Счетчик этапов
    acc_hist = [];          % 各阶段接受率 / Коэффициенты принятия на этапах

    % ================ 主循环 / Основной цикл ================
    while beta_cur < 1.0

        stage = stage + 1;
        if stage > opts.max_stages
            warning('Достигнут максимальный предел числа этапов %d; выполнение остановлено досрочно', opts.max_stages);
            break;
        end

        % (a) 自适应选择下一个 beta / Адаптивный выбор следующего beta
        beta_new = find_next_beta(log_L, beta_cur, cov_tgt);
        db       = beta_new - beta_cur;

        % (b) 计算重要性权重 / Вычисление весов важности
        log_w = db * log_L;

        % 数值稳定：对有限权重取最大值 / Численная устойчивость: берем максимум по конечным лог-весам
        finite_mask = isfinite(log_w);
        lw_max      = max(log_w(finite_mask));
        w           = exp(log_w - lw_max);
        w(~isfinite(w)) = 0;
        w_sum = sum(w);

        if w_sum <= 0 || ~isfinite(w_sum)
            warning('Сумма весов некорректна; выполнение остановлено досрочно');
            break;
        end

        % 归一化权重 / Нормализация весов
        w_n = w / w_sum;

        % (c) 更新对数模型证据 / Обновление логарифма свидетельства модели
        log_evid = log_evid + log(w_sum / N) + lw_max;

        % (d) 按权重重采样（多项式重采样） / Ресэмплинг по весам (мультиномиальный ресэмплинг)
        idx   = datasample(1:N, N, 'Weights', w_n);
        pts_r = particles(idx, :);      % 重采样粒子 / Ресэмплированные частицы
        lL_r  = log_L(idx);             % 对应似然 / Соответствующие правдоподобия

        % (e) 构造提议分布协方差 / Построение ковариации распределения предложения
        mu_w = w_n' * particles;                            % 加权均值 / Взвешенное среднее
        dif  = bsxfun(@minus, particles, mu_w);             % 中心化 / Центрирование
        Sig  = (bsxfun(@times, dif, w_n))' * dif;           % 加权协方差 / Взвешенная ковариация
        Sig_prop = sc^2 * Sig + 1e-10 * eye(n_p);           % 缩放 + 正则化 / Масштабирование + регуляризация

        % Cholesky 分解（若失败则退化为对角矩阵） / Разложение Холецкого (при неудаче используется диагональная матрица)
        [L_ch, flg] = chol(Sig_prop, 'lower');
        if flg ~= 0
            L_ch = diag(sqrt(max(diag(Sig_prop), 1e-12)));
        end

        % (f) MCMC 扰动（组件级 Metropolis-Hastings） / MCMC-возмущение (покомпонентный Metropolis-Hastings)
        pts_new = pts_r;
        lL_new  = lL_r;
        n_acc   = 0;

        for i = 1:N
            th_c = pts_r(i, :);             % 当前位置 / Текущее положение
            llc  = beta_new * lL_r(i);      % 当前加温似然 / Текущее подогретое правдоподобие

            for j = 1:nMH
                % 提议新样本：th_p = th_c + L_ch * z / Предлагаемый новый образец: th_p = th_c + L_ch * z
                z    = randn(n_p, 1);
                th_p = th_c + (L_ch * z)';

                % 检查是否在先验边界内 / Проверка нахождения в границах априора
                if any(th_p < lb) || any(th_p > ub)
                    continue;   % 越界则直接拒绝 / При выходе за границы предложение сразу отклоняется
                end

                % 计算提议样本的加温似然 / Вычисление подогретого правдоподобия для предложенного образца
                ll_raw = log_likelihood(th_p, theta_f, data, cond);
                llp    = beta_new * ll_raw;

                % Metropolis-Hastings 接受准则 / Критерий принятия Metropolis-Hastings
                if log(rand()) < llp - llc
                    th_c  = th_p;
                    llc   = llp;
                    n_acc = n_acc + 1;
                end
            end

            % 保存扰动后的粒子 / Сохранение возмущенной частицы
            pts_new(i, :) = th_c;
            lL_new(i)     = log_likelihood(th_c, theta_f, data, cond);
        end

        % ---- 统计与输出 / Статистика и вывод ----
        ar = n_acc / (N * nMH);
        acc_hist(end+1) = ar; %#ok<AGROW>

        fprintf('  Этап %2d: beta %.4f -> %.4f  |  коэффициент принятия MCMC = %.3f\n', ...
            stage, beta_cur, beta_new, ar);

        % 更新状态 / Обновление состояния
        particles = pts_new;
        log_L     = lL_new;
        beta_cur  = beta_new;
        betas(end+1) = beta_new; %#ok<AGROW>

        if beta_cur >= 1.0, break; end
    end

    % ================ 汇总结果 / Сводные результаты ================
    th_mean = mean(particles, 1);           % 后验均值 / Апостериорное среднее
    th_std  = std(particles, 0, 1);         % 后验标准差 / Апостериорное стандартное отклонение

    % MAP 估计（最大后验密度粒子） / MAP-оценка (частица с максимальной апостериорной плотностью)
    [~, bi] = max(log_L);
    th_map  = particles(bi, :);

    % 95% 置信区间 / 95%-й доверительный интервал
    ci95 = zeros(n_p, 2);
    for i = 1:n_p
        ci95(i, :) = quantile(particles(:, i), [0.025, 0.975]);
    end

    % 打包结果 / Упаковка результатов
    results.particles    = particles;
    results.log_L        = log_L;
    results.theta_mean   = th_mean;
    results.theta_std    = th_std;
    results.theta_map    = th_map;
    results.theta_ci95   = ci95;
    results.stage_betas  = betas;
    results.acc_hist     = acc_hist;
    results.log_evidence = log_evid;
    results.n_stages     = stage;
end

%% ---------- 绘图函数 / Функция построения графиков ----------
function plot_results(results, theta_true, lb, ub, param_labels, param_labels_latex)
    pts = results.particles;
    n_p = size(pts, 2);

    % 配色方案 / Цветовая схема
    col_hist    = [0.08, 0.38, 0.65];      % 直方图蓝色 / Синий для гистограммы
    col_prior   = [1.00, 0.75, 0.10];      % 先验黄色 / Желтый для априора
    col_kde     = [0.15, 0.75, 0.35];      % KDE 绿色 / Зеленый для KDE
    col_scatter = [0.10, 0.65, 0.60];      % 散点青色 / Бирюзовый для точек

    % ================================================================
    % 图1：Corner plot（联合后验） / Рисунок 1: corner plot (совместное апостериорное распределение)
    % ================================================================
    figure('Name', 'Совместное апостериорное распределение TMCMC (corner plot)', ...
           'Position', [60, 40, 1050, 1000]);

    ng = 45;        % 2D KDE 网格分辨率 / Разрешение сетки для двумерного KDE
    h_hist  = [];
    h_prior = [];
    h_kde   = [];

    for row = 1:n_p
        for col = 1:n_p
            ax = subplot(n_p, n_p, (row - 1) * n_p + col);

            if row == col
                % ---- 对角线：边缘后验直方图 / Диагональ: гистограмма маргинального апостериорного распределения ----
                h_h = histogram(pts(:, col), 50, 'Normalization', 'pdf', ...
                    'FaceColor', col_hist, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
                hold on;

                xlims = [lb(col), ub(col)];
                xlim(xlims);

                % 均匀先验 PDF / Плотность равномерного априорного распределения
                prior_pdf = 1 / (ub(col) - lb(col));
                xp_vec = linspace(lb(col), ub(col), 200);
                h_p = plot(xp_vec, prior_pdf * ones(1, 200), '-', ...
                    'Color', col_prior, 'LineWidth', 2.0);

                % KDE 平滑曲线 / Сглаженная кривая KDE
                x_eval = linspace(lb(col), ub(col), 200);
                try
                    [f_kde, xi_kde] = ksdensity(pts(:, col), x_eval);
                    h_k = plot(xi_kde, f_kde, '--', ...
                        'Color', col_kde, 'LineWidth', 2.0);
                catch
                    h_k = [];
                end

                % 保存图例句柄（仅第一个对角格） / Сохранение дескрипторов легенды (только для первой диагональной ячейки)
                if row == 1
                    h_hist  = h_h;
                    h_prior = h_p;
                    if ~isempty(h_k), h_kde = h_k; end
                end

                % 真值竖线 / Вертикальная линия истинного значения
                xline(theta_true(col), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
                set(ax, 'Box', 'on');
                ylabel('Плотность', 'FontSize', 8);

            elseif row > col
                % ---- 下三角：2D KDE 热力图 / Нижний треугольник: тепловая карта двумерного KDE ----
                xg = linspace(lb(col), ub(col), ng);
                yg = linspace(lb(row), ub(row), ng);
                [Xg, Yg] = meshgrid(xg, yg);
                grid_pts = [Xg(:), Yg(:)];

                try
                    f2d = ksdensity([pts(:, col), pts(:, row)], grid_pts);
                    F2d = reshape(f2d, ng, ng);
                catch
                    F2d = zeros(ng, ng);
                end

                imagesc(ax, xg, yg, F2d);
                colormap(ax, parula);
                set(ax, 'YDir', 'normal');
                hold on;

                % 真值十字标记 / Крестообразная отметка истинного значения
                plot(theta_true(col), theta_true(row), 'w+', ...
                    'MarkerSize', 9, 'LineWidth', 2.0);

                xlim([lb(col), ub(col)]);
                ylim([lb(row), ub(row)]);
                set(ax, 'Box', 'on');

            else
                % ---- 上三角：散点图 / Верхний треугольник: диаграмма рассеяния ----
                scatter(pts(:, col), pts(:, row), 3, 'filled', ...
                    'MarkerFaceColor', col_scatter, ...
                    'MarkerFaceAlpha', 0.18, ...
                    'MarkerEdgeAlpha', 0.18);
                hold on;

                plot(theta_true(col), theta_true(row), 'r+', ...
                    'MarkerSize', 9, 'LineWidth', 2.0);

                xlim([lb(col), ub(col)]);
                ylim([lb(row), ub(row)]);
                set(ax, 'Box', 'on');
            end

            % 坐标轴标签（仅边缘子图） / Подписи осей (только на крайних подграфиках)
            if row == n_p
                xlabel(param_labels_latex{col}, 'FontSize', 9, 'Interpreter', 'latex');
            end
            if col == 1 && row ~= col
                ylabel(param_labels_latex{row}, 'FontSize', 9, 'Interpreter', 'latex');
            end
            if col == 1 && row == 1
                ylabel(param_labels_latex{1}, 'FontSize', 9, 'Interpreter', 'latex');
            end
            set(ax, 'FontSize', 7, 'TickLength', [0.015 0.015]);
        end
    end

    % 统一图例 / Единая легенда
    leg_handles = [];
    leg_labels  = {};
    if ~isempty(h_prior), leg_handles(end+1) = h_prior; leg_labels{end+1} = 'Априорное распределение'; end
    if ~isempty(h_kde),   leg_handles(end+1) = h_kde;   leg_labels{end+1} = 'Апостериорное распределение'; end
    if ~isempty(h_hist),  leg_handles(end+1) = h_hist;  leg_labels{end+1} = 'TMCMC'; end

    if ~isempty(leg_handles)
        lgd = legend(leg_handles, leg_labels, ...
            'Orientation', 'horizontal', 'FontSize', 10, ...
            'Location', 'southoutside', 'Box', 'on');
        lgd.Position(2) = 0.01;
    end
    sgtitle('Совместные и маргинальные апостериорные распределения шести корректирующих параметров', 'FontSize', 12, 'FontWeight', 'bold');

    % ================================================================
    % 图2：单独边缘后验 / Рисунок 2: отдельные маргинальные апостериорные распределения
    % ================================================================
    figure('Name', 'Маргинальные апостериорные распределения TMCMC (детально)', ...
           'Position', [50, 50, 1200, 500]);

    nc = ceil(n_p / 2);
    for i = 1:n_p
        subplot(2, nc, i);

        histogram(pts(:, i), 60, 'Normalization', 'pdf', ...
            'FaceColor', col_hist, 'EdgeColor', 'none', 'FaceAlpha', 0.80);
        hold on;

        xline(theta_true(i), '-', 'LineWidth', 2.2, ...
            'Color', [0.1 0.7 0.2], 'DisplayName', 'Истинное значение');
        xline(results.theta_mean(i), '--', 'LineWidth', 1.8, ...
            'Color', [0.85 0.15 0.1], 'DisplayName', 'Апостериорное среднее');
        xline(results.theta_map(i), ':', 'LineWidth', 1.8, ...
            'Color', [0.7 0.1 0.8], 'DisplayName', 'MAP');

        x_eval = linspace(lb(i), ub(i), 200);
        try
            [fk, xk] = ksdensity(pts(:, i), x_eval);
            plot(xk, fk, '--', 'Color', col_kde, 'LineWidth', 1.5, ...
                'DisplayName', 'Апостериорное KDE');
        catch
        end

        xlim([lb(i), ub(i)]);
        xlabel(param_labels_latex{i}, 'FontSize', 11, 'Interpreter', 'latex');
        ylabel('Плотность', 'FontSize', 10);

        title(sprintf('%s  Истина=%.4f | Среднее=%.4f | MAP=%.4f', ...
            param_labels{i}, theta_true(i), ...
            results.theta_mean(i), results.theta_map(i)), ...
            'FontSize', 8.5, 'Interpreter', 'none');

        grid on; grid minor;
        if i == 1, legend('Location', 'best', 'FontSize', 8); end
    end
    sgtitle('Гистограммы маргинальных апостериорных распределений', 'FontSize', 12);

    % ================================================================
    % 图3：TMCMC 诊断（beta 演化 + 接受率） / Рисунок 3: диагностика TMCMC (эволюция beta + коэффициент принятия)
    % ================================================================
    figure('Name', 'Диагностика процесса TMCMC', ...
           'Position', [200, 50, 950, 380]);

    ns = results.n_stages;

    subplot(1, 2, 1);
    plot(0:ns, results.stage_betas, 'bo-', ...
        'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', 'b');
    xlabel('Номер этапа', 'FontSize', 11);
    ylabel('\beta', 'FontSize', 11);
    title('Эволюция температурного параметра \beta', 'FontSize', 12);
    ylim([0, 1.05]);
    grid on; grid minor;

    subplot(1, 2, 2);
    bar(1:ns, results.acc_hist, 0.6, 'FaceColor', [0.30, 0.68, 0.38]);
    xlabel('Номер этапа', 'FontSize', 11);
    ylabel('Коэффициент принятия MCMC', 'FontSize', 11);
    title('Коэффициент принятия на каждом этапе', 'FontSize', 12);
    ylim([0, 1]);
    grid on; grid minor;

    sgtitle('Диагностика процесса TMCMC', 'FontSize', 12);
end
