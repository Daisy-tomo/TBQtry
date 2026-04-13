function [y, aux] = engine_forward_12(theta, cfg)

if nargin < 2 || ~isfield(cfg, 'cond')
    error('engine_forward_12 需要 cfg.cond 工况输入。');
end

cond = cfg.cond;

% 12 个修正参数
eta_k      = theta(1);
eta_t      = theta(2);
eta_m      = theta(3);
eta_v      = theta(4);
eta_tv     = theta(5);
eta_c1     = theta(6);
eta_c2     = theta(7);
sigma_cc   = theta(8);
sigma_kan  = theta(9);
sigma_kask = theta(10);
sigma_ks   = theta(11);
eta_T      = theta(12);   

% ========== 固定参数 ==========
lambda = 1.03;

% ========== 工况参数 ==========
T_H      = cond.T_H;
M_flight = cond.M_flight;
m        = cond.m;
pi_k     = cond.pi_k;
T_g      = cond.T_g;

y = [NaN, NaN];
aux = struct();

try
    % (1) 基本常数与飞行速度
    k_air = 1.4;
    R_air = 287.3;
    a = sqrt(k_air * R_air * T_H);
    if ~isfinite(a) || a <= 0
        return;
    end
    V_flight = a * M_flight;

    % (2) 分段参数
    kT = piecewise_kT(T_g);
    RT = piecewise_RT(T_g);
    d = delta_cooling(T_g);

    % (3) 入口总压比与压气机入口温度
    inner = 1 + V_flight^2 / (2 * (k_air / (k_air - 1)) * R_air * T_H);
    if inner <= 0
        return;
    end
    tau_v = inner^(k_air / (k_air - 1));
    T_B = T_H * inner^k_air;
    if ~isfinite(T_B) || T_B <= 0
        return;
    end

    % (4) 压气机出口温度
    pi_k_ratio = pi_k^((k_air - 1) / k_air);
    if ~isfinite(pi_k_ratio) || pi_k_ratio < 1
        return;
    end
    T_k = T_B * (1 + (pi_k_ratio - 1) / eta_k);
    if ~isfinite(T_k) || T_k <= 0
        return;
    end

    % (5) 相对耗油量 — 燃烧室能量平衡 + 查表插值
    %     g_T = (CpT_g - CpT_k) / (H_u * eta_T - iT_g + CpT_k)
    %     表中和 H_u 统一使用 кДж/кг，结果无量纲
    H_u = 43000;                      % 燃料低热值, кДж/кг
    CpT_k  = interp_CpT(T_k);            % 压气机出口焓, кДж/кг
    CpT_g  = interp_CpT(T_g);            % 涡轮前燃气焓, кДж/кг
    iT_g   = interp_iT(T_g);             % 燃料焓, кДж/кг

    denom_gT = H_u * eta_T - iT_g + CpT_k;
    if abs(denom_gT) < 1e-6
        return;
    end
    g_T = (CpT_g - CpT_k) / denom_gT;
    if ~isfinite(g_T) || g_T <= 0
        return;
    end

   
    % (6) 进气总压恢复系数
    sigma_bx = sigma_cc * sigma_kan;

    % (7) 单位自由能
    exp_T = (kT - 1) / kT;
    expansion_pr_denom = tau_v * sigma_bx * pi_k * sigma_kask * sigma_ks;
    if expansion_pr_denom <= 0
        return;
    end

    expansion_term = (1.0 / expansion_pr_denom)^exp_T;
    if ~isfinite(expansion_term)
        return;
    end

    term1 = (kT / (kT - 1)) * RT * T_g * (1 - expansion_term);
    compress_work = (k_air / (k_air - 1)) * R_air * T_B ...
                  * (pi_k^((k_air - 1) / k_air) - 1);
    denom2 = (1 + g_T) * eta_k *  eta_t * eta_m * (1 - d);
    if abs(denom2) < 1e-10
        return;
    end

    term2 = compress_work / denom2;
    L_sv = lambda * (term1 - term2);
    if ~isfinite(L_sv) || L_sv <= 0
        return;
    end

    % (8) 最优自由能分配系数
    V2_term = m * V_flight^2;
    num_xpc = 1 + V2_term / (2 * L_sv * eta_tv * eta_v * eta_c2);
    den_xpc = 1 + (m * eta_tv * eta_v * eta_c2) / (eta_c1 * lambda);
    if abs(den_xpc) < 1e-10
        return;
    end

    x_pc = num_xpc / den_xpc;
    if ~isfinite(x_pc) || x_pc <= 0
        return;
    end

    % (9) 比推力
    inner_sq1 = 2 * eta_c1 * lambda * x_pc * L_sv;
    if inner_sq1 < 0
        return;
    end
    V_j1 = (1 + g_T) * sqrt(inner_sq1) - V_flight;

    inner_sq2 = 2 * (1 - x_pc) / m * L_sv * eta_tv * eta_v * eta_c2 ...
              + V_flight^2;
    if inner_sq2 < 0
        return;
    end
    V_j2 = sqrt(inner_sq2) - V_flight;

    R_ud = (1 / (1 + m)) * V_j1 + (m / (1 + m)) * V_j2;
    if ~isfinite(R_ud) || R_ud <= 0
        return;
    end

    % (10) 比油耗
    denom_C = R_ud * (1 + m);
    if abs(denom_C) < 1e-10
        return;
    end
    C_ud = 3600 * g_T  / denom_C;
    if ~isfinite(C_ud) || C_ud <= 0
        return;
    end

    y = [R_ud, C_ud];

    % 中间变量输出
    aux.T_B = T_B;
    aux.T_k = T_k;
    aux.tau_v = tau_v;
    aux.g_T = g_T;
    aux.CpT_k = CpT_k;
    aux.CpT_g = CpT_g;
    aux.iT_g = iT_g;
    aux.lambda_heat = lambda_heat;
    aux.sigma_bx = sigma_bx;
    aux.L_sv = L_sv;
    aux.x_pc = x_pc;
    aux.kT = kT;
    aux.RT = RT;
    aux.delta = d;
    aux.V_flight = V_flight;
    aux.lambda = lambda;
catch ME
    warning('engine_forward_demo_12p 捕获到异常: %s', ME.message);
end
end

% 热力学查表插值
function val = interp_CpT(T)
% CpT(T) — 线性插值; 超范围 clamp 到边界
persistent T_tab CpT_tab
if isempty(T_tab)
    [T_tab, CpT_tab, ~] = thermo_table_data();
end
T_c = max(T_tab(1), min(T_tab(end), T));
val = interp1(T_tab, CpT_tab, T_c, 'linear');
end

function val = interp_iT(T)
persistent T_tab iT_tab
if isempty(T_tab)
    [T_tab, ~, iT_tab] = thermo_table_data();
end
T_c = max(T_tab(1), min(T_tab(end), T));
val = interp1(T_tab, iT_tab, T_c, 'linear');
end

function [T_tab, CpT_tab, iT_tab] = thermo_table_data()

%  T [K]    CpT* [кДж/кг]   iT* [кДж/кг]
data = [ ...
   280     280.13     426.83
   290     290.18     446.26
   300     300.23     472.73
   310     310.28     485.50
   320     320.32     505.58
   330     330.41     525.89
   340     340.50     546.41
   350     350.50     567.17
   360     360.60     588.19
   370     370.69     609.71
   380     380.78     631.06
   390     390.91     653.04
   400     401.40     675.02
   410     411.13     697.50
   420     421.31     720.28
   430     431.52     743.02
   440     441.66     766.17
   450     451.87     789.30
   460     462.09     812.89
   470     472.30     836.46
   480     482.60     860.16
   490     492.86     884.02
   500     503.16     908.14
   510     513.46     932.34
   520     523.82     956.72
   530     534.10     981.40
   540     544.48    1006.20
   550     554.87    1031.20
   560     565.29    1056.30
   570     575.67    1081.50
   580     586.18    1107.00
   590     596.60    1132.50
   600     607.10    1158.40
   610     617.60    1184.20
   620     628.20    1210.30
   630     638.80    1236.60
   640     649.30    1263.10
   650     659.90    1289.80
   660     670.50    1316.60
   670     681.30    1343.50
   680     692.00    1370.70
   690     702.70    1398.10
   700     713.40    1425.70
   710     724.20    1453.20
   720     735.00    1480.80
   730     745.90    1508.60
   740     756.70    1536.70
   750     767.50    1564.90
   760     778.40    1592.90
   770     789.40    1621.00
   780     800.29    1649.30
   790     811.60    1677.80
   800     822.20    1706.50
   810     833.20    1735.14
   820     844.17    1763.99
   830     855.18    1793.04
   840     866.20    1822.20
   850     877.20    1851.60
   860     888.38    1881.10
   870     899.47    1910.90
   880     910.60    1940.70
   890     921.80    1970.80
   900     932.99    2001.08
   910     944.15    2031.40
   920     955.37    2061.90
   930     966.59    2092.60
   940     977.85    2123.50
   950     989.11    2154.50
   960    1000.42    2185.50
   970    1011.72    2216.60
   980    1023.10    2248.00
   990    1034.45    2279.50
  1000    1045.86    2311.19
  1010    1057.27    2342.63
  1020    1068.70    2374.30
  1030    1080.17    2406.10
  1040    1091.66    2438.10
  1050    1103.16    2470.33
  1060    1114.67    2502.41
  1070    1126.20    2534.50
  1080    1137.85    2567.01
  1090    1149.40    2599.58
  1100    1161.00    2632.28
  1110    1172.68    2664.86
  1120    1184.30    2697.60
  1130    1195.98    2730.46
  1140    1207.68    2763.49
  1150    1219.38    2796.69
  1160    1231.10    2830.15
  1170    1242.80    2863.77
  1180    1254.57    2897.50
  1190    1266.30    2931.40
  1200    1278.10    2965.51
  1210    1289.76    2999.59
  1220    1301.44    3033.79
  1230    1313.16    3068.21
  1240    1324.90    3102.75
  1250    1336.60    3137.46
  1260    1348.40    3172.04
  1270    1360.25    3206.76
  1280    1372.10    3241.60
  1290    1383.95    3276.67
  1300    1395.80    3311.80
  1310    1407.70    3346.89
  1320    1419.62    3382.10
  1330    1431.55    3417.47
  1340    1443.39    3452.98
  1350    1455.50    3488.65
  1360    1467.40    3524.20
  1370    1479.36    3559.99
  1380    1491.30    3595.87
  1390    1503.30    3631.92
  1400    1515.30    3668.09
  1410    1527.26    3704.23
  1420    1539.25    3740.53
  1430    1551.20    3776.95
  1440    1563.19    3813.45
  1450    1575.28    3850.22
  1460    1587.17    3886.81
  1470    1599.10    3923.53
  1480    1611.00    3960.80
  1490    1622.98    3997.43
  1500    1634.98    4034.52
  1510    1647.25    4071.29
  1520    1659.50    4108.21
  1530    1671.83    4145.22
  1540    1684.10    4182.40
  1550    1696.49    4219.71
  1560    1708.67    4256.93
  1570    1720.80    4294.32
  1580    1733.00    4331.83
  1590    1745.22    4369.43
  1600    1757.45    4407.19
  1610    1769.67    4444.92
  1620    1781.90    4482.76
  1630    1794.10    4520.74
  1640    1806.40    4558.84
  1650    1818.66    4597.06
  1660    1830.85    4635.29
  1670    1843.10    4673.68
  1680    1855.40    4712.16
  1690    1867.70    4750.76
  1700    1880.00    4789.49
  1710    1892.30    4827.79
  1720    1904.62    4866.71
];

T_tab   = data(:,1);
CpT_tab = data(:,2);
iT_tab  = data(:,3);
end

function kT = piecewise_kT(T_g)
if T_g > 800 && T_g <= 1400
    kT = 1.33;
elseif T_g > 1400 && T_g <= 1600
    kT = 1.30;
elseif T_g > 1600
    kT = 1.25;
else
    kT = 1.33;
end
end

function RT = piecewise_RT(T_g)
if T_g > 800 && T_g <= 1400
    RT = 287.6;
elseif T_g > 1400 && T_g <= 1600
    RT = 288.0;
elseif T_g > 1600
    RT = 288.6;
else
    RT = 287.6;
end
end

function d = delta_cooling(T_g)
d = 0.02 + (T_g - 1200) / 100 * 0.02;
d = max(0.0, min(d, 0.15));
end