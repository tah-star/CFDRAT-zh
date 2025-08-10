function precond = init_precond(threshold_increasment_ratio, threshold_maxmin_ratio, uplimit_recorder, A_dp)
%INIT_PRECOND 初始化线性求解器的预条件器配置
%   为动量和压力方程设置预条件器结构。对动量方程使用自适应ILU
%   预条件，对压力泊松方程使用静态不完全Cholesky分解。
%
%   输入参数:
%       threshold_increasment_ratio - 更新的求解器迭代增长阈值
%       threshold_maxmin_ratio      - 更新的最大/最小迭代比率阈值
%       uplimit_recorder           - 强制更新前的最大步数
%       A_dp                       - 压力泊松矩阵(对称正定)
%
%   输出参数:
%       precond - 包含预条件器配置的结构:
%                 .precond_u  - U动量方程的自适应ILU
%                 .precond_v  - V动量方程的自适应ILU
%                 .precond_dp - 压力方程的静态ICT
%
%   预条件策略:
%   - 动量: 基于性能的更新的自适应ILU
%   - 压力: 静态不完全Cholesky分解(一次计算)
%   - 丢弃容差根据网格大小缩放以增强鲁棒性
%
%   See also SOLVE_LINEAR_SYSTEM, GET_AB_U, GET_AB_V.

global Ny Nx 

precond = struct();

%% 动量方程的自适应ILU预条件器
% ILU因子在仿真过程中基于求解器性能指标动态计算和更新

% 根据网格大小缩放丢弃容差以获得最佳性能
grid_size = Ny * Nx;
if grid_size >= 1e5
    droptol = 5e-6;                              % 大网格的严格容差
elseif grid_size >= 1e4 
    droptol = 5e-5;                              % 中等容差
else
    droptol = 5e-4;                              % 小网格的宽松容差
end

% 两个速度分量的通用ILU配置
setup_uv.type = 'ilutp';                         % 带阈值和选主元的ILU
setup_uv.droptol = droptol;                      % 分解的丢弃容差
setup_uv.udiag = true;                           % 替换零对角元素

% U动量预条件器配置
precond_u.type = 'u';
precond_u.setup = setup_uv;
precond_u.threshold_increasment_ratio = threshold_increasment_ratio;
precond_u.threshold_maxmin_ratio = threshold_maxmin_ratio;
precond_u.failure_trigger = 200;                 % 强制更新前的迭代限制
precond_u.uplimit_recorder = uplimit_recorder;
precond_u.L = [];                                % 仿真过程中计算的ILU因子
precond_u.U = [];
precond_u.recorder = [];                         % 性能历史

% V动量预条件器(相同配置)
precond_v.type = 'v';
precond_v.setup = setup_uv;
precond_v.threshold_increasment_ratio = threshold_increasment_ratio;
precond_v.threshold_maxmin_ratio = threshold_maxmin_ratio;
precond_v.failure_trigger = 200;
precond_v.uplimit_recorder = uplimit_recorder;
precond_v.L = [];
precond_v.U = [];
precond_v.recorder = [];

%% 压力方程的静态不完全Cholesky预条件器
% 预计算一次分解(压力矩阵为常数)
setup_p.type = "ict";
setup_p.droptol = 1e-8;                          % 压力求解的高精度
setup_p.diagcomp = 1e-6;                         % 稳定性的对角正则化

% 预计算Cholesky分解: A_dp ≈ L*L'
L = ichol(A_dp, setup_p);

% 压力预条件器结构
precond_dp.L = L;                                % 下三角因子
precond_dp.U = L';                               % 上三角(转置)
precond_dp.type = 'dp';

precond.precond_u = precond_u;
precond.precond_v = precond_v;
precond.precond_dp = precond_dp;

end
