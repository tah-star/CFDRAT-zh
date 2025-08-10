function [u_all, v_all, p_all, t_solve] = solve_pde(t_simu, T_record, u_inlet_func, grid_u, grid_v, grid_p, A_dp, precond, solver_opt)
%SOLVE_PDE 2D不可压缩Navier-Stokes方程的时间步进循环
%   在交错网格上使用PISO算法执行主要的CFD仿真。处理动量预测和
%   两阶段压力修正以实现精确的压力-速度耦合。
%
%   输入参数:
%       t_simu        - 总仿真时间
%       T_record      - 数据记录间隔
%       u_inlet_func  - 入口速度函数句柄 @(t,y)
%       grid_u        - 带掩码和索引的U速度网格结构
%       grid_v        - 带掩码和索引的V速度网格结构
%       grid_p        - 带边界分类的压力网格结构
%       A_dp          - 压力泊松矩阵(稀疏，常数)
%       precond       - 预条件器配置结构
%       solver_opt    - 线性求解器参数(容差，最大迭代数)
%
%   输出参数:
%       u_all, v_all  - 速度场时间序列 [m/s]
%       p_all         - 压力场时间序列 [Pa]
%       t_solve       - 总计算时间 [s]
%
%   算法(PISO):
%   1. 动量预测: 显式求解动量方程
%   2. 第一次压力修正: 强制连续性方程
%   3. 第二次压力修正: 提高动量方程精度
%   4. 更新压力场并记录数据
%
%   注意:
%   - 使用全局变量(dt, Nx, Ny)存储网格参数
%   - 动量方程采用自适应ILU预条件
%   - 进度更新打印到命令窗口
%
%   See also PISO_PREDICT, PISO_CORRECT1, PISO_CORRECT2.

% 访问全局网格和时间步长参数
global dt Nx Ny simulation_control

% 提取每个场的预条件器配置
precond_u = precond.precond_u;                    % U动量预条件器
precond_v = precond.precond_v;                    % V动量预条件器
precond_dp = precond.precond_dp;                  % 压力修正预条件器

% 计算仿真时间参数
total_steps = ceil(t_simu / dt);                  % 执行的总时间步数
record_interval_steps = max(1, round(T_record / dt));  % 数据记录间隔步数
num_records = floor(total_steps / record_interval_steps) + 1;  % 总数据帧数

% 在交错网格上初始化场变量
u = zeros(Ny+2, Nx+1);                           % U速度(垂直面，包含虚拟层)
v = zeros(Ny+1, Nx+2);                           % V速度(水平面，包含虚拟层)
p = zeros(Ny, Nx);                               % 压力(单元中心)
dp1 = zeros(Ny, Nx);                             % 第一次压力修正
dp2 = zeros(Ny, Nx);                             % 第二次压力修正

% 为时间序列数据预分配存储数组
u_all = zeros(Ny+2, Nx+1, num_records);
v_all = zeros(Ny+1, Nx+2, num_records);
p_all = zeros(Ny, Nx, num_records);

% 存储初始条件(静止流动)
record_count = 1;
u_all(:,:,record_count) = u;
v_all(:,:,record_count) = v;
p_all(:,:,record_count) = p;
record_count = record_count + 1;

fprintf('--- 仿真开始 ---\n');
fprintf('总时间: %.2fs | 时间步长(dt): %.4fs | 总步数: %d\n', t_simu, dt, total_steps);

% 主PISO时间步进循环
tic;
for step = 1:total_steps
    pause(0.001);
    if strcmp(simulation_control, 'stopped')
        break;  % 跳出时间循环
    end

    t = step * dt;

    % PISO步骤1: 动量预测(求解不含压力梯度的动量方程)
    [u_star1, v_star1, precond_u, precond_v] = piso_predict(t, u, v, p, grid_u, grid_v, u_inlet_func, precond_u, precond_v, solver_opt);

    % PISO步骤2: 第一次压力修正(强制连续性)
    [u_star2, v_star2, dp1, precond_dp] = piso_correct1(u_star1, v_star1, grid_u, grid_v, grid_p, A_dp, dp1, precond_dp, solver_opt);

    % PISO步骤3: 第二次压力修正(提高动量精度)
    [u, v, dp2, precond_dp] = piso_correct2(u_star1, v_star1, u_star2, v_star2, grid_u, grid_v, grid_p, A_dp, dp2, precond_dp, solver_opt);

    % 用修正量更新压力场
    p = p + dp1 + dp2;

    % 按指定间隔记录仿真数据
    if mod(step, record_interval_steps) == 0 && record_count <= num_records
        u_all(:,:, record_count) = u;
        v_all(:,:, record_count) = v;
        p_all(:,:, record_count) = p;
        record_count = record_count + 1;

         % 向用户显示进度
        fprintf('时间: %.3fs (%.1f%%), 时间步: %d/%d, 已用时: %.2f 秒\n', ...
            t, 100*t/t_simu, step, total_steps, toc);
    end
end

t_solve = toc;                                   % 记录总计算时间

if strcmp(simulation_control, 'stopped') && step < total_steps
    fprintf('--- 仿真被用户停止，已完成 %d/%d 步，用时 %.2f 秒 ---\n', step, total_steps, t_solve);
else
    fprintf('--- 仿真完成，用时 %.2f 秒 ---\n', t_solve);
end

end
