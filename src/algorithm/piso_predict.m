function [u_star, v_star, precond_u, precond_v] = piso_predict(t, u, v, p, grid_u, grid_v, inlet_func, precond_u, precond_v, solver_opt)
%PISO_PREDICT PISO算法的动量预测步骤
%   求解动量方程以获得中间速度场(u*, v*)。
%   根据系统配置自动在串行和并行执行之间切换。
%
%   输入参数:
%       t          - 当前时间
%       u, v       - 前一时刻的速度场
%       p          - 前一时刻的压力场
%       grid_u/v   - 速度网格结构
%       inlet_func - 入口速度边界条件函数
%       precond_u/v- 预条件子
%       solver_opt - 线性求解器配置
%
%   输出参数:
%       u_star/v_star - 中间速度场
%       precond_u/v   - 更新的预条件器结构
%
%   算法:
%   1. 组装动量矩阵: Au* = b (排除压力梯度)
%   2. 使用预条件迭代求解器求解线性系统
%   3. 基于求解器性能更新预条件器
%
%
%   See also PISO_CORRECT, GET_AB_U, GET_AB_V, SOLVE_LEQ.

global dt mu rho 

% 根据并行池可用性选择执行策略
if ~isempty(gcp('nocreate'))
    % 并行执行: 异步组装矩阵
    f_Au = parfeval(@get_Ab_u, 2, t, u, v, p, inlet_func, dt, mu, rho, grid_u);
    f_Av = parfeval(@get_Ab_v, 2, u, v, p, dt, mu, rho, grid_v);
    
    % 准备就绪时获取组装的矩阵
    [A_u, b_u] = fetchOutputs(f_Au);
    [A_v, b_v] = fetchOutputs(f_Av);
    
    % 并行求解动量系统
    f_sol_u = parfeval(@solve_leq, 2, A_u, b_u, u, precond_u, solver_opt);
    f_sol_v = parfeval(@solve_leq, 2, A_v, b_v, v, precond_v, solver_opt);
    
    % 收集解和更新的预条件器
    [u_star, precond_u] = fetchOutputs(f_sol_u);
    [v_star, precond_v] = fetchOutputs(f_sol_v);
    
else
    % 串行执行: 顺序组装和求解
    [A_u, b_u] = get_Ab_u(t, u, v, p, inlet_func, dt, mu, rho, grid_u);
    [A_v, b_v] = get_Ab_v(u, v, p, dt, mu, rho, grid_v);
    
    % 顺序求解动量系统
    [u_star, precond_u] = solve_leq(A_u, b_u, u, precond_u, solver_opt);
    [v_star, precond_v] = solve_leq(A_v, b_v, v, precond_v, solver_opt);
    
end

end
