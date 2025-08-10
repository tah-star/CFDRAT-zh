function [u_star_star, v_star_star, dp, precond_dp] = ...
         piso_correct1(u_star, v_star, grid_u, grid_v, grid_p, A_dp, dp_old, precond_dp, solver_opt)
%PISO_CORRECT1 PISO算法的第一次修正步骤
%   计算第一次压力修正并修正中间速度场以更好地满足连续性方程。
%
%   输入参数:
%       u_star, v_star - 预测步骤的中间速度场 [m/s]
%       grid_u/v/p     - 交错网格结构
%       A_dp           - 压力泊松矩阵(常数)
%       dp_old         - 前一次压力修正(初始猜测) [Pa]
%       precond_dp     - 压力方程预条件器
%       solver_opt     - 线性求解器配置
%
%   输出参数:
%       u_star_star, v_star_star - 第一次修正的速度场 [m/s]
%       dp             - 第一次压力修正 [Pa]
%       precond_dp     - 更新的预条件器(如果适用)
%
%   算法步骤:
%   1. 组装右端项: b = (ρ/Δt) * ∇·u*
%   2. 求解: ∇²p' = b
%   3. 修正: u** = u* - (Δt/ρ) * ∇p'
%
%   See also PISO_PREDICT, PISO_CORRECT2, GET_B_DP1, PRESSURE_CORRECTION.

% 从速度散度组装压力泊松方程的右端项
b_dp = get_b_dp1(u_star, v_star, grid_p);

% 求解第一次压力修正
[dp, precond_dp] = solve_leq(A_dp, b_dp, dp_old, precond_dp, solver_opt);

% 将压力修正应用到速度场
[u_star_star, v_star_star] = pressure_correction(u_star, v_star, dp, grid_u, grid_v);

end
