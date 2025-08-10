function [u_new, v_new, dp_corr, precond_dp] =...
         piso_correct2(u_star, v_star, u_star_star, v_star_star, grid_u, grid_v, grid_p, A_dp, dp_corr_old, precond_dp, solver_opt)
%PISO_CORRECT2 PISO算法的第二次修正步骤
%   计算第二次压力修正并产生当前时间步的最终速度场。
%
%   输入参数:
%       u_star, v_star         - 预测器的中间速度 [m/s]
%       u_star_star, v_star_star - 第一次修正的速度场 [m/s]
%       grid_u/v/p             - 交错网格结构
%       A_dp                   - 压力泊松矩阵(常数)
%       dp_corr_old           - 前一次第二修正(初始猜测) [Pa]
%       precond_dp            - 压力方程预条件器
%       solver_opt            - 线性求解器配置
%
%   输出参数:
%       u_new, v_new  - 时间步的最终修正速度场 [m/s]
%       dp_corr       - 第二次压力修正 [Pa]
%       precond_dp    - 更新的预条件器(如果适用)
%
%   算法步骤:
%   1. 组装右端项: b ≈ ρ * ∇·(L(u**) - L(u*))
%   2. 求解: ∇²p'' = b  
%   3. 修正: u_new = u** - (Δt/ρ) * ∇p''
%
%   See also PISO_CORRECT1, GET_B_DP2, PRESSURE_CORRECTION.

% 基于修正速度和中间速度差异组装右端项
b_dp_corr = get_b_dp2(u_star, v_star, u_star_star, v_star_star, grid_p);

% 求解第二次压力修正
[dp_corr, precond_dp] = solve_leq(A_dp, b_dp_corr, dp_corr_old, precond_dp, solver_opt);

% 应用第二次修正获得最终速度场
[u_new, v_new] = pressure_correction(u_star_star, v_star_star, dp_corr, grid_u, grid_v);

end
