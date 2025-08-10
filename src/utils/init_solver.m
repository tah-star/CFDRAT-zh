% configure_solver.m
function solver_opt = init_solver(max_iter, residual_tol)
%INIT_SOLVER 配置迭代线性求解器参数
%   为PISO算法中使用的基于GMRES的线性系统求解器
%   设置收敛准则。
%
%   输入参数:
%       max_iter     - 失败前的最大求解器迭代数
%       residual_tol - 收敛的相对残差容差
%
%   输出参数:
%       solver_opt - 求解器配置结构
%
%   用法:
%       solver_opt = init_solver(1000, 1e-6);
%
%   See also SOLVE_LINEAR_SYSTEM, GMRES.

solver_opt.max_iter = max_iter;                  % 迭代限制
solver_opt.tol = residual_tol;                   % 收敛容差

end
