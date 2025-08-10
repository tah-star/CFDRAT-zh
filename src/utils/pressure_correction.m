function [u_new, v_new] = pressure_correction(u_old, v_old, dp, grid_u, grid_v)
%PRESSURE_CORRECTION 对速度场应用压力梯度修正
%   通过用压力梯度修正速度场来强制质量守恒，实现投影方法的核心步骤。
%
%   输入参数:
%       u_old, v_old - 未修正的速度场
%       dp           - 压力修正场
%       grid_u/v     - 速度网格结构
%
%   输出参数:
%       u_new, v_new - 压力修正后的速度场
%
%   数学操作:
%   u_new = u_old - (Δt/ρ)∇dp
%   
%   修正将速度场投影到无散度空间，同时保持边界条件的符合性。
%
%   算法:
%   1. 在交错网格位置计算压力梯度
%   2. 仅对内部流体节点应用修正
%   3. 重新强制所有物理边界条件
%
%   See also PISO_CORRECT1, PISO_CORRECT2, GET_B_DP1, GET_B_DP2.

global rho dt h slip_opt

%% 在速度位置计算压力梯度
% U速度节点的X梯度(垂直面)
grad_dp_x = (dp(:, 2:end) - dp(:, 1:end-1)) / h;
% V速度节点的Y梯度(水平面)
grad_dp_y = (dp(2:end, :) - dp(1:end-1, :)) / h;

%% 对内部节点应用压力修正
% 初始化修正场
u_new = u_old;
v_new = v_old;

% 修正公式: u_new = u_old - (dt/rho)*grad_p
u_new(2:end-1, 2:end-1) = u_old(2:end-1, 2:end-1) - (dt / rho) * grad_dp_x;
v_new(2:end-1, 2:end-1) = v_old(2:end-1, 2:end-1) - (dt / rho) * grad_dp_y;

%% 重新强制边界条件
% 压力修正仅更新内部节点；必须重新应用边界条件以保持物理一致性

% 固体边界: 零速度
u_new(grid_u.is_solid) = 0;
v_new(grid_v.is_solid) = 0;
u_new(grid_u.is_solid_boundary) = 0;
v_new(grid_v.is_solid_boundary) = 0;

% U速度的壁面边界条件(上下)
if strcmp(slip_opt, '无滑移')
   % 无滑移: u_ghost = -u_interior (在壁面强制u=0)
    u_new([1, end], :) = -u_new([2, end-1], :);
elseif strcmp(slip_opt, '滑移')
    % 滑移: u_ghost = u_interior (在壁面强制du/dy=0)
    u_new([1, end], :) = u_new([2, end-1], :);
else
    error('无效的壁面条件选项！');
end

% V速度边界条件
% 左边界(入口): 无穿透v=0
v_new(2:end-1, 1) = -v_new(2:end-1, 2);

% 右边界(出口): 零法向梯度
v_new(2:end-1, end) = v_new(2:end-1, end-1);

end
