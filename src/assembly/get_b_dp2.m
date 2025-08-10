function b = get_b_dp2(u_star, v_star, u_star_star, v_star_star, grid_p)
%GET_B_DP2 组装第二次压力修正方程的右端项（PISO）
%   计算PISO算法中第二次压力修正的右端项，以更好地满足
%   动量方程残差。
%
%   输入参数:
%       u_star, v_star           - 预测器的中间速度 [m/s]
%       u_star_star, v_star_star - 第一次修正后的速度 [m/s]  
%       grid_p                   - 压力网格结构
%
%   输出参数:
%       b - 第二次修正的右端项向量: ∇²p'' = b
%
%   数学背景:
%   第二次修正通过计算以下内容来处理动量方程残差：
%   ∇²p'' ≈ ρ ∇·[L(u**) - L(u*)]
%   其中L(u) = -对流 + 扩散是空间动量算子。
%   此修正改善了超越第一次修正的压力-速度耦合。
%
%   See also GET_B_DP1, PISO_CORRECT2, GET_CONV, GET_DIFF, GET_DIV.

global rho  

%% 计算两个速度场的空间动量算子
% L(u) = -对流 + 扩散项来自动量方程
[conv_u_s, conv_v_s] = get_conv(u_star, v_star);
[diff_u_s, diff_v_s] = get_diff(u_star, v_star);
L_u_star = -conv_u_s + diff_u_s;
L_v_star = -conv_v_s + diff_v_s;

[conv_u_ss, conv_v_ss] = get_conv(u_star_star, v_star_star);
[diff_u_ss, diff_v_ss] = get_diff(u_star_star, v_star_star);
L_u_star_star = -conv_u_ss + diff_u_ss;
L_v_star_star = -conv_v_ss + diff_v_ss;

%% 计算动量算子差的散度
% ∇·[L(u**) - L(u*)] 表示动量残差的空间分布
delta_L_u = L_u_star_star - L_u_star;
delta_L_v = L_v_star_star - L_v_star;

div_delta_L = get_div(delta_L_u, delta_L_v);

%% 组装右端项向量
b_matrix = rho * div_delta_L;

% 固体节点的右端项为零
b_matrix(grid_p.is_solid) = 0;

% 负Laplacian矩阵约定的符号修正
b_matrix = -b_matrix;

% 转换为线性求解器的列向量
b_transposed = b_matrix';
b = b_transposed(:);

end

function [conv_u, conv_v] = get_conv(u, v)
%GET_CONV 计算压力节点处的对流项
%   使用中心差分计算-(u·∇)u和-(u·∇)v。

global  h
    % 将速度插值到压力节点位置（单元中心）
    u_on_p = (u(2:end-1, 2:end) + u(2:end-1, 1:end-1)) / 2;
    v_on_p = (v(2:end, 2:end-1) + v(1:end-1, 2:end-1)) / 2;

    % 压力节点处的速度梯度（2阶中心差分）
    grad_u_x_on_p = (u(2:end-1, 2:end) - u(2:end-1, 1:end-1)) / h;
    grad_u_y = (u(3:end, :) - u(1:end-2, :)) / (2*h);
    grad_u_y_on_p = (grad_u_y(:, 1:end-1) + grad_u_y(:, 2:end)) / 2;

    grad_v_y_on_p = (v(2:end, 2:end-1) - v(1:end-1, 2:end-1)) / h;
    grad_v_x = (v(:, 3:end) - v(:, 1:end-2)) / (2*h);
    grad_v_x_on_p = (grad_v_x(1:end-1, :) + grad_v_x(2:end, :)) / 2;

    % 对流项: u·∇u 和 u·∇v
    conv_u = u_on_p .* grad_u_x_on_p + v_on_p .* grad_u_y_on_p;
    conv_v = u_on_p .* grad_v_x_on_p + v_on_p .* grad_v_y_on_p;
end


function [diff_u, diff_v] = get_diff(u, v)
%GET_DIFF 计算压力节点处的粘性扩散项
%   使用5点Laplacian模板计算μ∇²u和μ∇²v。

global mu h slip_opt
    % 将速度插值到压力节点位置
    u_on_p = (u(2:end-1, 2:end) + u(2:end-1, 1:end-1)) / 2;
    v_on_p = (v(2:end, 2:end-1) + v(1:end-1, 2:end-1)) / 2;

    % 使用适当边界条件构建填充的U速度矩阵
    u_inlet_values = u(2:end-1, 1);
    u_p_ghost_left = 2 * u_inlet_values - u_on_p(:, 1);  % 入口 (Dirichlet)
    u_p_ghost_right = u_on_p(:, end);                    % 出口 (Neumann)
    u_on_p_padded_lr = [u_p_ghost_left, u_on_p, u_p_ghost_right];

    % U速度的壁面边界条件
    if strcmp(slip_opt, '无滑移')
        u_p_ghost_bottom = -u_on_p_padded_lr(1, :);      % 无滑移: u_ghost = -u_interior
        u_p_ghost_top = -u_on_p_padded_lr(end, :);
    elseif strcmp(slip_opt, '滑移')
        u_p_ghost_bottom = u_on_p_padded_lr(1, :);       % 滑移: u_ghost = u_interior
        u_p_ghost_top = u_on_p_padded_lr(end, :);
    else
        error('无效的滑移选项！');
    end

    u_on_p_padded = [u_p_ghost_bottom; u_on_p_padded_lr; u_p_ghost_top];

    % 使用边界条件构建填充的V速度矩阵
    % 壁面: 无穿透条件 (v=0 → v_ghost = -v_interior)
    v_p_ghost_bottom = -v_on_p(1, :);
    v_p_ghost_top = -v_on_p(end, :);
    v_on_p_padded_tb = [v_p_ghost_bottom; v_on_p; v_p_ghost_top];

    % V速度的入口/出口条件
    v_p_ghost_left = -v_on_p_padded_tb(:, 1);           % 入口 (Dirichlet v=0)
    v_p_ghost_right = v_on_p_padded_tb(:, end);         % 出口 (Neumann)
    
    v_on_p_padded = [v_p_ghost_left, v_on_p_padded_tb, v_p_ghost_right];

    % 应用5点Laplacian模板计算扩散项
    laplacian_u = (u_on_p_padded(2:end-1, 1:end-2) + u_on_p_padded(2:end-1, 3:end) + ...
                   u_on_p_padded(1:end-2, 2:end-1) + u_on_p_padded(3:end, 2:end-1) - ...
                   4 * u_on_p_padded(2:end-1, 2:end-1)) / h^2;
               
    laplacian_v = (v_on_p_padded(2:end-1, 1:end-2) + v_on_p_padded(2:end-1, 3:end) + ...
                   v_on_p_padded(1:end-2, 2:end-1) + v_on_p_padded(3:end, 2:end-1) - ...
                   4 * v_on_p_padded(2:end-1, 2:end-1)) / h^2;

    diff_u = mu * laplacian_u;
    diff_v = mu * laplacian_v;
end


function div_L = get_div(L_u, L_v)
%GET_DIV 计算压力节点处矢量场的散度
%   使用2阶差分和适当的边界处理计算∇·L。

    global h 
    [Ny, Nx] = size(L_u);
    grad_L_u_x = zeros(Ny, Nx);
    grad_L_v_y = zeros(Ny, Nx);

    % 内部节点：2阶中心差分
    grad_L_u_x(:, 2:end-1) = (L_u(:, 3:end) - L_u(:, 1:end-2)) / (2 * h);
    grad_L_v_y(2:end-1, :) = (L_v(3:end, :) - L_v(1:end-2, :)) / (2 * h);
 
    % 边界节点：2阶单侧差分
    % 假设边界处散度的法向梯度为零
    
    % 左/右边界 (x方向)
    grad_L_u_x(:, 1)   = (-3*L_u(:, 1) + 4*L_u(:, 2) - L_u(:, 3)) / (2 * h);
    grad_L_u_x(:, end) = (3*L_u(:, end) - 4*L_u(:, end-1) + L_u(:, end-2)) / (2 * h);
    
    % 下/上边界 (y方向)
    grad_L_v_y(1, :)   = (-3*L_v(1, :) + 4*L_v(2, :) - L_v(3, :)) / (2 * h);
    grad_L_v_y(end, :) = (3*L_v(end, :) - 4*L_v(end-1, :) + L_v(end-2, :)) / (2 * h);

    div_L = grad_L_u_x + grad_L_v_y;
end
