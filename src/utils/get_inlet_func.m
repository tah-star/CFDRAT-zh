function inlet_func = get_inlet_func(expr_str, ramp_opt, t_ramp, H, slip_opt)
%GET_INLET_FUNC 创建时间相关的入口速度函数句柄
%   将用户定义的速度表达式转换为鲁棒的矢量化函数句柄，
%   具有可选的时间渐变和无滑移壁面的自动边界条件强制。
%
%   输入参数:
%       expr_str  - 使用变量t, y, H的速度剖面字符串表达式
%                   示例: '0.15*4*y*(H-y)/H^2' (抛物线剖面)
%       ramp_opt  - 时间渐变类型: '线性上升'|'S型上升'|'无'
%       t_ramp    - 渐变持续时间，设为0表示无渐变
%       H         - 计算域高度
%       slip_opt  - 壁面边界条件: '无滑移'|'滑移'
%
%   输出参数:
%       inlet_func - 函数句柄@(t,y)返回速度
%                   其中t是时间，y是垂直坐标
%
%   功能特性:
%   - 用户表达式的自动矢量化
%   - 平滑的时间启动以避免冲击性流动启动
%   - 无滑移边界一致性的Tukey窗口强制
%   - 无效表达式的鲁棒错误处理
%
%   示例:
%       % 带线性渐变的抛物线剖面
%       inlet_func = get_inlet_func('0.1*4*y*(H-y)/H^2', 'linear', 0.5, 0.3, 'no-slip');
%       u_inlet = inlet_func(1.0, [0:0.01:0.3]'); % 在t=1s时求值
%
%   See also STR2FUNC, SMOOTHSTEP.

% 用于无滑移边界强制的Tukey窗口参数
tukey_alpha = 0.2;                               % 过渡区域分数 (20%)

% 将用户表达式转换为逐元素运算以实现向量兼容性
safe_expr_str = strrep(expr_str, '*', '.*');    % 启用逐元素乘法
safe_expr_str = strrep(safe_expr_str, '/', './'); % 启用逐元素除法
safe_expr_str = strrep(safe_expr_str, '^', '.^'); % 启用逐元素幂运算

% 将用户表达式解析为可执行的函数句柄
try
    profile_handle = str2func(['@(t,y,H)' safe_expr_str]);
catch ME
    error('无效的入口速度表达式: "%s".\n错误: %s\n请确保仅使用变量t, y和H.', safe_expr_str, ME.message);
end

% 配置时间渐变函数
if t_ramp <= 0
    ramp_factor_handle = @(t) 1.0;               % 无渐变
else
    switch ramp_opt
        case '线性上升' 
            ramp_factor_handle = @(t) min(1.0, t / t_ramp);        % 线性渐变
        case 'S型上升' 
            ramp_factor_handle = @(t) smoothstep_internal(t / t_ramp); % 平滑渐变
        otherwise
            ramp_factor_handle = @(t) 1.0;       % 默认：无渐变
    end
end

% 将空间剖面与时间渐变结合
% 逐元素乘法处理标量渐变与矢量剖面
base_inlet_func = @(t, y) profile_handle(t, y, H) .* ramp_factor_handle(t);

% 对无滑移壁面应用边界条件强制
if strcmp(slip_opt, '无滑移')
    % 在时间范围内评估用户函数在壁面边界处的值
    vel_y0 = max(abs(base_inlet_func(0:0.1:10, 0)));       % 底壁速度
    vel_yH = max(abs(base_inlet_func(0:0.1:10, H)));       % 顶壁速度
    
    % 检查边界条件一致性
    if abs(vel_y0) > 1e-3 || abs(vel_yH) > 1e-3
        % 用户剖面违反无滑移条件 - 应用修正窗口
        fprintf('\n[警告] 在"无滑移"情况下，入口剖面在壁面处非零。\n');
        fprintf('        求解器将应用平顶窗口以强制壁面处速度为零，\n');
        fprintf('        保留剖面的中心80%%部分。\n');
        fprintf('        用户定义的y=0处速度为%.4f，y=H处速度为%.4f。\n\n', vel_y0, vel_yH);
        
        % 应用Tukey窗口强制壁面处速度为零
        inlet_func = @(t, y) base_inlet_func(t, y) .* tukey_window_vectorized(y, H, tukey_alpha);
    else
        % 用户剖面已满足无滑移条件
        inlet_func = base_inlet_func;
    end
else
    % 对于滑移条件，直接使用用户剖面
    inlet_func = base_inlet_func;
end

end


function W = tukey_window_vectorized(y, H, alpha)
%TUKEY_WINDOW_VECTORIZED 用于边界强制的平顶余弦窗
%   应用Tukey窗口，从1（内部）平滑过渡到0（壁面）
%   以强制无滑移边界条件。
    y_norm = y / H;                              % 将坐标标准化到[0,1]
    W = ones(size(y_norm));                      % 初始化窗口值
    transition_length = alpha / 2;               % 过渡区域的半宽度
    
    % 在底壁附近应用余弦锥化 (y=0)
    left_mask = y_norm <= transition_length;
    if any(left_mask)
        W(left_mask) = 0.5 * (1 + cos(pi * (y_norm(left_mask) / transition_length - 1)));
    end
    
    % 在顶壁附近应用余弦锥化 (y=H)
    right_mask = y_norm >= (1 - transition_length);
    if any(right_mask)
        W(right_mask) = 0.5 * (1 + cos(pi * (y_norm(right_mask) - (1 - transition_length)) / transition_length));
    end
end


function val = smoothstep_internal(x)
%SMOOTHSTEP_INTERNAL 用于渐变的平滑Hermite插值
%   使用三次Hermite多项式提供从0到1的C2连续过渡。
    x_clamped = max(0, min(1, x));               % 将输入限制到[0,1]范围
    val = x_clamped.^2 .* (3.0 - 2.0 * x_clamped); % 三次Hermite多项式
end
