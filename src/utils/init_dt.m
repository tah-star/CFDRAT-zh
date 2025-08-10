function dt = init_dt(speed_option, U_max)
%INIT_DT 基于CFL条件计算稳定时间步长
%   使用Courant-Friedrichs-Lewy(CFL)准则确定数值稳定性的
%   适当时间步长。将结果舍入为便于一致仿真输出的方便值。
%
%   输入参数:
%       speed_option - 速度偏好: 'fast'|'medium'
%       U_max        - 计算域中的最大速度
%
%   输出参数:
%       dt           - 稳定时间步长
%
%   CFL数:
%   - 'fast': CFL = 1.5 (较大时间步长，更快仿真)
%   - 'medium': CFL = 0.75 (较小时间步长，更保守)
%
%   See also GET_U_MAX, SOLVE_PDE.

global  h

%% 基于用户偏好的CFL数选择
CFL_FAST = 1.5;                                  % 激进的时间步进
CFL_MEDIUM = 0.75;                               % 保守的时间步进

switch lower(speed_option)
    case '快'
        cfl = CFL_FAST;
    case '中'
        cfl = CFL_MEDIUM;
    otherwise
        warning('INIT_DT:UnknownOption', ...
                '未知速度选项"%s"。使用默认"medium"(CFL=%.1f)。', ...
                speed_option, CFL_MEDIUM);
        cfl = CFL_MEDIUM;
end

% CFL约束: dt <= cfl * h / U_max
dt_max = cfl * h / U_max;

%% 向下舍入为可重现仿真的便利值
% 舍入为格式: 1×10^n, 2×10^n, 或 5×10^n (例如, 0.001, 0.002, 0.005)
power_of_10 = 10^floor(log10(dt_max));           % 找到适当的10的幂
first_digit = floor(dt_max / power_of_10);       % 提取首位数字

if first_digit >= 5
    dt = 5 * power_of_10;
elseif first_digit >= 2
    dt = 2 * power_of_10;
else
    dt = 1 * power_of_10;
end

end
