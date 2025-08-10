function [solution_matrix, preconditioner] = solve_leq(A, b, x0_matrix, preconditioner, solver_setup)
%SOLVE_LEQ 使用自适应预条件求解线性系统
%   带智能预条件器管理的迭代线性求解器包装器。
%   动量方程使用BiCGSTAB，压力方程使用PCG。
%
%   输入参数:
%       A              - 系统矩阵
%       b              - 右端项向量
%       x0_matrix      - 初始猜测(2D矩阵形式)
%       preconditioner - 带自适应参数的预条件器结构
%       solver_setup   - 求解器配置(容差，最大迭代次数)
%
%   输出参数:
%       solution_matrix - 2D矩阵形式的解
%       preconditioner  - 更新的预条件器结构
%
%   自适应策略:
%   - 动量(u,v): 基于性能更新的ILU预条件器
%   - 压力(dp): 固定Cholesky预条件器(一次计算)
%
%   预条件器更新条件:
%   - 首次使用或超过使用限制
%   - 求解器迭代次数显著增加
%   - 性能下降超过阈值

%% 设置并向量化初始猜测
type = preconditioner.type;
tol = solver_setup.tol;
max_iter = solver_setup.max_iter;

[Ny, Nx] = size(x0_matrix);
x0_vec = reshape(x0_matrix', [], 1); 

%% 根据方程类型求解
if strcmp(type, 'u') || strcmp(type, 'v')
    % 动量方程: 带BiCGSTAB的自适应ILU预条件
    
    % 检查预条件器是否需要更新
    should_renew = should_renew_preconditioner(preconditioner);
    
    % 必要时更新ILU因子
    if should_renew
        [L, U] = ilu(A, preconditioner.setup);
        preconditioner.L = L;
        preconditioner.U = U;
    else
        L = preconditioner.L;
        U = preconditioner.U;
    end

    % 使用BiCGSTAB求解
    [solution_vec, flag, ~, iter, ~] = bicgstab(A, b, tol, max_iter, L, U, x0_vec);

    % 检查收敛性并在有问题时警告
    if flag ~= 0 || iter > 0.8*max_iter
        fprintf('警告: %s求解器 - 标志:%d, 迭代:%d/%d (%.1f%%)\n', type, flag, iter, max_iter, 100*iter/max_iter);
    end

    % 更新性能历史
    if should_renew
        preconditioner.recorder = iter;           % 更新后重置历史
    else
        preconditioner.recorder(end + 1) = iter; % 追加到历史
    end
    
elseif strcmp(type, 'dp')
    % 压力方程: 带PCG的固定Cholesky预条件
    L = preconditioner.L;
    U = preconditioner.U;                        % U = L' for Cholesky
    
    [solution_vec, flag, ~, iter, ~] = pcg(A, b, tol, max_iter, L, U, x0_vec);

    % 检查收敛性
    if flag ~= 0 || iter > 0.8*max_iter
        fprintf('警告: %s求解器 - 标志:%d, 迭代:%d/%d (%.1f%%)\n', type, flag, iter, max_iter, 100*iter/max_iter);
    end

else
    error('未知的求解器类型: %s。必须是 ''u'', ''v'', 或 ''dp''。', type);
end

%% 将解转换回矩阵形式
solution_matrix = reshape(solution_vec, Nx, Ny)';

end


function renew_flag = should_renew_preconditioner(preconditioner)
%SHOULD_RENEW_PRECONDITIONER 决定是否更新ILU预条件器
%   分析求解器性能历史以确定预条件器更新是否会提高效率。

recorder = preconditioner.recorder;
failure_trigger = preconditioner.failure_trigger;
uplimit_recorder = preconditioner.uplimit_recorder;
threshold_increasment_ratio = preconditioner.threshold_increasment_ratio;
threshold_maxmin_ratio = preconditioner.threshold_maxmin_ratio;

% 首次使用: 预条件器尚未计算
if isempty(preconditioner.L)
    renew_flag = true;
    return;
end

% 老化: 预条件器使用了太多时间步
if length(recorder) >= uplimit_recorder
    renew_flag = true;
    return;
end

% 性能下降分析(需要历史数据)
if length(recorder) >= 2
    last_iter = recorder(end);
    prev_iter = recorder(end - 1);
    min_iter = min(recorder);
    
    % 绝对失败: 迭代次数过多
    if last_iter >= failure_trigger
        renew_flag = true;
        return;
    end
    
    % 急剧增长: 突然的性能下降
    if prev_iter > 0 && (last_iter / prev_iter) > threshold_increasment_ratio
        renew_flag = true;
        return;
    end
    
    % 渐进下降: 性能偏离历史最佳
    if min_iter > 0 && (last_iter / min_iter) > threshold_maxmin_ratio
        renew_flag = true;
        return;
    end
end

% 默认: 保持现有预条件器
renew_flag = false;

end
