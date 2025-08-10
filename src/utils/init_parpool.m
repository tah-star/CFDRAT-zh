% manage_parallel_pool.m
function init_parpool(enable_parallel)
%INIT_PARPOOL 管理并行池的启动和停止
%   根据用户选择启动或停止并行池
if strcmp(enable_parallel, '开')
    if isempty(gcp('nocreate'))
        fprintf('正在启动并行池，请稍等...\n');
        parpool;
    end
    if ~isempty(gcp('nocreate'))
        fprintf('并行池启动成功。\n');
    else
        fprintf('启动并行池失败，将以串行模式运行。\n');
    end
else
    if ~isempty(gcp('nocreate'))
        fprintf('关闭并行池中...\n');
        delete(gcp('nocreate'));
    else
        fprintf('并行池已关闭\n');
    end
end
end
