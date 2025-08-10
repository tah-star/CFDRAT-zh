function cfdrat()
%CFDRAT 二维CFD仿真与可视化图形用户界面

global simulation_control  % 全局控制变量
simulation_control = 'stopped';  % 'running',  'stopped'

root_path = fileparts(mfilename('fullpath'));

% Define all source subdirectories based on your structure
src_folders = {
    fullfile(root_path, 'src', 'algorithm');
    fullfile(root_path, 'src', 'assembly');
    fullfile(root_path, 'src', 'driver');
    fullfile(root_path, 'src', 'grid');
    fullfile(root_path, 'src', 'utils');
    fullfile(root_path, 'src', 'visualization')
};
 
% Add all source folders to the MATLAB path
for i = 1:length(src_folders)
    if exist(src_folders{i}, 'dir')
        addpath(src_folders{i});
    else
        fprintf('Warning: Folder not found: %s\n', src_folders{i});
    end
end

try
    %% GUI界面布局

    % 主窗口
    fig = figure('Name', 'CFDRAT: 二维流体仿真平台', ...
        'Position', [100, 100, 900, 700], ...
        'MenuBar', 'none', 'ToolBar', 'none', 'Resize', 'off', ...
        'NumberTitle', 'off', 'Color', [0.95, 0.95, 0.95], ...
        'CloseRequestFcn', @closeGUI);

    % 主要面板
    simPanel = uipanel(fig, 'Title', '仿真参数设置', ...
        'FontSize', 14, 'FontWeight', 'bold', 'Units', 'pixels', ...
        'Position', [20, 20, 420, 660], 'BackgroundColor', [0.98, 0.98, 0.98]);
    visPanel = uipanel(fig, 'Title', '结果可视化', ...
        'FontSize', 14, 'FontWeight', 'bold', 'Units', 'pixels', ...
        'Position', [460, 20, 420, 660], 'BackgroundColor', [0.98, 0.98, 0.98]);
    
    % 仿真参数面板组件（左侧）
    y_start = 580;
    dy = 38;
    x_offset = -30;
    y_pos = y_start;

    % 计算域高度
    uicontrol(simPanel, 'Style', 'text', 'String', '计算域高度 H (m):', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.H = uicontrol(simPanel, 'Style', 'edit', 'String', '0.3', 'Position', [150+x_offset, y_pos, 260, 25], 'BackgroundColor', 'white');
    y_pos = y_pos - dy;

    % 入口速度函数
    uicontrol(simPanel, 'Style', 'text', 'String', '左入口法向速度 (m/s):', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.U_inlet_expr = uicontrol(simPanel, 'Style', 'edit', 'String', '0.15*4*y*(H-y)/H^2', 'Position', [150+x_offset, y_pos, 260, 25], 'BackgroundColor', 'white', 'TooltipString', '可使用变量 t, y, H');
    y_pos = y_pos - dy;

    % 入口条件渐变
    uicontrol(simPanel, 'Style', 'text', 'String', '入口渐变缓冲:', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.ramp_type = uicontrol(simPanel, 'Style', 'popupmenu', 'String', {'S型上升', '线性上升', '无'}, 'Position', [150+x_offset, y_pos, 100, 25]);
    uicontrol(simPanel, 'Style', 'text', 'String', '渐变时间 (s):', 'Position', [255+x_offset, y_pos, 85, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.t_ramp = uicontrol(simPanel, 'Style', 'edit', 'String', '0.5', 'Position', [350+x_offset, y_pos+7, 50, 19], 'BackgroundColor', 'white');
    y_pos = y_pos - dy;

    % 壁面条件
    uicontrol(simPanel, 'Style', 'text', 'String', '壁面边界条件:', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.wall_cond = uicontrol(simPanel, 'Style', 'popupmenu', 'String', {'无滑移', '滑移'}, 'Position', [150+x_offset, y_pos, 260, 25]);
    y_pos = y_pos - dy;

    % 出口条件
    uicontrol(simPanel, 'Style', 'text', 'String', '出口边界条件:', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.outlet_cond = uicontrol(simPanel, 'Style', 'popupmenu', 'String', {'定压力'}, 'Position', [150+x_offset, y_pos, 260, 25]);
    y_pos = y_pos - dy;

    % 流体密度
    uicontrol(simPanel, 'Style', 'text', 'String', '流体密度 (kg/m³):', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.rho = uicontrol(simPanel, 'Style', 'edit', 'String', '1', 'Position', [150+x_offset, y_pos, 260, 25], 'BackgroundColor', 'white');
    y_pos = y_pos - dy;

    % 动力粘度
    uicontrol(simPanel, 'Style', 'text', 'String', '动力粘度 (Pa·s):', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.mu = uicontrol(simPanel, 'Style', 'edit', 'String', '0.0001', 'Position', [150+x_offset, y_pos, 260, 25], 'BackgroundColor', 'white');
    y_pos = y_pos - dy;

    % 网格规模
    uicontrol(simPanel, 'Style', 'text', 'String', '节点数:', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.scale = uicontrol(simPanel, 'Style', 'edit', 'String', '50000', 'Position', [150+x_offset, y_pos, 260, 25], 'BackgroundColor', 'white');
    y_pos = y_pos - dy;

    % 记录时间间隔
    uicontrol(simPanel, 'Style', 'text', 'String', '记录时间间隔 (s):', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.T_record = uicontrol(simPanel, 'Style', 'edit', 'String', '0.05', 'Position', [150+x_offset, y_pos, 260, 25], 'BackgroundColor', 'white');
    y_pos = y_pos - dy;

    % 仿真总时长
    uicontrol(simPanel, 'Style', 'text', 'String', '仿真总时长 (s):', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.t_end = uicontrol(simPanel, 'Style', 'edit', 'String', '5', 'Position', [150+x_offset, y_pos, 260, 25], 'BackgroundColor', 'white');
    y_pos = y_pos - dy;

    % 求解器速度
    uicontrol(simPanel, 'Style', 'text', 'String', '求解器速度:', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.speed_opt = uicontrol(simPanel, 'Style', 'popupmenu', 'String', {'快', '中'}, 'Position', [150+x_offset, y_pos, 260, 25]);
    y_pos = y_pos - dy;

    % 并行计算
    uicontrol(simPanel, 'Style', 'text', 'String', '并行计算:', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.parallel_opt = uicontrol(simPanel, 'Style', 'popupmenu', 'String', {'关', '开'}, 'Position', [150+x_offset, y_pos, 260, 25]);
    y_pos = y_pos - dy; 

    % 障碍物图像
    default_image_path = fullfile(root_path, 'examples', 'hamster_demo.png');
    uicontrol(simPanel, 'Style', 'text', 'String', '障碍物图像:', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.image_path = uicontrol(simPanel, 'Style', 'edit', 'String', default_image_path, 'Position', [150+x_offset, y_pos, 155, 25], 'BackgroundColor', 'white');
    uicontrol(simPanel, 'Style', 'pushbutton', 'String', '浏览...', 'Position', [315+x_offset, y_pos, 95, 25], 'Callback', @browseImage);
    y_pos = y_pos - dy;

    % 保存路径
    default_save_path = fullfile(root_path, 'results', 'hamster_demo_result.mat');
    uicontrol(simPanel, 'Style', 'text', 'String', '结果保存位置:', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.save_path = uicontrol(simPanel, 'Style', 'edit', 'String', default_save_path, 'Position', [150+x_offset, y_pos, 155, 25], 'BackgroundColor', 'white');
    uicontrol(simPanel, 'Style', 'pushbutton', 'String', '选择...', 'Position', [315+x_offset, y_pos, 95, 25], 'Callback', @browseSaveFile);

    % === 修改：按钮布局改为两个按钮并排 ===
    % 运行仿真按钮（左）
    h.runBtn = uicontrol(simPanel, 'Style', 'pushbutton', 'String', '开始仿真', 'Position', [20, 30, 180, 45], 'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.2, 0.7, 0.2], 'ForegroundColor', 'white', 'Callback', @runSimulation);
    
    % === 新增：停止仿真按钮（右）===
    h.stopBtn = uicontrol(simPanel, 'Style', 'pushbutton', 'String', '停止仿真', 'Position', [220, 30, 180, 45], 'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.8, 0.2, 0.2], 'ForegroundColor', 'white', 'Callback', @stopSimulation, 'Enable', 'off');

    % 可视化面板组件（右侧）
    y_pos_vis = 580;
    x_offset_vis = -40;
    uicontrol(visPanel, 'Style', 'text', 'String', '数据文件:', 'Position', [20+x_offset_vis, y_pos_vis, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.load_path = uicontrol(visPanel, 'Style', 'edit', 'String', '', 'Position', [150+x_offset_vis, y_pos_vis, 155, 25], 'BackgroundColor', 'white');
    uicontrol(visPanel, 'Style', 'pushbutton', 'String', '浏览...', 'Position', [315+x_offset_vis, y_pos_vis, 95, 25], 'Callback', @browseLoadFile);
    y_pos_vis = y_pos_vis - dy;
    uicontrol(visPanel, 'Style', 'text', 'String', '显示模式:', 'Position', [20+x_offset_vis, y_pos_vis, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.vis_display_mode = uicontrol(visPanel, 'Style', 'popupmenu', 'String', {'内部色块', '仅边界'}, 'Position', [150+x_offset_vis, y_pos_vis, 260, 25]);
    h.visBtn = uicontrol(visPanel, 'Style', 'pushbutton', 'String', '播放动画', 'Position', [20, 30, 380, 45], 'FontSize', 16, 'FontWeight', 'bold', 'BackgroundColor', [0.2, 0.4, 0.8], 'ForegroundColor', 'white', 'Callback', @runVisualization);

    % 状态栏
    h.status_label = uicontrol(fig, 'Style', 'text', 'String', '状态: 就绪', 'Position', [20, 5, 400, 20], 'FontSize', 12, 'ForegroundColor', [0.2, 0.6, 0.2], 'BackgroundColor', [0.95, 0.95, 0.95], 'HorizontalAlignment', 'left');

    % 存储句柄
    guidata(fig, h);

catch ME
    fprintf('GUI创建失败: %s\n', ME.message);
    if exist('fig', 'var') && ishandle(fig), delete(fig); end
    rethrow(ME);
end

%% 回调函数和辅助函数

    function browseImage(~, ~)
        h = guidata(gcf);
        [file, path] = uigetfile({'*.png;*.jpg;*.jpeg;*.bmp;*.gif;*.tiff', '图像文件'; '*.*', '所有文件'}, '选择障碍物图像文件');
        if file ~= 0, set(h.image_path, 'String', fullfile(path, file)); end
    end

    function browseSaveFile(~, ~)
        h = guidata(gcf);
        [file, path] = uiputfile('*.mat', '选择保存位置', 'sim_result.mat');
        if file ~= 0
            full_path = fullfile(path, file);
            set(h.save_path, 'String', full_path);
            set(h.load_path, 'String', full_path);
        end
    end

    function browseLoadFile(~, ~)
        h = guidata(gcf);
        [file, path] = uigetfile('*.mat', '选择仿真数据文件');
        if file ~= 0, set(h.load_path, 'String', fullfile(path, file)); end
    end

    function runSimulation(~, ~)
        h = guidata(gcf);
        simulation_control = 'running';
        setButtonsState(h, 'running', '状态: 准备仿真中...');
        drawnow;

        try
            % 从GUI收集所有输入参数并打包成结构体
            params = struct();
            params.H = str2double(get(h.H, 'String'));
            params.u_inlet_func_str = get(h.U_inlet_expr, 'String');

            ramp_opts = get(h.ramp_type, 'String');
            params.ramp_opt = ramp_opts{get(h.ramp_type, 'Value')};
            params.t_ramp = str2double(get(h.t_ramp, 'String'));

            wall_opts = get(h.wall_cond, 'String');
            params.slip_opt = wall_opts{get(h.wall_cond, 'Value')};

            params.rho = str2double(get(h.rho, 'String'));
            params.mu = str2double(get(h.mu, 'String'));
            params.node_scale = str2double(get(h.scale, 'String'));
            params.t_simu = str2double(get(h.t_end, 'String'));
            params.T_record = str2double(get(h.T_record, 'String'));

            speed_opts = get(h.speed_opt, 'String');
            params.speed_opt = speed_opts{get(h.speed_opt, 'Value')};

            parallel_opts = get(h.parallel_opt, 'String');
            params.enable_parallel = parallel_opts{get(h.parallel_opt, 'Value')};

            params.image_filename = get(h.image_path, 'String');

            display_modes = get(h.vis_display_mode, 'String');
            params.image_display_opt = display_modes{get(h.vis_display_mode, 'Value')};

            params.save_filename = get(h.save_path, 'String');

            % 输入验证
            if any(isnan([params.H, params.t_ramp, params.rho, params.mu, params.node_scale, params.t_simu, params.T_record]))
                error('所有数值输入必须为有效数字。');
            end
            if isempty(params.image_filename) || ~exist(params.image_filename, 'file')
                error('请选择有效的障碍物图像文件。');
            end

            fprintf('\n 仿真启动！\n');

            % 调用主运行函数
            run_simu(params);
            
            % 检查仿真是否被用户停止 
            if strcmp(simulation_control, 'stopped')
                set(h.status_label, 'String', '状态: 仿真已被用户停止');
            else
                set(h.load_path, 'String', params.save_filename);
            end

        catch ME
            msgbox(['错误: ' ME.message], '错误', 'error');
            fprintf('仿真错误详情: %s\n', getReport(ME, 'extended'));
        end

        simulation_control = 'stopped';
        setButtonsState(h, 'stopped', '状态: 就绪');
    end

    function stopSimulation(~, ~)
        h = guidata(gcf);
        simulation_control = 'stopped';
        setButtonsState(h, 'stopped', '状态: 正在停止仿真...');
        fprintf('*** 用户请求停止仿真 ***\n');
        % 稍后状态会在 runSimulation 函数结束时更新为"就绪"
    end

    function runVisualization(~, ~)
        h = guidata(gcf);
        setButtonsState(h, 'visualizing', '状态: 加载动画中...');
        drawnow;
        try
            load_filename = get(h.load_path, 'String');
            if isempty(load_filename) || ~exist(load_filename, 'file'), error('请选择有效的数据文件。'); end
            display_modes = {'内部色块', '仅边界'};
            display_mode = display_modes{get(h.vis_display_mode, 'Value')};
            draw_uv(load_filename, display_mode);

        catch ME
            msgbox(['可视化错误: ' ME.message], '错误', 'error');
            fprintf('可视化错误详情: %s\n', getReport(ME, 'extended'));
        end
        setButtonsState(h, 'stopped', '状态: 就绪');
    end

% 按钮状态管理函数
    function setButtonsState(h, state, statusText)
        switch state
            case 'running'
                set(h.runBtn, 'Enable', 'off');
                set(h.stopBtn, 'Enable', 'on');
                set(h.visBtn, 'Enable', 'off');
            case 'stopped'
                set(h.runBtn, 'Enable', 'on');
                set(h.stopBtn, 'Enable', 'off');
                set(h.visBtn, 'Enable', 'on');
            case 'visualizing'
                set([h.runBtn, h.stopBtn, h.visBtn], 'Enable', 'off');
            otherwise  % 默认为 stopped 状态
                set(h.runBtn, 'Enable', 'on');
                set(h.stopBtn, 'Enable', 'off');
                set(h.visBtn, 'Enable', 'on');
        end
        set(h.status_label, 'String', statusText);
    end

    function closeGUI(~, ~)
        simulation_control = 'stopped';
        fprintf('GUI已关闭。\n');
        delete(gcf);
    end

end
