function draw_uv(data_source, image_display_mode)
%DRAW_UV 交互式速度场动画带GUI控制
%   创建CFD仿真时间相关速度数据的动画可视化，具有播放控制、
%   时间滑块和速度调节功能。
%
%   用法:
%       draw_uv('path/to/simulation_data.mat')
%       draw_uv('path/to/simulation_data.mat', 'bound only')
%
%   输入参数:
%       data_source        - 包含仿真结果的MAT文件路径
%       image_display_mode - (可选) 障碍物显示模式:
%                           'show image': 纹理映射 (默认)
%                           'bound only': 仅边界渲染
%
%   功能特性:
%   - 可变速度控制的实时播放
%   - 通过滑块进行交互式时间导航
%   - 自动障碍物遮罩和边界可视化
%   - 使用jet色彩方案的速度幅值颜色映射
%
%   GUI控制:
%   - 播放/暂停按钮用于动画控制
%   - 时间滑块用于直接帧导航
%   - 速度选择器用于播放率调节
%   - 时间显示显示当前仿真时间
%
%   See also SAVE_DATA, GET_BLOCKED_MASK, PCOLOR, COLORBAR.

%% 数据加载和初始化
fprintf('--- 速度场动画播放 ---\n');
try
    s = load(data_source, 'info');
    info = s.info;
catch ME
    error('加载文件失败: %s.\n请确保路径正确且文件包含"info"结构体.\nMATLAB错误: %s', data_source, ME.message);
end

% 提取仿真数据和网格参数
u_all = info.u_all;
v_all = info.v_all;
p_all = info.p_all;
grid_p = info.grid_p;
h = grid_p.h;
[Ny, Nx] = size(p_all(:, :, 1));
T_record = info.T_record;
block_info = info.block_info;
num_frames = size(u_all, 3);
if num_frames <= 1, warning('仅找到一个时间帧。无法播放动画。'), return; end

% 配置障碍物可视化模式
show_image = strcmpi(image_display_mode, '内部色块') && isfield(block_info{1}, 'image');
if strcmpi(image_display_mode, '内部色块') && ~isfield(block_info{1}, 'image')
    fprintf('注意: 请求显示图像，但未找到图像数据。对障碍物使用实体填充。\n');
end

% 设置计算网格和障碍物遮罩
L = grid_p.L; H = grid_p.H;
x_vec = h*(0:1:Nx-1);
y_vec = h*(0:1:Ny-1);
[X, Y] = meshgrid(x_vec, y_vec);
mask_block = get_blocked_mask(block_info, X, Y);
dt_record = T_record;

% 为所有时间帧预计算速度幅值
speed_all = zeros(size(u_all), 'single');
for i = 1:num_frames
    speed = sqrt(u_all(:,:,i).^2 + v_all(:,:,i).^2);
    speed(mask_block) = NaN; % 遮罩障碍物内部以进行适当可视化
    speed_all(:,:,i) = speed;
end

%% GUI布局和控制设置
% 创建具有适当宽高比的主图形
fig_width = 900;
fig_height = fig_width * (H/L) + 150;
fig = figure('Name', sprintf('速度动画播放器: %s', data_source), ...
    'Position', [150, 150, fig_width, fig_height], ...
    'NumberTitle', 'off', 'Visible', 'off');
ax = axes('Parent', fig, 'Position', [0.1, 0.25, 0.85, 0.7]);

% 创建GUI控制元素
btn_play = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', '▶ 播放', 'Position', [50, 20, 100, 40], 'FontSize', 12, 'Callback', @play_pause_callback);
slider_time = uicontrol('Parent', fig, 'Style', 'slider', 'Position', [170, 25, fig_width-350, 30], 'Value', 1, 'Min', 1, 'Max', num_frames, 'SliderStep', [1/(num_frames-1), 10/(num_frames-1)], 'Callback', @slider_callback);
uicontrol('Parent', fig, 'Style', 'text', 'String', '播放速度:', 'Position', [fig_width-170, 45, 80, 20], 'HorizontalAlignment', 'center');
popup_speed = uicontrol('Parent', fig, 'Style', 'popupmenu', 'String', {'0.1倍', '0.2倍', '0.5倍', '1倍 (实时)', '2倍', '5倍'}, 'Value', 4, 'Position', [fig_width-170, 20, 80, 30], 'Callback', @speed_change_callback);
lbl_time = uicontrol('Parent', fig, 'Style', 'text', 'Position', [(fig_width-150)/2, 70, 150, 20],'String', '时间: 0.000 s', 'FontSize', 11, 'HorizontalAlignment', 'center');

%% 动画引擎和回调函数
% 应用程序状态变量
current_frame = 1; 
is_playing = false; 
playback_speed_multiplier = 1;
h_pcolor = []; % 用于高效绘图更新的句柄

    % 具有静态/动态元素分离的绘图更新函数
    function update_plot(frame_idx)
        frame_idx = round(frame_idx);
        if frame_idx > num_frames, frame_idx = num_frames; end
        if frame_idx < 1, frame_idx = 1; end
        
        if isempty(h_pcolor) || ~isvalid(h_pcolor)
            % 初始绘图：创建静态元素（障碍物、坐标轴、颜色条）
            axes(ax);
            h_pcolor = pcolor(X, Y, speed_all(:,:,frame_idx));
            shading interp; hold on;
            
            % 使用纹理或仅边界显示渲染障碍物
            for k = 1:length(block_info)
                b = block_info{k};
                if show_image
                    h_img = image(b.x_coords, b.y_coords, flipud(b.image));
                    if isfield(b, 'mask') && ~isempty(b.mask)
                        set(h_img, 'AlphaData', flipud(b.mask));
                    end
                    plot([b.points(:,1); b.points(1,1)], [b.points(:,2); b.points(1,2)], 'k-', 'LineWidth', 1.0);
                else
                    fill(b.points(:,1), b.points(:,2), [0.98 0.98 0.98], 'EdgeColor', 'k', 'LineWidth', 1.0);
                end
            end
            
            hold off; axis equal; axis([0 L 0 H]);
            colormap(ax, 'jet'); c = colorbar; c.Label.String = '速度 (m/s)';
            xlabel('X (m)'); ylabel('Y (m)');
        else
            % 后续更新：仅刷新速度数据
            set(h_pcolor, 'CData', speed_all(:,:,frame_idx));
        end

        % 更新时间相关的显示元素
        physical_time = (frame_idx - 1) * dt_record;
        title(ax, sprintf('速度幅值 (时间: %.3f s)', physical_time));
        set(lbl_time, 'String', sprintf('时间: %.3f s', physical_time));
        set(slider_time, 'Value', frame_idx);
        drawnow('limitrate');
    end

    % 具有实时同步的播放/暂停按钮回调
    function play_pause_callback(~, ~)
        is_playing = ~is_playing;
        if is_playing
            set(btn_play, 'String', '❚❚ 暂停');
            real_time_start = tic;
            physical_time_at_start = (current_frame - 1) * dt_record;
            while is_playing && isvalid(fig)
                real_time_elapsed = toc(real_time_start);
                target_physical_time = physical_time_at_start + real_time_elapsed * playback_speed_multiplier;
                target_frame = floor(target_physical_time / dt_record) + 1;
                if target_frame > num_frames
                    % 无缝循环动画
                    current_frame = 1; physical_time_at_start = 0; real_time_start = tic;
                    update_plot(current_frame); continue;
                end
                if target_frame > current_frame
                    update_plot(target_frame); current_frame = target_frame;
                end
                pause(0.01); % 将控制权交给GUI和系统
            end
            if isvalid(btn_play), set(btn_play, 'String', '▶ 播放'); end
            is_playing = false;
        else
            set(btn_play, 'String', '▶ 播放');
        end
    end

    % 用于直接导航的时间滑块回调
    function slider_callback(source, ~)
        if is_playing, play_pause_callback(); end % 手动控制期间自动暂停
        new_frame = round(get(source, 'Value'));
        if new_frame ~= current_frame
            update_plot(new_frame); current_frame = new_frame;
        end
    end

    % 播放速度选择器回调
    function speed_change_callback(source, ~)
        speed_str = source.String{source.Value};
        speed_val_str = regexp(speed_str, '[\d\.]+', 'match');
        playback_speed_multiplier = str2double(speed_val_str{1});
        if is_playing, play_pause_callback(); play_pause_callback(); end % 使用新速度重启
    end

%% 初始化和显示
% 渲染初始帧并激活GUI
update_plot(1);
set(fig, 'Visible', 'on');

end
