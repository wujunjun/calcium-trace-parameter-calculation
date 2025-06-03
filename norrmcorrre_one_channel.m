%%% xujun  %%%%%%%%%%%

clear

data_path='F:\xujun\20250517\t';
file_expr='*.tif';
frame_batch = 100;  % 根据内存调整分块大小 (建议32-100)
batch_size=12;


data_num = dir(fullfile(data_path,file_expr));
data_num = {data_num.name}';

warning('off','all') % 关闭所有警告
if isempty(data_num)
    error('No channel1 files found');
end

for i=1:length(data_num)
    % ================= 文件路径处理 =================
    [~,base_name] = fileparts(data_num{i});
    base_name = base_name(1:end); % 移除'_chan1'后缀
    
    chan1_file = fullfile(data_path, data_num{i});
    % chan2_file = fullfile(data_path, [base_name '_chan2.tif']);
    tifinfo=imfinfo(chan1_file);
    d1=tifinfo(1).Height;
    d2=tifinfo(1).Width;




    % ================= 参数初始化 =================
    output_suffix = '_corrected.tif';
    opts = NoRMCorreSetParms('d1',d1,'d2',d2,...  % 根据实际尺寸修改
        'grid_size',[batch_size,batch_size],'mot_uf',4,'max_dev',3,'bin_width',100,...      % 保持与原参数一致
        'max_shift',10,'us_fac',10,'iter',3);                % 重要参数勿修改

    % ================= 通道1处理 =================
    fprintf('Processing %s (%d/%d)\n',base_name,i,length(data_num));
    [~,shifts_total] = process_single_channel(chan1_file, opts,...
        fullfile(data_path,[base_name '_red' output_suffix]),d1,d2,batch_size);

    % ================= 通道2应用位移 =================
    % apply_shifts_to_channel(chan2_file, shifts_total, opts,...
    %     fullfile(data_path,[base_name '_green' output_suffix]),...
    %     frame_batch);
end

function [M_corrected, shifts] = process_single_channel(input_file, opts, output_path,d1,d2,batch_size)
% 处理参考通道（chan1）
    data = read_file(input_file);
    data = single(data); % 转换为单精度节省内存

    opts2 = NoRMCorreSetParms('d1',d1,'d2',d2,...  % 根据实际尺寸修改
        'grid_size',[batch_size,batch_size],'mot_uf',4,'max_dev',3,'bin_width',100,...      % 保持与原参数一致
        'max_shift',20,'us_fac',10,'init_batch',100,'iter',1);      
    
    % 执行配准
    [~,~,template] = normcorre_batch(data, opts2);
    [M_corrected,shifts] = normcorre_batch(data, opts,template);
    
    % 保存结果并添加平均帧
    save_with_mean_frame(M_corrected, output_path);
    
end
function apply_shifts_to_channel(input_file, shifts, opts, output_path, frame_batch)
% 应用位移到目标通道（chan2）并分块保存
    tiff_info = imfinfo(input_file);
    total_frames = length(tiff_info);
    
    % 初始化输出文件
%     tiff_obj = Tiff(output_path, 'w');
%     setTag(tiff_obj, tiff_info(1));
%     
%     % 预写平均帧位置（保持原逻辑）
%     avg_frame = get_channel_average(input_file);
%     write(tiff_obj, uint16(avg_frame));
    
    % 分块处理
    for bidx = 1:ceil(total_frames/frame_batch)
        batch_range = (bidx-1)*frame_batch+1 : min(bidx*frame_batch, total_frames);
        data_batch = single(read_file(input_file, batch_range(1),length(batch_range)));
        
        try
            % 应用位移时添加边界处理
            corrected_batch = apply_shifts(data_batch,...
                shifts(batch_range), opts);
        catch ME
            fprintf('Batch %d error: %s\n', bidx, ME.message);
            corrected_batch = data_batch; % 失败时保留原数据
        end
        
        % 逐帧写入并清理内存
        for fidx = 1:size(corrected_batch,3)
            frame = corrected_batch(:,:,fidx); % 填充缺失值
            imwrite(uint16(frame),output_path, 'writemode', 'append');
        end
        clear data_batch corrected_batch frame
    end
    
    M2=read_file(output_path);
    writeTiff(uint16(cat(3,mean(M2,3),M2)),output_path);



end
function save_with_mean_frame(data, output_path)
% 保持原cat(3,mean...)逻辑的写入方式
    avg_frame = mean(data,3);
    data = cat(3, avg_frame, data);
    writeTiff(uint16(data), output_path);
end
