%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   By Xujun Wu                      %
%   Email:xujunwu@hotmail.com                % 
%   2021.04.29                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

cd('C:\Users\Miao Li\Desktop\neuropile-camk2\#011 remd\dff0')
baseline_start=1;
baseline=60;


par_start=61;
par_end=120;

threshold=3;
baseline_window=10;



data_num = dir('*.xlsx');
mkdir result
for j=1:length(data_num)
    data=xlsread(data_num(j).name); %%data=xlsread(data_num(j).name,'Sheet1');
    str=data_num(j).name;
    title_name=[str(1:end-5) '_result'];
    pat='_.*.xlsx';
    m = regexpi(str, pat, 'match');
    t_per_frame=str2num(m{1}(2:end-5));
    numFrames=size(data,1);
    
    %%%baseline_window=round(numFrames/10);
    
    
    spine_number=size(data,2);
    threshold_all=zeros(1,spine_number);
    for i=1:spine_number
        
        for k=baseline_start:(baseline-baseline_window+1)
            if k==baseline_start
                min_win=mean(data(k:(k+baseline_window-1),i));
                tem_data=data(k:(k+baseline_window-1),i);
            else
                min_win_now=mean(data(k:(k+baseline_window-1),i));
                if min_win_now<min_win
                    min_win=min_win_now;
                    tem_data=data(k:(k+baseline_window-1),i);
            
                end
            end
        end
        
        threshold_all(i)=min_win+threshold*std(tem_data);
        


        
        
        
    end
    
    
   data=data(par_start:par_end,:);
    [average_amplitude,frequency,duration,total_activity]=extract_calcium_transient_parameter(data,threshold_all,t_per_frame);
    title_x=[];
    for i=1:spine_number
        title_xplus={['cell' num2str(i)]};
        title_x=[title_x ,title_xplus];
    end
    output=[average_amplitude;frequency;duration;total_activity];
    
    %%%%%%%%%%%%%%??????????????????excel????%%%%%%
    [m, n] = size(output);
    output_cell = mat2cell(output, ones(m,1), ones(n,1));  % ??data??????m*n??cell????
    result = [title_x; output_cell];   % ??????????????????????result
    title_y=[{'Begin'},'average_amplitude','frequency_number','duration_seconds','total_activity_Height_times_time']';
    cd result
    xlswrite([title_name,'.xls'],result,1,'B1');
    xlswrite([title_name,'.xls'],title_y,1,'A1');
    
    
    [thr_m, thr_n] = size(threshold_all);
    thr_output_cell = mat2cell(threshold_all, ones(thr_m,1), ones(thr_n,1)); 
    xlswrite([title_name,'threshold.xlsx'],thr_output_cell);
    
    
    
    cd ..
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function [startinx,endinx]=find_trasient_interval(raster)
% % startinx=[];
% % endinx=[];
% % b=raster';
% % x=find(b);
% % n=size(raster,1);
% % %find start end inx
% % for i=x
% %     switch i
% %         case 1
% %             if b(2)==1
% %                 startinx=[startinx 1];
% %             elseif b(2)==0
% %                 startinx=[startinx 1];
% %                 endinx=[endinx 1];
% %             end
% %         case n
% %             if b(n-1)==1
% %                 endinx=[endinx n];
% %             elseif b(n-1)==0
% %                 startinx=[startinx n];
% %                 endinx=[endinx n];
% %             end
% %         otherwise
% %             if b(i+1)==1 && b(i-1)==0
% %                 startinx=[startinx i];
% %                 
% %             elseif b(i+1)==0 && b(i-1)==0
% %                 startinx=[startinx i];
% %                 endinx=[endinx i];
% %             elseif b(i+1)==0 && b(i-1)==1
% %                 endinx=[endinx i];
% %             end
% %     end
% % end
% % startinx=sort(startinx);
% % endinx=sort(endinx);
% % endinx=endinx+1;%%??????????????????????threshold????????????????????????????
% % if length(endinx)>0
% %     if endinx(end)>n
% %         endinx(end)=endinx(end)-1;
% %     end
% % end
% % end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [average_amplitude,frequency,duration,total_activity]=extract_calcium_transient_parameter(deltaF_F0,threshold_data,t_per_frame)
% raster=double(bsxfun(@gt,deltaF_F0,threshold_data));
% for i=1:size(raster,2)  %find parameters for each spine
%     if sum(raster(:,i))>0
%         % transient_inx=find(raster(:,i));
%         [startinx,endinx]=find_trasient_interval(raster(:,i));
%         n=size(startinx,2); % number of transient intervals for each spine
%         timepoint=size(deltaF_F0,1);
%         amplitude=zeros(n,1);
%         for j=1:n   %find peak amplitude in each interval for each spine
%             amplitude(j,1)=max(deltaF_F0(startinx(j):endinx(j),i));
%         end
%         frequency(1,i)=n;  %%%% one minute
%         total_activity(1,i)=sum(deltaF_F0(raster(:,i)>0,i))*t_per_frame; % calculate total_activity
%         duration(1,i)=sum(raster(:,i))*t_per_frame;
%         average_amplitude(1,i)=mean(amplitude,1);
%     else
%         average_amplitude(1,i)=0;
%         total_activity(1,i)=0;
%         frequency(1,i)=0;
%         duration(1,i)=0;
%         
%     end
% end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
