%%%%%%   xujunwu@hotmail.com   %%%%%%%%%%%%
 cd('C:\Users\Miao Li\Desktop\neuropile-camk2\#011 remd')
start_baseline=1;
resting=60;
baseline_window=10;


data_num = dir('*.xlsx');
mkdir result
for j=1:length(data_num)
    raw=xlsread(data_num(j).name);
    data=xlsread(data_num(j).name);
    
    str=data_num(j).name;
    title_name=str(1:end-5);
    numFrames=size(data,1);
    
    
    %baseline_window=round(numFrames/10);
    
    
    soma_number=size(data,2);
    threshold=zeros(1,soma_number);
    for i=1:soma_number
        
        for k=start_baseline:(resting-baseline_window+1)
            if k==start_baseline
                min_win=mean(data(k:(k+baseline_window-1),i));
            else
                min_win_now=mean(data(k:(k+baseline_window-1),i));
                if min_win_now<min_win
                    min_win=min_win_now;
                end
            end
        end
        threshold(i)=min_win;
    end
    data_dff=(data-repmat(threshold,numFrames,1))./repmat(threshold,numFrames,1);
    
    cd result
    xlswrite([title_name,'.xlsx'],data_dff);
    cd ..
end

clear
