%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Written by Xujun Wu;Avital;Rongwei Si    %
%   Email:xujunwu@hotmail.com                % 
%   2019.12.20                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cd('G:\two-phton\qq\8#_original_intensity_trace\shaftA')
cd('C:\Users\wxj\Desktop\zhongyuan_spineordendriteactivity\spine_dff\thr3per10')
warning('off','MATLAB:polyfit:PolyNotUnique');


figure;
set(gcf,'unit','centimeters','position',[14 1 15 26]);
set(gcf,'WindowStyle','normal');

% Plot1 = questdlg('Plot results of robust regression for each spine?', 'Plot?', 'Yes','No','Yes');



% t_win=7;  %%%  t_win * 2
Plot1='Yes';
threshold=3;
prctile_num=8;
thr_length=4;
corr_method=3;% % 1 for raw mean; 2 for spine above throshold; 3 for shaft above throshold
shaft_contribution_method=2; % 1 for remove all always; 2 for informative shaft and P<0.05
corr_time=60;  %%%%  length of time(S) for calculating the correlation between spines  %%%  for corr_method == 1
min_threshold=0.1;  %  the min threshold for DF/F
max_cc_shaft_spine=0.6;  %%%%  if the correlation coefficient between shaft and spine is great than this value,don't remove backgroud noise from shaft.
threshold_for_negetive_value=-2;
wdw=10;  %%%  the number of frames for prctile

paramater=['thres_',num2str(threshold),'_','prctile_',num2str(prctile_num)];


data_num = dir('*.txt');

for j=1:length(data_num)
    data=load(data_num(j).name);
    str=data_num(j).name;
    pat='_.*txt';
    m = regexpi(str, pat, 'match');
    t_per_frame=str2double(m{1}(2:end-4));
    % wdw=round(t_win/t_per_frame);
    
    numFrames=size(data,1);
    spine_number=size(data,2)-1;
    smoothBaseline=zeros(size(data));
    if numFrames > 2*wdw
        for k=1:size(data,2)
            dataSlice=data(:,k);
            temp=zeros(numFrames-2*wdw,1);
            for i=(wdw+1):(numFrames-wdw)
                temp(i-wdw)=prctile(dataSlice(i-wdw:i+wdw),prctile_num);
            end
            smoothBaseline(:,k)=[prctile(dataSlice(1:wdw),prctile_num)*ones(wdw,1) ; temp; prctile(dataSlice(end-wdw:end),prctile_num)*ones(wdw,1)];
            smoothBaseline(:,k)=runline(smoothBaseline(:,k),wdw,1);
        end
    else
        for k=1:size(data,2)
            smoothBaseline(:,k)=[ones(numFrames,1)*prctile(data(:,k),prctile_num)];
        end
    end
    
    

    deltaF_F0=(data-smoothBaseline)./smoothBaseline;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    data_raw_dff=deltaF_F0;
    
    corr_time_numFrames=round(corr_time/t_per_frame);
    
    [timepoint,C]=size(deltaF_F0);

    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    spine_specific_calcium=zeros(timepoint,spine_number);
    shaft=deltaF_F0(:,C);
    regression_slope=zeros(1,spine_number);
    regression_slope_pvalue=zeros(1,spine_number);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%-%%%%%%%%%%%%%%spine specific signal %%%%%%%%
    [smoothDist,x] = ksdensity(shaft);
    [~,indPeak]=max(smoothDist);
    
        xFit=[x(1:indPeak) (x(indPeak)+(x(2)-x(1))*[1:(indPeak-1)])];
    dataToFit=[smoothDist(1:indPeak) fliplr(smoothDist(1:(indPeak-1)))];

    
    [sigma_shaft,mean_shaft,~]=mygaussfit(xFit',dataToFit,0.1,1);
    if mean_shaft==-inf
        dataToFit=interp1(xFit,dataToFit,xFit(1):(xFit(end)-xFit(1))/10:xFit(end),'spline');
        [sigma_shaft,mean_shaft,~]=mygaussfit((xFit(1):(xFit(end)-xFit(1))/10:xFit(end))',dataToFit);
    end
    
    if sum((shaft-mean_shaft)>(3*sigma_shaft))==0
        regression_slope=7*ones(1,spine_number);
        regression_slope_pvalue=7*ones(1,spine_number);
        spine_specific_calcium=deltaF_F0(:,1:(end-1));
        data_no_remove=spine_specific_calcium;
    elseif sum((shaft-mean_shaft)>(3*sigma_shaft))>0
        data_no_remove=deltaF_F0(:,1:(end-1));
        if shaft_contribution_method==2
            for k=1:spine_number
                curr_spine = deltaF_F0(:,k);
                [Coeff,stats] = robustfit(shaft(:), curr_spine(:)); %robustfit for each spine in a dendrite
                regression_slope(1,k)= Coeff(2);
                regression_slope_pvalue(1,k)= stats.p(2);
                if stats.p(2)<0.05
                    tem_value=corrcoef(deltaF_F0(:,k),deltaF_F0(:,end));
                    
                    [smoothDist,x] = ksdensity(curr_spine);
                    [~,indPeak]=max(smoothDist);
                    
                    
                    xFit=[x(1:indPeak) (x(indPeak)+(x(2)-x(1))*[1:(indPeak-1)])];
                    dataToFit=[smoothDist(1:indPeak) fliplr(smoothDist(1:(indPeak-1)))];

                    [sigma_curr_spine,mean_curr_spine,~]=mygaussfit(xFit',dataToFit,0.1,1);
                    if mean_curr_spine==-inf
                        dataToFit=interp1(xFit,dataToFit,xFit(1):(xFit(end)-xFit(1))/10:xFit(end),'spline');
                        [sigma_curr_spine,mean_curr_spine,~]=mygaussfit((xFit(1):(xFit(end)-xFit(1))/10:xFit(end))',dataToFit);
                    end
                    
                    if tem_value(2)>max_cc_shaft_spine
                        spine_specific_calcium(:,k) = curr_spine-((shaft-mean_shaft)>(3*sigma_shaft)).*(regression_slope(1,k)*shaft);
                        tem_value2=spine_specific_calcium(:,k);
                        tem_value2(tem_value2<threshold_for_negetive_value)=mean_curr_spine +sigma_curr_spine*randn(sum(tem_value2<threshold_for_negetive_value),1);
                        spine_specific_calcium(:,k)=tem_value2;
                    else
                        spine_specific_calcium(:,k) = curr_spine-(regression_slope(1,k)*shaft);
                        tem_value2=spine_specific_calcium(:,k);
                        tem_value2(tem_value2<threshold_for_negetive_value)=mean_curr_spine +sigma_curr_spine*randn(sum(tem_value2<threshold_for_negetive_value),1);
                        spine_specific_calcium(:,k)=tem_value2;
                    end
                else
                    spine_specific_calcium(:,k) = curr_spine;
                    
                end
                
            end
        elseif shaft_contribution_method==1
            for k=1:spine_number
                curr_spine = deltaF_F0(:,k);
                
                [smoothDist,x] = ksdensity(curr_spine);
                [~,indPeak]=max(smoothDist);
                
                
                xFit=[x(1:indPeak) (x(indPeak)+(x(2)-x(1))*[1:(indPeak-1)])];
                dataToFit=[smoothDist(1:indPeak) fliplr(smoothDist(1:(indPeak-1)))];
                
                
                
                [sigma_curr_spine,mean_curr_spine,~]=mygaussfit(xFit',dataToFit,0.1,1);
                if mean_curr_spine==-inf
                    dataToFit=interp1(xFit,dataToFit,xFit(1):(xFit(end)-xFit(1))/10:xFit(end),'spline');
                    [sigma_curr_spine,mean_curr_spine,~]=mygaussfit((xFit(1):(xFit(end)-xFit(1))/10:xFit(end))',dataToFit);
                end
                
                
                [Coeff,stats] = robustfit(shaft(:), curr_spine(:)); %robustfit for each spine in a dendrite
                regression_slope(1,k)= Coeff(2);
                regression_slope_pvalue(1,k)= stats.p(2);
                tem_value=corrcoef(deltaF_F0(:,k),deltaF_F0(:,end));
                if tem_value(2)>max_cc_shaft_spine
                    spine_specific_calcium(:,k) = curr_spine-((shaft-mean_shaft)>(3*sigma_shaft)).*(regression_slope(1,k)*shaft);
                    tem_value2=spine_specific_calcium(:,k);
                    tem_value2(tem_value2<threshold_for_negetive_value)=mean_curr_spine +sigma_curr_spine*randn(sum(tem_value2<threshold_for_negetive_value),1);
                    spine_specific_calcium(:,k)=tem_value2;
                    
                    
                else
                    spine_specific_calcium(:,k) = curr_spine-(regression_slope(1,k)*shaft);
                    tem_value2=spine_specific_calcium(:,k);
                    tem_value2(tem_value2<threshold_for_negetive_value)=mean_curr_spine +sigma_curr_spine*randn(sum(tem_value2<threshold_for_negetive_value),1);
                    spine_specific_calcium(:,k)=tem_value2;
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%robust regression???? %%%%%%%%%%%%%
    if strcmp(Plot1,'Yes')
        figure
        set(gcf,'unit','centimeters','position',[14 1 15 40]);
        pat='.*_';
        m=regexpi(str,pat,'match');
        title_name=m{1}(1:end-1);
        
%         figure;
        
        for i=1:spine_number
            curr_spine = deltaF_F0(:,i);
            if sum(shaft>(3*sigma_shaft))>0
                [Coeff,stats] = robustfit(shaft(:), curr_spine(:)); %robustfit for each spine in a dendrite
                if 1 %stats.p(2)<0.05
                    ax(i)=subplot(ceil(spine_number/3),3,i);
                    %scatter diagram(plot)
                    plot(shaft(:), deltaF_F0(:,i), '.k');
                    hold on;
                    %robust regression plot
                    y = (shaft(:)*regression_slope(1,i))+Coeff(1);
                    plot(shaft(:), y, 'r');
                    xlabel('shaft');
                    ylabel(['spine' num2str(i)]);
                    title(['p=' num2str(stats.p(2))])
                    axis square

                end
            end
        end
        
        if (exist('ax','var'))            
            linkaxes(ax,'xy')
            %         set(gca,'LooseInset',get(gca,'TightInset'))
            print(gcf,[title_name,'_contribution.pdf'],'-dpdf','-bestfit')
        end
    end
    set(gcf,'unit','centimeters','position',[14 1 15 26]);
    clear ax
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    
    for i=1:size(spine_specific_calcium,2)
        [smoothDist,x]=ksdensity(spine_specific_calcium(:,i));
        [~,indPeak]=max(smoothDist);
        
         xFit=[x(1:indPeak) (x(indPeak)+(x(2)-x(1))*[1:(indPeak-1)])];
        dataToFit=[smoothDist(1:indPeak) fliplr(smoothDist(1:(indPeak-1)))];
        [sigma,mu,~]=mygaussfit(xFit',dataToFit,0.1,1);
        

        if mu==-inf
            dataToFit=interp1(xFit,dataToFit,xFit(1):(xFit(end)-xFit(1))/10:xFit(end),'spline');
            [sigma,mu,~]=mygaussfit((xFit(1):(xFit(end)-xFit(1))/10:xFit(end))',dataToFit);
        end
        
        
        spine_specific_calcium(:,i)=spine_specific_calcium(:,i)-mu;
    end
    
    
    
    for i=1:size(data_no_remove,2)
        [smoothDist,x]=ksdensity(data_no_remove(:,i));
        [valuePeak,indPeak]=max(smoothDist);
        
        xFit=[x(1:indPeak) (x(indPeak)+(x(2)-x(1))*[1:(indPeak-1)])];
        dataToFit=[smoothDist(1:indPeak) fliplr(smoothDist(1:(indPeak-1)))];
        [sigma,mu,~]=mygaussfit(xFit',dataToFit,0.1,1);
        
        
        if mu==-inf
            dataToFit=interp1(xFit,dataToFit,xFit(1):(xFit(end)-xFit(1))/10:xFit(end),'spline');
            [sigma,mu,~]=mygaussfit((xFit(1):(xFit(end)-xFit(1))/10:xFit(end))',dataToFit);
        end
        data_no_remove(:,i)=data_no_remove(:,i)-mu;
    end
    
    
    
    
    
   
    [smoothDist,x]=ksdensity(shaft);
    [valuePeak,indPeak]=max(smoothDist);
    
    xFit=[x(1:indPeak) (x(indPeak)+(x(2)-x(1))*[1:(indPeak-1)])];
    dataToFit=[smoothDist(1:indPeak) fliplr(smoothDist(1:(indPeak-1)))];
    [sigma,mu,~]=mygaussfit(xFit',dataToFit,0.1,1);
    
    
    if mu==-inf
        dataToFit=interp1(xFit,dataToFit,xFit(1):(xFit(end)-xFit(1))/10:xFit(end),'spline');
        [sigma,mu,~]=mygaussfit((xFit(1):(xFit(end)-xFit(1))/10:xFit(end))',dataToFit);
    end
    shaft=shaft-mu;
    %%%%
    
    
    %%%%%%%%%%%%%%%  ????threshold_data  %%%%%%%%%%%
    threshold_spine=[];
    for i=1:size(spine_specific_calcium,2)
        [smoothDist,x]=ksdensity(spine_specific_calcium(:,i));
        [valuePeak,indPeak]=max(smoothDist);
        xFit=[x(1:indPeak) (x(indPeak)+(x(2)-x(1))*[1:(indPeak-1)])];
        dataToFit=[smoothDist(1:indPeak) fliplr(smoothDist(1:(indPeak-1)))];
        [sigma,mu,~]=mygaussfit(xFit',dataToFit,0.1,1);
        if mu==-inf
            dataToFit=interp1(xFit,dataToFit,xFit(1):(xFit(end)-xFit(1))/10:xFit(end),'spline');
            [sigma,mu,~]=mygaussfit((xFit(1):(xFit(end)-xFit(1))/10:xFit(end))',dataToFit);
        end
        threshold_spine(i)=max(sigma*threshold,min_threshold);
    end
    %%%% shaft threshold  %%%%
    [smoothDist,x]=ksdensity(shaft);
    [valuePeak,indPeak]=max(smoothDist);
    xFit=[x(1:indPeak) (x(indPeak)+(x(2)-x(1))*[1:(indPeak-1)])];
    dataToFit=[smoothDist(1:indPeak) fliplr(smoothDist(1:(indPeak-1)))];
    [sigma,mu,~]=mygaussfit(xFit',dataToFit,0.1,1);
    if mu==-inf
        dataToFit=interp1(xFit,dataToFit,xFit(1):(xFit(end)-xFit(1))/10:xFit(end),'spline');
        [sigma,mu,~]=mygaussfit((xFit(1):(xFit(end)-xFit(1))/10:xFit(end))',dataToFit);
    end
    threshold_shaft=max(sigma*threshold,min_threshold);
    threshold_all=[threshold_spine threshold_shaft];
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%
    data_all=[spine_specific_calcium shaft];
    %%%%%%%%%%%%%%%%%%%%%%%%%
    data_all_no_remove=[data_no_remove shaft];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    threshold_for_negetive=threshold;
    for i=1:size(data_all,2)
        curr_col=data_all(:,i);
        [smoothDist,x] = ksdensity(curr_col);
        [~,indPeak]=max(smoothDist);
        
        
        xFit=[x(1:indPeak) (x(indPeak)+(x(2)-x(1))*[1:(indPeak-1)])];
        dataToFit=[smoothDist(1:indPeak) fliplr(smoothDist(1:(indPeak-1)))];
  
        
        [sigma_col,mean_col,~]=mygaussfit(xFit',dataToFit,0.1,1);
        if mean_col==-inf
            dataToFit=interp1(xFit,dataToFit,xFit(1):(xFit(end)-xFit(1))/10:xFit(end),'spline');
            [sigma_col,mean_col,~]=mygaussfit((xFit(1):(xFit(end)-xFit(1))/10:xFit(end))',dataToFit);
        end
        curr_col(curr_col<(mean_col-threshold_for_negetive*sigma_col))=mean_col+sigma_col*randn(sum(curr_col<(mean_col-threshold_for_negetive*sigma_col)),1);
        data_all(:,i)=curr_col;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%
    corr_all=zeros(size(data_all,2),size(data_all,2));
    p_all=ones(size(data_all,2),size(data_all,2));
    
    if corr_method==2
        corr_spine_line=0;
        corr_shaft_line=0;
        for p=1:size(data_all,2)
            for q=1:size(data_all,2)
                if sum(data_all(:,p)>threshold_all(p))> thr_length
                    if sum(data_all(:,q)>threshold_all(q))> thr_length
                        [corr_all(p,q) p_all(p,q)]=corr(data_all(((data_all(:,p)>threshold_all(p))+(data_all(:,q)>threshold_all(q)))>0,p),data_all(((data_all(:,p)>threshold_all(p))+(data_all(:,q)>threshold_all(q)))>0,q));
                        %[corr_all(q,p) p_all(q,p)]=corr(data_all(((data_all(:,p)>threshold_all(p))+(data_all(:,q)>threshold_all(q)))>0,p),data_all(((data_all(:,p)>threshold_all(p))+(data_all(:,q)>threshold_all(q)))>0,q));
                        if p==q && p<C
                            corr_spine_line=corr_spine_line+1;
                        end
                        if p==q && p==C
                            corr_shaft_line=1;
                        end
                    end
                end
            end
        end
        
        corr_spine=corr_all(1:(end-1),1:(end-1));
        corr_spine_one=sum(sum(abs(tril(corr_spine))))-corr_spine_line;
        corr_shaft_one=sum(abs(corr_all(end,:)))-corr_shaft_line;
        pat='.*_';
        m=regexpi(str,pat,'match');
        title_name=[m{1}(1:end-1) '_spine_' num2str(corr_spine_one) '_shaft_' num2str(corr_shaft_one) '_' paramater];
        var_name_spine=strcat('spine',string(mat2cell(1:spine_number,1,ones(spine_number,1))));
        var_name_all=[var_name_spine 'shaft'];
        
    elseif corr_method==1
        frames_for_corr=min(corr_time_numFrames,size(spine_specific_calcium,1));
        corr_spine=corr(spine_specific_calcium(1:frames_for_corr,:));
        [corr_all,p_all]=corr(data_all(1:frames_for_corr,:));
        corr_spine_one=sum(sum(abs(tril(corr_spine))))-size(corr_spine,1);
        corr_shaft_one=sum(abs(corr_all(end,:)))-1;
        pat='.*_';
        m=regexpi(str,pat,'match');
        title_name=[m{1}(1:end-1) '_spine_' num2str(corr_spine_one) '_shaft_' num2str(corr_shaft_one) '_' paramater];
        var_name_spine=strcat('spine',string(mat2cell(1:spine_number,1,ones(spine_number,1))));
        var_name_all=[var_name_spine 'shaft'];
        
    elseif corr_method==3
        corr_spine_line=0;
        corr_shaft_line=0;
        if sum(shaft>threshold_shaft)>0
            for p=1:size(data_all,2)
                for q=1:size(data_all,2)
                    
                    [corr_all(p,q) p_all(p,q)]=corr(data_all((shaft>threshold_shaft),p),data_all((shaft>threshold_shaft),q));
                    %[corr_all(q,p) p_all(q,p)]=corr(data_all(((data_all(:,p)>threshold_all(p))+(data_all(:,q)>threshold_all(q)))>0,p),data_all(((data_all(:,p)>threshold_all(p))+(data_all(:,q)>threshold_all(q)))>0,q));
                    if p==q && p<C
                        corr_spine_line=corr_spine_line+1;
                    end
                    if p==q && p==C
                        corr_shaft_line=1;
                    end
                    
                    
                end
            end
        end
        
        corr_spine=corr_all(1:(end-1),1:(end-1));
        corr_spine_one=sum(sum(abs(tril(corr_spine))))-corr_spine_line;
        corr_shaft_one=sum(abs(corr_all(end,:)))-corr_shaft_line;
        pat='.*_';
        m=regexpi(str,pat,'match');
        title_name=[m{1}(1:end-1) '_spine_' num2str(corr_spine_one) '_shaft_' num2str(corr_shaft_one) '_' paramater];
        var_name_spine=strcat('spine',string(mat2cell(1:spine_number,1,ones(spine_number,1))));
        var_name_all=[var_name_spine 'shaft'];
        
        
    end
    
    
    
    
    %%%%%%%%%%  spearman correlation  %%%%%%%%%%%%%%%%
    
    
    corr_all_spearman=zeros(size(data_all,2),size(data_all,2));
    p_all_spearman=ones(size(data_all,2),size(data_all,2));
    
    if corr_method==2
        corr_spine_line=0;
        corr_shaft_line=0;
        for p=1:size(data_all,2)
            for q=1:size(data_all,2)
                if sum(data_all(:,p)>threshold_all(p))> thr_length
                    if sum(data_all(:,q)>threshold_all(q))> thr_length
                        [corr_all_spearman(p,q),p_all_spearman(p,q)]=corr(data_all(((data_all(:,p)>threshold_all(p))+(data_all(:,q)>threshold_all(q)))>0,p),data_all(((data_all(:,p)>threshold_all(p))+(data_all(:,q)>threshold_all(q)))>0,q),'Type','Spearman');
                        %[corr_all(q,p) p_all(q,p)]=corr(data_all(((data_all(:,p)>threshold_all(p))+(data_all(:,q)>threshold_all(q)))>0,p),data_all(((data_all(:,p)>threshold_all(p))+(data_all(:,q)>threshold_all(q)))>0,q));
                        if p==q && p<C
                            corr_spine_line=corr_spine_line+1;
                        end
                        if p==q && p==C
                            corr_shaft_line=1;
                        end
                    end
                end
            end
        end
        

        
    elseif corr_method==1
        frames_for_corr=min(corr_time_numFrames,size(spine_specific_calcium,1));
%         corr_spine=corr(spine_specific_calcium(1:frames_for_corr,:));
        [corr_all_spearman,p_all_spearman]=corr(data_all(1:frames_for_corr,:),'Type','Spearman');

    elseif corr_method==3
        corr_spine_line=0;
        corr_shaft_line=0;
        for p=1:size(data_all,2)
            for q=1:size(data_all,2)
                if sum(data_all(:,p)>threshold_all(p))> thr_length
                    if sum(data_all(:,q)>threshold_all(q))> thr_length
                        [corr_all_spearman(p,q),p_all_spearman(p,q)]=corr(data_all((shaft>threshold_shaft),p),data_all((shaft>threshold_shaft),q),'Type','Spearman');
                        %[corr_all(q,p) p_all(q,p)]=corr(data_all(((data_all(:,p)>threshold_all(p))+(data_all(:,q)>threshold_all(q)))>0,p),data_all(((data_all(:,p)>threshold_all(p))+(data_all(:,q)>threshold_all(q)))>0,q));
                        if p==q && p<C
                            corr_spine_line=corr_spine_line+1;
                        end
                        if p==q && p==C
                            corr_shaft_line=1;
                        end
                    end
                end
            end
        end



    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clf
    for r=1:C
        ax(r)=subplot(C,1,r);
        plot(data(:,r));
        hold on
        plot(smoothBaseline(:,r),'red');
        if r<C
            title(['spine' num2str(r)]);
        else
            title('shaft');
        end
        hold off
    end
    linkaxes(ax,'xy')
    print(gcf,[title_name,'_baseline.pdf'],'-dpdf','-bestfit')
    clf
    clear ax
    for r=1:C
        ax(r)=subplot(C,1,r);
        plot(data_all(:,r));
        hold on
        plot(ones(1,timepoint)*threshold_all(r),'red');
        if r<C
            title(['spine' num2str(r)]);
        else
            title('shaft');
        end
        hold off
    end
    linkaxes(ax,'xy')
    print(gcf,[title_name,'_dff.pdf'],'-dpdf','-bestfit')
    
    
    clf
    clear ax

    for r=1:C
        ax(r)=subplot(C,1,r);
        plot(data_all_no_remove(:,r));
        if r<C
            title(['spine' num2str(r)]);
        else
            title('shaft');
        end
    end
    linkaxes(ax,'xy')
    print(gcf,[title_name,'_no_remove.pdf'],'-dpdf','-bestfit')
    
    
    clf
    clear ax
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    [m, n] = size([threshold_all;data_all]);
    output_cell = mat2cell([threshold_all;data_all], ones(m,1), ones(n,1));  
    
    xlswrite(['Trace',title_name,'.xls'],output_cell,1,'A1');
    
    
    %%%%
    [average_amplitude,frequency,duration,total_activity]=extract_calcium_transient_parameter(data_all,threshold_all,t_per_frame);
    
    data_active_time=ones(1,size(data_all,2));
    data_active_num=ones(1,size(data_all,2));
    for a1=1:(size(data_all,2)-1)
        if sum(data_all(:,a1)>threshold_all(a1))>0
            spine_signal=(data_all(:,a1)>=threshold_all(a1));
            shaft_signal=(data_all(:,end)>=threshold_all(end));
            if sum(shaft_signal)>0
                data_active_time(a1)=sum(spine_signal & shaft_signal)/sum(spine_signal);
                
                active_num_spine=IDXtoINT_ss(double(spine_signal));
                active_num_spine=size(active_num_spine{1},1);
                
                if sum(spine_signal & shaft_signal)>0
                    
                    active_num_spine_shaft=IDXtoINT_ss(double(spine_signal & shaft_signal));
                    active_num_spine_shaft=size(active_num_spine_shaft{1},1);
                    
                    data_active_num(a1)=active_num_spine_shaft/active_num_spine;
                else
                    data_active_num(a1)=0;
                end
                
            else
                data_active_time(a1)=-1;
                data_active_num(a1)=-1;
            end
            
        else
            shaft_signal=(data_all(:,end)>=threshold_all(end));
            if sum(shaft_signal)>0
                data_active_time(a1)=-3;
                data_active_num(a1)=-3;
            else
                data_active_time(a1)=-2;
                data_active_num(a1)=-2;
            end
        end
    end
    
    data_active=[data_active_num;data_active_time];
    
    
    
    

    
    %%%%
    %%%%%%%%%%%%%%%plot calcium transients  %%%%%%%
    title_x=[];
    % Plot2 = questdlg('Plot calcium transients in the trace?', 'Plot?', 'Yes','No','Yes');
    Plot2='No';
    
    if strcmp(Plot2,'Yes')
        raster=double(bsxfun(@gt,spine_specific_calcium,threshold_data));
        frame=1:size(spine_specific_calcium,1);
        for i=1:spine_number
            %subplot(1,spine_number,i);
            figure(i+1);
            plot(frame(raster(:,i)>0),spine_specific_calcium(raster(:,i)>0,i), '.r');
            hold on;
            plot(frame(raster(:,i)==0),spine_specific_calcium(raster(:,i)==0,i), '.k');
            xlabel('frame');
            ylabel('DelataF/F0_spine_specific' );
            
            title_xplus={['spine' num2str(i)]};
            title_x=[title_x ,title_xplus];
        end
    else
        for i=1:spine_number
            title_xplus={['spine' num2str(i)]};
            title_x=[title_x ,title_xplus];
        end
    end
    
    regression_slope=[regression_slope 0];
    regression_slope_pvalue=[regression_slope_pvalue 0];
    output=[average_amplitude;frequency;duration;total_activity;regression_slope;regression_slope_pvalue];
    %output_tem=[output zeros(4,1)];
    output_all=[output;corr_all;p_all;corr_all_spearman;p_all_spearman];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    title_x=[title_x,'shaft'];
    
    %%%%%%%%%%%%%%??????????????????excel????%%%%%%
    [m, n] = size(output_all);
    output_cell = mat2cell(output_all, ones(m,1), ones(n,1));  % ??data??????m*n??cell????
    result = [title_x; output_cell];   % ??????????????????????result
    title_y=[{'Begin'},'average_amplitude','frequency','duration','total_activity','regression_slope','regression_slope_pvalue',title_x,title_x,title_x,title_x]';
    
  
    xlswrite([title_name,'.xls'],result,1,'B1');
    xlswrite([title_name,'.xls'],title_y,1,'A1');
    
    [m_active, n_active] = size(data_active);
    output_cell_active = mat2cell(data_active, ones(m_active,1), ones(n_active,1));  % ??data??????m*n??cell????
    result_active = [title_x; output_cell_active];
    data_active_title=['active' '_' title_name];
    title_y_active=[{'Begin'},'num_percent','time_percent']';
    xlswrite([data_active_title,'.xls'],result_active,1,'B1');
    xlswrite([data_active_title,'.xls'],title_y_active,1,'A1');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    data_all=data_raw_dff;
    data_active_time=ones(1,size(data_all,2));
    data_active_num=ones(1,size(data_all,2));
    for a1=1:(size(data_all,2)-1)
        if sum(data_all(:,a1)>threshold_all(a1))>0
            spine_signal=(data_all(:,a1)>=threshold_all(a1));
            shaft_signal=(data_all(:,end)>=threshold_all(end));
            if sum(shaft_signal)>0
                data_active_time(a1)=sum(spine_signal & shaft_signal)/sum(spine_signal);
                
                active_num_spine=IDXtoINT_ss(double(spine_signal));
                active_num_spine=size(active_num_spine{1},1);
                
                if sum(spine_signal & shaft_signal)>0
                    
                    active_num_spine_shaft=IDXtoINT_ss(double(spine_signal & shaft_signal));
                    active_num_spine_shaft=size(active_num_spine_shaft{1},1);
                    
                    data_active_num(a1)=active_num_spine_shaft/active_num_spine;
                else
                    data_active_num(a1)=0;
                end
                
            else
                data_active_time(a1)=-1;
                data_active_num(a1)=-1;
            end
            
        else
            shaft_signal=(data_all(:,end)>=threshold_all(end));
            if sum(shaft_signal)>0
                data_active_time(a1)=-3;
                data_active_num(a1)=-3;
            else
                data_active_time(a1)=-2;
                data_active_num(a1)=-2;
            end
        end
    end
    
    data_active=[data_active_num;data_active_time];
    
    [m_active, n_active] = size(data_active);
    output_cell_active = mat2cell(data_active, ones(m_active,1), ones(n_active,1));  % ??data??????m*n??cell????
    result_active = [title_x; output_cell_active];
    data_active_title=['active' '_' 'raw' title_name];
    title_y_active=[{'Begin'},'num_percent','time_percent']';
    xlswrite([data_active_title,'.xls'],result_active,1,'B1');
    xlswrite([data_active_title,'.xls'],title_y_active,1,'A1');
    
    clearvars -except Plot1 threshold prctile_num thr_length wdw corr_method shaft_contribution_method corr_time min_threshold max_cc_shaft_spine threshold_for_negetive_value paramater data_num j
    
    
end
clear all
warning('on','MATLAB:polyfit:PolyNotUnique');





function [ INT ] = IDXtoINT_ss( IDX ,numstates)
%IDXtoINT(IDX) Converts state indices to state on/offsets
%
%INPUT
%   IDX:    [t x 1] vector of state indices, where states are identified by
%           integers starting from 1. Times with IDX 0 will not be counted
%           in any interval INT
%   numstates (optional)  number of interval types (for use
%
%OUTPUT
%   INT:    {nstates} cell array of intervals - start and end times
%%

if ~exist('numstates','var')
    numstates = max(IDX);
end

states = 1:numstates;

if isrow(IDX)
    IDX = IDX';
end

IDX = [0; IDX; 0];
for ss = 1:numstates
    statetimes = IDX==states(ss);
    INT{ss} = [find(diff(statetimes)==1) find(diff(statetimes)==-1)-1];
end
end



