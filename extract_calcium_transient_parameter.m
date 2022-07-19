function [average_amplitude,frequency,duration,total_activity]=extract_calcium_transient_parameter(deltaF_F0,threshold_data,t_per_frame)


raster=double(bsxfun(@gt,deltaF_F0,threshold_data)); 

  %%  amplitude && frequency && duration && accumulative activity
   
  %find peak
    
   for i=1:size(raster,2)  %find parameters for each spine
     if sum(raster(:,i))>0
       
       % transient_inx=find(raster(:,i));
       [startinx,endinx]=find_trasient_interval(raster(:,i));
       n=size(startinx,2); % number of transient intervals for each spine
       timepoint=size(deltaF_F0,1);  
       amplitude=zeros(n,1);
       for j=1:n   %find peak amplitude in each interval for each spine  
          amplitude(j,1)=max(deltaF_F0(startinx(j):endinx(j),i));
       end
          frequency(1,i)=n/(timepoint*t_per_frame)*60;  %%%% one minute
          total_activity(1,i)=sum(deltaF_F0(raster(:,i)>0,i))*(60/size(deltaF_F0,1)); % calculate total_activity
          duration(1,i)=sum(raster(:,i))*t_per_frame/n;
          average_amplitude(1,i)=mean(amplitude,1);
     else
       average_amplitude(1,i)=0;
       total_activity(1,i)=0;
       frequency(1,i)=0;
       duration(1,i)=0;
       
     end
   end
   

 
     
end
  
       

