%==========================================================================
%              demo for decompositon of multi-channel sEMG                *
%                                                                         *
% Read sEMG Data                                                          *
% Format: Multi - Channel sEMG Data should be organized into a M * N      *
% matrix, which contain N channels and M points in each channel.          *
% x = [x_ch1_time1, ... , x_chN_time1;                                    *
%      x_ch1_time2, ... , x_chN_time2;                                    *
%              ..., ... ,         ...;                                    *
%      x_ch1_timeM, ... , x_chN_timeM]                                    *
%                                                                         *
%                                                                         *
%  WARNINGS:                                                              *
%    None                                                                 *
%                                                                         *
%  HISTORY:                                                               *
%    07/08/2020  : XuY create.                                            *
%==========================================================================
clear all;clc;
%close all;
%%
%load('R05M1T25L25_trial1');
load('data_tringle.mat');
%load('decomp_lpy_cycle.mat');
%%
data = cell(40,10);

for i = 1:1:40
    for j = 1:2:19
        %data{i,(j-1)/2+1}=filter(filterBandpass2K,EMG{i,j});
        data{i,(j-1)/2+1}=EMG{i,j};
    end
end
for j=1:10
    for i=1:39
        data{i,j}=data{i,j}-data{i+1,j};
    end
end
for i = 20
    for j = 1:10
        sEMG.data((i-19-1)*10+j,:)=data{i,j};
    end
end

sEMG.dt = 1/2048;
sEMG.t0 = 0;
sEMG.ch = [1,10];
[~,N] = size(sEMG.data);
t = (sEMG.t0:N-1)' * sEMG.dt;

figure
for j=1
    for i=1:sEMG.ch(2)
        select_ch = [j,i];
        plot(t,sEMG.data((j-1)*sEMG.ch(2)+i,:)+300*(i-1))
        hold on
    end
end

figure
for i=1:sEMG.ch(1)
    select_ch = [i,9];
    plot(t,sEMG.data((i-1)*sEMG.ch(2)+select_ch(2),:)+300*i+100*40);
    hold on
end
% for i=1:39
%     select_ch = [i,10];
%     plot(t,sEMG.data{select_ch(1),select_ch(2)}+100*(i+41))
%     hold on
% end
%%
    decomp_result = dMulti(sEMG);
%%
figure
hold on
for i=1:64
    plot(data1(i,:)+5*(i-1))
end
%%
data1 = data(2:65,:);
elearr = [4 5 11 10 24 32 34 39 40 49 50 62 61;
          3 6 12  9 23 31 33 38 48 41 51 63 60;
          2 7 13 17 22 30 27 37 47 42 52 64 59;
          1 8 14 18 21 29 26 36 46 43 53 56 58;
       NaN 16 15 19 20 28 25 35 45 44 54 55 57];
for i=1:5
    for j=1:11
        data_d(i,j,:) = data1(elearr(i,j+2),:)-data1(elearr(i,j+1),:);
    end
end
data_plot = permute(data_d,[2 3 1]);
figure
hold on
for i=1:11
    plot(data_plot(i,:,1)+1*(i-1))
end
%plot(data(10,:))
%% 8*8 array,inter-electrode distance 10mm
elearr = [8 16 24 32 40 48 56 64;
          7 15 23 31 39 47 55 63;
          6 14 22 30 38 46 54 62;
          5 13 21 29 37 45 53 61;
          4 12 20 28 36 44 52 60;
          3 11 19 27 35 43 51 59;
          2 10 18 26 34 42 50 58;
          1 9  17 25 33 41 49 57];
for i=1:7
    for j=1:8
        data_d(i,j,:) = data1(elearr(i,j),:)-data1(elearr(i+1,j),:);
    end
end
data_plot = permute(data_d,[1 2 3]);
figure
hold on
for i=1:8
    plot(data_plot(:,i,4)+1*(i-1))
end
% for j=1:8
%     data_raw_test(j,:) = data1(elearr(4,j),:);
% end
%%
clc;clear all;close all;
%load('F:\CC\Data\CC181026_Pro2_Task3_Trial1_sEMG.mat')
load('F:\CC\Data\CC181026_M2_Trial1_sEMG.mat')
load('F:\CC\Data\CC181026_M2_Trial1_decomp_grid_elec5_epoch30_diff1_bip0.mat')
%data1 = Data_sEMG(1:64,:);
data1 = Data_sEMG(65:128,:);
%%
sEMG.dt = 1/2048/2;
sEMG.t0 = 0;
sEMG.ch = [1,8];
%sEMG.data = data_plot(:,:,:)';
data_f_d(1,5,:) = 0.01;
data_f_d_u = upsampleSig(data_f_d);
sEMG.data = data_f_d_u;
%sEMG.data = data_f;

%sEMG.data=data1(elearr(4,:),:);
%data_re = reshape(data_plot,size(data_plot,1),[]);

% for i =1:7
%     data_re(:,(i-1)*8+1:i*8) = data_plot(:,:,i);
% end
%sEMG.data = data_re(:,17:40)';
%sEMG.data = data_raw_test;
%a = reshape(data_plot,length(data_plot(:,1,1)),[]);
%sEMG.data = a';
figure
hold on
for i=1:8
    plot(sEMG.data(i,:)+1*(i-1))
end
clear filter_test;
filter_test = filter(filter2k,data1(:,:));
filter_test = filter(notching50,filter_test);
for i=1:8
    smoothed(i,:) = smooth(1:61440,sEMG.data(i,:),0.00015,'loess');
end

%figure
hold on
for i=1:8
    plot(smoothed(i,:)+1*(i-1))
end
sEMG.data = smoothed;

%%

filter_test = filter(filter2k,data1')';
figure
hold on
for i=1:8
    for j=1:8
        
        plot(data1((i-1)*8+j,:) + (i-1)*8+j);
        plot(filter_test((i-1)*8+j,:) + (i-1)*8+j,'r');
    end
end
    
clear data_d_f
for i=1:7
    for j=1:8
        data_d_f(i,j,:) = filter_test(elearr(i,j),:)-filter_test(elearr(i+1,j),:);
        %data_d(i,j,:) = data1(elearr(i,j),:)-data1(elearr(i+1,j),:);
    end
end

for i=1:8
    for j=1:8
        data_f(i,j,:) = filter_test(elearr(i,j),:);
        data_raw(i,j,:) = data1(elearr(i,j),:);
    end
end


figure
hold on
for i=1:7
    for j=1:8
        %plot(reshape(data_d_f(i,j,:),[],1) + (i-1)*8+j);
        %plot(reshape(data_d(i,j,:),[],1) + (i-1)*8+j);
        %plot((j-1)*70000:(j-1)*70000+61439,reshape(data_d(i,j,:),[],1)+1.5*i)
        %plot((j-1)*70000:(j-1)*70000+61439,reshape(data_d_f(i,j,:),[],1)+1.5*i)
        
        %plot((i-1)*70000:(i-1)*70000+61439,reshape(filter_test((i-1)*8+j,:),[],1)+1.5*j)
        
         %plot((j-1)*70000:(j-1)*70000+61439,reshape(data_raw(i,j,:),[],1)+1.5*i,'r')
         %plot((j-1)*70000:(j-1)*70000+61439,reshape(data_f(i,j,:),[],1)+1.5*i,'b')
         %plot((j-1)*70000:(j-1)*70000+61439,reshape(data_raw_d(i,j,:),[],1)+1.5*i,'b')
         %plot((j-1)*70000:(j-1)*70000+61439,reshape(data_f_d(i,j,:),[],1)+1.5*i,'r')
         plot((j-1)*130000:(j-1)*130000+122878,reshape(data_f_d_u(i,j,:),[],1)+1.5*i,'r')
         %plot((j-1)*70000:(j-1)*70000+61439,reshape(smoothed(i,j,:),[],1)+1.5*i,'r')
        %plot((j-1)*70000:(j-1)*70000+61439,reshape(data_smooth_d(i,j,:),[],1)+1.5*i,'r')
    end
end

for i=1:8
    for j=1:8
        data_f(i,j,:) = filter_test((i-1)*8+j,:);
    end
end

clear smoothed
for i=1:8
    for j=1:8
        smoothed(i,j,:) = smooth(1:61440, data_f(i,j,:),0.00015,'loess');
    end
end

for i=1:7
    for j=1:8
        %data_smooth_d(i,j,:) = smoothed(i+1,j,:) - smoothed(i,j,:);
        data_f_d(i,j,:) = data_f(i+1,j,:) - data_f(i,j,:);
        data_raw_d(i,j,:)=data_raw(i+1,j,:) - data_raw(i,j,:);
    end
end


sEMG.data = data_f(2:7,2:7,:);

%%

    width = 0.015;
    threshold = 0.1;
    dr_coe = 3;
    [spikes, templates_result] = getTemplates(sEMG,width,threshold,dr_coe);
    
    upsample_temp = upsampleTemp(templates_result.templates,0.2);
    for i=1:7
        warped_temp(:,:,:,:,i) = warp_template(upsample_temp, templates_result.width, 0.2+0.2*i, 1);
        %plot(1:templates_result.width,warped_temp(:,1,4,4,i))
        %hold on
    end
    templates_result.warped = warped_temp;   
    
    plotTemplates(spikes, templates_result);
   
%%
    A=ones(2,2,5,2);
    for i=1:5
        for j=1:2
            A(:,:,i,j)=j*i*A(:,:,i,j);
        end
    end
    A = permute(A,[3,4,1,2]);
    S2 = A.^2;
    S = sqrt(sum(S2,1));
    W = A./S;
    W=permute(W,[3,4,1,2]);





















