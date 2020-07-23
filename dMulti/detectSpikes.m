function [spikes_con,timings, spike_num] = detectSpikes(sEMG,width,threshold)
%==========================================================================
%                        detect spikes                                    *
%                                                                         *
% INPUT:                                                                  *
%    sEMG            -- filtered multi-channel sEMG data                  *
%    threshold       -- threshold of peak detection                       *
%                                                                         *
% OUTPUT:                                                                 *
%    spikes          -- spikes for clustering                             *
%    timings         -- spike locations                                   *
%                                                                         *
%                                                                         *
%                                                                         *
%  WARNINGS:                                                              *
%    None                                                                 *
%                                                                         *
%  HISTORY:                                                               *
%    07/08/2020  : XuY create.                                            *
%==========================================================================
    %%

    data = sEMG.data;
    [ch_num,data_len]=size(data);
    ch_r = sEMG.ch(1);
    ch_c = sEMG.ch(2);
    dt = sEMG.dt;
    
    
    %=======test parameter================
    
    width = goodwidth(0.02/dt);%spike width
    temp_width = goodwidth(0.02/dt);
    threshold = 0.07;
    %=================================
    [max_vector(1,:),max_vector(2,:)]=max(abs(data));
    [pks,pks_locs] = findpeaks(max_vector(1,:),'MinPeakHeight', threshold,'MinPeakDistance',floor(temp_width/2));
    
    %peak_vector = [max_vector(:,pks_locs);pks_locs];
    spikes = sigsegment(data',pks_locs,width);%width * peaks * channels
    spikes_temp = sigsegment(data',pks_locs,temp_width);
    %spikes = sigsegment(data',pks_locs,15);
    spikes = spikes.*hanning(width);
    timings = pks_locs;
    [~,spike_num]=size(timings);
    
    spikes_con = reshape(permute(spikes,[1 3 2]),[],spike_num);%concatenate spikes by channels 
    spikes_con_temp = reshape(permute(spikes_temp,[1 3 2]),[],spike_num);
    %spikes_hanning = reshape(permute(spikes.*hanning(width),[1 3 2]),[],spike_num);
    %spikes_con = spikes_hanning;
    %%test pca
    %[feature_P,feature_Z] = featurePCA(spikes_con,width,5);
    
   % dist_test = zeros(spike_num*(spike_num-1)/2,3);
%     k=1;
%     for i =1:(spike_num-1)
%         for j = (i+1):spike_num
%             dist(k,1) = i;
%             dist(k,2) = j;
%             dist(k,3) = pdist2(feature_P(i,:),feature_P(j,:));
%             k = k+1;
%         end
%     end
    S = floor(0.01*spike_num);
%     dist = pdist2(feature_P,feature_P,'euclidean','Smallest',S);
%     rho = mean(dist);
%    dist = pdist2(feature_P,feature_P);
%     dist = zeros(spike_num,spike_num);
%     spikes = permute(spikes,[1 3 2]);
%     for i=1:spike_num
%         for j=1:spike_num
%             %[dist(i,j),~,~] = myDTW_new(spikes_con(:,i),spikes_con(:,j),width,ch_num);
%             %[dist(i,j),~,~] = dtw(spikes_con(:,i),spikes_con(:,j));
%             [dist(i,j),~,~] = dtw(spikes(:,:,i),spikes(:,:,j));
%         end
%         num = i
%     end
%%
    for i =1:spike_num
        for j=1:spike_num
            [~,~,dist(i,j),~]=caldiff(spikes_con(:,i),spikes_con(:,j));
        end
    end
    %%
    for i = 1:spike_num
        for j = 1:spike_num
            dist_test(i,j) = max(dist(i,j),dist(j,i));
        end
    end
    dist = dist_test;
    %dist = pdist2(spikes_con',spikes_con');
    %dist = dist./A;
    %%
    S=floor(0.01*spike_num);
    dist_sort = sort(dist);
    rho = mean(dist_sort(1:S,:));
    delta = zeros(1,spike_num);

    for i = 1:spike_num
        try
            delta(i) = min(dist(i,find(rho<rho(i))));
        catch
            delta(i) = max(dist(:,i));
        end
    end
    dr = delta./rho;%cal delta/rho rate
    dr_sort = sort(dr);
    dr_threshold = mean(dr_sort(1:floor(length(dr)/2)))*3;
    
    dr_threshold = 1.6;
    
    cluster_ind = find(delta./rho>dr_threshold);
    template_num = length(cluster_ind);
    spike_cluster = zeros(1,spike_num);
    
    templates = zeros(width*ch_num,template_num);
    [~,nearest] = sort(dist(:,cluster_ind));
    for i = 1:length(cluster_ind)
        %templates(:,i) = mean(spikes_con(:,nearest(1:10,i)),2);
        templates(:,i) = mean(spikes_con_temp(:,nearest(1:10,i)),2);
    end
    for i=1:spike_num
        [~,spike_cluster(i)] = min(dist(i,cluster_ind));
    end
    
    %spike_feature = [timings ; feature_P'];
%     clear sig
%     for j=1:ch_r
%         for i=1:ch_c
%             select_ch = [j,i];
%             sig.data(i,:) = data((j-1)*ch_c+i,:);
%         end
%     end
%     sig.data = sig.data';
%     %sig.data = max_vector(1,:)';
%     sig.dt = sEMG.dt;
%     sig.t0 = sEMG.t0;
%     figure();
%     hold on
%     h = markpeaks(sig, pks_locs);
%     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   