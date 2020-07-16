function [spikes,timings] = detectSpikes(sEMG,threshold)
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
    threshold = 25;
    width = goodwidth(0.020/dt);
    %=================================
    [max_vector(1,:),max_vector(2,:)]=max(abs(data));
    [pks,pks_locs] = findpeaks(max_vector(1,:),'MinPeakHeight', threshold,'MinPeakDistance',floor(width/2));
    
    %peak_vector = [max_vector(:,pks_locs);pks_locs];
    spikes = sigsegment(data',pks_locs,width);%width * peaks * channels
    timings = pks_locs;
    [~,spike_num]=size(timings);
    spikes_con = reshape(permute(spikes,[1 3 2]),[],spike_num);%concatenate spikes by channels 
    %%test pca
    [feature_P,feature_Z] = featurePCA(spikes_con,width,10);
    
%     dist_test = zeros(spike_num*(spike_num-1)/2,3);
%     k=1;
%     for i =1:(spike_num-1)
%         for j = (i+1):spike_num
%             dist(k,1) = i;
%             dist(k,2) = j;
%             dist(k,3) = pdist2(feature_P(i,:),feature_P(j,:));
%             k = k+1;
%         end
%     end
    S = floor(0.05*spike_num);
%     dist = pdist2(feature_P,feature_P,'euclidean','Smallest',S);
%     rho = mean(dist);
    dist = pdist2(feature_P,feature_P);
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
    cluster_ind = find(delta./rho>5);
    template_num = length(cluster_ind);
    spike_cluster = zeros(1,spike_num);
    
    templates = zeros(width*ch_num,template_num);
    [~,nearest] = sort(dist(:,cluster_ind));
    for i = 1:length(cluster_ind)
        templates(:,i) = mean(spikes_con(:,nearest(1:10,i)),2);
    end
    for i=1:spike_num
        [~,spike_cluster(i)] = min(dist(i,cluster_ind));
    end
    
    %spike_feature = [timings ; feature_P'];
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   