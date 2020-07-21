function templates = genTemplates(spikes_con, ch_num, width, num_compoents)
%==========================================================================
%                        generate templates                               *
%                                                                         *
% INPUT:                                                                  *
%    spikes_con      -- concatenation spikes                              *
%    width           -- spike width                                       *
%    num_compoents   -- num for ICA compoents                             *
%                                                                         *
% OUTPUT:                                                                 *
%    templates       -- MU templates                                      *
%                                                                         *
%                                                                         *
%                                                                         *
%  WARNINGS:                                                              *
%    None                                                                 *
%                                                                         *
%  HISTORY:                                                               *
%    07/08/2020  : XuY create.                                            *
%==========================================================================

    [~,spike_num] = size(spikes_con);
    [feature_P,feature_Z] = featurePCA(spikes_con,width,num_compoents);

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
    dr = delta./rho;%cal delta/rho rate
    dr_sort = sort(dr);
    dr_threshold = mean(dr_sort(1:floor(length(dr)/2)))*4;
    cluster_ind = find(delta./rho>dr_threshold);
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