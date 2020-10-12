function [spikes, results] = getTemplates(sEMG,width,threshold,dr_coe)
%==========================================================================
%                        get templates                                    *
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
%    10/10/2020  : XuY create.                                            *
%==========================================================================
    %%

    data = sEMG.data;
    [ch_num_r, ch_num_c, data_len] = size(data);
    dt = sEMG.dt;
    
    
    %=======test parameter================
    
    width = goodwidth(width/dt);%spike width
    temp_width = goodwidth(0.015/dt);
    %threshold = .1;
    
    
    %=================================
    %[max_c, index_c] = max(abs(data),[],2);
    [max_c, index_c] = max(data,[],2);
    [max_r, index_r] = max(max_c,[],1);

    max_index = zeros(2,data_len);
    max_index(1,:) = index_r;
    for i=1:data_len
        max_index(2,i) = index_c(index_r(:,:,i),:,i);
    end
    
    [pks,pks_locs] = findpeaks(reshape(max_r(1,1,:),[],1),'MinPeakHeight', threshold,'MinPeakDistance',floor(temp_width/2));
    
    
    spikes_all_test = sigsegment(data,pks_locs,width);
    
    
    spikes_temp_all_test = sigsegment(data,pks_locs,temp_width);
    timings_all = pks_locs;
    [spike_num_all,~]=size(timings_all);
    
    spikes.spikes = spikes_temp_all_test;
    spikes.timings = timings_all;

    %%
    result_cluster_ind = zeros(1,0);
    result_templates = zeros(temp_width,0,ch_num_r,ch_num_c);
    result_nearest = zeros(20,0);
    for c = 1:ch_num_r
        k=0;
        pks_ch = zeros(1,0);
        locs_ch = zeros(0,1);
        for i=1:length(pks_locs)
            %if(max_index(2,pks_locs(i)) == c && max_index(1,pks_locs(i))==sel_test)
            if(max_index(2,pks_locs(i)) == c )
            %if(i<=length(pks_locs))
                k=k+1;
                pks_ch(k) = pks_locs(i);
                locs_ch(k) = i;
            end
        end
        if(length(pks_ch)<20)
            continue;
        end

        
        %spikes_test = sigsegment(data(:,select{c},:),pks_ch,width);
        spikes_test = sigsegment(data(:,:,:),pks_ch,width);
        spikes_test = spikes_test .*hanning(width);
        spikes_temp = sigsegment(data,pks_ch,temp_width);
        timings = pks_ch;
        [~,spike_num]=size(timings);
        
        
%         spikes_hanning = spikes.*hanning(width);
%         spikes_temp_hanning = spikes_temp.*hanning(temp_width);
        
    %%
        dist = zeros(spike_num,spike_num);
        for i =1:spike_num
            for j=1:spike_num
                clear sel_row sel_col;
                %[~,~,dist(i,j),~]=caldiff(spikes_con(:,i),spikes_con(:,j));
                [sel_row,sel_col] = chooseChannel(max_index(1,i), max_index(2,i), max_index(1,j), max_index(2,j), ch_num_r, ch_num_c);
                %[~,~,dist(i,j),~]=caldiff(spikes_test(:,:,:,i),spikes_test(:,:,:,j));
                %[~,~,dist(i,j),~]=caldiff(spikes_test(:,sel_row,sel_col,i),spikes_test(:,sel_row,sel_col,j));
                [~,~,dist(i,j),~]=caldiff(spikes_test(:,:,sel_col,i),spikes_test(:,:,sel_col,j));
                %[~,~,dist(i,j),~]=caldiff(spikes_test(:,sel_test,:,i),spikes_test(:,sel_test,:,j));
                %dist(i,j) = dtw(spikes_con(:,i)',spikes_con(:,j)',3);
                %[~,~,dist(i,j),~]=caldiff(spikes_con_hanning(:,i),spikes_con_hanning(:,j));
                %dist(i,j) = norm(spikes_con(:,i)-spikes_con(:,j));
            end
%             clc
%             i
            
        end
        %%
        dist_test = zeros(spike_num,spike_num);
        for i = 1:spike_num
            for j = 1:spike_num
                dist_test(i,j) = max(dist(i,j),dist(j,i));
            end
        end
        dist = dist_test;
        %dist = pdist2(spikes_con',spikes_con');
        %dist = dist./A;
        %%
        %S=floor(0.02*spike_num_all);
        S=8;
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
        dr_threshold = mean(dr_sort(1:floor(length(dr)/2)))*2;

        %dr_threshold = 1.4;

        cluster_ind = find(delta./rho>dr_threshold);
        template_num = length(cluster_ind);
        spike_cluster = zeros(1,spike_num);

        templates = zeros(temp_width,ch_num_r,ch_num_c,template_num);
        [~,nearest] = sort(dist(:,cluster_ind));
        
        for i = 1:length(cluster_ind)
            %templates(:,:,:,i) = mean(spikes_temp_all_test(:,:,:,locs_ch(nearest(1:S,i))).*hanning(temp_width),4);%hanning
            templates(:,:,:,i)=mean(spikes_temp_all_test(:,:,:,locs_ch(nearest(1:S,i))),4);%without hanning
 
        end
        for i=1:spike_num
            [~,spike_cluster(i)] = min(dist(i,cluster_ind));
        end

        %spike_feature = [timings ; feature_P'];
%         clear sig
%         for j=1:ch_r
%             for i=1:ch_c
%                 sig.data(i,:) = data((j-1)*ch_c+i,:);
%             end
%         end
%         sig.data = sig.data';
%         %sig.data = max_vector(1,:)';
%         sig.dt = sEMG.dt;
%         sig.t0 = sEMG.t0;
%         figure();
%         hold on
%         h = markpeaks(sig, pks_locs);
%     %     
    
        result_cluster_ind = [result_cluster_ind,cluster_ind];
        result_templates = [result_templates,permute(templates,[1 4 2 3])];
        result_nearest = [result_nearest,reshape(locs_ch(nearest(1:20,:)),20,[])];
    end
    result_template_num = length(result_cluster_ind);
%     [energy,en_locs] = sort(sum(result_templates.^2),'descend');
%     result_templates = result_templates(:,en_locs);
%     result_nearest = result_nearest(:,en_locs);
%     result_cluster_ind = result_cluster_ind(:,en_locs);
	results.templates = result_templates;
    results.a(1,:) = sqrt(sum(result_templates.^2,[1,3,4]));
    results.a(2,:) = 1.1*results.a(1,:);
    results.a(3,:) = 0.9*results.a(1,:);
    results.templates_w = result_templates./results.a;
    results.template_num = result_template_num;
    results.nearest = result_nearest;
    results.cluster_ind = result_cluster_ind;
    results.width = temp_width;
    
    
    
    
    
    
    
    
    
    
    
   