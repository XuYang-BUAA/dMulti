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
    width = goodwidth(0.015/dt);
    %=================================
    [max_vector(1,:),max_vector(2,:)]=max(abs(data));
    [pks,pks_locs] = findpeaks(max_vector(1,:),'MinPeakHeight', threshold,'MinPeakDistance',floor(width/2));
    
    %peak_vector = [max_vector(:,pks_locs);pks_locs];
    spikes = sigsegment(data,pks_locs,width);
    timings = pks_locs;
    
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   