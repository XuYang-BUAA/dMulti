function [decomp_result] = dMulti(sEMG)
%==========================================================================
%              top of decompositon of multi-channel sEMG                  *
%                                                                         *
% INPUT:                                                                  *
%    sEMG            -- filtered multi-channel sEMG data                  *
%                                                                         *
% OUTPUT:                                                                 *
%    decomp_result   -- cell for the decompsition result                  *
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
% ========== Data information ========================================-----
     data = sEMG.data;
%     % fdata = doFilter10KHz(data);
%     fdata = data;
     t0 = sEMG.t0; dt = sEMG.dt;
% %     %ch_num = sEMG.ch_num;
% %     [~,data_len] = size(sEMG.data{1,1});
% %     % Time series
% %     t = t0 + dt * (0:data_len-1);
%     data = sEMG.data;
    [ch_num,data_len]=size(data);
%     ch_r = sEMG.ch(1);
%     ch_c = sEMG.ch(2);
    width = goodwidth(0.020/dt);
    threshold = 5;
    num_compoents = 10;
    
    [spikes, timings, spike_num]= detectSpikes(sEMG,threshold,width);
    templates = genTemplates(spikes, ch_num, width, num_compoents);
    
    decomp_result = [];
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    