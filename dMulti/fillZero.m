function H = fillZero(spikes, width, ch_num, varargin)
%==========================================================================
%                        fill templates with zero-padding                 *
%   Gligorijevi I , Johannes P. van Dijk¡­.                               *
%   A new and fast approach towards sEMG decomposition[J].                *
%   Medical & Biological Engineering & Computing, 2013, 51(5):593-605.    *
%                                                                         *
% INPUT:                                                                  *
%    templates       -- concatenation templates                           *
%    width           -- templates width                                   *
%    ch_num          -- channel numbers                                   *
%                                                                         *
% OUTPUT:                                                                 *
%    H               -- templates with shifted zero padding               *
%                                                                         *
%  WARNINGS:                                                              *
%    None                                                                 *
%                                                                         *
%  HISTORY:                                                               *
%    07/18/2020  : XuY created.                                           *
%==========================================================================
    if nargin == 4  
        op_fill_mode = varargin{1};
    else 
        of_fill_mode = 'Spike';
    end
    %%
    [spike_len,spike_num] = size(spikes);
    switch op_fill_mode
        case 'Template'
            
            H = zeros(((2*width-1)*ch_num+width-1),2*width-1,spike_num);
            for i=1:spike_num
                for j=1:ch_num
                    H(((2*width-1)*(j-1)+1):((2*width-1)*(j-1)+width),1,i) ...
                        = spikes((width*(j-1)+1):width*j,i);
                end
            end
            for i = 2:2*width-1
                H(:,i,:) = circshift(H(:,1,:),i-1);
            end
        case 'Spike'
            H = zeros(((2*width-1)*ch_num+width-1),spike_num);
            for i = 1:spike_num
                for j=1:ch_num
                    H(((2*width-1)*(j-1)+width):((2*width-1)*(j-1)+2*width-1),i) ...
                        = spikes((width*(j-1)+1):width*j,i);
                end
            end
    end
    
    


            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            