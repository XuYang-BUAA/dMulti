function [row, col] = chooseChannel(r1, c1, r2, c2, ch_num_r, ch_num_c)

    row = [];
    col = [];
    %row = channelIndex(r1, ch_num_r, channelIndex(r2, ch_num_r, row));
    row = channelIndex(r1, ch_num_r, row);
    col = channelIndex(c1, ch_num_c, channelIndex(c2, ch_num_c, col));
    row = sort(unique(row));
    col = sort(unique(col));
    
    
    
    
    function index = channelIndex(x, ch_num, index)
        if(x>1)
            index = [index, x-1, x];
        else
            index = [index, x];
        end
        if(x<ch_num)
            index = [index, x+1];
        end

        
    