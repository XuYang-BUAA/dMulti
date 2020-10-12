function result = templateMatching(sig, templates)
    data = sig.data ;
    [ch_num_r, ch_num_c, data_len] = size(data);
    dt = sEMG.dt;
    width = sig.width;
    
    [max_c, index_c] = max(data,[],2);
    [max_r, index_r] = max(max_c,[],1);
    max_index = zeros(2,data_len);
    max_index(1,:) = index_r;
    for i=1:data_len
        max_index(2,i) = index_c(index_r(:,:,i),:,i);
    end
    [pks,pks_locs] = findpeaks(reshape(max_r(1,1,:),[],1),'MinPeakHeight', threshold,'MinPeakDistance',floor(width/2));
    spikes = sigsegment(data,pks_locs,width);
    timings = pks_locs;
    [spike_num,~]=size(timings);
    
    for n=1:spike_num %matching each spike
        
    end
    