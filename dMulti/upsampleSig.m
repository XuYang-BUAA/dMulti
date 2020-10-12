function upsample_sig = upsampleSig( sig )
    [r, c, l] = size(sig);
    step = 0.5;
    upsample_sig = zeros( r, c, (l-1)/step + 1);
    x = 1:l;
    xi = 1:step:l;
    for i=1:r
        for j=1:c
                upsample_sig(i,j,:) = interp1(x,reshape(sig(i,j,:),[],1),xi,'spline');
        end
    end