function upsample_temp = upsampleTemp( templates, step )
    [l, r, c, n] = size(templates);
    upsample_temp = zeros( (l-1)/step + 1, r, c, n);
    x = 1:l;
    xi = 1:step:l;
    for i=1:r
        for j=1:c
            for k=1:n
                upsample_temp(:,i,j,k) = interp1(x,templates(:,i,j,k),xi,'spline');
            end
        end
    end