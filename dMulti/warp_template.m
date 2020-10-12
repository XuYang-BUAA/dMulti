function warped_temp = warp_template(template, len, k, a)
    
    [l, num,r, c] = size(template);
    warped_temp = zeros(len,num,r,c);
    v_s = mean(template(1:10,:,:,:));
    v_e = mean(template(l-9:l,:,:,:));
    step = (len-1)/(l-1);
    for t_n = 1:num
        for i=1:r
            for j=1:c
                for m=1:len
                    mid = (len+1)/2;
                    n = round( ( ( (m-mid)/k+mid )-1)/step )+1;
                    try
                        warped_temp(m,t_n,i,j) = a * template(n, t_n,i, j);
                    catch
                        if(n<1)
                            warped_temp(m,t_n,i,j) = v_s(1,t_n,i,j);
                        else 
                            warped_temp(m,t_n,i,j) = v_e(1,t_n,i,j);
                        end
                    end
                end
            end
        end
    end
    