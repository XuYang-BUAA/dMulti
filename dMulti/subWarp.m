function [residual,c_warp] = subWarp(q,c,ix,iy,width,ch_num)

    y=c;
    trace_len = length(ix);
	c_warp = zeros(width,ch_num);
    c_temp = c(iy,:);
    cx = 1;
    i=1;
    while(i<trace_len)
        k=0;
        try
            while(ix(i+k+1) == ix(i))
                k=k+1;
            end
        catch
        end
        c_warp(cx,:) = mean(c_temp(i:i+k,:),1);
        i=i+k+1;
        cx=cx+1;
    end
        
        
    
    residual = q - c_warp;
%     figure
%     hold on
%     for tc=1:8
%         plot(q(:,tc)+0.5*(tc-1),'b');
%         plot(c(:,tc)+0.5*(tc-1),'k');
%         plot(c_warp(:,tc)+0.5*(tc-1),'r');
%     end
    stop = 1;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    