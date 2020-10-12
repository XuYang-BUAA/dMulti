%function [dist,ix,iy,sub] = myDTW_new(q,c,width,ch_num) 
function [dist, c_wrap] = myDTW_new(q,c,width,ch_num,select_ch) 
%   References: 
%   Matlab dtw.m

x = reshape(q,width,ch_num);
y = reshape(c,width,ch_num);

x_c = x(:,select_ch);%.*hanning(width);
y_c = y(:,select_ch);
width_h = floor(width/5);
% x_c = x;
% y_c = y;

validateattributes(x,{'single','double'},{'nonnan','2d','finite'},'dtw','X',1);
validateattributes(y,{'single','double'},{'nonnan','2d','finite'},'dtw','Y',2);


% if ~isempty(x) && ~isempty(y)
%   C = myCumulativeDistance(x, y,width,ch_num);
%   dist=C(size(x,1),size(y,1));
%   [ix,iy] = traceback(C);
if ~isempty(x) && ~isempty(y)
  C = myCumulativeDistance(x_c, y_c,width,ch_num);
  dist=C(size(x,1),size(y,1));
  [ix,iy] = traceback(C);

  %test
  [x_len,~] = size(ix);
%   figure
% 
%   plot(x(:,2),'linewidth',5)
%   hold on
%   plot(y(:,2)+0.5,'linewidth',5)
%   for i=1:x_len
%      line([ix(i),iy(i)],[x(ix(i),2),y(iy(i),2)+0.5])
%      hold on
%   end
   trace_len = length(ix);
  c_wrap = zeros(width,ch_num);
  c_wrap(1,:) = y(1,:);
  i=2;
  while(i<=trace_len)
      k = 0;
      try
          while(ix(i+k+1)==ix(i)&&i+k<trace_len-1)
              k=k+1;
          end
      catch
      end
%       if(k>1)
%           stop=0;
%       end
    if(ix(i)>1)
      	c_wrap(ix(i),:) = mean(y(iy(i):iy(i+k),:),1);
    end
      i=i+k+1;
  end
  i=0;
  while(i<=trace_len)
      k = 0;
      try
          if(iy(i+1)==iy(i)&&i<trace_len-1)
              c_wrap(iy(i),:) = mean(c_wrap(iy(i-1),:),c_wrap(iy(i+1),:),1);
          end
      catch
      end
%       if(k>1)
%           stop=0;
%       end
      
      i=i+1;
  end 
  i=2;
  while(i<length(c_wrap)-1)
      if(c_wrap(i)==c_wrap(i-1))
          c_wrap(i) = (c_wrap(i-1)+c_wrap(i+1))/2;
      end
      i=i+1;
  end
%   
  
  %c_wrap(width,:) = y(width,:);
  %y(:,select_ch) = c_wrap(:,select_ch);
  c_test = reshape(c_wrap,[],1);
  %c_test = reshape(y,[],1);
%   figure
%       hold on
%  plot(q,'k')
%  plot(c)
%  plot(c_test,'color','r')
%     for j=1:ch_num
%         line([width*j,width*j],[-0.5,.5],'linestyle','--');
%         hold on
%     end
%sub = x(:,center-1:center+1).*hanning(width) - c_wrap(:,center-1:center+1);
sub = x - c_wrap;
%sub = x - c_wrap;
%dist = sum(sum(sub(:,select_ch).^2)./sum(c_wrap(:,select_ch).^2));

%dist = mean(sum(x_c.*hanning(width) .* c_wrap(:,select_ch)./sqrt(sum(c_wrap(:,select_ch).^2))));
stop = 1;
%plot(sub,'linewidth',1)
% figure
% plot(q);
% hold on
% plot(c+0.5);
% x_len = length(ix);
% for i=1:x_len
%      line([ix(i),iy(i)],[q(ix(i)),c(iy(i))+0.5])
%      hold on
% end

   
%   for i =1:ch_num-1
%         ix = [ix;ix(1:trace_len)+i*width];
%         iy = [iy;iy(1:trace_len)+i*width];
%   end
    
    
else
  dist = NaN;
  ix = zeros(1,0);
  iy = zeros(1,0);
end





metric = 'euclidean';
if nargout==0
  dtwplot(q, c, ix, iy, dist, metric)
% elseif needsTranspose
%   ix = ix';
%   iy = iy';
end

%--------------------------------------------------------------------------
function C = myCumulativeDistance(x, y, width, ch_num)
    m = width;
    n = width;
    width_h = floor(width/2);
    width_h =1;
    C = zeros(m,n); %m*n 
    

    d = d_ddtw(x,y,m,n);
    %d = d_sdtw(x,y,m,n);
    %d = d_dtw(x,y,m,n);
%--------------------------------------------------------------------------
    %  symmetic1
    %  o       o
    %  |     / 
    %  1   1   
    %  | /     
    %  * - 1 - o     
%     C(1,1) = d(1,1);
%     for i = 2:m
%         C(i,1) = C(i-1,1) + d(i,1);
%     end
%     for j = 2:n
%         C(1,j) = C(1,j-1) + d(1,j);
%     end
%     for i =2:m
%         for j=2:n
%                 if(j<=i+width_h && j>=i-width_h)
%                     C1 = C(i-1,j) + d(i,j);
%                     C2 = C(i-1,j-1) + d(i,j);
%                     C3 = C(i,j-1) + d(i,j);
%                     C(i,j) = min([C1,C2,C3]);
%                 end
%         end
%     end
%--------------------------------------------------------------------------
    %       o
    %     / o o
    %    | / /
    %    * -
    C = 1./zeros(m,n);
    C(1,1) = d(1,1);
    C(2,2) = C(1,1)+d(2,2);
    C(3,2) = C(1,1)+d(3,2);
    C(2,3) = C(1,1)+d(2,3);
    for i = 3:m
        for j = 3:n
            if(j<=i+width_h && j>=i-width_h)
                C1 = C(i-1,j-2) + d(i,j);
                C2 = C(i-2,j-1) + d(i,j);
                C3 = C(i-1,j-1) + d(i,j);
                C(i,j) = min([C1,C2,C3]);
            end
        end
    end
    a=1;
%--------------------------------------------------------------------------
%dtw euclidean distance
function d = d_dtw(x,y,m,n)
    d = zeros(m,n);
    for i = 1:m
        for j=1:n
            d(i,j) = mean( abs(x(i,:)-y(j,:)) );
        end
    end

%--------------------------------------------------------------------------
%%ddtw D = {[q(i) - q(i-1)]+[q(i+1)-q(i-1)/2]}/2
function d = d_ddtw(x,y,m,n)
    d = zeros(m,n);
    Dx = zeros(m,1);
    Dy = zeros(n,1);
    Dx(1) = mean(x(2,:)-x(1,:));
    Dx(m) = mean(x(m,:)-x(m-1,:));
    Dy(1) = mean(y(2,:)-y(1,:));
    Dy(n) = mean(y(n,:)-y(n-1,:));
    for i = 2:m-1
        Dx(i) = mean( ( ( x(i,:)-x(i-1,:) )+( x(i+1,:)-x(i-1,:)/2 ) )/2 );
    end
    for i = 2:n-1
        Dy(i) = mean( ( ( y(i,:)-y(i-1,:) )+( y(i+1,:)-y(i-1,:)/2 ) )/2 );
    end
    for i = 1:m
        for j=1:n
            d(i,j) = abs(Dx(i)-Dy(j));
        end
    end
%--------------------------------------------------------------------------
%%sdtw  sL = x_len / k;
function d = d_sdtw(x,y,m,n)
    d = zeros(m,n);
    k=floor(m/5);
    k=1;
    sx = zeros(m,2*k+1);
    sy = zeros(n,2*k+1);
    for i=1:m
        sx(i,k+1) = x(i);
        sy(i,k+1) = y(i);
        for j=1:k
            if(i-j<=0)
                sx(i,k+1-j)=sx(i,k+2-j); 
                sy(i,k+1-j)=sy(i,k+2-j);
            else
                sx(i,k+1-j) = x(i-j);
                sy(i,k+1-j) = y(i-j);
            end
            if(i+j>m)
                sx(i,k+1+j)=sx(i,k+j);
                sy(i,k+1+j)=sy(i,k+j);
            else
                sx(i,k+1+j)=x(i+j);
                sy(i,k+1+j)=y(i+j);
            end
        end
    end
%     for i=1:n
%         sy(i,k+1) = y(i);
%         for j=1:k
%             if(i-j<=0)
%                 sy(i,k+1-j)=sy(i,k+2-j);
%             else
%                 sy(i,k+1-j) = y(i-j);
%             end
%             if(i+j>n)
%                 sy(i,k+1+j)=sy(i,k+j);
%             else
%                 sy(i,k+1+j)=y(i+j);
%             end
%         end
%     end
    for i = 1:m
        for j=1:n
            d(i,j) = pdist2(sx(i,:),sy(j,:));
        end
    end
    
%-------------------------------------------------------------------------
function [ix,iy] = traceback(C)
m = size(C,1);
n = size(C,2);

% pre-allocate to the maximum warping path size.
ix = zeros(m+n,1);
iy = zeros(m+n,1);

ix(1) = m;
iy(1) = n;

i = m;
j = n;
k = 1;
%--------------------------------------------------------------------------
%  o       o
%  |     / 
%  1   1   
%  | /     
%  * - 1 - o  
% while i>1 || j>1
%   if j == 1
%     i = i-1;
%   elseif i == 1
%     j = j-1;
%   else
%     % trace back to the origin, ignoring any NaN value
%     % prefer i in a tie between i and j
%     cij = C(i-1,j-1);
%     ci = C(i-1,j);
%     cj = C(i,j-1);
%     i = i - (ci<=cj | cij<=cj | cj~=cj);
%     j = j - (cj<ci | cij<=ci | ci~=ci);
%   end
%   k = k+1;
%   ix(k) = i;
%   iy(k) = j;
% end

%--------------------------------------------------------------------------
%       o
%     / o o
%    | / /
%    * -
while i>2 && j>2
    i=i-1;
    j=j-1;
    k=k+1;
    ix(k) = i;
    iy(k) = j;
    
    c1 = C(i,j);
    c2 = C(i-1,j);
    c3 = C(i,j-1);
    index1 = c2<=c1&&c2<=c3;
    index2 = c3<=c1&&c3<=c2; 
    if(index1 || index2)
        i = i - index1;
        j = j - index2;
        k = k+1;
        ix(k) = i;
        iy(k) = j;
    end
end
if(i==2&&j==2)
    k=k+1;
    ix(k) = 1;
    iy(k) = 1;
else
    k=k+1;
    ix(k) = i-1;
    iy(k) = j-1;
    k=k+1;
    ix(k) = 1;
    iy(k) = 1;
end
%--------------------------------------------------------------------------
ix = ix(k:-1:1);
iy = iy(k:-1:1);












