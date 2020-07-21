function [dist,ix,iy] = myDTW(x,y,varargin) 
%DTW    Distance between signals via Dynamic Time Warping
%   DIST = DTW(X,Y) computes the total distance, DIST, as the minimum sum
%   of Euclidean distances between the samples of vectors X and Y, where
%   samples in either X or Y may repeat consecutively any number of times.
%   If X and Y are matrices, then X and Y must have the same number of rows
%   and DTW minimizes the total Euclidean distances between the column
%   vectors of X and Y, allowing columns of X and Y to consecutively
%   repeat.
%
%   [DIST,IX,IY] = DTW(X,Y) additionally returns the warping path, IX and
%   IY, that minimizes the total Euclidean distance between X(IX) and Y(IY)
%   when X and Y are vectors and between X(:,IX) and Y(:,IY) when X and Y
%   are matrices.
%
%   [DIST,IX,IY] = DTW(X,Y,MAXSAMP) additionally restricts IX and IY so
%   that they must be within MAXSAMP samples of a straight-line fit between
%   X and Y.  If MAXSAMP is unspecified, then no restriction will be placed
%   upon IX or IY.
%
%   [DIST,IX,IY] = DTW(X,Y,...,METRIC) will return in DIST the summed
%   distance between the corresponding entries of X and Y according to the
%   distance metric defined by METRIC.  The default is 'euclidean'.
%      'absolute'  the sum of the absolute (manhattan) differences
%      'euclidean' the root sum squared differences
%      'squared'   the squared Euclidean distance
%      'symmkl'    the symmetric Kullback-Leibler distance 
%   When data and signal are matrices, the distances are taken between
%   corresponding column vectors; when data and signal are vectors,
%   distances are taken between the corresponding elements.  Note that when
%   data and signal are vectors, 'absolute' and 'euclidean' are equivalent.
%
%   DTW(...) without output arguments plots the original and aligned
%   signals.  DTW displays the alignment of X and Y via a line plot when
%   they are vectors, and via horizontally aligned images when they are
%   matrices.  If the matrices are complex, the real and imaginary
%   portions appear in the top and bottom half of each image, respectively.
%   
%   % Example 1:
%   %   Compute and plot the best Euclidean distance between real
%   %   chirp and sinusoidal signals using dynamic time warping.
%   x = chirp(0:999,0,1000,1/100);
%   y = cos(2*pi*5*(0:199)/200);
%   dtw(x,y)
%
%   % Example 2:
%   %   Compute and plot the best Euclidean distance between complex
%   %   chirp and sinusoidal signals using dynamic time warping.
%   x = exp(2i*pi*(3*(1:1000)/1000).^2);
%   y = exp(2i*pi*9*(1:399)/400);
%   dtw(x,y)
%
%   % Example 3:
%   %   Align handwriting samples along the x-axis.
%   load blockletterex
%   dtw(MATLAB1,MATLAB2);
%
%   See also EDR, ALIGNSIGNALS, FINDSIGNAL, FINDDELAY, XCORR.

%   References: 
%   * H. Sakoe and S. Chiba, "Dynamic Programming Algorithm Optimization
%     for Spoken Word Recognition" IEEE Transactions on Acoustics, Speech
%     and Signal Processing, Vol. ASSP-26, No. 1, Feb 1978, pp. 43-49.
%   * K.K. Paliwal, A. Agarwal and S.S. Sinha, "A Modification over Sakoe
%     and Chiba's dynamic time warping algorithm for isolated word
%     recognition" IEEE International Conference on ICASSP 1982., Vol. 7,
%     pp. 1259-1261

%   Copyright 2015 The MathWorks, Inc.

needsTranspose = iscolumn(x);

if iscolumn(x) && isvector(y)
  x = x.';
end

if iscolumn(y) && isvector(x)
  y = y.';
end

validateattributes(x,{'single','double'},{'nonnan','2d','finite'},'dtw','X',1);
validateattributes(y,{'single','double'},{'nonnan','2d','finite'},'dtw','Y',2);

if size(x,1) ~= size(y,1)
  error(message('signal:dtw:RowMismatch'));
end

if ~isempty(x) && ~isempty(y)
  C = myCumulativeDistance(x, y);
  dist=C(size(x,2),size(y,2));
  [ix,iy] = traceback(C);
else
  dist = NaN;
  ix = zeros(1,0);
  iy = zeros(1,0);
end
metric = 'euclidean';
if nargout==0
  dtwplot(x, y, ix, iy, dist, metric)
elseif needsTranspose
  ix = ix';
  iy = iy';
end

%-------------------------------------------------------------------------
function C = myCumulativeDistance(x, y)
    m = length(x);
    n = length(y);
    C = zeros(m,n); %m*n 
    

    d = d_ddtw(x,y,m,n);
    %d = d_sdtw(x,y,m,n);
    %d = d_dtw(x,y,m,n);
    
    %  symmetic1
    %  o - 1 - *
    %        / |
    %      1   1
    %    /     |
    %  o       o
    C(1,1) = d(1,1);
    for i = 2:m
        C(i,1) = C(i-1,1) + d(i,1);
    end
    for j = 2:n
        C(1,j) = C(1,j-1) + d(1,j);
    end
    for i =2:m
        for j=2:n
            C1 = C(i-1,j) + d(i,j);
            C2 = C(i-1,j-1) + d(i,j);
            C3 = C(i,j-1) + d(i,j);
            C(i,j) = min([C1,C2,C3]);
            C(i,j) = min([C2,C3]);
        end
    end
            
    a=1;
%--------------------------------------------------------------------------
%dtw euclidean distance
function d = d_dtw(x,y,m,n)
    d = zeros(m,n);
    for i = 1:m
        for j=1:n
            d(i,j) = abs(x(i)-y(j));
        end
    end

%--------------------------------------------------------------------------
%%ddtw D = {[q(i) - q(i-1)]+[q(i+1)-q(i-1)/2]}/2
function d = d_ddtw(x,y,m,n)
    d = zeros(m,n);
    Dx = zeros(m,1);
    Dy = zeros(n,1);
    Dx(1) = x(2)-x(1);
    Dx(m) = x(m)-x(m-1);
    Dy(1) = y(2)-y(1);
    Dy(n) = y(n)-y(n-1);
    for i = 2:m-1
        Dx(i) = ( ( x(i)-x(i-1) )+( x(i+1)-x(i-1)/2 ) )/2;
    end
    for i = 2:n-1
        Dy(i) = ( ( y(i)-y(i-1) )+( y(i+1)-y(i-1)/2 ) )/2;
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
    k=floor(m/10);
    sx = zeros(m,2*k+1);
    sy = zeros(n,2*k+1);
    for i=1:m
        sx(i,k+1) = x(i);
        for j=1:k
            if(i-j<=0)
                sx(i,k+1-j)=sx(i,k+2-j); 
            else
                sx(i,k+1-j) = x(i-j);
            end
            if(i+j>m)
                sx(i,k+1+j)=sx(i,k+j);
            else
                sx(i,k+1+j)=x(i+j);
            end
        end
    end
    for i=1:n
        sy(i,k+1) = y(i);
        for j=1:k
            if(i-j<=0)
                sy(i,k+1-j)=sy(i,k+2-j);
            else
                sy(i,k+1-j) = y(i-j);
            end
            if(i+j>n)
                sy(i,k+1+j)=sy(i,k+j);
            else
                sy(i,k+1+j)=y(i+j);
            end
        end
    end
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

while i>1 || j>1
  if j == 1
    i = i-1;
  elseif i == 1
    j = j-1;
  else
    % trace back to the origin, ignoring any NaN value
    % prefer i in a tie between i and j
    cij = C(i-1,j-1);
    ci = C(i-1,j);
    cj = C(i,j-1);
    i = i - (ci<=cj | cij<=cj | cj~=cj);
    j = j - (cj<ci | cij<=ci | ci~=ci);
  end
  k = k+1;
  ix(k) = i;
  iy(k) = j;
end

ix = ix(k:-1:1);
iy = iy(k:-1:1);