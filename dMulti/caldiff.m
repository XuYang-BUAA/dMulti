function [K1sqr, K2sqr, K3sqr, A] = caldiff(x, y)
% Calculate the degree of difference between candi & template.
% If candi(Mc x Nc) or template(Mt x Nt) is matrix, calculate on column, 
% and the result is Nc x Nt.

% Denote candi, template as vector x, y.
% A = x' * y/(y' * y);
% K1 = |A - 1|
% K2 = sqrt(|x|^2 / |y|^2 - A^2);
% K3 = sqrt(K2^2 + (A - 1)^2).
% Return value K1^2, K2^2, K3^2.
%==========================================================================
%                        cal distance                                     *
%                                                                         *
% INPUT:                                                                  *
%    x              -- length * row * column                              *
%    y              -- length * row * column                              *
%                                                                         *
% OUTPUT:                                                                 *
%    K1sqr          -- K1^2                                               *
%    K2sqr          -- K2^2                                               *
%    K3sqr          -- K3^2                                               *
%    A              -- A                                                  *
%                                                                         *
%  WARNINGS:                                                              *
%    None                                                                 *
%                                                                         *
%  HISTORY:                                                               *
%    07/08/2020  : XuY create.                                            *
%==========================================================================
    candi = reshape(x,[],1);
    template = reshape(y,[],1);
    [cnd_len, cnd_num] = size(candi);
    [tmp_len, ~] = size(template);
    %[cnd_len, cnd_r, cnd_c] = size(candi);
    %[tmp_len, ~, ~] = size(template);
    
    if (cnd_len ~= tmp_len)
        error('Vectors should have the same length!');
    end
    tmp2 = sum(template.^2);
    cnd2 = sum(candi.^2);
    A = candi' * template ./ repmat(tmp2, cnd_num, 1);
    K1sqr = abs(A - 1).^2;
    K2sqr = cnd2' * (1./tmp2) - A.^2;
    %K2sqr = abs(cnd2'- A.^2);
    K3sqr = K2sqr + K1sqr;
    
end

%------------------------------EOF-----------------------------------------

