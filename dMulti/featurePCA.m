function [P,Z] = featurePCA(X,width,num_compoents)
% load hald
% X = ingredients;

% 计算原始信号的特征均值
% mu = mean(X, 1);
% 
% %PCA降维，维度为2
% [P, Z]  = pca(X, 'NumComponents', num_compoents);
%  %[P, Z]  = pca(X);
% % 信号重构
% try
%     X_re = Z*P'+mu;
% catch
%     X_re = Z*P'+repmat(mu, size(X, 1), 1);
% end

    [seg_len,spike_num] = size(X);
    ch_num = seg_len / width;
    P = zeros(spike_num,num_compoents,ch_num);
    Z = zeros(seg_len,num_compoents);
    mu = zeros(ch_num,spike_num);
    for i = 1:ch_num
        [P(:,:,i),Z((i-1)*width+1:i*width,:)]= pca(X((i-1)*width+1:i*width,:), 'NumComponents', num_compoents);
        mu(i,:) = mean(X((i-1)*width+1:i*width,:), 1);
    end
    % 信号重构
    X_re = zeros(seg_len,spike_num);
    for i = 1:ch_num
        X_re((i-1)*width+1:i*width,:) = Z((i-1)*width+1:i*width,:)*P(:,:,i)'+repmat(mu(i,:),width,1);
    end
    
    %P_test = reshape(permute(P,[3 1 2]),[],spike_num);
    P = reshape(P,spike_num,[]);

%  可视化结果
figure
%set(gcf, 'position', [100 100 600 600])
for i = 11:35
    subplot(5, 5, i-10)
    hold on
    grid on
    plot(X(:, i), '-', 'LineWidth', 1.5, 'color', [254, 67, 101]/255)
    plot(X_re(:, i), '-', 'LineWidth', 1.5, 'color', [0, 114, 189]/255)
    for j=1:ch_num
        line([width*j,width*j],[-40,40],'linestyle','--');
    end
    %legend('原始信号', '重构信号')
    title(['第 ',num2str(i),' 个特征'])

end