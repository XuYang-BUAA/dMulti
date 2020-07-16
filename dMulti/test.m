clc;
M_test = zeros(20,5);
for i=1:20
    M_test(i,:)=i;
end
dist = pdist2(M_test(1,:),M_test);
test_dist = pdist2(M_test,M_test);
test_min = pdist2(M_test,M_test,'euclidean','Smallest',5);
test_mean = mean(test_min);
%%
figure
for i =1:spike_num
    subplot(5, 1, spike_cluster(i))
    hold on
    grid on
    %plot(spikes_con(:,i), '-', 'LineWidth', 1.5, 'color', [254, 67, 101]/255)
    plot(spikes_con(:,i),'.-' , 'LineWidth', 0.5);
    %plot(X_re(:, i), '-', 'LineWidth', 1.5, 'color', [0, 114, 189]/255)

    %legend('原始信号', '重构信号')
    %title(['第 ',num2str(x(i)),' 个特征'])

end
for i=1:5
    subplot(5, 1, i)
    for j=1:ch_num
        line([width*j,width*j],[-200,200],'linestyle','--');
        hold on
    end
end
%%
for i=1:5
    subplot(5, 1, i)
    plot(templates(:,i), '-', 'LineWidth', 1.5, 'color', [0, 0, 0]/255)
    hold on 
end
    
    
    