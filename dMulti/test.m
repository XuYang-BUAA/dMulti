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
    plot(spikes_con(:,i))%,'color','b');
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
%% plot templates
%figure
[~,temp_len] = size(templates);
for i=1:temp_len
    %subplot(temp_len, 1, i)
    plot(templates(:,i)+i-1, '-', 'LineWidth', 1.5, 'color', [0, 0, 0]/255)
    hold on 
end
for i=1:temp_len
    %subplot(temp_len, 1, i)
    for j=1:ch_num
        line([width*j,width*j],[-1,temp_len],'linestyle','--');
        hold on
    end
end
%%
n_r = 15*8;
figure
hold on
for i=1:1

    %subplot(temp_len, 1, i)
    for j = 1:7
       
        plot(templates(((j-1)*n_r+1):j*n_r,i)+j-1,  'color', [0, 0, i*30]/255)
        hold on 
    end
    for j=1:8
        line([width*j,width*j],[-1,7],'linestyle','--');
        hold on
    end
end
%% test H
figure
inter = 10;
sub_num = 8;
for i =1:sub_num
    subplot(sub_num,1,i)
    plot(test_H(:,i*inter));
end

%% test fillZero
H_test = fillZero(spikes_con, width, ch_num, 'Spike');
H_temp = fillZero(templates,width, ch_num,'Template');
%%
for i =1:temp_len
    value_test(:,i) = H_temp(:,:,i)'*H_test(:,50)/(H_temp(:,:,i)'*H_temp(:,:,i));
end
[value_max,delay]=max(value_test);
figure

plot(H_temp(:,delay(6),6),'r')
hold on
plot(H_test(:,50));

%plot(H_test(:,50) - H_temp(:,delay,8),'k')
    
%% test dtw
%close all
clear ix iy

%test_seg = width*3+1:width*4;
test_seg = 1:size(templates,1);
figure
test_spike = circshift(spikes_con_temp(:,1057),0);
test_spike_han = spikes_con(:,1057);
%plot(test_spike)
hold on
plot(test_spike(test_seg));
for j=1:ch_num
    line([width*j,width*j],[-0.5,1],'linestyle','--');
    hold on
end
ch_num = 1;
% for j=1:ch_num
%     line([width*j,width*j],[-0.5,0.5],'linestyle','--');
%     hold on
% end
stop =1;
for i = 3
    plot(templates(test_seg,i),'k')
    for j=1:8
        %figure
        %[dist(j,:),ix(j,:),iy(j,:)]=myDTW(test_spike(width*(j-1)+1:width*j),templates(width*(j-1)+1:width*j,i));
        %[dist(j,:),ix(j,:),iy(j,:)]=myDTW_new(test_spike,templates(:,i),width,ch_num);
        %[dist(i),~,~,sub]=myDTW_new(test_spike,templates(:,i),width,ch_num);
        %myDTW_new(test_spike,templates(:,i),width,ch_num);
        figure
        test_seg = width*(j-1)+1:width*j;
        test_dtw(j) = myDTW_new(test_spike(test_seg),templates(test_seg,i),width,ch_num);
        
        %[dist(j,i),~,~]=myDTW(test_spike(width*(j-1)+1:width*j),templates(width*(j-1)+1:width*j,i));
        %[dist(j,:,i),~,~]=dtw(test_spike((2*width-1)*(ch_num-1)+1:(2*width-1)*(ch_num-1)+3*width-2),templates(width*(ch_num-1)+1:width*ch_num,i));
        %myDTW(test_spike(width*(j-1)+1:width*j),templates(width*(j-1)+1:width*j,i));
        %dtw(test_spike,templates(:,i),width,ch_num);
    end
end

%dist = sum(dist);
%%
figure
ch_sel = 1;
plot(test_spike(width*(ch_sel-1)+1:width*ch_sel));
hold on
plot(templates(width*(ch_sel-1)+1:width*ch_sel,3)+0.5);
[~,x_len] = size(ix);
for i=1:x_len
     line([ix(ch_sel,i),iy(ch_sel,i)],[test_spike(ix(ch_sel,i)+width*(ch_sel-1)),templates(iy(ch_sel,i)+width*(ch_sel-1),3)+0.5])
     hold on
end
% ali_temp = 
    
%%
close all
figure
q = spikes_con(:,40);
c = spikes_con(:,606);
plot(q)
hold on
plot(c)
for j=1:ch_num
    line([width*j,width*j],[-0.5,1],'linestyle','--');
    hold on
end

%[dist,~,~] = 
myDTW_new(q,c,width,ch_num);
%%
figure
for j = 1:template_num
    for i = 1:10
        %plot(spikes_con(:,nearest(i,j))+j-1)
        plot(spikes_con_temp(:,nearest(i,j))+j-1)
        hold on
    end
end
% for j=1:ch_num
%     line([width*j,width*j],[-1,template_num],'linestyle','--');
%     hold on
% end  
for j=1:ch_num
    line([width*j,width*j],[-1,template_num],'linestyle','--');
    hold on
end  
    
    
    
%%
figure
for j = 1:8
	plot(data(j,:)+j-1)
	hold on
end
%% test laplace

%test_laplace = data1(37,:)*4 - data1(36,:) - data1(38,:)- data1(29,:)- data1(45,:) ;
test_d = data1(37,:)-data1(38,:);
figure
plot(test_laplace)
hold on
%figure
plot(data1(37,:))



















    