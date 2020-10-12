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

%% test H
figure
inter = 10;
sub_num = 8;
for i =1:sub_num
    subplot(sub_num,1,i)
    plot(test_H(:,i*inter));
end

%% test fillZero
% H_test = fillZero(spikes_con_temp_all, temp_width, ch_num, 'Spike');
H_temp = fillZero(result_templates(temp_width+1:temp_width*4,:),temp_width, 3,'Template');

H_test_seg = fillZero(test_spike(temp_width+1:temp_width*4,1), temp_width,3,'Spike');
%H_test_seg = fillZero(test_spike, temp_width,ch_num,'Spike');

%%
for i =1:template_num
    %value_test(:,i) = H_temp(:,:,i)'*H_test(:,500);%/(H_temp(:,:,i)'*H_temp(:,:,i));
    value_test(:,i) = H_temp(:,:,i)'*H_test_seg;
end
[value_max,delay]=max(value_test);
figure

plot(H_temp(:,delay(4),4),'r')
hold on
%plot(H_test(:,500));
plot(H_test_seg);
    
%% test dtw
%close all
clear ix iy

%test_seg = width*3+1:width*4;
test_seg = 1:size(templates,1);
figure
test_spike = circshift(spikes_con_temp_all(:,501),0);
test_spike_han = spikes_con_all(:,501);
%plot(test_spike)
hold on
plot(test_spike(test_seg));
for j=1:ch_num
    line([temp_width*j,temp_width*j],[-0.5,1],'linestyle','--');
    hold on
end
ch_num = 8;
% for j=1:ch_num
%     line([width*j,width*j],[-0.5,0.5],'linestyle','--');
%     hold on
% end
plot(result_templates(:,3))
%plot(test_spike - result_templates(:,3))
stop =1;
warp = zeros(temp_width*ch_num,0);
for i = 1:result_template_num
%for i = 6
    %plot(result_templates(test_seg,i),'k')
    
    template_ch = reshape(result_templates(:,i),[],ch_num);
    [~,select_ch] = sort(sum(template_ch.^2),'descend');
    for j=5
        %myDtwCenter(test_spike(temp_width*(j-1)+1:temp_width*j), result_templates(temp_width*(j-1)+1:temp_width*j,i),temp_width,1,center);
        [dist_temp(i),warp(:,i)] = myDtwCenter(test_spike, result_templates(:,i),temp_width,ch_num,6);%select_ch(1:5));
        %P(i) = calProb(result_templates(:,i),test_spike);
        %figure
        %[dist(j,:),ix(j,:),iy(j,:)]=myDTW(test_spike(width*(j-1)+1:width*j),templates(width*(j-1)+1:width*j,i));
        %[dist(j,:),ix(j,:),iy(j,:)]=myDTW_new(test_spike,templates(:,i),width,ch_num);
        %[dist(i),~,~,sub]=myDTW_new(test_spike,templates(:,i),width,ch_num);
        %myDTW_new(test_spike,templates(:,i),width,ch_num);
%         figure
%         test_seg = width*(j-1)+1:width*j;
%         test_dtw(j) = myDTW_new(test_spike(test_seg),templates(test_seg,i),width,ch_num);
        
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
    
%% test nearest
figure
for j = 1:result_template_num
    for i = 1:S
        %plot(spikes_con(:,nearest(i,j))+j-1)
        plot(spikes_con_temp_all(:,result_nearest(i,j))+j-1)
        %plot(spikes_con_hanning_all(:,result_nearest(i,j))+j-1)
        hold on
    end
end
for j=1:ch_num
    line([temp_width*j,temp_width*j],[-1,result_template_num],'linestyle','--');
    hold on
end  
%% plot templates
%figure
%[~,temp_len] = size(templates);
for i=1:result_template_num
    %subplot(temp_len, 1, i)
    plot(result_templates(:,i)+i-1, '-', 'LineWidth', 1.5, 'color', [0, 0, 0]/255)
    %plot(H_temp(:,i)+i-1, '-', 'LineWidth', 1.5, 'color', [0, 0, 0]/255)
    hold on 
end
for i=1:result_template_num
    %subplot(temp_len, 1, i)
    for j=1:ch_num
        line([temp_width*j,temp_width*j],[-1,result_template_num],'linestyle','--');
        hold on
    end
end    
    
%% test laplace

%test_laplace = data1(37,:)*4 - data1(36,:) - data1(38,:)- data1(29,:)- data1(45,:) ;
test_d = data1(37,:)-data1(38,:);
figure
plot(test_laplace)
hold on
%figure
plot(data1(37,:))

%% plot firing
figure
hold on
for i=1:result_template_num
    i
    for j=1:data_len
        if(result_timings(i,j) ==1)
            line([j,j],[(i-1),(i-1)+0.5])
        end
    end
end

%% test gaussian kernel

%% test all channels
channels = data_plot(:,:,:);
clear spikes_all_ch templates_ch
for i=1:7
    spikes_ch_temp = sigsegment(channels(:,:,i),pks_locs,temp_width);
    spikes_ch_con = reshape(permute(spikes_ch_temp,[1 3 2]),[],spike_num_all);
    spikes_all_ch(:,:,i) = spikes_ch_con;
end
for j=1:7
    for i = 1:13
        templates_ch(:,i,j) = mean(spikes_all_ch(:,result_nearest(1:S,i),j),2);
    end
end
ch_row = 7;

    for i=1:result_template_num
    %for i =4
        figure
        for j=1:ch_row
            for k =1 : 8
                plot([(k-1)*40+3:k*40-3],templates_ch((k-1)*35+1:k*35 ,i,j)+j-1, '-', 'LineWidth', 1.5, 'color', [0, 0, 0]/255)
                hold on
            end
        end
    end  


%% test peaks
    figure
    hold on
    for i=1:7
        for j=1:8
            plot((j-1)*40:(j-1)*40+34,spikes_temp_all_test(:,i,j,602)+i)
        end
    end
    figure 
    hold on
    for i=1:4
        for j=1:8
            plot((j-1)*70000:(j-1)*70000+61439,reshape(sEMG.data(i,j,:),[],1)+1.5*i)
        end
    end

%%
    for n=1:result_template_num
        figure
        hold on
        for m=1:8
            for i=1:7
                for j=1:8
                    %plot((j-1)*40:(j-1)*40+34,spikes_temp(:,i+1,j,nearest(m,n))-spikes_temp(:,i,j,nearest(m,n))+i,'b')
                    %plot((j-1)*40:(j-1)*40+34,spikes_temp(:,i,j,nearest(m,n))+i,'b')
                    plot((j-1)*40:(j-1)*40+34,spikes_temp_all_test(:,i,j,result_nearest(m,n))+0.5*i,'b')
                    if(m==8)
                        plot((j-1)*40:(j-1)*40+34,result_templates(:,n,i,j)+0.5*i,'k')
                    end
                end
            end
        end
    end
%%
    sel_piece = 2;
    [~,pulses_num] = size(DecompResult.Pulses{1,sel_piece});
    figure
    hold on
    for i = 1:pulses_num
        for j=1:length(DecompResult.Pulses{1,sel_piece}{1,i})
            line([DecompResult.Pulses{1,sel_piece}{1,i}(1,j) DecompResult.Pulses{1,sel_piece}{1,i}(1,j)],[0+i-1,0.5+i-1])
        end
    end

%% 
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        




    