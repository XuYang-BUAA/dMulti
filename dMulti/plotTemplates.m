function plotTemplates(spikes, templates)
    [width, num, ch_num_r, ch_num_c] = size(templates.templates);
    dc=hsv;
%     for n=1:num
%         figure
%         hold on
%         for m=1:20
%             for i=1:ch_num_r
%                 for j=1:ch_num_c
%                     %plot((j-1)*40:(j-1)*40+34,spikes_temp(:,i+1,j,nearest(m,n))-spikes_temp(:,i,j,nearest(m,n))+i,'b')
%                     %plot((j-1)*40:(j-1)*40+34,spikes_temp(:,i,j,nearest(m,n))+i,'b')
%                     plot((j-1)*(width+5):(j-1)*(width+5)+(width-1),spikes.spikes(:,i,j,templates.nearest(m,n))+0.5*i,'color',dc(20+2*m,:))
%                     if(m==8)
%                         plot((j-1)*(width+5):(j-1)*(width+5)+(width-1),templates.templates(:,n,i,j)+0.5*i,'k')
%                         plot((j-1)*(width+5):(j-1)*(width+5)+(width-1),templates.templates_w(:,n,i,j)+0.5*i,'r')
%                     end
%                 end
%             end
%         end
%     end
    for n=1:num
        figure
        hold on
        for m=1:7
            for i=1:ch_num_r
                    for j=1:ch_num_c
                        plot((j-1)*(width+5):(j-1)*(width+5)+(width-1),templates.warped(:,n,i,j,m)+0.5*i,'color',dc(9*m,:))
                    end
            end
        end
    end