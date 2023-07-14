Rsq_sum = 0;
idx = 0;
figure(201)
set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[17 1 17 12]);
for i = 1:n_pools
    
    for k = 1:n_draws
        for j = 4:size(data.tInt,3)
%             plot(time_soln(k, :, i, j),[limb_angles(k*10-9:k*10-5,:,i, j)],'color',[.17 .17 .17],'LineWidth',1.5)
%             h1 = plot(data.tInt(1, :, j),data.xInt(1:5,:, j),'r--','LineWidth',1.1);
            subplot(5,n_pools,i)
            hold on
            title(['Simulation ' num2str(i)])
            plot(time_soln(k, :, i, j),limb_angles(k*10-7,:,i, j),'LineWidth',1+.1,'Color',[.7 .7 .7]-.05*j);
%             plot(time_soln(k, :, i, j),limb_angles(k*10-7,:,i, j),'LineWidth',1+.1*i_v,'Color',[.7 .7 .7]);
%             plot(data.tInt(1, :, j),data.xInt(3,:, j),'LineWidth',1,'Color',[.7 .7 .7]-.05*j);
            plot(data.tInt(1, :, j),data.xInt(3,:, j),'--','LineWidth',1.3,'Color',[.2 .2 .2]);
%             
            
            Rsq_tmp = fitlm(data.xInt(1,:, j), limb_angles(k*10-9,:,i, j)).Rsquared.Adjusted + ...
                        fitlm(data.xInt(2,:, j), limb_angles(k*10-8,:,i, j)).Rsquared.Adjusted + ...
                        fitlm(data.xInt(3,:, j), limb_angles(k*10-7,:,i, j)).Rsquared.Adjusted + ...
                        fitlm(data.xInt(4,:, j), limb_angles(k*10-6,:,i, j)).Rsquared.Adjusted + ...
                        fitlm(data.xInt(5,:, j), limb_angles(k*10-5,:,i, j)).Rsquared.Adjusted ;
            
            idx = idx + 1;
            Rsq_sum = (Rsq_sum + Rsq_tmp ) ;
        
        
        subplot(5,n_pools,[3+i,6+i])
        hold on
%         plot(time_soln(k, :, i, j),limb_angles(k*10-9,:,i, j),'LineWidth',1+.1*i_v,'Color',[1-.1*i_v 0 1-.1*i_v]);
%         plot(time_soln(k, :, i, j),limb_angles(k*10-6,:,i, j),'LineWidth',1+.1*i_v,'Color',[0 0 1-.1*i_v]);
        plot(time_soln(k, :, i, j),limb_angles(k*10-8,:,i, j),'LineWidth',1.1,'Color',[1-.1*j 0 0]);
        plot(time_soln(k, :, i, j),limb_angles(k*10-6,:,i, j),'LineWidth',1.1,'Color',[0 0 1-.1*j]);
        plot(data.tInt(1, :, j),data.xInt(4,:, j),'--','LineWidth',1.3,'Color',[.2 .2 .2]);
        plot(data.tInt(1, :, j),data.xInt(2,:, j),'--','LineWidth',1.3,'Color',[.2 .2 .2]);

        
        
        subplot(5,n_pools,[9+i,12+i])
        hold on
        plot(time_soln(k, :, i, j),limb_angles(k*10-9,:,i, j),'LineWidth',1.1,'Color',[1-.1*j 0 1-.1*j]);
        plot(time_soln(k, :, i, j),limb_angles(k*10-5,:,i, j),'LineWidth',1.1,'Color',[0 1-.1*j 1-.1*j]);
%         plot(time_soln(k, :, i, j),limb_angles(k*10-8,:,i, j),'LineWidth',1+.1*i_v,'Color',[1 0 1]);
%         plot(time_soln(k, :, i, j),limb_angles(k*10-5,:,i, j),'LineWidth',1+.1*i_v,'Color',[0 1 1]);
        plot(data.tInt(1, :, j),data.xInt(1,:, j),'--','LineWidth',1.3,'Color',[.2 .2 .2]);
        plot(data.tInt(1, :, j),data.xInt(5,:, j),'--','LineWidth',1.3,'Color',[.2 .2 .2]);

    end
%     for k = 2:4
%          plot(time,position_noise(:,k),'r--','LineWidth',2);
%     end

%     set(gca,'fontsize',10)
    xlabel('Time (s)')
%     ylabel('Angle (rad)')
%     box off    
    

    end
    subplot(5,n_pools,1)
    ylabel('Torso Angle (rad)')
    subplot(5,n_pools,3)
    legend('Torso Angles', 'Actual Trajectories')
    legend boxoff

    subplot(5,n_pools,[3+1,6+1])
    ylabel('Femur Angle (rad)')
    subplot(5,n_pools,[3+3,6+3])
    legend('Stance Leg', 'Swing Leg')
    legend boxoff

    subplot(5,n_pools,[9+1,12+1])
    ylabel('Tibia Angle (rad)')
    subplot(5,n_pools,[9+3,12+3])
    legend('Stance Leg', 'Swing Leg')
    legend boxoff

   
end
set(findobj(gcf,'type','axes'),'FontSize',14, ...
'LineWidth', 1.2);
set(findobj(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'in'     , ...
            'TickLength'  , [.01 .01] , ...
            'YGrid'       , 'on'      , ...
            'XColor'      , [.3 .3 .3], ...
            'YColor'      , [.3 .3 .3], ...
            'LineWidth'   , 2.0       ,...
            'xtick'       , (0:.2:1)    ,...
            'ytick'       , [-.5, .25 0, .25 .5]    ,...
            'XLim'          , [0;.8]     ,...
            'YLim'          , [-.5; .5], ...
            'FontSize'    , 14));
Rsq_mean = Rsq_sum / (5*idx)