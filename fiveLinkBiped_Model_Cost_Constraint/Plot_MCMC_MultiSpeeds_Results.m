close all
clear all

addpath ../MCMC_IOC/
addpath ../DirectCollocation_OC/
addpath ../MCMC_IOC/Violinplot_Matlab/

filename = 'Results_800iters.mat';
% Resultsresults_20221027T150543.mat
% Resultschain_20221027T150543.mat
% results = load('Resultsresults_20221112T175333.mat').results;
% chain = load('Resultschain_20221112T175333.mat').chain;


synthetic = load('Simulated_Speeds_5Bases.mat');
W_true = synthetic.W;
% results = load(filename).results;
% chain = load(filename).chain;


load('Cost_Bases_Structure_Results.mat');
sPCR_loadings   = sPCR_Object.loadings;
% sPCR_mean       = sPCR_Object.center;
% sPCR_dev        = sPCR_Object.scale;
% W = [.25, .3, .15, .1, .2];
sPCR_loadings = abs(sPCR_loadings);





for ii=1:3
%     filename = strcat('Results_', num2str(ii),'_1500iter.mat');
    filename = strcat('Results_', num2str(ii),'_2500iter.mat');
%     filename = strcat('Results_', num2str(ii),'_2500iter_Bases.mat');
    results(1,1,ii) = load(filename).results;
    chain(:,:,ii)  = load(filename).chain;
%     chain(:,:,ii)  = load(filename).chain * pinv(sPCR_loadings);
    W_true = load(filename).W ;
%     chain(:,:,ii)  = load(filename).chain * sPCR_loadings';
%     W_true = load(filename).W * sPCR_loadings';
%     W_true = load(filename).W * pinv(sPCR_loadings);
    options.nsimu = size(chain,1);
    n_params = size(chain,2);
    n_pools = size(chain,3);
    options.stats = 1;
    options.stats2 = 1;
    options.waitbar = 0;

% define burn-in as a percentage of the number of simulations... 
    burn_in = options.nsimu *0.80;
%     chain_(:,:,ii) = load(filename).chain(burn_in:end,:,:);
    
    indices = find(results(1,1,ii).sss2 <= .05);
%     indices = find(results(1,1,ii).accechain == 1 );
    indices = find(results(1,1,ii).accechain == 1 & results(1,1,ii).sss2 <= .05);
%     indices = find(results(1,1,ii).accechain(burn_in:end) == 1 & results(1,1,ii).sss2(burn_in:end) <= .3);
%     chain_(:,:,ii) = load(filename).chain(indices(end-200:end),:,:);
    chain_(:,:,ii) = chain(indices(end-200:end),:, ii);
    weights_selected(:,:,ii) = chain_(:,:,ii);
    for jj=1:size(chain,2)
        [y,x]=density(chain_(:,jj,ii)./sum(chain_(:,:,ii),2),[], 2);
        xx = linspace(min(x), max(x), size(chain_,1));
        yy = interp1(x, y, xx);
        x_density(:,jj,ii) = xx';
        y_density(:,jj,ii) = yy';
    end
end
chain = chain./sum(chain,2);
chain_ = chain_./sum(chain_,2)*sum(W_true);
% chain_ = chain_./sqrt(sum(chain_.^2,2))*sqrt(sum(W_true.^2));

% chain_ = chain_./sqrt(sum(chain_.^2,2));

% [h,p,k] = kstest2(.2*ones(size(chain_,1),1), chain_(:,1,3))
% [h,p,k] = kstest2(.5*ones(size(chain_,1),1), chain_(:,2,3))
% [h,p,k] = kstest2(.3*ones(size(chain_,1),1), chain_(:,3,3))

for iii = 1:3
    for ik = 1:16
        [h,p,kstat(ik)] = kstest2(W_true(ik), chain_(:,ik,iii));
        Median_Dist(ik) = norm( median(x_density(:,ik,iii))-W_true(ik) );
%         disp(['W_' num2str(ik) ': KS-Stat=' num2str(kstat(ik), '%.2f') '  |  P=' num2str(p, '%.2f') '  |  h=' num2str(h, '%.2f')])
        disp(['W_' num2str(ik) ': Median Distance=' num2str( Median_Dist(ik) , '%.2f') ])
    end
    disp(['Simulation ' num2str(iii) ':  Mean Median Distance = ' num2str(mean(Median_Dist))])
%     sqrt(sum(kstat.^2)/16)
end 

[h,p,k] = kstest2(.2, chain_(:,1,1));
disp(['W_1: KS-Stat=' num2str(k, '%.2f') '|P=' num2str(p, '%.2f')])
[h,p,k] = kstest2(.5, chain_(:,2,1));
disp(['W_2: KS-Stat=' num2str(k, '%.2f') '|P=' num2str(p, '%.2f')])
[h,p,k] = kstest2(.3, chain_(:,3,1));
disp(['W_3: KS-Stat=' num2str(k, '%.2f') '|P=' num2str(p, '%.2f')])
% save('ResultsChains1500iter_20230116_BayesianIOC.mat', 'chain', 'chain_', 'x_density', 'y_density', 'results');
% save('ResultsChains3000iter_20230207_BayesianIOC.mat', 'chain', 'chain_', 'x_density', 'y_density', 'results', 'W_true');

% save('ResultsChains2500iter_20230222_BayesianIOC.mat', 'chain', 'chain_', 'x_density', 'y_density', 'results', 'W_true', 'weights_selected');
% save('ResultsChains2500iter_20230222_BayesianIOC_Bases.mat', 'chain', 'chain_', 'x_density', 'y_density', 'results', 'W_true', 'weights_selected');
% 


Posterior=[];
for ii=1:3
    filename = strcat('Results_', num2str(ii),'_2500iter.mat');
%     filename = strcat('Results_', num2str(ii),'_2500iter_Bases.mat');
    indices = find(load(filename).results.accechain(1:500) == 1);
    ss = load(filename).results.sss2;
    ss(ss>=1e3) = 30;
%     filename = strcat('Results_', num2str(ii),'_2500iter_Bases.mat');
    Posterior(:,ii) = ss;
   plot(indices,ss(indices,1))
   hold on
   ylim([-1 25])

end

save('Posterior2500iter_20230222_BayesianIOC.mat', 'Posterior');
save('Posterior2500iter_20230222_BayesianIOC_Bases.mat', 'Posterior');


plot([Posterior(1:2500,:)])
ylim([-.1 5])


titles = ["w1", "w2", "w3"];

options.nsimu = size(chain,1);
n_params = size(chain,2);
n_pools = size(chain,3);
options.stats = 1; 
options.stats2 = 1; 
options.waitbar = 0;

% define burn-in as a percentage of the number of simulations... 
burn_in = options.nsimu *0.90;

%%
data_sol = load('Simulated_Speeds_.2dq.u+.5du^2+.3d3q^2.mat');
data_sol = load('Simulated_Speeds_5Bases.mat');
V_Treadmill_Array = cell2mat(data_sol.soln_map.keys);

for i_v = 1:length(V_Treadmill_Array)
    t = data_sol.soln_map(V_Treadmill_Array(i_v)).grid.time;
    data.tInt(:,:,i_v) = linspace(t(1),t(end),10*length(t)+1);
    data.xInt(:,:,i_v) = data_sol.soln_map(V_Treadmill_Array(i_v)).interp.state(data.tInt(:,:,i_v));
end
data.V_Treadmill_Array = V_Treadmill_Array;
%%


% assemble the results
results_one = results(:,:,1);
results_two = results(:,:,2);
% % results_thr = results(:,:,3); 
% results_for = results(:,:,4);
% results_fiv = results(:,:,5); 

% plot the chains
figure(9);
set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 12])
set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[17 1 17 12]);


for i =1:size(chain,2)
    subplot(3,1,i)
%     plot(chain(:,i,1),'k')
    plot(chain(:,i,1)./sum(chain(:,:,1),2),'k')
    hold on 
%     plot(chain(:,i,2),'b')
    plot(chain(:,i,2)./sum(chain(:,:,2),2),'b')
    hold on 
%     plot(chain(:,i,3),'r')
    plot(chain(:,i,3)./sum(chain(:,:,3),2),'r')
    hold on 
%     plot(chain(:,i,4),'g')
%     plot(chain(:,i,4)./sum(chain(:,:,4),2),'g')
%     hold on 
%     plot(chain(:,i,5),'c')
    title(titles{i})
    xlabel('Iteration No.')
    ylabel('Weight Value')
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'in'     , ...
        'TickLength'  , [.01 .01] , ...
        'YGrid'       , 'off'      , ...
        'XColor'      , [.3 .3 .3], ...
        'YColor'      , [.3 .3 .3], ...
        'LineWidth'   , 2.0       ,...
        'FontSize'    , 14);
end



w1_real = .2;
w2_real = .5;
w3_real = .3;

w_true = [w1_real w2_real w3_real];

weightnames = ["Work" "Torque" "Jerk"];
% weightnames = [1 2 3];

% plot more results from MCMC toolbox - this creates a bunch of plots,
% uncomment if you want to look at each of the chains
figure(); clf    
set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[17 1 17 12]);

for i = 1:size(chain,3)
%     figure(); clf
%     set(gcf,'units','centimeters','Position',[7.5935 4.2863 10 13])
%     mcmcplot(chain(burn_in:end,:,i),[],results(:,:,i),'pairs');

    
%     set(gcf,'units','centimeters','Position',[7.5935 4.2863 10 13])

%     mcmcplot(chain(burn_in:end,:,i),[],results(:,:,i),'denspanel-ranges',2);
%     subplot(1,size(chain,3),i)
    subplot(2,2,i)
    mcmcplot(chain_(:,:,i),[],results(:,:,i),'denspanel-ranges',2);
    title(['Simulation ' num2str(i)])
    ylim([0, 4.5])
    hold on
    for ii = 1:size(chain,2)
%         subplot(2,2,ii)
        hold on
%         xline(w_true(ii),':r', 'LineWidth', 1.5)
        plot([w_true(ii) w_true(ii)], [4-ii-.1 4-ii+.3],'r', 'LineWidth', 1.5)
        hold on
        xlim([0,1])
        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'in'     , ...
            'TickLength'  , [.01 .01] , ...
            'YGrid'       , 'on'      , ...
            'XColor'      , [.3 .3 .3], ...
            'YColor'      , [.3 .3 .3], ...
            'LineWidth'   , 2.0       ,...
            'xtick'       , (0:.25:1)    ,...
            'ytick'       , [1, 2, 3]    ,...
            'yticklabels'       , {'W_3', 'W_2', 'W_1'}    ,...
            'FontSize'    , 14);
    end
%     figure(101)
%     subplot(2,2,i)
%     ch = reshape([chain_(:,1)' chain_(:,2)' chain_(:,3)'], length(chain_)*length(weightnames),1);
%     we = categorical( reshape(repmat(weightnames, length(chain_),1), length(chain_)*length(weightnames),1) , weightnames);
%     swarmchart( we , ch , 'filled', 'k')
%     hold on
%     plot(['Work'-.1 .2], ['work'+.1 .2])
%     plot([categorical([1-.1]) .2], [categorical([1+.1]) .2])
%     for ii = 1:size(chain,2)
%         subplot(2,2,ii)
%         hold on
%         xline(w_true(ii),':r', 'LineWidth', 1.5)
%         set(gca, ...
%             'Box'         , 'off'     , ...
%             'TickDir'     , 'in'     , ...
%             'TickLength'  , [.01 .01] , ...
%             'YGrid'       , 'off'      , ...
%             'XColor'      , [.3 .3 .3], ...
%             'YColor'      , [.3 .3 .3], ...
%             'LineWidth'   , 2.0       ,...
%             'xtick'       , (0:.25:1)    ,...
%             'ytick'       , [1, 2, 3]    ,...
%             'yticklabels'       , {'W_1', 'W_2', 'W_3'}    ,...
%             'FontSize'    , 14);
%     end
%     figure
%     vs = violinplot(chain(burn_in:end,:,i), ones(101,1));
%     ylabel('Fuel Economy in MPG');
%     xlim([0.5, 7.5]);
    
    
    
    

%     figure(); clf
%     set(gcf,'units','centimeters','Position',[7.5935 4.2863 10 13])
%     mcmcplot(chain(burn_in:end,:,i),[],results(:,:,i),'denspanel2',2);
end
legend('Posterior Density', 'True Values')
sgtitle('Cost Weights - Bayesian Estimates vs. True Values')

% do 10 random draws from each chain, then plot the results to see how
% the results fit the data
n_draws = 5; 
% Draw = zeros(n_draws); 
Draw_Results = zeros(n_draws,n_params); 
% y0 = zeros(n_draws,2,n_pools); 
% t = zeros(length(time),n_draws); 
% % % time = opt_tInt;
% oscillator = zeros(n_draws*2,length(time),n_pools); 

n_pools = size(chain_, 3);
draw_pool = size(chain_, 1);

for i = 1:n_pools
    indices = find(results(1,1,i).accechain == 1);
    for k = 1:n_draws
%        Draw(k) = randi([burn_in options.nsimu]);
       Draw = randi([int64(.5*draw_pool) draw_pool]);
       Draw_Results(k,:,i) = chain_(Draw,:,i);
%        y0(k,:,i) = [Draw_Results(k,6,i),Draw_Results(k,7,i)];
%        [t(:,k),oscillator(:,k*2-1:k*2,i)] = ode15s(@MSD_sys,time,y0(k,:,i),[],Draw_Results(k,:,i));
%        [xInt,~] = ForwardSimulation_simple(Draw_Results(k,:,i));
%        limb_angles(k*10-9:k*10,:,i) = xInt;
       [xInt,tInt] = ForwardSimulation_Speeds_simple(Draw_Results(k,:,i), data.V_Treadmill_Array);
       limb_angles(k*10-9:k*10,:,i, :) = xInt;
       time_soln (k, :, i, :) = tInt; 
    end
end

% plot the random draws 
% figure()
set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 13])

Rsq_sum = 0;
idx = 0;

for i = 1:n_pools
    figure(101)
    subplot(2,2,i)
    for k = 1:n_draws
        hold on
        for j = 1:size(tInt,3)
%             plot(time_soln(k, :, i, j),[limb_angles(k*10-9:k*10-5,:,i, j)],'color',[.17 .17 .17],'LineWidth',1.5)
%             h1 = plot(data.tInt(1, :, j),data.xInt(1:5,:, j),'r--','LineWidth',1.1);
            plot(time_soln(k, :, i, j),limb_angles(k*10-9,:,i, j),'LineWidth',1+.1*i_v,'Color',[1-.1*i_v 0 1-.1*i_v]);
            plot(time_soln(k, :, i, j),limb_angles(k*10-8,:,i, j),'LineWidth',1+.1*i_v,'Color',[1-.1*i_v 0 0]);
            plot(time_soln(k, :, i, j),limb_angles(k*10-7,:,i, j),'LineWidth',1+.1*i_v,'Color',[.5 .5 .5]-.05*i_v);
            plot(time_soln(k, :, i, j),limb_angles(k*10-6,:,i, j),'LineWidth',1+.1*i_v,'Color',[0 0 1-.1*i_v]);
            plot(time_soln(k, :, i, j),limb_angles(k*10-5,:,i, j),'LineWidth',1+.1*i_v,'Color',[0 1-.1*i_v 1-.1*i_v]);
            plot(data.tInt(1, :, j),data.xInt(1,:, j),'LineWidth',1+.1*i_v,'Color',[1-.1*i_v 0 1-.1*i_v]);
            plot(data.tInt(1, :, j),data.xInt(2,:, j),'LineWidth',1+.1*i_v,'Color',[1-.1*i_v 0 0]);
            plot(data.tInt(1, :, j),data.xInt(3,:, j),'LineWidth',1+.1*i_v,'Color',[.5 .5 .5]-.05*i_v);
            plot(data.tInt(1, :, j),data.xInt(4,:, j),'LineWidth',1+.1*i_v,'Color',[0 0 1-.1*i_v]);
            plot(data.tInt(1, :, j),data.xInt(5,:, j),'LineWidth',1+.1*i_v,'Color',[0 1-.1*i_v 1-.1*i_v]);
            
            Rsq_tmp = fitlm(data.xInt(1,:, j), limb_angles(k*10-9,:,i, j)).Rsquared.Adjusted + ...
                        fitlm(data.xInt(2,:, j), limb_angles(k*10-8,:,i, j)).Rsquared.Adjusted + ...
                        fitlm(data.xInt(3,:, j), limb_angles(k*10-7,:,i, j)).Rsquared.Adjusted + ...
                        fitlm(data.xInt(4,:, j), limb_angles(k*10-6,:,i, j)).Rsquared.Adjusted + ...
                        fitlm(data.xInt(5,:, j), limb_angles(k*10-5,:,i, j)).Rsquared.Adjusted ;
            
            idx = idx + 1;
            Rsq_sum = (Rsq_sum + Rsq_tmp ) ;
        end
    end
    hold on
%     for k = 2:4
%          plot(time,position_noise(:,k),'r--','LineWidth',2);
%     end
    if i == 1
%         legend(h1,'Ref. Trajectory')
%         legend('boxoff')
    end
    set(gca,'fontsize',10)
    xlabel('Time (s)')
    ylabel('Angle (rad)')
    box off
    title(['Result from worker No. ', num2str(i)])
    
    
%     figure(102)
%     subplot(2,2,i)
%     for k = 1:n_draws
%         hold on
%         for j = 1:size(tInt,3)
%             plot(time_soln(k, :, i, j),[limb_angles(k*10-4:k*10,:,i, j)],'color',[.17 .17 .17],'LineWidth',1.5)
%             plot(data.tInt(1, :, j),data.xInt(6:10,:, j),'r--','LineWidth',1.1);
%         end
%     end
%     hold on
% %     yline
% %     for k = 2:4
% %          plot(time,position_noise(:,k),'r--','LineWidth',2);
% %     end
%     if i == 1
%         legend(h1,'Ref. Trajectory')
%         legend('boxoff')
%     end
%     set(gca,'fontsize',10)
%     xlabel('Time (s)')
%     ylabel('Velocity (rad/s)')
%     box off
%     title(['Result from worker No. ', num2str(i)])
end 
Rsq_mean = Rsq_sum / (5*idx)








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
    legend('Torso Angles', 'Actual Trajectories')
    ylabel('Torso Angle (rad)')
    legend boxoff

    subplot(5,n_pools,[3+3,6+1])
    ylabel('Femur Angle (rad)')
    legend('Stance Leg', 'Swing Leg')
    legend boxoff

    subplot(5,n_pools,[9+3,12+1])
    ylabel('Tibia Angle (rad)')
    legend('Stance Leg', 'Swing Leg')
    legend boxoff

   
end
set(findobj(gcf,'type','axes'),'FontSize',14, ...
'LineWidth', 1.2,'layer','top');grid on
set(findobj(gca), ...
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
            'FontSize'    , 14);
Rsq_mean = Rsq_sum / (5*idx)






End_results_worker1 = chain(end,:,1)';
End_results_worker2 = chain(end,:,2)';
End_results_worker3 = chain(end,:,3)';
End_results_worker4 = chain(end,:,4)';
True_Weights = [.2, .5, .3]';
table(End_results_worker1, End_results_worker2, End_results_worker3, End_results_worker4, True_Weights, 'RowNames', {'w1', 'w2' , 'w3'})

% Plot the rank without the first burn in 
rank_plot_MSD(chain(burn_in:end,:,:),options,n_pools);

% Plot the density of 
results_one_dist = mvnrnd(results_one.mean, results_one.cov, 2000);

figure(20)
histogram(results_one_dist(:,1))