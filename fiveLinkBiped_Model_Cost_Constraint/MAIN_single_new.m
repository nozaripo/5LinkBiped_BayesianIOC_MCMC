% global scaling_vector
% function costFcn = MAIN(W, opt_soln)
% scaling_data = load('Cost_Scaling1.mat');
% scaling_vector = scaling_data.cost_bases;
% scaling_vector = scaling_vector';
% scaling_vector(scaling_vector==0) = 1;
% scaling_vector(1,3) = 1;
% numBasis = 15;
% for i = 1:numBasis
%     W = zeros(1,numBasis);
%     W(1,i) = 1;


% [u^2  du^2  work  ddq^2  t]
scaling_vector = [1  78  2.3  94  2.6  1  10];
numBasis = 7;
scaling_vector = ones(1,15);
numBasis = 15;





for i = 1:1
    
%     i=k-1;
    
%     W = rand(1,15);
% % % %     W = randsample(0:.01:1,15,true)
% %     W = [1 0 1 0 4 0 0 0 0 4 2 0 0 1 2]./scaling_vector;
% % %         W = zeros(1,numBasis);
% % %     W(1,9) = 1;
%     W = W/norm(W);


% W = zeros(1,numBasis);
% W(1,3) = 1;


%%
W = zeros(1,numBasis);
W(1,i) = 1/scaling_vector(1,i);
% W = rand(1,numBasis);
% r = randi([1 numBasis-1]);
% spInd = randperm(15,r);
% 
% W(spInd) = 0;
% W = W/norm(W);
% 
%     W = W./scaling_vector;



    
% W = [1; 0; 0; 0; 0]';
% W = [1, 3, 2, 4];
% W = W/norm(W);
% W = [.1826 .5477 .3651 .7303]

% W = [2 3 1 1 4];
% W = [2 3 1 4 2.5];

% W = [0 0 2 0 0];
%MAIN.m  --  simple walker trajectory optimization
%
% This script sets up a trajectory optimization problem for a simple model
% of walking, and solves it using OptimTraj. The walking model is a double
% pendulum, with point feet, no ankle torques, impulsive heel-strike (but
% not push-off), and continuous hip torque. Both legs have inertia. Cost
% function is minimize integral of torque-squared.
%
%

% clc; clear;
addpath ../../

% global param
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                  Parameters for the dynamics function                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% param.dyn.m = 10;  %leg mass
% param.dyn.I = 1;  %leg inertia about CoM
% param.dyn.g = 9.81;  %gravity
% param.dyn.l = 1;  %leg length
% param.dyn.d = 0.2;  %Leg CoM distance from hip

param.dyn.m = 10;  %leg mass
param.dyn.I = .06;  %leg inertia about CoM
param.dyn.g = 9.81;  %gravity
param.dyn.l = .76;  %leg length
param.dyn.d = 0.2;  %Leg CoM distance from hip


param.dyn.m = 10;  %leg mass
param.dyn.I = .01;  %leg inertia about CoM
param.dyn.g = 9.81;  %gravity
param.dyn.l = .76;  %leg length
param.dyn.d = 0.2;  %Leg CoM distance from hip

maxTorque = 20;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                            Weights Vector                               %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% global W
% 
% W = W';


% for j = 1:4
%     W = zeros(1,4);
%     W(1,j) = 1;
    
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up function handles                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% a handle to give you the double dot terms
problem.func.dynamics = @(t,x,u)( dynamics(x,u,param.dyn) );

problem.func.kinematics = @(t,x,f)( kinematics(x,f,param.dyn) );

%% ,diff(f)/(t(end)/(length(t)-1))

problem.func.pathObj = @(t,x,kin,dyn,u)(  costFun(t,  kin  ,  dyn , u , param.dyn , W) );

problem.func.bndCst = @(t0,x0,tF,xF)( periodicGait(xF,x0,param.dyn) );

% problem.func.pathCst= @(t,u)(torqueRates(t,u,maxTorque,1));
problem.func.pathCst= @(t,x)(pathVel(x));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%               Set up bounds on time, state, and control                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
t0 = 0;  tF = 1;

% tF = opt_soln.soln.grid.time(end);

problem.bounds.initialTime.low = t0;
problem.bounds.initialTime.upp = t0;
problem.bounds.finalTime.low = .2;
problem.bounds.finalTime.upp = 1.5;
problem.bounds.control.low = -maxTorque;
problem.bounds.control.upp = maxTorque;

% State: [q1;q2;dq1;dq2];

problem.bounds.state.low = [-pi/3; -pi/3; -inf(2,1)];
problem.bounds.state.upp = [ pi/3;  pi/3;  inf(2,1)];

% problem.bounds.state.low = [-pi/2; -pi/2; -inf(2,1)];
% problem.bounds.state.upp = [ pi/2;  pi/2;  inf(2,1)];

stepAngle = 0.3;
% problem.bounds.initialState.low = [stepAngle-.1; -stepAngle-.1; -inf(2,1)];
% problem.bounds.initialState.upp = [stepAngle+.1; -stepAngle+.1;  inf(2,1)];
% problem.bounds.initialControl.low = -maxTorque;
% problem.bounds.initialControl.upp = maxTorque;
% problem.bounds.finalControl.low = -maxTorque;
% problem.bounds.finalControl.upp = maxTorque;

problem.bounds.initialState.low = [stepAngle-.2; -stepAngle+.2; -inf(2,1)];
problem.bounds.initialState.upp = [stepAngle-.2; -stepAngle+.2;  inf(2,1)];
problem.bounds.initialControl.low = -maxTorque;
problem.bounds.initialControl.upp = maxTorque;
problem.bounds.finalControl.low = -maxTorque;
problem.bounds.finalControl.upp = maxTorque;



% % problem.bounds.initialState.low = [stepAngle-.1; 0; -inf(2,1)];
% % problem.bounds.initialState.upp = [stepAngle-.1; 0;  inf(2,1)];
% % problem.bounds.initialControl.low = -maxTorque;
% % problem.bounds.initialControl.upp = maxTorque;
% % problem.bounds.finalControl.low = -maxTorque;
% % problem.bounds.finalControl.upp = maxTorque;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%              Create an initial guess for the trajectory                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% For now, just assume a linear trajectory between boundary values

problem.guess.time = [t0, tF-.2];

stepRate = 1*(2*stepAngle)/(tF-t0);
x0 = [stepAngle; -stepAngle; stepRate; stepRate];
xF = [-stepAngle; stepAngle; stepRate; stepRate];
problem.guess.state = [x0, xF];



% % stepRate = 1*(2*stepAngle)/(tF-t0);
% % x0 = [stepAngle-.1; 0; stepRate; stepRate];
% % xF = [0; stepAngle-.1; stepRate; stepRate];
% % problem.guess.state = [x0, xF];


problem.guess.control = [0, maxTorque];




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Options:                                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


%NOTE:  Here I choose to run the optimization twice, mostly to demonstrate
%   functionality, although this can be important on harder problems. I've
%   explicitly written out many options below, but the solver will fill in
%   almost all defaults for you if they are ommitted.

% method = 'trapezoid';
method = 'hermiteSimpson';
% method = 'chebyshev';
% method = 'rungeKutta';
% method = 'gpops';

switch method
    case 'trapezoid'
        
        % First iteration: get a more reasonable guess
        problem.options(1).nlpOpt = optimset(...
            'Display','iter',...   % {'iter','final','off'}
            'TolFun',1e-3,...
            'MaxFunEvals',1e4);   %options for fmincon
        problem.options(1).verbose = 3; % How much to print out?
        problem.options(1).method = 'trapezoid'; % Select the transcription method
        problem.options(1).trapezoid.nGrid = 10;  %method-specific options
        
        
        % Second iteration: refine guess to get precise soln
        problem.options(2).nlpOpt = optimset(...
            'Display','iter',...   % {'iter','final','off'}
            'TolFun',1e-6,...
            'MaxFunEvals',5e4);   %options for fmincon
        problem.options(2).verbose = 3; % How much to print out?
        problem.options(2).method = 'trapezoid'; % Select the transcription method
        problem.options(2).trapezoid.nGrid = 25;  %method-specific options
        
    case 'hermiteSimpson'
        
        % First iteration: get a more reasonable guess
        problem.options(1).nlpOpt = optimset(...
            'Display','iter',...   % {'iter','final','off'}
            'TolFun',1e-4,...
            'TolX',1e-30,...
            'TolCon',1e-8,...
            'MaxIter',inf,...
            'MaxFunEvals',8e4);   %options for fmincon
        problem.options(1).verbose = 3; % How much to print out?
        problem.options(1).method = 'hermiteSimpson'; % Select the transcription method
%         problem.options(1).hermiteSimpson.nSegment = 6;  %method-specific options
        problem.options(1).hermiteSimpson.nSegment = 15;  %method-specific options

        
% %         % Second iteration: refine guess to get precise soln
% %         problem.options(2).nlpOpt = optimset(...
% %             'Display','iter',...   % {'iter','final','off'}
% %             'TolFun',1e-6,...
% %             'MaxFunEvals',5e4);   %options for fmincon
% %         problem.options(2).verbose = 3; % How much to print out?
% %         problem.options(2).method = 'hermiteSimpson'; % Select the transcription method
% %         problem.options(2).hermiteSimpson.nSegment = 15;  %method-specific options
        
        
    case 'chebyshev'
        
        % First iteration: get a more reasonable guess
        problem.options(1).nlpOpt = optimset(...
            'Display','iter',...   % {'iter','final','off'}
            'TolFun',1e-3,...
            'MaxFunEvals',1e4);   %options for fmincon
        problem.options(1).verbose = 3; % How much to print out?
        problem.options(1).method = 'chebyshev'; % Select the transcription method
        problem.options(1).chebyshev.nColPts = 9;  %method-specific options
        
        
        % Second iteration: refine guess to get precise soln
        problem.options(2).nlpOpt = optimset(...
            'Display','iter',...   % {'iter','final','off'}
            'TolFun',1e-8,...
            'MaxFunEvals',5e4);   %options for fmincon
        problem.options(2).verbose = 3; % How much to print out?
        problem.options(2).method = 'chebyshev'; % Select the transcription method
        problem.options(2).chebyshev.nColPts = 15;  %method-specific options
     
    case 'rungeKutta'
        problem.options(1).method = 'rungeKutta'; % Select the transcription method
        problem.options(1).defaultAccuracy = 'low';
        problem.options(2).method = 'rungeKutta'; % Select the transcription method
        problem.options(2).defaultAccuracy = 'medium';
        
    case 'gpops'
        problem.options.method = 'gpops';
        problem.options.defaultAccuracy = 'medium';
        
    otherwise
        error('Invalid method!');
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Solve!                                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


objVal_th = inf;

for j=1:15


%%
i
j

%%%%% THE KEY LINE:
soln = optimTraj(problem);

constraintViolation = soln.info.constrviolation;
costEval    = soln.info.objVal;






%     Transcription Grid points:
% t(j,:,i)    = soln(end).grid.time;
% q1(j,:,i)   = soln(end).grid.state(1,:);
% q2(j,:,i)   = soln(end).grid.state(2,:);
% dq1(j,:,i)  = soln(end).grid.state(3,:);
% dq2(j,:,i)  = soln(end).grid.state(4,:);
% u(j,:,i)    = soln(end).grid.control;
% 
% constviol(j,i) = soln.info.constrviolation;
% objectVal(j,i)    = soln.info.objVal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if costEval < objVal_th  &&  constraintViolation < 1e-4

%     Transcription Grid points:
t(i,:) = soln(end).grid.time;
x(:,:,i)= soln(end).grid.state;
q1(i,:) = soln(end).grid.state(1,:);
q2(i,:) = soln(end).grid.state(2,:);
dq1(i,:) = soln(end).grid.state(3,:);
dq2(i,:) = soln(end).grid.state(4,:);
u(i,:) = soln(end).grid.control;

dyn = dynamics(soln(end).grid.state,u(i,:),param.dyn);
kin = kinematics(soln(end).grid.state,dyn,param.dyn);
[cost_eval(i,:) cost_vector(:,:,i)] = costFun(t(i,:),  kin  ,  dyn , u(i,:) , param.dyn , W);

% rand_cost_bases_eval_scaled(i,:) = rand_cost_bases_eval(i,:)./scaling_vector;

% weight_vector(i,:) = W;
% 
% % rand_cost_bases_eval_scaled(i,:) = rand_cost_bases_eval(i,:)./scaling_vector;
% % 
% weight_vector(i,:) = W.*scaling_vector;

constviol(i,1) = soln.info.constrviolation;
objectVal(i,1)    = soln.info.objVal;

%% Regular optim vs. reproducibility
%%% regular
% tInt(i,:)   = linspace(t(1),t(end),10*length(t)+1);
% xInt        = soln(end).interp.state(tInt(i,:));
% q1Int(i,:)  = xInt(1,:);
% q2Int(i,:)  = xInt(2,:);
% dq1Int(i,:) = xInt(3,:);
% dq2Int(i,:) = xInt(4,:);
% uInt(i,:)   = soln(end).interp.control(tInt(i,:));

% % % % % tEnd(i,1) = soln.grid.time(end);
% % % % % cost_bases(i,1) = soln.info.objVal;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%% reproducbility
% tInt(i,:) = linspace(t(1),t(end),10*length(t)+1);
% xInt = soln(end).interp.state(tInt(i,:));
% q1Int(i,:) = xInt(1,:);
% q2Int(i,:) = xInt(2,:);
% dq1Int(i,:) = xInt(3,:);
% dq2Int(i,:) = xInt(4,:);
% uInt(i,:) = soln(end).interp.control(tInt(i,:));
%%

% % figure(5)
% % plot(tInt,q1Int,'LineWidth',1.5)
% % xlabel('Time(s)')
% % ylabel('Hip Angle(rad)')
% % hold on
% % plot(tInt,q2Int,'-.','LineWidth',1.5)
% % grid on

% figure(6)
% plot(tInt,q2Int)
% xlabel('Time(s)')
% ylabel('Swing Leg Angle(rad)')
% hold on

% % figure(7)
% % plot(tInt,dq1Int,'LineWidth',1.5)
% % xlabel('Time(s)')
% % ylabel('Angular Velocity(rad/s)')
% % hold on
% % plot(tInt,dq2Int,'-.','LineWidth',1.5)
% % grid on


% figure(8)
% plot(tInt,dq2Int)
% xlabel('Time(s)')
% ylabel('Swing Leg Angular Velocity(rad/s)')
% hold on

% % figure(9)
% % plot(tInt,uInt,'LineWidth',1.5)
% % xlabel('Time(s)')
% % ylabel('Hip Relative Torque(N.m)')
% % hold on
% % grid on

% cost_terms_val = costTerms(tInt, xInt, uInt, param.dyn);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                Define Cost Function with optimal soln                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% save('opt_soln_2-3-1-4-2.5_tFvar.mat','soln');
% opt_soln = load('opt_soln_2-3-1-4-2.5_tfVar.mat', 'soln');

% % Interpolated optimal solution:
% tInt_opt = linspace(opt_soln.soln.grid.time(1),opt_soln.soln.grid.time(end),10*length(t)+1);
% xInt_opt = opt_soln.soln(end).interp.state(tInt);
% q1Int_opt = xInt_opt(1,:);
% q2Int_opt = xInt_opt(2,:);
% dq1Int_opt = xInt_opt(3,:);
% dq2Int_opt = xInt_opt(4,:);
% uInt_opt = opt_soln.soln(end).interp.control(tInt);



% costFcn = 5 * norm(q1Int_opt-q1Int) + ...
%             5 * norm(q2Int_opt-q2Int) + ...
%                 norm(uInt_opt-uInt);

% costFcn = norm(uInt_opt-uInt);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Cost Function %%%%%%%%%%'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if tInt_opt>=tInt
%     t_min = tInt;
% else
%     t_min = tInt_opt;
% end
% xInt_opt_cst = spline(tInt_opt,xInt_opt,t_min);
% xInt_cst     = spline(tInt,xInt,t_min);

% costFcn = 5*norm(tInt(end)-tInt_opt(end))+sum(sum((xInt_opt_cst(1:2,:)-xInt_cst(1:2,:)).^2,2).^.5);

% costFcn = sum(sum((xInt_opt(1:2,:)-xInt(1:2,:)).^2,2).^.5);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Plot the solution                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% % % figure(100); clf;
% % % 
% % % subplot(3,1,1); hold on;
% % % plot(tInt,q1Int,'r-'); plot(tInt,q2Int,'b-');
% % % plot([t(1),t(end)],[0,0],'k--','LineWidth',1);
% % % plot(t,q1,'ro'); plot(t,q2,'bo');
% % % legend('leg one','leg two')
% % % xlabel('time (sec)')
% % % ylabel('angle (rad)')
% % % title('Leg Angles')
% % % 
% % % subplot(3,1,2); hold on;
% % % plot(tInt,dq1Int,'r-'); plot(tInt,dq2Int,'b-');
% % % plot(t,dq1,'ro'); plot(t,dq2,'bo');
% % % legend('leg one','leg two')
% % % xlabel('time (sec)')
% % % ylabel('rate (rad/sec)')
% % % title('Leg Angle Rates')
% % % 
% % % subplot(3,1,3); hold on;
% % % plot(t,u,'mo'); plot(tInt,uInt,'m-');
% % % xlabel('time (sec)')
% % % ylabel('torque (Nm)')
% % % title('Hip Torque')


% A.plotFunc = @(t,z)( drawWalker(t,z,param.dyn) );
% A.speed = 0.25;
% A.figNum = 1;
% animate(tInt,xInt,A)
% hold on
% patch(p1(1,:),p1(2,:),0,'EdgeColor','r','FaceColor',...
%     'none','EdgeAlpha','interp','LineWidth',2)

% [p1,p2] = kinematicsAnimate(xInt,param.dyn);
% figure(1);plot(p1(1,:),p1(2,:));
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


for k=1:4
    A.plotFunc = @(t,z)( drawWalker(t,z,param.dyn) );
    A.speed = 0.3;
    A.figNum = 10;
%     animate(t(k,:),x(:,:,k),A)    
    animate(t(k,:),[q1(k,:);q2(k,:);dq1(k,:);dq2(k,:)],A)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save('cost_range_scaled_new.mat' , 'soln' , 'cost_eval',...
%     'rand_cost_bases_eval' , 'rand_cost_bases_eval_scaled',...
%     'scaling_vector','weight_vector','constviol');
% save('cost_range_scaled_init.mat' , 'soln' , 'cost_eval',...
%     'rand_cost_bases_eval' ,'rand_cost_bases_eval_scaled',...
%     'scaling_vector','weight_vector','constviol');

% save('Cost_Scaling1.mat','cost_bases');



% 
% 
% figure(5)
% legend('Torque Squared Integrated - Stance Leg','Torque Squared Integrated - Swing Leg','|Work| - Stance Leg','|Work| - Swing Leg','1/MoS Integrated - Stance Leg','1/MoS Integrated - Swing Leg','Accel Squared Integrated - Stance Leg','Accel Squared Integrated - Swing Leg')
% 
% % figure(6)
% % legend('Torque Squared Integrated','Torque Squared Integrated','|Work|','|Work|','1/MoS Integrated','1/MoS Integrated','Accel Squared Integrated','Accel Squared Integrated')
% 
% figure(7)
% legend('Torque Squared Integrated - Stance Leg','Torque Squared Integrated - Swing Leg','|Work| - Stance Leg','|Work| - Swing Leg','1/MoS Integrated - Stance Leg','1/MoS Integrated - Swing Leg','Accel Squared Integrated - Stance Leg','Accel Squared Integrated - Swing Leg')
% 
% % figure(8)
% % legend('Torque Squared Integrated','|Work|','1/MoS Integrated','Accel Squared Integrated')
% 
% figure(9)
% legend('Torque Squared Integrated','|Work|','1/MoS Integrated','Accel Squared Integrated')
% end

% % % % for i = 1:30
% % % %    y(i,:) = q1Int(1,:,i);
% % % % end
% % % % figure(1)
% % % % shadedErrorBar(tInt, q1Int, {@mean,@std}, 'lineprops', '-r')
% % % 
% % % % for i = 1:30
% % % %    y(i,:) = q2Int(1,:,i);
% % % % end
% % % % figure(2)
% % % % shadedErrorBar(tInt, q2Int, {@mean,@std}, 'lineprops', '-r')



% figure(1)
% plot(tInt,q1Int_opt,'LineWidth',1)
% xlabel('Time (s)','FontSize',14)
% ylabel('Hip Angle (rad)','FontSize',14)
% hold on
% grid on
% for i = 1:nIOC
%    y(i,:) = soln.xInt(1,:,i);
% end
% shadedErrorBar(tInt, y, {@mean,@std}, 'lineprops', '-c')
% 
% 
% 
% figure(1)
% plot(tInt,q2Int_opt,'r','LineWidth',1)
% for i = 1:nIOC
%    y(i,:) = soln.xInt(2,:,i);
% end
% shadedErrorBar(tInt, y, {@mean,@std}, 'lineprops', '-m')
% legend ('Stance Leg Demonstration','Stance Leg IOC Optimized','Swing Leg Demonstration','Swing Leg IOC Optimized')
% %%%%%%%%%%%%%%%%%%
% 
% 
% figure(2)
% plot(tInt,xInt_opt(3,:),'LineWidth',1)
% xlabel('Time (s)','FontSize',14)
% ylabel('Hip Angular Velocity (rad/s)','FontSize',14)
% hold on
% grid on
% for i = 1:nIOC
%    y(i,:) = soln.xInt(3,:,i);
% end
% shadedErrorBar(tInt, y, {@mean,@std}, 'lineprops', 'c-')
% 
% 
% figure(2)
% plot(tInt,xInt_opt(4,:),'r','LineWidth',1)
% for i = 1:nIOC
%    y(i,:) = soln.xInt(4,:,i);
% end
% shadedErrorBar(tInt, y, {@mean,@std}, 'lineprops', 'm-')
% legend ('Stance Leg Demonstration','Stance Leg IOC Optimized','Swing Leg Demonstration','Swing Leg IOC Optimized')
% 
% % figure(8)
% % plot(tInt,dq2Int)
% % xlabel('Time(s)')
% % ylabel('Swing Leg Angular Velocity(rad/s)')
% % hold on
% 
% figure(3)
% plot(tInt,uInt_opt,'LineWidth',2)
% xlabel('Time (s)','FontSize',14)
% ylabel('Hip Torque (N.m)','FontSize',14)
% hold on
% grid on
% 
% for i = 1:nIOC
%    y(i,:) = soln.uInt(1,:,i);
% end
% shadedErrorBar(tInt, y, {@mean,@std}, 'lineprops', '-r')
% 
% legend('Demonstration', 'IOC Optimized')
% 
% figure(4)
% boxplot(soln.Dec_Var./Dec_Var(:,1))
% set( gca, 'XTickLabel', {'Torques Squared Int ((N.m)^2.s)' 'Absolute Work (J)' '1/MOS Int (s/m)' 'Sum of Accel Int (m/s)'}','FontSize',12 )
% xtickangle(50)
% xlabel('Cost Components','FontSize',14)
% ylabel('Weights Normalized to 1st','FontSize',14)
% hold on
% grid on
% 
% figure(5)
% boxplot(soln.cost_terms_val)
% set( gca, 'XTickLabel', {'Torques Squared Int ((N.m)^2.s)' 'Absolute Work (J)' '1/MOS Int (s/m)' 'Sum of Accel Int (m/s)'}','FontSize',12 )
% xtickangle(50)
% xlabel('Cost Components','FontSize',14)
% ylabel('IOC Optimal Values','FontSize',14)
% hold on
% grid on
