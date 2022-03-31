% Main script for optimizing the gait of a walking robot model
% Copyright 2017 The MathWorks, Inc.

%% Set initial parameters
clc; clear; 
addpath ../../

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up parameters and options                     %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
num_forward_sim = 4;

param = getPhysicalParameters();

param.stepLength = 0.35;
param.stepTime = 0.7;


q0 = [...
    0.3; % stance leg tibia angle
    0.3; % stance leg femur angle
    0.0; % torso angle
    -0.3; % swing leg femur angle
    -0.3]; % swing leg tibia angle


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%           Load the walking data and set initial config                  %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
walk_data = importdata('WBDS23walkT03ang.txt');
l_hip  = walk_data.data(:,13) * pi/180;
l_knee = walk_data.data(:,19) * pi/180;
% l_Ankle = walk_data.data(:,25);

rep = 2;
shift1 = -25;
shift0 = 50;

r_hip  = circshift(l_hip,shift0);
r_knee = circshift(l_knee,shift0);
% r_ankle = walk_data.data(:,25);

angles = [l_hip , r_hip , l_knee , r_knee];


GRF_data = importdata('WBDS23walkT03grf.txt');

GRF_r_N = GRF_data.data(:,3);
GRF_l_N = GRF_data.data(:,10);
GRF_r_f = GRF_data.data(:,2);
GRF_l_f = GRF_data.data(:,9);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = (walk_data.data(:,1)-1)*.01;

Time = time;

t = time;


l_hip   = lowpass(l_hip,6,100);
r_hip   = lowpass(r_hip,6,100);
l_knee  = lowpass(l_knee,6,100);
r_knee  = lowpass(r_knee,6,100);

l_knee(l_knee<0)=0;
r_knee(r_knee<0)=0;

% joint kinematics
q = [l_hip-l_knee , l_hip , zeros(length(l_hip),1) , r_hip , r_hip-r_knee]';

q0 = q(:,1);

q0(5) = -acos((param.l1*cos(q0(1)) + param.l2*cos(q0(2)) - param.l4*cos(q0(4)))/param.l5);


q0 = [...
    0.3; % stance leg tibia angle
    0.3; % stance leg femur angle
    0.0; % torso angle
    -0.3; % swing leg femur angle
    -0.3]; % swing leg tibia angle



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up function handles                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics =  @(t,x,u)( dynamics(t,x,u,param) );

% problem.func.pathObj = @(t,x,u)( obj_torqueSquared(u) );
% problem.func.pathObj = @(t,x,u)( composite_objective(t,x,u,param) );

problem.func.bndCst = @(t0,x0,tF,xF)( stepConstraint(x0,xF,param) );

problem.func.pathCst = @(t,x,u)( pathConstraint(x,u,q0,param) );


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%               Set up bounds on time, state, and control                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = param.stepTime-.3;
problem.bounds.finalTime.upp = param.stepTime-.0;

problem.bounds.finalTime.low = param.stepTime-.4;
problem.bounds.finalTime.upp = param.stepTime-.1;
% State: (absolute reference frames)
%   1 = stance leg tibia angle
%   2 = stance leg femur angle
%   3 = torso angle
%   4 = swing leg femur angle
%   5 = swing leg tibia angle

qLow = (-pi/2)*ones(5,1);
qUpp = (pi/2)*ones(5,1);
dqLow = -10*ones(5,1);
dqUpp = 10*ones(5,1);
problem.bounds.state.low = [qLow; dqLow];
problem.bounds.state.upp = [qUpp; dqUpp];
problem.bounds.initialstate.low = [qLow; dqLow];
problem.bounds.initialstate.upp = [qUpp; dqUpp];
problem.bounds.finalstate.low = [qLow; dqLow];
problem.bounds.finalstate.upp = [qUpp; dqUpp];

% 
% % % problem.bounds.state.low(3) = -.15;
% % % problem.bounds.state.upp(3) = +.15;


uMax = 150;  %Nm
problem.bounds.control.low = -uMax*ones(5,1);
problem.bounds.control.upp = uMax*ones(5,1);

% Disable the stance ankle motor:
problem.bounds.control.low(1) = 0;
problem.bounds.control.upp(1) = 0;
problem.bounds.control.low(3) = -150;
problem.bounds.control.upp(3) = 150;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%              Create an initial guess for the trajectory                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% For now, just assume a linear trajectory between boundary values

problem.guess.time = [0, problem.bounds.finalTime.low + (problem.bounds.finalTime.upp-problem.bounds.finalTime.low)*rand(1)];

% q0 = [...
%     -0.3; % stance leg tibia angle
%     0.7; % stance leg femur angle
%     0.0; % torso angle
%     -0.5; % swing leg femur angle
%     -0.6]; % swing leg tibia angle
qF = q0([5;4;3;2;1]);   %Flip left-right

dq0 = (qF-q0)/problem.guess.time(2);
dqF = dq0;

problem.guess.state = [q0, qF; dq0, dqF];

problem.guess.control = zeros(5,2);  %Start with passive trajectory


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Options:                                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


%NOTE:  Here I choose to run the optimization twice, mostly to demonstrate
%   functionality, although this can be important on harder problems. I've
%   explicitly written out many options below, but the solver will fill in
%   almost all defaults for you if they are ommitted.

method = 'trapezoid';
% method = 'trapGrad';   % This one is also good
% method = 'hermiteSimpson';
% method = 'hermiteSimpsonGrad';   % Suggested method
% method = 'chebyshev';   
% method = 'rungeKutta';  %slow!
% method = 'rungeKuttaGrad';
% method = 'gpops';

%%%% Method-independent options:
% problem.options(1).nlpOpt = optimset(...
%     'Display','iter',...   % {'iter','final','off'}
%     'TolFun',1e-3,...
%     'MaxFunEvals',1e4);   %options for fmincon
% problem.options(2).nlpOpt = optimset(...
%     'Display','iter',...   % {'iter','final','off'}
%     'TolFun',1e-6,...
%     'MaxFunEvals',1e6);   %options for fmincon

problem.options(1).nlpOpt = optimset(...
    'Display','iter',...   % {'iter','final','off'}
    'TolFun',1e-6,...
    'MaxIter',1000,...
    'MaxFunEvals',5e5);   %options for fmincon

switch method
    
    case 'trapezoid'
        problem.options(1).method = 'trapezoid'; % Select the transcription method
        problem.options(1).trapezoid.nGrid = 10;  %method-specific options
%         
%         problem.options(2).method = 'trapezoid'; % Select the transcription method
%         problem.options(2).trapezoid.nGrid = 20;  %method-specific options
        
    case 'trapGrad'  %trapezoid with analytic gradients
        
        problem.options(1).method = 'trapezoid'; % Select the transcription method
        problem.options(1).trapezoid.nGrid = 10;  %method-specific options
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'off';
        
        problem.options(2).method = 'trapezoid'; % Select the transcription method
        problem.options(2).trapezoid.nGrid = 45;  %method-specific options
        problem.options(2).nlpOpt.GradConstr = 'on';
        problem.options(2).nlpOpt.GradObj = 'on';
        
    case 'hermiteSimpson'
        
        % First iteration: get a more reasonable guess
        problem.options(1).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(1).hermiteSimpson.nSegment = 10;  %method-specific options
        
        % Second iteration: refine guess to get precise soln
%         problem.options(2).method = 'hermiteSimpson'; % Select the transcription method
%         problem.options(2).hermiteSimpson.nSegment = 15;  %method-specific options
        
    case 'hermiteSimpsonGrad'  %hermite simpson with analytic gradients
        
        problem.options(1).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(1).hermiteSimpson.nSegment = 6;  %method-specific options
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'off';
        
        problem.options(2).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(2).hermiteSimpson.nSegment = 15;  %method-specific options
        problem.options(2).nlpOpt.GradConstr = 'on';
        problem.options(2).nlpOpt.GradObj = 'on';
        
        
    case 'chebyshev'
        
        % First iteration: get a more reasonable guess
        problem.options(1).method = 'chebyshev'; % Select the transcription method
        problem.options(1).chebyshev.nColPts = 9;  %method-specific options
        
        % Second iteration: refine guess to get precise soln
        problem.options(2).method = 'chebyshev'; % Select the transcription method
        problem.options(2).chebyshev.nColPts = 15;  %method-specific options
        
    case 'multiCheb'
        
        % First iteration: get a more reasonable guess
        problem.options(1).method = 'multiCheb'; % Select the transcription method
        problem.options(1).multiCheb.nColPts = 6;  %method-specific options
        problem.options(1).multiCheb.nSegment = 4;  %method-specific options
        
        % Second iteration: refine guess to get precise soln
        problem.options(2).method = 'multiCheb'; % Select the transcription method
        problem.options(2).multiCheb.nColPts = 9;  %method-specific options
        problem.options(2).multiCheb.nSegment = 4;  %method-specific options
        
    case 'rungeKutta'
        problem.options(1).method = 'rungeKutta'; % Select the transcription method
        problem.options(1).defaultAccuracy = 'low';
        problem.options(2).method = 'rungeKutta'; % Select the transcription method
        problem.options(2).defaultAccuracy = 'medium';
    
    case 'rungeKuttaGrad'
      
        problem.options(1).method = 'rungeKutta'; % Select the transcription method
        problem.options(1).defaultAccuracy = 'low';
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'off';
        
        problem.options(2).method = 'rungeKutta'; % Select the transcription method
        problem.options(2).defaultAccuracy = 'medium';
        problem.options(2).nlpOpt.GradConstr = 'on';
        problem.options(2).nlpOpt.GradObj = 'on';
        
    case 'gpops'
        problem.options = [];
        problem.options.method = 'gpops';
        problem.options.defaultAccuracy = 'high';
        problem.options.gpops.nlp.solver = 'snopt';  %Set to 'ipopt' if you have GPOPS but not SNOPT
        
    otherwise
        error('Invalid method!');
end



% opt_soln = load('opt_soln_1324_Tvar_TorqueRateCst.mat', 'soln');
% optimal_solution = load('Single_Cost_Optim_11.26.mat');
% % % optimal_solution = load('Single_Cost_Optim_12.02.19.mat');

optimal_solution = load('Simulated_Single_Component_Reference_Solutions.mat');
optimal_solution = optimal_solution.solution;


for i=1:length(optimal_solution)
scaling_vector_simul(1,i) = optimal_solution{i,1}.info.objVal;
end



scaling_vector = max(importdata('Cost_Comp_Eval_Traj.txt'));


scaling_vector = ones(1, length(optimal_solution));

% for robustness: will run GA 3 times for each single criterion
nIOC = 1;



Criterion_Index = [4 5 6 10 13 14];

Criterion_Index = [4 6 8 10 13 14 17];


for i_rob=1:nIOC
 
    
% Flags to speed up simulation
accelFlag = true;
parallelFlag = false;



% Uncomment to create initial conditions from scratch
numPoints = length(Criterion_Index);              % Number of joint angle points

upperBnd = ones(1,numPoints);  
lowerBnd = zeros(1,numPoints); 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p0 = lowerBnd+(upperBnd-lowerBnd).*rand(1,numPoints);  % Create zero motion initial conditions
% p0 = [.25 .047 .72];  % Create zero motion initial conditions
% p0 = [.25 .1];  % Create zero motion initial conditions
% p0 = .25;  % Create zero motion initial conditions



% ind = randperm(randi([numPoints-4 numPoints-1]));

% Do the IOC for recovery of each single cost optimization case
% numPoints
for i_criteria = 1:numPoints

%% Set optimization options
opts = optimoptions('ga');
opts.Display = 'iter';
opts.FitnessLimit = 1e-5;
opts.FunctionTolerance = 1e-12;
opts.ConstraintTolerance = 1e-8;
opts.MaxGenerations = 2;
opts.PopulationSize = 5;
pop1 = ceil(.70*opts.PopulationSize);
spPop = opts.PopulationSize-pop1;

map_mat = ones(opts.PopulationSize,numPoints);
for m = pop1+1:opts.PopulationSize
    
    ind =randperm(numPoints, randi([2 numPoints-1]));
    map_mat(m,ind) = 0;
end

% opts.InitialPopulationMatrix = repmat(p0,[opts.PopulationSize 1]); % Add copies of initial gait
opts.InitialPopulationMatrix = lowerBnd+(upperBnd-lowerBnd).*rand(opts.PopulationSize,numPoints).*map_mat; % Add copies of initial gait
% opts.InitialPopulationMatrix((ind-1)*po
opts.PlotFcn = {@gaplotbestf, @gaplotbestindiv, @gaplotmaxconstr}; % Add progress plot of fitness function
opts.UseParallel = parallelFlag;
opts.UseVectorized = false;

        

%% Run commands to set up parallel/accelerated simulation
% doSpeedupTasks;

% scaling_vector = [1  78  2.3  94  2.6  1  10];
% scaling_vector = [1  78  2.3  94];





% opt_soln.tInt   = optimal_solution.tInt(k,:);
% opt_soln.q1Int  = optimal_solution.q1Int(k,:);    
% opt_soln.q2Int  = optimal_solution.q2Int(k,:);
% opt_soln.dq1Int = optimal_solution.dq1Int(k,:);
% opt_soln.dq2Int = optimal_solution.dq2Int(k,:);
% opt_soln.uInt   = optimal_solution.uInt(k,:);

% % % opt_soln.tInt   = optimal_solution{Criterion_Index(i_criteria)}.t(k,:);
% % % opt_soln.q1Int  = optimal_solution.q1(k,:);    
% % % opt_soln.q2Int  = optimal_solution.q2(k,:);
% % % opt_soln.dq1Int = optimal_solution.dq1(k,:);
% % % opt_soln.dq2Int = optimal_solution.dq2(k,:);
% % % opt_soln.uInt   = optimal_solution.u(k,:);


    t  = optimal_solution{Criterion_Index(i_criteria),1}.grid.time;
    x  = optimal_solution{Criterion_Index(i_criteria),1}.grid.state;
    u  = optimal_solution{Criterion_Index(i_criteria),1}.grid.control;
    
    opt_soln.tInt = linspace(0,t(end),200);
    opt_soln.xInt = optimal_solution{Criterion_Index(i_criteria),1}.interp.state(opt_soln.tInt);
        opt_soln.qInt   = opt_soln.xInt(1:5,:);
        opt_soln.dqInt  = opt_soln.xInt(6:10,:);
    opt_soln.uInt = optimal_solution{Criterion_Index(i_criteria),1}.interp.control(opt_soln.tInt);

%% Run optimization
% costFcn = @(p)MAIN_composite_new(problem, p, scaling_vector, num_forward_sim, opt_soln{i_criteria,1});
costFcn = @(variables)MAIN_ForwardOptimization_Weighted(problem, param, variables, Criterion_Index, scaling_vector, num_forward_sim, opt_soln);

% nonlcon = @(p)weightcons(p); % Constraint for height
% tic
disp(['Running optimization! Population Size: ' num2str(opts.PopulationSize) ...
      ', Fitness Limit: ' num2str(opts.FitnessLimit) ...
      ', Generations No: ' num2str(opts.MaxGenerations)])
   
% ', max Generations No.: ' num2str(opts.maxGenerations) ...  
  
% [pFinal,reward] = ga(costFcn,numPoints,[],[],[],[], ... 
%                      lowerBnd,upperBnd, ...
%                      nonlcon,[],opts);
 tic                
 [pFinal,reward] = ga(costFcn,numPoints,[],[],ones(1,numPoints),1, ... 
                     lowerBnd,upperBnd, ...
                     [],[],opts);
                 
disp(['Final reward function value: ' num2str(-reward)])
% toc

IOC_Results.Weights_row{i_rob,i_criteria} = pFinal;
% [IOC_err(k,1), soln, constviol(k,1), objectVal(k,1)] = MAIN_composite_new(pFinal./scaling_vector,opt_soln);


% [IOC_err(k,1), soln, constviol(k,1), objectVal(k,1)] = 

% [Behavior_Prediction_Error, solution] = ...
[IOC_Results.IOC_err(i_rob,i_criteria), IOC_Results.IOC_Behavior_Prediction{i_rob,i_criteria}] = ...
    MAIN_ForwardOptimization_Weighted(problem, param, pFinal, Criterion_Index, scaling_vector, num_forward_sim , opt_soln);


soln = IOC_Results.IOC_Behavior_Prediction{i_rob,i_criteria};

t = soln(end).grid.time;
% Interpolated solution:
OC_Results.tInt{i_rob,i_criteria}     = linspace(t(1),t(end),200);
xInt   = soln(end).interp.state(OC_Results.tInt{i_rob,i_criteria});
OC_Results.qInt{i_rob,i_criteria}    = xInt(1:5,:);
OC_Results.dqInt{i_rob,i_criteria}    = xInt(6:10,:);
OC_Results.uInt{i_rob,i_criteria}     = soln(end).interp.control(OC_Results.tInt{i_rob,i_criteria});


% % % dyn = dynamics(xInt(:,:,k),uInt(k,:),param.dyn);
% % % kin = kinematics(xInt(:,:,k),dyn,param.dyn);
% % % [composite_cost_integrand(k,:) cost_integrands(:,:,k)] = costFun(tInt,  kin  ,  dyn , uInt(k,:) , param.dyn , pFinal./scaling_vector);

% % % [composite_cost_integrand(k,:) cost_integrands(:,:,k)] = costFun(tInt,  kin  ,  dyn , uInt(k,:) , param.dyn , pFinal./scaling_vector);



W = zeros(1,length(scaling_vector));

W(1,Criterion_Index) = pFinal;

IOC_Results.Weight_full = W;

[IOC_Results.Composite_Objective{i_rob,i_criteria}, IOC_Results.Criteria_TimeSeries{i_rob,i_criteria}] = ...
    composite_objective_weighted(soln(end).grid.time,soln(end).grid.state,soln(end).grid.control,param, W, scaling_vector);



% % % composite_cost_eval(k,1) = simps(tInt(k,:)',composite_cost_integrand(k,:)');
IOC_Results.Criteria_Objective_eval{i_rob,i_criteria}     = simps(soln(end).grid.time',IOC_Results.Criteria_TimeSeries{i_rob,i_criteria}');
% IOC_Results.comp_func_estim{i_rob,i_criteria}     = sum(Weights_row(k,:).*cost_terms_eval(k,:)./scaling_vector);

% scaling_vector = [1 78 2.3 94];
IOC_Results.contribution{i_rob,i_criteria} = W./scaling_vector.*IOC_Results.Criteria_Objective_eval{i_rob,i_criteria}/...
                                                  ((W./scaling_vector)*IOC_Results.Criteria_Objective_eval{i_rob,i_criteria}');

% cost_total(i,1) = -reward;
% cost_term_weighted(i,:)  = -reward*pFinal/sum(pFinal);
% cost_term_raw(i,:)  = -reward*pFinal/sum(pFinal);

IOC_Results.el_time(i_rob,i_criteria) = toc/3600;

end

end


save('IOC_Prediction_Results.mat','OC_Results','IOC_Results');

% el_time = toc/3600;
% el_time = [num2str(floor(el_time)), ' hrs & ', num2str(floor(60*(el_time-floor(el_time)))), ' mins']


% save('IOC_data_new_2-3-1-4-2.5_6.14.mat','Dec_Var','cost_terms_val', 'tInt', 'xInt', 'uInt','integrands');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)
% boxplot(Dec_Var)
% xlabel('Cost Terms')
% ylabel('Weight')
% 
% 
%% Animate walking 



% Anim.figNum = 1; clf(Anim.figNum);
figure(1000)
Anim.speed = 0.1;
Anim.plotFunc = @(t,q)( drawRobot(q,param) );
Anim.verbose = true;
animate(t,q,Anim);



% % % for k=1:1
% % %     A.plotFunc = @(t,z)( drawWalker(t,z,param.dyn) );
% % %     A.speed = 0.3;
% % %     A.figNum = 10;
% % %     animate(tInt(k,:),xInt(:,:,k),A)
% % % %     animate(time',[qL';qR';dqL';dqR'],A)
% % % end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%
% K_v   = pFinal(1);
% Ki_v  = pFinal(2);% Ktau = pFinal(2);
% Kd_v  = pFinal(3); % Kd = pFinal(2);
% K_t   = pFinal(4);
% Ki_t  = pFinal(5);% Ktau = pFinal(2);
% Kd_t  = pFinal(6); % Kd = pFinal(2);
% K_t1  = pFinal(7);
% Ki_t1 = pFinal(8);% Ktau = pFinal(2);
% Kd_t1 = pFinal(9); % Kd = pFinal(2);
% K_t2  = pFinal(10);
% Ki_t2 = pFinal(11);% Ktau = pFinal(2);
% Kd_t2 = pFinal(12); % Kd = pFinal(2);
% K_t3  = pFinal(13);
% Ki_t3 = pFinal(14);% Ktau = pFinal(2);
% Kd_t3 = pFinal(15); % Kd = pFinal(2);
% K_t4  = pFinal(16);
% Ki_t4 = pFinal(17);% Ktau = pFinal(2);
% Kd_t4 = pFinal(18); % Kd = pFinal(2);
% K_t5  = pFinal(19);
% Ki_t5 = pFinal(20);% Ktau = pFinal(2);
% Kd_t5 = pFinal(21); % Kd = pFinal(2);
%%

% el_time = toc/3600
% el_time = [num2str(floor(el_time)), ' hrs & ', num2str(floor(60*(el_time-floor(el_time)))), ' mins']


%% Save results to MAT-file
% Convert from optimization integer search space to trajectories in radians
% pScaled = scalingFactor*pFinal;
% time = linspace(0,gaitPeriod,numPoints+1)';
% ankle_motion = deg2rad([pScaled(1:numPoints) pScaled(1)]');
% knee_motion = deg2rad([pScaled(numPoints+1:2*numPoints) pScaled(numPoints+1)]');
% hip_motion = deg2rad([pScaled(2*numPoints+1:3*numPoints) pScaled(2*numPoints+1)]');
% curveData = createSmoothTrajectory(ankle_motion,knee_motion,hip_motion,gaitPeriod);
% outFileName = ['optimizedData_' datestr(now,'ddmmmyy_HHMM')];
% save(outFileName,'curveData','reward','gaitPeriod','time', ... 
%                  'hip_motion','knee_motion','ankle_motion');
% plotSmoothTrajectory;
%% Cleanup
% bdclose(mdlName);

% [Traj, time, ta, init_height] = run_rcp(pFinal(1,1),pFinal(1,2:end));

% Data=[time, Traj];
% save('Dataaaa.mat','Data','ta','init_height')

% if parallelFlag
%    delete(gcp('nocreate')); 
% end
% 
% save('pFinal.mat','K_t', 'Ki_t', 'Kd_t','K_v', 'Ki_v', 'Kd_v')
% 
% savefig('FitnessPlot.fig')
% robotParameters;
% 
% youBot_PARAM;
% 
% 
% function state = plotfun(options,state,flag)
% 
% 
% 
% end