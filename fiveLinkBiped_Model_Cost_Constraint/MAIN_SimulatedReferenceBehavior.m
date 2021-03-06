% MAIN.m  --  Five Link Biped trajectory optimization
%
% This script sets up and then solves the optimal trajectory for the five
% link biped, assuming that the walking gait is compused of single-stance
% phases of motion connected by impulsive heel-strike (no double-stance or
% flight phases).
%
% The equations of motion and gradients are all derived by:
%   --> Derive_Equations.m 
%

% clc; clear; 
addpath ../../



% % % Cost_Components = importdata('Cost_Comp_Eval.txt');
% % % St_Dev = std(Cost_Components);

Cost_Components = importdata('Cost_Components_Eval_Fakuchi_AllData_NoDynamics.txt');
St_Dev = std(Cost_Components);


numBasis = size(Cost_Components,2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up parameters and options                     %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
num_simulation = 5;

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
%                       Set up function handles                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% for i = 1:numBasis
for i = [4 7 9]

W = zeros(1,numBasis);
W(1,i) = 1;


problem.func.dynamics =  @(t,x,u)( dynamics(t,x,u,param) );

% problem.func.pathObj = @(t,x,u)( obj_torqueSquared(u) );
% problem.func.pathObj = @(t,x,u)( composite_objective(t,x,u,param) );
problem.func.pathObj = @(t,x,u)( composite_objective_weighted(t,x,u,param, W, ones(1,numBasis)) );

problem.func.bndCst = @(t0,x0,tF,xF)( stepConstraint(x0,xF,param) );

problem.func.pathCst = @(t,x,u)( pathConstraint(x,u,q0,param) );


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%               Set up bounds on time, state, and control                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = param.stepTime-.3;
problem.bounds.finalTime.upp = param.stepTime-.3;

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
problem.bounds.state.low(3) = -.005;
problem.bounds.state.upp(3) = +.005;


uMax = 100;  %Nm
problem.bounds.control.low = -uMax*ones(5,1);
problem.bounds.control.upp = uMax*ones(5,1);

% Disable the stance ankle motor:
problem.bounds.control.low(1) = 0;
problem.bounds.control.upp(1) = 0;
problem.bounds.control.low(3) = -0;
problem.bounds.control.upp(3) = 0;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%              Create an initial guess for the trajectory                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% For now, just assume a linear trajectory between boundary values

problem.guess.time = [0, param.stepTime];
problem.guess.time = [0, problem.bounds.finalTime.low + (problem.bounds.finalTime.upp-problem.bounds.finalTime.low)*rand(1)];

% q0 = [...
%     -0.3; % stance leg tibia angle
%     0.7; % stance leg femur angle
%     0.0; % torso angle
%     -0.5; % swing leg femur angle
%     -0.6]; % swing leg tibia angle
qF = q0([5;4;3;2;1]);   %Flip left-right

dq0 = (qF-q0)/param.stepTime;
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
    'MaxFunEvals',3e5);   %options for fmincon

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



objVal_th = inf;

% for j=1:15
for j=1:num_simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Solve!                                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

i
j

%%%%% THE KEY LINE:
soln = optimTraj(problem);



% if strcmp(soln.info.message(1,2:6),'Local') || strcmp(soln.info.message(1,2:6),'Solve')
if strcmp(soln.info.message(1,2:6),'Local') || strcmp(soln.info.message(1,2:6),'Feasi') || soln.info.constrviolation<1e-3

    objVal_th(soln.info.objVal<objVal_th) = soln.info.objVal;
    solution{i,1} = soln;
        

end
solutions_all{i,j} = soln;

objective_values (i,j) = soln.info.objVal;
constviol_values (i,j) = soln.info.constrviolation;
messages {i,j} = soln.info.message(1,2:6)

end
end

% dlmwrite('Simulated_Single_Component_Reference_Solutions.txt',solution)
save('Simulated_Single_Component_Reference_Solutions_0.4sec.mat','solution')
save('Simulated_Single_Component_Reference_Solutions_0.4sec.mat','solution', 'objective_values', 'constviol_values')

% % % for i = 1:numBasis
% % %     i
% % %     solution{i,1}.info.message(1,2:6)
% % %     solution{i,1}.info.constrviolation
% % %     
% % % end

% Transcription Grid points:
% % % t = soln(end).grid.time;
% % % q = soln(end).grid.state(1:5,:);
% % % dq = soln(end).grid.state(6:10,:);
% % % u = soln(end).grid.control;



% % % cost_index = inf*ones(length(solutions_eval),1);
% % % for i=1:num_simulation
% % % if solutions_eval{i,1}.constrviolation<1e-5 && solutions_eval{i,1}.objVal<300
% % % cost_index(i,1)=solutions_eval{i,1}.objVal;
% % % end
% % % end
% % % index = find(cost_index==min(cost_index));
% % % 
% % % soln = solution{index,1};
% % % % Transcription Grid points:
% % % t = soln(end).grid.time;
% % % q = soln(end).grid.state(1:5,:);
% % % dq = soln(end).grid.state(6:10,:);
% % % u = soln(end).grid.control;







%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Plot the solution                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%Anim.figNum = 1; clf(Anim.figNum);
% % % Anim.speed = 0.2;
% % % Anim.plotFunc = @(t,q)( drawRobot(q,param) );
% % % Anim.verbose = true;
% % % animate(t,q,Anim);
% % % 
% % % figure(2); clf;
% % % subplot(1,2,1);
% % % plot(t,q);
% % % legend('q1','q2','q3','q4','q5');
% % % xlabel('time')
% % % ylabel('link angles')
% % % subplot(1,2,2);
% % % plot(t,u);
% % % legend('u1','u2','u3','u4','u5');
% % % xlabel('time')
% % % ylabel('joint torques')
% % % 
% % % if isfield(soln(1).info,'sparsityPattern')
% % %    figure(3); clf;
% % %    spy(soln(1).info.sparsityPattern.equalityConstraint);
% % %    axis equal
% % %    title('Sparsity pattern in equality constraints')
% % % end
% % % 
% % % 
% % % 
% % % end

