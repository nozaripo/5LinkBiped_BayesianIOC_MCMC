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

clc; clear; 
addpath ../../

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up parameters and options                     %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
num_simulation = 1;

param = getPhysicalParameters();

param.stepLength = 0.35;
param.stepTime = 0.7;


q0 = [...
    0.3; % stance leg tibia angle
    0.3; % stance leg femur angle
    0.0; % torso angle
    -0.3; % swing leg femur angle
    -0.3]; % swing leg tibia angle





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

% q0 = [...
%     0.3; % stance leg tibia angle
%     0.3; % stance leg femur angle
%     0.0; % torso angle
%     -0.3; % swing leg femur angle
%     -0.3]; % swing leg tibia angle

% q0 = q(:,1);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up function handles                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics =  @(t,x,u)( dynamics(t,x,u,param) );

% problem.func.pathObj = @(t,x,u)( obj_torqueSquared(u) );
problem.func.pathObj = @(t,x,u)( composite_objective(t,x,u,param) );

problem.func.bndCst = @(t0,x0,tF,xF)( stepConstraint(x0,xF,param) );
% problem.func.bndCst = @(t0,x,tF)( stepConstraint(x(1:10,1),x(1:10,end),param) );
% problem.func.bndCst = @(t0,x,tF)( stepConstraint(x(1:10,1),x(1:10,end),param) );

problem.func.pathCst = @(t,x,u)( pathConstraint(x,u,q0,param) );


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%               Set up bounds on time, state, and control                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = param.stepTime-.3;
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


problem.bounds.initialstate.low = [[.25  ; .3 ; -.01 ; -.40 ; -.4] ; dqLow];
problem.bounds.initialstate.upp = [[.35  ;.33 ; +.01 ; -.20 ; -.1] ; dqUpp];
problem.bounds.finalstate.low   = [[-.4 ; -.5 ; -.01 ; .4  ; .10] ; dqLow];
problem.bounds.finalstate.upp   = [[-.1 ; -.2 ; +.01 ; .2  ; .40] ; dqUpp];

% % % % problem.bounds.initialstate.low = [[.25 ; .25 ; -.01 ; -.35 ; -.35] ; dqLow];
% % % % problem.bounds.initialstate.upp = [[.35 ; .35 ; +.01 ; -.25 ; -.25] ; dqUpp];
% % % % problem.bounds.finalstate.low   = [[-.35 ; -.35 ; -.01 ; .25  ; .25] ; dqLow];
% % % % problem.bounds.finalstate.upp   = [[-.25 ; -.25 ; +.01 ; .35  ; .35] ; dqUpp];

problem.bounds.state.low(3) = -.1;
problem.bounds.state.upp(3) = +.1;


uMax = 50;  %Nm
problem.bounds.control.low = -uMax*ones(5,1);
problem.bounds.control.upp = uMax*ones(5,1);

% Disable the stance ankle motor:
% problem.bounds.control.low(1) = 0;
% problem.bounds.control.upp(1) = 0;
% problem.bounds.control.low(3) = -0;
% problem.bounds.control.upp(3) = 0;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%              Create an initial guess for the trajectory                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% For now, just assume a linear trajectory between boundary values

% problem.guess.time = param.stepTime;
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
    'Algorithm','interior-point',...
    'Display','iter',...   % {'iter','final','off'}
    'TolFun',1e-3,... %TolFun is a lower bound on the change in the value of the objective function during a step
    'TolCon',1e-6,...
    'DiffMinChange',1e-2,...
    'DiffMaxChange',1e-1,...
    'FinDiffType','Central',...
    'MaxFunEvals',1e6);   %options for fmincon

switch method
    
    case 'trapezoid'
        problem.options(1).method = 'trapezoid'; % Select the transcription method
        problem.options(1).trapezoid.nGrid = 10;  %method-specific options
%         
%         problem.options(2).method = 'trapezoid'; % Select the transcription method
%         problem.options(2).trapezoid.nGrid = 20;  %method-specific options
        
    case 'trapGrad'  %trapezoid with analytic gradients
        
        problem.options(1).method = 'trapezoid'; % Select the transcription method
        problem.options(1).trapezoid.nGrid = 30;  %method-specific options
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
        problem.options(1).hermiteSimpson.nSegment = 6;  %method-specific options
%         problem.options(1).hermiteSimpson.nSegment = 15;  %method-specific options
        
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


for i=1:num_simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Solve!                                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%%% THE KEY LINE:
soln = optimTraj(problem);

% Transcription Grid points:
t = soln(end).grid.time;
q = soln(end).grid.state(1:5,:);
dq = soln(end).grid.state(6:10,:);
u = soln(end).grid.control;


solution{i,1} = soln;
solutions_eval{i,1} = soln.info;


end


cost_index = inf*ones(length(solutions_eval),1);
for i=1:num_simulation
if solutions_eval{i,1}.constrviolation<1e-5 && solutions_eval{i,1}.objVal<300
cost_index(i,1)=solutions_eval{i,1}.objVal;
end
end
index = find(cost_index==min(cost_index));

soln = solution{index,1};
% Transcription Grid points:
t = soln(end).grid.time;
q = soln(end).grid.state(1:5,:);
dq = soln(end).grid.state(6:10,:);
u = soln(end).grid.control;





for i=1:num_simulation
i
solution{i,1}.info.message


end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Plot the solution                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Anim.figNum = 1; clf(Anim.figNum);
figure(1000)
set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[17 1 17 12]);

Anim.speed = 0.1;
Anim.plotFunc = @(t,q)( drawRobot(q,param) );
Anim.verbose = true;
animate(t,q,Anim);


figure(2);
set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[17 1 17 12]);

hold on
for i=1:num_simulation
    t = solution{i,1}.grid.time;
    q = solution{i,1}.grid.state(1:5,:);
    dq = solution{i,1}.grid.state(6:10,:);
    u = solution{i,1}.grid.control;
    
    t = linspace(0,t(end),200);
    x = solution{i,1}.interp.state(t);
        q  = x(1:5,:);
        dq = x(6:10,:);
    u = solution{i,1}.interp.control(t);
% 
%     u = lowpass(u',.01,100)';
%     figure(99)
%     plot(t', [u' u1])
%     du = diff(u1)'/(t(end)/(length(t)-1));
%     du = [du(:,1) du];
% %     [u2 , du2] = Cubic_Bspline(t' , u(2,:)')
%     figure(999)
%     plot(t , [du])

    
%     if strcmp(soln.info.message(1,2:6),'Feasi') || strcmp(soln.info.message(1,2:6),'Local')
subplot(2,1,1);
plot(t,q);
ax = gca;
ax.ColorOrder = [1 0 1; 1 0 0; 0 0 0; 0 0 1; 0 1 1];
hold on

subplot(2,1,2);
plot(t,u);
ax = gca;
ax.ColorOrder = [1 0 1; 1 0 0; 0 0 0; 0 0 1; 0 1 1];
hold on

%     end
end
subplot(2,1,1)
legend('Stance Tibia','Stance Femur','Torso','Swing Femur','Swing Tibia');
legend boxoff
xlabel('time')
ylabel('link angles')
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.2        ,...
  'FontSize'    , 17);
subplot(2,1,2)
legend('u1','u2','u3','u4','u5');
legend('Stance Ankle','Stance Knee','Stance Hip','Swing Hip','Swing Knee');
legend boxoff

xlabel('time')
ylabel('joint torques')
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.2        ,...
  'FontSize'    , 17);




if isfield(soln(1).info,'sparsityPattern')
   figure(3); clf;
   spy(soln(1).info.sparsityPattern.equalityConstraint);
   axis equal
   title('Sparsity pattern in equality constraints')
end



[c_step, ceq_step] = stepConstraint([q(:,1);dq(:,1)],[q(:,end);dq(:,end)],param)
[c_path, ceq_path] = pathConstraint([q;dq],u,q(:,1),param)

% if max(c_step)<=0 && max(c_path)<=0



