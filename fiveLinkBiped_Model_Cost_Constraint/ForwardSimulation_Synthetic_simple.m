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
addpath ../DirectCollocation_OC/

clear all
V_treadmill = 1.0;

V_Treadmill_Array = [.4, .7, 1.0, 1.2, 1.4, 1.6];
soln_map = containers.Map('KeyType','double','ValueType','any');


for i_v=1:length(V_Treadmill_Array)% comment when you would want to run a single optim
    V_treadmill = V_Treadmill_Array(i_v);% comment when you would want to run a single optim
    
    
    % % % Cost_Components = importdata('Cost_Comp_Eval.txt');
    % % % St_Dev = std(Cost_Components);
    
    Cost_Components = importdata('Cost_Components_Eval_Fakuchi_AllData_NoDynamics.txt');
    St_Dev = std(Cost_Components);
    
    
    numBasis = size(Cost_Components,2);
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %                       Set up parameters and options                     %
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    % num_simulation = 5;
    
    param = getPhysicalParameters();
    
    param = getPhysicalParameters_Anthropometric(mass, height);
    
    param.stepLength = 0.35;
    param.stepTime = 0.7;
    
    
    q0 = [...
        0.2; % stance leg tibia angle
        0.3; % stance leg femur angle
        0.0; % torso angle
        -0.2; % swing leg femur angle
        -0.3]; % swing leg tibia angle
    
    q0 = [...
        0.25; % stance leg tibia angle
        0.35; % stance leg femur angle
        0.0; % torso angle
        -0.20; % swing leg femur angle
        -0.3]; % swing leg tibia angle
    
    qF = q0([5;4;3;2;1]);   %Flip left-right
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %                       Set up function handles                           %
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    
    % for i = 1:numBasis
    for i = 1
        
        W = zeros(1,numBasis);
        W(1,i) = 1;
        
        W = [.2, .5, .3];
        
        problem.func.dynamics =  @(t,x,u)( dynamics(t,x,u,param) );
        
        % problem.func.pathObj = @(t,x,u)( obj_torqueSquared(u) );
        % problem.func.pathObj = @(t,x,u)( composite_objective(t,x,u,param) );
        problem.func.pathObj = @(t,x,u)( composite_objective_weighted_simple(t,x,u,param, W, ones(1,numBasis)) );
        
        problem.func.bndCst = @(t0,x0,tF,xF)( stepConstraint(x0,xF,param) );
        
        problem.func.pathCst = @(t,x,u)( pathConstraint(t,x,u,q0,param, V_treadmill) );
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        %               Set up bounds on time, state, and control                 %
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        problem.bounds.initialTime.low = 0;
        problem.bounds.initialTime.upp = 0;
        problem.bounds.finalTime.low = param.stepTime-.35;
        problem.bounds.finalTime.upp = param.stepTime;
        
        % State: (absolute reference frames)
        %   1 = stance leg tibia angle
        %   2 = stance leg femur angle
        %   3 = torso angle
        %   4 = swing leg femur angle
        %   5 = swing leg tibia angle
        
        qLow = (-pi/4)*ones(5,1);
        qUpp = (pi/3)*ones(5,1);
        dqLow = -10*ones(5,1);
        dqUpp = 10*ones(5,1);
        problem.bounds.state.low = [qLow; dqLow];
        problem.bounds.state.upp = [qUpp; dqUpp];
        
        problem.bounds.initialState.low = [qLow; dqLow];
        problem.bounds.initialState.upp = [qUpp; dqUpp];
        % problem.bounds.initialState.low = [q0; dqLow];
        % problem.bounds.initialState.upp = [q0; dqUpp];
        
        problem.bounds.finalState.low = [qLow; dqLow];
        problem.bounds.finalState.upp = [qUpp; dqUpp];
        % problem.bounds.finalState.low = [qF; dqLow];
        % problem.bounds.finalState.upp = [qF; dqUpp];
        
        %
        % % problem.bounds.state.low(3) = -.1;
        % % problem.bounds.state.upp(3) = +.1;
        
        
        uMax = 30;  %Nm
        problem.bounds.control.low = -uMax*ones(5,1);
        problem.bounds.control.upp = uMax*ones(5,1);
        
        % Disable the stance ankle motor:
        % % problem.bounds.control.low(1) = 0;
        % % problem.bounds.control.upp(1) = 0;
        % problem.bounds.control.low(3) = -0;
        % problem.bounds.control.upp(3) = 0;
        
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
        method = 'hermiteSimpson';
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
            'TolFun',1e-10,...%TolFun is a lower bound on the change in the value of the objective function during a step
            'MaxFunEvals',3e5,...
            'TolCon', 1e-6, ...
            'Tolx',1e-6)   %size smallest step. Smaller step causes optimizer to stop
        %     'TolCon',1e-3,...
        %     'DiffMinChange',1e-2,...
        %     'DiffMaxChange',1e-1)
        % %     'FinDiffType','Central');   %options for fmincon
        
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
                problem.options(1).hermiteSimpson.nSegment = 5;  %method-specific options
                
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
        % for j=1:num_simulation
        for j=1:1
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
            %                           Solve!                                        %
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
            
            i
            j
            
            %%%%% THE KEY LINE:
            tic;
            % soln = optimTraj(problem);
            % el_t(j)=toc;
            
            
            
            
            % t = soln(end).grid.time;
            % q = soln(end).grid.state;
            % u = soln(end).grid.control;
            
            
            
            
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
            %                     Plot the solution                                   %
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
            
            % % Anim.figNum = 1; clf(Anim.figNum);
            % figure(1); clf;
            % Anim.speed = 0.2;
            % Anim.plotFunc = @(t,q)( drawRobot(q,param) );
            % Anim.verbose = true;
            % animate(t,q,Anim);
            
            % % figure(2); clf;
            % % subplot(1,2,1);
            % % plot(t,q);
            % % legend('q1','q2','q3','q4','q5');
            % % xlabel('time')
            % % ylabel('link angles')
            % % subplot(1,2,2);
            % % plot(t,u);
            % % legend('u1','u2','u3','u4','u5');
            % % xlabel('time')
            % % ylabel('joint torques')
            % %
            % % if isfield(soln(1).info,'sparsityPattern')
            % %    figure(3); clf;
            % %    spy(soln(1).info.sparsityPattern.equalityConstraint);
            % %    axis equal
            % %    title('Sparsity pattern in equality constraints')
            % % end
            
            % figure(20);
            % subplot(1,3,1);
            % plot(t,q(1:5,:));
            % hold on
            % legend('q1','q2','q3','q4','q5');
            % xlabel('time')
            % ylabel('link angles')
            % title('angles')
            % subplot(1,3,2);
            % plot(t,q(6:10,:));
            % hold on
            % legend('dq1','dq2','dq3','dq4','dq5');
            % xlabel('time')
            % ylabel('link angle velocities')
            % title('angular velocity')
            % subplot(1,3,3);
            % plot(t,u);
            % legend('u1','u2','u3','u4','u5');
            % xlabel('time')
            % ylabel('joint torques')
            % title('torques')
            % hold on
            
        end
        
        
        
    end
    
    
    
    % Anim.figNum = 1; clf(Anim.figNum);
    % Anim.speed = 0.2;
    % Anim.plotFunc = @(t,q)( drawRobot(t,q,param, V_treadmill) );
    % Anim.verbose = true;
    % animate(t,q,Anim);
    
    
    
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %                           Solve!                                        %
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    
    
    objVal_th = inf;
    constraintViolation = .1;
    firstorderopt = inf;
    % j=0;
    while constraintViolation>=1e-7 || firstorderopt>.1
        
        
        %%
        % i
        % j = j+1
        
        %%%%% THE KEY LINE:
        soln = optimTraj(problem);
        
        constraintViolation = soln.info.constrviolation;
        firstorderopt = soln.info.firstorderopt;
        costEval    = soln.info.objVal;
        
        t = soln(end).grid.time;
        tInt   = linspace(t(1),t(end),10*length(t)+1);
        xInt   = soln(end).interp.state(tInt);
        
        q = soln(end).grid.state;
        % if constraintViolation < 1e-6
        
        %     Transcription Grid points:
        % t = soln(end).grid.time;
        % % q1 = soln(end).grid.state(1,:);
        % % q2 = soln(end).grid.state(2,:);
        % % dq1 = soln(end).grid.state(3,:);
        % % dq2 = soln(end).grid.state(4,:);
        u = soln(end).grid.control;
        
        
        
        
        %
        %
        %
        % tInt   = linspace(t(1),t(end),10*length(t)+1);
        % xInt   = soln(end).interp.state(tInt);
        
        % break
        
        % q1Int  = xInt(1,:);
        % q2Int  = xInt(2,:);
        % dq1Int = xInt(3,:);
        % dq2Int = xInt(4,:);
        % uInt   = soln(end).interp.control(tInt);
        %
        % dyn = dynamics(xInt,uInt,param.dyn);
        % kin = kinematics(xInt,dyn,param.dyn);
        
        
        % % cost_eval = costFun(tInt,  kin  ,  dyn , uInt , param.dyn , W);
        
        % rand_cost_bases_eval_scaled(i,:) = rand_cost_bases_eval(i,:)./scaling_vector;
        
        % weight_vector(i,:) = W;
        %
        % % rand_cost_bases_eval_scaled(i,:) = rand_cost_bases_eval(i,:)./scaling_vector;
        % %
        % weight_vector(i,:) = W.*scaling_vector;
        
        % constviol    = soln.info.constrviolation;
        % objectVal    = soln.info.objVal;
        
        %% Regular optim vs. reproducibility
        %%% regular
        
        
        % % % % % tEnd(i,1) = soln.grid.time(end);
        % % % % % cost_bases(i,1) = soln.info.objVal;
        
        % end
        
    end
    
%     plot(x(1:2),y(1:2),'LineWidth',6,'Color',[1 0 1]);
% plot(x(2:3),y(2:3),'LineWidth',6,'Color',[1 0 0]);
% plot(x(3:4),y(3:4),'LineWidth',6,'Color',[.3 .3 .3]);
% plot(x([3,5]),y([3,5]),'LineWidth',6,'Color',[0 0 1]);
% plot(x(5:6),y(5:6),'LineWidth',6,'Color',[0 1 1]);
    
    figure(20);
    subplot(3,1,1);
%     plot(t,q(1:5,:), 'LineWidth', 1.2);
        hold on
    plot(t,q(1,:),'LineWidth',1+.1*i_v,'Color',[1-.1*i_v 0 1-.1*i_v]);
    plot(t,q(2,:),'LineWidth',1+.1*i_v,'Color',[1-.1*i_v 0 0]);
    plot(t,q(3,:),'LineWidth',1+.1*i_v,'Color',[.5 .5 .5]-.05*i_v);
    plot(t,q(4,:),'LineWidth',1+.1*i_v,'Color',[0 0 1-.1*i_v]);
    plot(t,q(5,:),'LineWidth',1+.1*i_v,'Color',[0 1-.1*i_v 1-.1*i_v]);
    % patchline(t',q(1:5,:)', 'LineWidth', 1.2, 'edgealpha',.5 + i_v/2/length(V_Treadmill_Array));
    legend('q1','q2','q3','q4','q5');
    xlabel('Time (s)')
    ylabel('Limb Angles (rad)')
    title('Angles')
    subplot(3,1,2);
%     plot(t,q(1:5,:), 'LineWidth', 1.2);
    % patchline(t',q(6:10,:)', 'LineWidth', 1.2, 'edgealpha',.5 + i_v/2/length(V_Treadmill_Array));
    hold on
    plot(t,q(6,:),'LineWidth',1+.1*i_v,'Color',[1-.1*i_v 0 1-.1*i_v]);
    plot(t,q(7,:),'LineWidth',1+.1*i_v,'Color',[1-.1*i_v 0 0]);
    plot(t,q(8,:),'LineWidth',1+.1*i_v,'Color',[.5 .5 .5]-.05*i_v);
    plot(t,q(9,:),'LineWidth',1+.1*i_v,'Color',[0 0 1-.1*i_v]);
    plot(t,q(10,:),'LineWidth',1+.1*i_v,'Color',[0 1-.1*i_v 1-.1*i_v]);
%     legend('dq1','dq2','dq3','dq4','dq5');
    xlabel('Time (s)')
    ylabel('link angle velocities')
    title('Angular Velocity (rad/s)')
    subplot(3,1,3);
%     plot(t,u, 'LineWidth', 1.2);
    hold on
    plot(t,u(1,:),'LineWidth',1+.1*i_v,'Color',[1-.1*i_v 0 1-.1*i_v]);
    plot(t,u(2,:),'LineWidth',1+.1*i_v,'Color',[1-.1*i_v 0 0]);
    plot(t,u(3,:),'LineWidth',1+.1*i_v,'Color',[.5 .5 .5]-.05*i_v);
    plot(t,u(4,:),'LineWidth',1+.1*i_v,'Color',[0 0 1-.1*i_v]);
    plot(t,u(5,:),'LineWidth',1+.1*i_v,'Color',[0 1-.1*i_v 1-.1*i_v]);
    % patchline(t',u', 'LineWidth', 1.2, 'edgealpha',.5 + i_v/2/length(V_Treadmill_Array));
    legend('u1','u2','u3','u4','u5');
    xlabel('Time (s)')
    ylabel('Joint Torques (N.m)')
    title('Torques')
    
    
    t_anim = [t t(2:end)+t(end)];
    q_anim = [q  [q(5:-1:1,2:end);q(10:-1:6,2:end)]];
    Anim.speed = 0.25;
%     Anim.plotFunc = @(t,q)( drawRobot(t,q,param, V_treadmill) );
    Anim.plotFunc = @(t, q)( drawRobot(t, q, param, V_treadmill) );
    Anim.verbose = true;
    Anim1.plotFunc = @(t, q)( drawRobot_2ndStep(t, q, param, V_treadmill) );
    Anim1.verbose = true;
    Anim1.speed= .25;
    animate(t,q,Anim);
    hold on
    animate(t+t(end),q ,Anim1)
    animate(t+2*t(end),q,Anim);
    hold on
    animate(t+3*t(end),q ,Anim1)
    
    constviol    = soln.info.constrviolation;
    objectVal    = soln.info.objVal;
    
    
    
    soln_map(V_treadmill) = soln; % comment when you would want to run a single optim
    
%     Anim.speed = 0.2;
%     Anim.plotFunc = @(t,q)( drawRobot(t,q,param, V_treadmill) );
%     Anim.verbose = true;
%     animate(t,q,Anim);
    
end % comment when you would want to run a single optim

save('Simulated_Speeds_.2dq.u+.5du^2+.3d3q^2.mat','soln_map')


% Anim.speed = 0.2;
% Anim.plotFunc = @(t,q)( drawRobot(t,q,param, V_treadmill) );
% Anim.verbose = true;
% animate(t,q,Anim);

% dis_cost = norm(q1-opt_soln.q1Int) + norm(q2-opt_soln.q2Int)+...
%             norm(dq1-opt_soln.dq1Int) + norm(dq2-opt_soln.dq2Int);

% % % dis_cost = rms(rms(xInt-opt_soln.xInt,2));


% save('Simulated_.2u^2+.5du^2+.3d3q^2.mat','soln')


% % save('Simulated_.2dq*u+.5du^2+.3d3q^2.mat','soln')
%%
% save('Simulated_.2du^2+.5u+.3AM.mat','soln')





% dlmwrite('Simulated_Single_Component_Reference_Solutions.txt',solution)
% % % save('Simulated_Single_Component_Reference_Solutions_0.4sec.mat','solution')
% % % save('Simulated_Single_Component_Reference_Solutions_0.4sec.mat','solution', 'objective_values', 'constviol_values')









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

