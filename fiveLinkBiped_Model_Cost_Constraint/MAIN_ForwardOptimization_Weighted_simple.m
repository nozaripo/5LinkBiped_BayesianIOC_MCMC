function [Behavior_Prediction_Error, solution] = MAIN_ForwardOptimization_Weighted_simple(problem, param, p, Criterion_Index, scaling, num_forward_sim , opt_soln)
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

W = zeros(1,length(scaling));

W(1,Criterion_Index) = p;


problem.func.pathObj = @(t,x,u)( composite_objective_weighted_simple(t,x,u,param, W, scaling ) );


objVal_th = inf;
solution = [];
% for j=1:15
for j=1:num_forward_sim
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Solve!                                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%



%%%%% THE KEY LINE:
soln = optimTraj(problem);

% if strcmp(soln.info.message(1,2:6),'Local') || strcmp(soln.info.message(1,2:6),'Solve')
if strcmp(soln.info.message(1,2:6),'Local') || strcmp(soln.info.message(1,2:6),'Feasi') || soln.info.constrviolation < 1e-5

    if soln.info.objVal<objVal_th
    solution = soln;
    objVal_th = soln.info.objVal
    end
%     solution{i,1} = soln;
        
 
    
    
end
% solutions_all{i,j} = soln;
end

if isempty(solution)==0
    
    t  = solution.grid.time;
    x  = solution.grid.state;
    u  = solution.grid.control;
    
    t_pred = linspace(0,t(end),200);
    x_pred = solution.interp.state(t_pred);
        q_pred   = x_pred(1:5,:);
        dq_pred  = x_pred(6:10,:);
    u_pred = solution.interp.control(t_pred);


Behavior_Prediction_Error = norm ( q_pred - opt_soln.qInt );

else
    Behavior_Prediction_Error = inf;
end
    

end

% dlmwrite('Simulated_Single_Component_Reference_Solutions.txt',solution)
% save('Simulated_Single_Component_Reference_Solutions.mat','solution')

% for i = 1:numBasis
%     i
%     solution{i,1}.info.message(1,2:6)
%     solution{i,1}.info.constrviolation
% end

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

