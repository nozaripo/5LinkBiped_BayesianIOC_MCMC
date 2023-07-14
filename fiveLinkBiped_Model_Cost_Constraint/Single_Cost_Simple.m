function [dObj, dObj_TimeSeries] = Single_Cost_Simple(t,x,u,p, W, scaling, V_tr, cost_index)
% [dObj, dObjGrad] = obj_torqueSquared(u)
%
% This function computes the torque-squared objective function and its
% gradients.
%

%% Test cubic bspline method for differentiation: it is great
% ddq = diff(dq(1,:)')'/(t(end)/(length(t)-1));
% ddq1 = [ddq(:,1) ddq];
% ddq2 = filter(-smooth_diff(4),1,dq(1,:));
% [~, dq3, ddq3] = Cubic_Bspline(t' , q(1,:)');
% figure(3)
% plot([ddq1', ddq3])
% hold on
% plot([dq(1,:)', dq3], ":")



%% KINEMATIC COSTS
% Joints
q  = x(1:5,:);
dq = x(6:10,:);

ddq = diff(dq')'/(t(end)/(length(t)-1));
ddq = [ddq(:,1) ddq];

dddq = diff(ddq')'/(t(end)/(length(t)-1));
dddq = [dddq(:,1) dddq];

% % [~, ddq, dddq] = Cubic_Bspline(t' , dq');
% % ddq = ddq';
% % dddq= dddq';

% ddq = filter(-smooth_diff(10),1,dq);


% Task
[P, G] = getPoints(t,q,p, V_tr);

dP = diff(P')'/(t(end)/(length(t)-1));
dP = [dP(:,1) dP];
dG = diff(G')'/(t(end)/(length(t)-1));
dG = [dG(:,1) dG];

ddP = diff(dP')'/(t(end)/(length(t)-1));
ddP = [ddP(:,1) ddP];
ddG = diff(dG')'/(t(end)/(length(t)-1));
ddG = [ddG(:,1) ddG];

dddP = diff(ddP')'/(t(end)/(length(t)-1));
dddP = [dddP(:,1) dddP];
dddG = diff(ddG')'/(t(end)/(length(t)-1));
dddG = [dddG(:,1) dddG];



% integrands

% Vel_Joint  = sum((dq.^2));
Accel_Joint = (sum(ddq.^2));
Jerk_Joint = (sum(dddq.^2));
% 
Accel_Task  = (sum(ddG.^2));
Jerk_Task  = (sum(dddG.^2));
% 
% Time = 1;


%% DYNAMIC COSTS
du = diff(u')'/(t(end)/(length(t)-1));
du = [du(:,1) du];

% [~, du, ddu] = Cubic_Bspline(t' , u');
% du = du';
% ddu= ddu';

ddu= diff(du')'/(t(end)/(length(t)-1));
ddu= [ddu(:,1) ddu];
% 
[Fx, Fy] = contactForces(q,dq,ddq,p);
dFx = diff(Fx')'/(t(end)/(length(t)-1));
dFx = [dFx(:,1) dFx];
dFy = diff(Fy')'/(t(end)/(length(t)-1));
dFy = [dFy(:,1) dFy];


% integrands
Torque_Squared           = (sum(u.^2));
Torque_Absolute          = sum(abs(u));
TorqueRate_Squared       = (sum(du.^2));
TorqueRate_Absolute      = sum(abs(du));
TorqueRateChange_Squared = (sum(ddu.^2));
TorqueRateChange_Absolute= sum(abs(ddu));

yGRF_Rate    =  (dFy.^2).^.5;
xGRF_Rate    =  (dFx.^2).^.5;


%% KINETIC/ENERGETIC COSTS

% integrands
Work_Positive = sum(max(u.*dq,0));
Work_Absolute = sum(abs(u.*dq));

% I = [p.I1; p.I2; p.I3; p.I4; p.I5];
% m = [p.m1; p.m2; p.m3; p.m4; p.m5];
% l = [p.l1; p.l2; p.l3; p.l4; p.l5];
% c = [p.c1; p.c2; p.c3; p.c4; p.c5];

Kinetic_Energy = 1/2*(p.I1+p.m1*(p.l1-p.c1)^2)*dq(1,:).^2 + ...
                  1/2*(p.I2+p.m2*(p.l2-p.c2)^2)*dq(2,:).^2 + ...
                  1/2*(p.I3+p.m3*(p.l3-p.c3)^2)*dq(3,:).^2 + ... ;   % Kinetic Energy
                  1/2*(p.I4+p.m4*(p.c4)^2)*dq(4,:).^2 + ...
                  1/2*(p.I5+p.m5*(p.c5)^2)*dq(5,:).^2 ;
Kinetic_Energy_Peak = max(Kinetic_Energy);
% 
% 
Angular_Momentum = (p.I1+p.m1*(p.l1-p.c1)^2)*dq(1,:) + ...
                  (p.I2+p.m2*(p.l2-p.c2)^2)*dq(2,:) + ...
                  (p.I3+p.m3*(p.l3-p.c3)^2)*dq(3,:) + ... ;   % Kinetic Energy
                  (p.I4+p.m4*(p.c4)^2)*dq(4,:) + ...
                  (p.I5+p.m5*(p.c5)^2)*dq(5,:) ;
Angular_Momentum_pk2pk = ones(1,size(q,2))*(max(Angular_Momentum)-min(Angular_Momentum))/t(end);
%% Composite Objective              
% dObj = Torque_Squared + TorqueRate_Squared ;
% 
% 
% 
% 
% % COST COMPONENTS EVALS
% kinematic
C(1,:) = Accel_Joint;
C(2,:) = Jerk_Joint;
C(3,:) = Accel_Task;
C(4,:) = Jerk_Task;
% C(5,:) = 1;

% dynamic
C(5,:) = Torque_Squared;
C(6,:) = Torque_Absolute;
C(7,:) = TorqueRate_Squared;
C(8,:) = TorqueRate_Absolute;
C(9,:) = TorqueRateChange_Squared;
C(10,:) = TorqueRateChange_Absolute;
C(11,:) = yGRF_Rate;
C(12,:) = xGRF_Rate;

% energetic
C(13,:) = Work_Positive;
C(14,:) = Work_Absolute;
C(15,:) = Kinetic_Energy_Peak/t(end);

% balance / stability
C(16,:) = (max(Angular_Momentum)-min(Angular_Momentum))/t(end);
% 
% 
% dObj = W*(C./scaling');
% 
% dObj_TimeSeries = C;


%% for 3 costs
cost_vector = [50* Torque_Squared  ;...                       % torques squared
                1* TorqueRate_Squared ;...                       % torque rate squared
                .01* Jerk_Joint];            % Jerk
% % cost_vector = [50* Torque_Squared  ;...                       % torques squared
% %                 1* TorqueRate_Squared ;...                       % torque rate squared
% %                 5* Vel_Joint];            % Vel                       

cost_vector = [1 * Work_Positive  ;...                       % torques squared
                1* TorqueRate_Squared ;...                       % torque rate squared
                .1* Jerk_Joint];            % Jerk
            
% cost_vector = [10 * TorqueRate_Squared  ;...                       % torques squared
%                 100* Torque_Absolute ;...                       % torque rate squared
%                 1* Angular_Momentum_pk2pk];            % Jerk






%% Load Bases Structure

Cost_Vector = [C];

% Cost_Vector_scaled = (Cost_Vector - center) ./ scale  +  1 ;
% Cost_Vector_scaled = (Cost_Vector - center) ./ scale  + 10;

% Bases = loadings' * ( Cost_Vector ) ;

% Bases = loadings' * ( Cost_Vector_scaled ) ;



dObj = Cost_Vector(cost_index, :) / t(end);

% % dObj = Torque_Squared + .02*TorqueRate_Squared;






%%
% if nargout == 1 % numerical gradients
%     
%     dObj = autoGen_obj_torqueSquared(u(1,:),u(2,:),u(3,:),u(4,:),u(5,:));
%     
% else  %Analytic gradients
%     
%     [dObj,fz,fzi] = autoGen_obj_torqueSquared(u(1,:),u(2,:),u(3,:),u(4,:),u(5,:));
%     dObjGrad = zeros(16,length(dObj));  % 16 = 1 + 5 + 5 + 5 = time + angle + rate + control
%     dObjGrad(fzi,:) = fz;
%     
% end

end