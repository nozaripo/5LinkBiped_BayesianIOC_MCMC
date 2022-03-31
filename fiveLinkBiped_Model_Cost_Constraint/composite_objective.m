function [dObj, dObjGrad] = composite_objective(t,x,u,p)
% [dObj, dObjGrad] = obj_torqueSquared(u)
%
% This function computes the torque-squared objective function and its
% gradients.
%

%% KINEMATIC COSTS
% Joints
q  = x(1:5,:);
dq = x(6:10,:);

ddq = diff(dq')'/(t(end)/(length(t)-1));
ddq = [ddq(:,1) ddq];

dddq = diff(ddq')'/(t(end)/(length(t)-1));
dddq = [dddq(:,1) dddq];



% Task
[P, G] = getPoints(q,p);

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
Accel_Joint = (sum(ddq.^2));
Jerk_Joint = (sum(dddq.^2)).^.5;

Accel_Task  = (sum(ddG.^2)).^.5;
Jerk_Task  = (sum(dddG.^2)).^.5;

Time = ones(1,length(t));


%% DYNAMIC COSTS
du = diff(u')'/(t(end)/(length(t)-1));
du = [du(:,1) du];
ddu= diff(du')'/(t(end)/(length(t)-1));
ddu= [ddu(:,1) ddu];

[Fx, Fy] = contactForces(q,dq,ddq,p);
dFx = diff(Fx')'/(t(end)/(length(t)-1));
dFx = [dFx(:,1) dFx];
dFy = diff(Fy')'/(t(end)/(length(t)-1));
dFy = [dFy(:,1) dFy];


% integrands
Torque_Squared           = (sum(u.^2)).^.5;
Torque_Absolute          = sum(abs(u));
TorqueRate_Squared       = (sum(du.^2)).^.5;
TorqueRate_Absolute      = sum(abs(du));
TorqueRateChange_Squared = (sum(ddu.^2)).^.5;
TorqueRateChange_Absolute= sum(abs(ddu));

yGRF_Rate    =  (dFy.^2).^.5;
xGRF_Rate    =  (dFx.^2).^.5;


%% KINETIC/ENERGETIC COSTS

% integrands
Work_Absolute = sum(abs(u.*dq));
Work_Positive = sum(max(u.*dq,0));

% I = [p.I1; p.I2; p.I3; p.I4; p.I5];
% m = [p.m1; p.m2; p.m3; p.m4; p.m5];
% l = [p.l1; p.l2; p.l3; p.l4; p.l5];
% c = [p.c1; p.c2; p.c3; p.c4; p.c5];

Kinetic_Energy = 1/2*(p.I1+p.m1*(p.l1-p.c1)^2)*dq(1,:).^2 + ...
                  1/2*(p.I2+p.m2*(p.l2-p.c2)^2)*dq(2,:).^2 + ...
                  1/2*(p.I3+p.m3*(p.l3-p.c3)^2)*dq(3,:).^2 + ... ;   % Kinetic Energy
                  1/2*(p.I4+p.m4*(p.c4)^2)*dq(4,:).^2 + ...
                  1/2*(p.I5+p.m5*(p.c5)^2)*dq(5,:).^2 ;


%% Composite Objective              

dObj = Torque_Squared + TorqueRate_Squared ;
dObj = Time;
dObj = Work_Positive;
dObj = Torque_Squared ;

% dObj = Work_Positive + 100*Torque_Squared + TorqueRate_Squared;
              
% dObj = Accel_Joint;
% 
% dObj = max(Kinetic_Energy)/t(end)*ones(1,length(t));

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