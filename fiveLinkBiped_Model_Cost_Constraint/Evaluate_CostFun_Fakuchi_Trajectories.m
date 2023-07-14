% function [cost cost_bases_eval] = Evaluate_CostFun(t,kin,dyn,u,p,W)
% cost = costFun(u)
% scaling_vector(scaling_vector==0)=1;
%%%%%%%%%%%%%%%%%%%%  Cost Function Terms  %%%%%%%%%%%%%%%%%%%
% [torques squared integrated;
%  accelerations squared integrated;
%  net work;
%  1/margin of stability]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cost is the integral of torque-squared.

clear all
close all


addpath ../data/

dt = .01;


%% Parameters
param.dyn.m = 10;  %leg mass
param.dyn.I = .06;  %leg inertia about CoM
param.dyn.g = 9.81;  %gravity
param.dyn.l = .76;  %leg length
param.dyn.d = 0.3;  %Leg CoM distance from hip


param = getPhysicalParameters();

p.m = 10;  %leg mass
p.I = .06;  %leg inertia about CoM
p.g = 9.81;  %gravity
p.l = .76;  %leg length
p.d = 0.3;  %Leg CoM distance from hip

param = getPhysicalParameters();
p=param;

%% Load Demonstrations Data
Demo = load('Demos_for_IOC.mat');
% Demo = load('Baseline.mat');

Time = Demo.Baseline.time;
t = Time;

% t = linspace()

% GRFs
GRF_r_N = Demo.grf_r_N;
GRF_l_N = Demo.grf_l_N;
GRF_r_f = Demo.grf_r_f;
GRF_l_f = Demo.grf_l_f;

GRF_r_N = lowpass(GRF_r_N,50,100);
GRF_l_N = lowpass(GRF_l_N,50,100);
GRF_r_f = lowpass(GRF_r_f,50,100);
GRF_l_f = lowpass(GRF_l_f,50,100);

GRF = [GRF_r_N GRF_l_N GRF_r_f GRF_l_f];
dGRF = diff(GRF)/dt;
dGRF = [dGRF(1,:);dGRF];

l_hip   = -Demo.Baseline.L_Hip_Angle(:,1)*pi/180;
r_hip   = -Demo.Baseline.R_Hip_Angle(:,1)*pi/180;
l_knee  = Demo.Baseline.L_Knee_Angle(:,1)*pi/180;
r_knee  = Demo.Baseline.R_Knee_Angle(:,1)*pi/180;

l_hip   = -Demo.Baseline.L_Hip_Angle(:,1);
r_hip   = -Demo.Baseline.R_Hip_Angle(:,1);
l_knee  = Demo.Baseline.L_Knee_Angle(:,1);
r_knee  = Demo.Baseline.R_Knee_Angle(:,1);

angles = [l_hip , r_hip , l_knee , r_knee];













% Hip Angles/Torques
% QL = Demo.Baseline.L_Hip_Angle(:,1)*pi/180;
% QR = Demo.Baseline.R_Hip_Angle(:,1)*pi/180;
% 
% QL = lowpass(QL,6,100);
% QR = lowpass(QR,6,100);

l_hip   = lowpass(l_hip,6,100);
r_hip   = lowpass(r_hip,6,100);
l_knee  = lowpass(l_knee,6,100);
r_knee  = lowpass(r_knee,6,100);

l_knee(l_knee<0)=0;
r_knee(r_knee<0)=0;

% joint kinematics
q = [l_hip-l_knee , l_hip , zeros(length(l_hip),1) , r_hip , r_hip-r_knee]';

dq = diff(q')'/(Time(end)/(length(Time)-1));
dq = [dq(:,1) dq];

ddq = diff(dq')'/(Time(end)/(length(Time)-1));
ddq = [ddq(:,1) ddq];

dddq = diff(ddq')'/(Time(end)/(length(Time)-1));
dddq = [dddq(:,1) dddq];


% task kinematics
[P, G] = getPoints(t',q,p,1);

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


%% center as [0 ; 0]
% % r1 = p.l/2*[sin(r_hip),-cos(r_hip)];
% % l1 = p.l/2*[sin(l_hip),-cos(l_hip)];
% % r2 = r1 + p.l/2*[sin(r_hip-r_knee),-cos(r_hip-r_knee)];
% % l2 = l1 + p.l/2*[sin(l_hip-l_knee),-cos(l_hip-l_knee)];
% % 
% % QL = atan2(l2(:,2),l2(:,1))-3*pi/2;
% % QR = atan2(r2(:,2),r2(:,1))-3*pi/2;
%%%
% Fs = 100;
% FF = fft(y);
% 
% P2 = abs(FF/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% f = Fs*(0:(L/2))/L;
% figure (1)
% plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')

%%%%
x = [q ; dq];
u = dynamics_forward(t,x,ddq,p);


du = diff(u')'/(t(end)/(length(t)-1));
du = [du(:,1) du];
ddu= diff(du')'/(t(end)/(length(t)-1));
ddu= [ddu(:,1) ddu];

[Fx, Fy] = contactForces(q,dq,ddq,p);
dFx = diff(Fx')'/(t(end)/(length(t)-1));
dFx = [dFx(:,1) dFx];
dFy = diff(Fy')'/(t(end)/(length(t)-1));
dFy = [dFy(:,1) dFy];







% U = 


%% Detect Cycles
j = 1; k = 1; m = 1;
len = length(GRF_l_f);
for i = 1:len-1
    if GRF_r_N(i)>15 && GRF_r_N(i+1)<15
        r_TO(j,1)=i;
        j = j+1;
    end
    if GRF_r_N(i)<15 && GRF_r_N(i+1)>15
        r_HS(k,1)=i;
        k = k+1;
    end
    
    if GRF_l_N(i)>15 && GRF_l_N(i+1)<15
        l_TO(m,1)=i;
        m = m+1;
    end
end

if r_TO(1)>r_HS(1)
    r_HS(1) = [];
end

if length(r_TO)>length(r_HS)
    r_TO(end)=[];
end



if r_TO(1)>l_TO(1)
    l_TO(1) = [];
end

if length(r_TO)>length(l_TO)
    r_TO(end)=[];
end

step_time = (l_TO-r_TO);


% st_err = zeros(length(l_TO),1);
% for num = 1:length(l_TO)
%     num
%     for g=1:6
%         [~, ~, ~, y] = Trajectory(g,r_TO(num),l_TO(num));
%         st_err(num,1) = st_err(num,1) + norm(y(1)-y(end));
%     end
% end


%% 



% ddu = [zeros(1,2) diff(u',2)'/(t(end)/(length(t)-1)).^2];

% % % % % for i=1:length(step_time)
% % % % %     
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % time = Time(r_TO(i):r_HS(i)-1,1)-Time(r_TO(i),1);
% % % % % 
% % % % % 
% % % % % cycl_err(i,1) = norm(q(:,r_TO(i))-q(:,r_HS(i)));
% % % % % 
% % % % % 
% % % % %         
% % % % %         
% % % % %         % .^.5
% % % % %         
% % % % % Accel_Joint = (sum(ddq(:,r_TO(i):r_HS(i)-1).^2));
% % % % % Jerk_Joint = (sum(dddq(:,r_TO(i):r_HS(i)-1).^2));
% % % % % 
% % % % % Accel_Task  = (sum(ddG(:,r_TO(i):r_HS(i)-1).^2));
% % % % % Jerk_Task  = (sum(dddG(:,r_TO(i):r_HS(i)-1).^2));
% % % % % 
% % % % % % Time = 1;
% % % % % 
% % % % % % Velocity_Deviation = dP(2,:)
% % % % % 
% % % % % 
% % % % % Torque_Squared           = (sum(u(:,r_TO(i):r_HS(i)-1).^2));
% % % % % Torque_Absolute          = sum(abs(u(:,r_TO(i):r_HS(i)-1)));
% % % % % TorqueRate_Squared       = (sum(du(:,r_TO(i):r_HS(i)-1).^2));
% % % % % TorqueRate_Absolute      = sum(abs(du(:,r_TO(i):r_HS(i)-1)));
% % % % % TorqueRateChange_Squared = (sum(ddu(:,r_TO(i):r_HS(i)-1).^2));
% % % % % TorqueRateChange_Absolute= sum(abs(ddu(:,r_TO(i):r_HS(i)-1)));
% % % % % 
% % % % % yGRF_Rate    =  (dFy(:,r_TO(i):r_HS(i)-1).^2);
% % % % % xGRF_Rate    =  (dFx(:,r_TO(i):r_HS(i)-1).^2);
% % % % % 
% % % % % 
% % % % % Work_Absolute = sum(abs(u(:,r_TO(i):r_HS(i)-1).*dq(:,r_TO(i):r_HS(i)-1)));
% % % % % Work_Positive = sum(max(u(:,r_TO(i):r_HS(i)-1).*dq(:,r_TO(i):r_HS(i)-1),0));
% % % % % 
% % % % % Kinetic_Energy = 1/2*(p.I1+p.m1*(p.l1-p.c1)^2)*dq(1,r_TO(i):r_HS(i)-1).^2 + ...
% % % % %                   1/2*(p.I2+p.m2*(p.l2-p.c2)^2)*dq(2,r_TO(i):r_HS(i)-1).^2 + ...
% % % % %                   1/2*(p.I3+p.m3*(p.l3-p.c3)^2)*dq(3,r_TO(i):r_HS(i)-1).^2 + ... ;   % Kinetic Energy
% % % % %                   1/2*(p.I4+p.m4*(p.c4)^2)*dq(4,r_TO(i):r_HS(i)-1).^2 + ...
% % % % %                   1/2*(p.I5+p.m5*(p.c5)^2)*dq(5,r_TO(i):r_HS(i)-1).^2 ;
% % % % % 
% % % % % Angular_Momentum = (p.I1+p.m1*(p.l1-p.c1)^2)*dq(1,r_TO(i):r_HS(i)-1) + ...
% % % % %                   (p.I2+p.m2*(p.l2-p.c2)^2)*dq(2,r_TO(i):r_HS(i)-1) + ...
% % % % %                   (p.I3+p.m3*(p.l3-p.c3)^2)*dq(3,r_TO(i):r_HS(i)-1) + ... ;   % Kinetic Energy
% % % % %                   (p.I4+p.m4*(p.c4)^2)*dq(4,r_TO(i):r_HS(i)-1) + ...
% % % % %                   (p.I5+p.m5*(p.c5)^2)*dq(5,r_TO(i):r_HS(i)-1) ;
% % % % % 
% % % % %               
% % % % % % COST COMPONENTS EVALS
% % % % % 
% % % % % % kinematic
% % % % % cost_components_eval(i,1) = simps(time, Accel_Joint');
% % % % % cost_components_eval(i,2) = simps(time, Jerk_Joint');
% % % % % cost_components_eval(i,3) = simps(time, Accel_Task');
% % % % % cost_components_eval(i,4) = simps(time, Jerk_Task');
% % % % % cost_components_eval(i,5) = time(end);
% % % % % 
% % % % % % dynamic
% % % % % cost_components_eval(i,6) = simps(time, Torque_Squared');
% % % % % cost_components_eval(i,7) = simps(time, Torque_Absolute');
% % % % % cost_components_eval(i,8) = simps(time, TorqueRate_Squared');
% % % % % cost_components_eval(i,9) = simps(time, TorqueRate_Absolute');
% % % % % cost_components_eval(i,10) = simps(time, TorqueRateChange_Squared');
% % % % % cost_components_eval(i,11) = simps(time, TorqueRateChange_Absolute');
% % % % % cost_components_eval(i,12) = simps(time, yGRF_Rate');
% % % % % cost_components_eval(i,13) = simps(time, xGRF_Rate');
% % % % % 
% % % % % % energetic
% % % % % cost_components_eval(i,14) = simps(time, Work_Positive');
% % % % % cost_components_eval(i,15) = simps(time, Work_Absolute');
% % % % % cost_components_eval(i,16) = max(Kinetic_Energy)/time(end);
% % % % % 
% % % % % % balance / stability
% % % % % cost_components_eval(i,17) = (max(Angular_Momentum)-min(Angular_Momentum))/time(end);       
% % % % %         
% % % % %         
% % % % % % cost_bases_eval(i,:) = simps(time,C');
% % % % % %        
% % % % % % % cost_vector = C./scaling_vector;
% % % % % % clear C
% % % % % clear time
% % % % % end
% % % % % 
% % % % % dlmwrite('Cost_Comp_Eval.txt',cost_components_eval)
% % % % % % cost_vector = [.15* u.^2  ;...                          % torques squared
% % % % % %                 .01* du.^2 ;...                         % Torque Change squared
% % % % % %                     .4* abs(u.*(dq1-dq2))  ;...         % absolute work
% % % % % %                         1/40* (ddp1.^2+ddp2.^2);...     % Accel squared
% % % % % %                            .01* (d3p1.^2+d3p2.^2)];     % Jerk
% % % % % 
% % % % % % cost_vector = [50* u.^2  ;...                          % torques squared
% % % % % %                 1* du.^2 ;...                         % Torque Change squared
% % % % % %                     10* abs(u.*(dq1-dq2))  ;...         % absolute work
% % % % % %                         .25* (ddp1.^2+ddp2.^2);...     % Accel squared
% % % % % %                            .01* (d3p1.^2+d3p2.^2)];     % Jerk
% % % % % 
% % % % % len_cost = size(cost_components_eval,2)
% % % % % for i = 1:len_cost
% % % % %     for j = i:len_cost
% % % % %         
% % % % % 
% % % % %         Rsq(i,j) = corr2(cost_components_eval(:,i),cost_components_eval(:,j));
% % % % %         
% % % % %         
% % % % %     end
% % % % % end
% % % % % 
% % % % % 
% % % % % 
% % % % % figure(9)
% % % % % imagesc(abs(Rsq))
% % % % % % impixelregion(imagesc(Rsq))
% % % % % set(gca, 'XTIck', [1:len_cost])
% % % % % set(gca, 'YTIck', [1:len_cost])
% % % % % set(gca, 'XAxisLocation', 'top')
% % % % % set(gca, 'YAxisLocation', 'right')
% % % % % colorbar('Location','eastoutside')
% % % % % colorbar('Location','southoutside')
% % % % % colormap Bone
% % % % % 
% % % % % textStrings = num2str(abs(Rsq), '%0.2f');       % Create strings from the matrix values
% % % % % % textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
% % % % % text_modif = strings(len_cost,len_cost);
% % % % % 
% % % % % for k = 1:len_cost
% % % % %     for i = 1:len_cost
% % % % %     if Rsq(i,k)<0
% % % % %         txt_sgn = '-';
% % % % %     else
% % % % %         txt_sgn = '';
% % % % %     end
% % % % %     text_modif (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
% % % % %     end
% % % % % end
% % % % % 
% % % % % [xs, ys] = meshgrid(1:len_cost);  % Create x and y coordinates for the strings
% % % % % text(xs(:), ys(:), text_modif(:), ...  % Plot the strings
% % % % %                 'HorizontalAlignment', 'center');
% % % % %             
% % % % % 
% % % % %  set( gca, 'XTickLabel', {'1.angl accel sqrd' '2.angl jerk sqrd' '3.cartes accel sqrd' '4.cartes jerk sqrd' '5.time' '6.torqs sqrd' '7.torqs abs' '8.torq rate sqrd' '9.torq rate abs' '10.torq rate change sqrd' '11.torq rate change abs' '12.yGRF rate sqrd' '13.xGRF rate sqrd' '14.positive work' '15.absolute work' '16.kinetic energy' '17.angular momentum'}','FontSize',12 )
% % % % % xtickangle(70)           
% % % % % 
% % % % %  set( gca, 'YTickLabel', {'angl accel sqrd' 'angl jerk sqrd' 'cartes accel sqrd' 'cartes jerk sqrd' 'time' 'torqs sqrd' 'torqs abs' 'torq rate sqrd' 'torq rate abs' 'torq rate change sqrd' 'torq rate change abs' 'yGRF rate sqrd' 'xGRF rate sqrd' 'positive work' 'absolute work' 'kinetic energy' 'angular momentum'}','FontSize',12 )
% % % % % 
% % % % %             
% % % % %             
% % % % % % % %  set( gca, 'XTickLabel', {'torqs sqrd' 'torqs abs' 'torq chng sqrd' 'torq chng abs' 'torq 2der sqrd' 'torq 2der abs' 'abs work' 'positv work' 'cartes accel sqrd' 'cartes jerk sqrd' 'angl accel sqrd' 'angl jerk sqrd' 'time' 'dev from vel_{ref}' 'kinetic energy' 'kinetic energy2'}','FontSize',12 )
% % % % % % % % xtickangle(60)           
% % % % % % % % 
% % % % % % % %  set( gca, 'YTickLabel', {'torqs sqrd' 'torqs abs' 'torq chng sqrd' 'torq chng abs' 'torq 2der sqrd' 'torq 2der abs' 'abs work' 'positv work' 'cartes accel sqrd' 'cartes jerk sqrd' 'angl accel sqrd' 'angl jerk sqrd' 'time' 'dev from vel_{ref}' 'kinetic energy' 'kinetic energy2'}','FontSize',12 )
% % % % % 
% % % % % 
% % % % % % figure(10)
% % % % % % corrplot(cost_components_eval)
% % % % % %                     
% % % % % % cost = W*cost_bases_eval'; 
% % % % % 
% % % % % % end
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % %% correlation based on trajecotry
% % % % % for i = 1:length(r_TO)
% % % % % % cycl_err(i,1) = norm(q(:,r_TO(i))-q(5:-1:1,r_HS(i)));
% % % % % cycl_err(i,1) = norm([l_hip(r_TO(i))-r_hip(r_HS(i)); l_knee(r_TO(i))-r_knee(r_HS(i))]);
% % % % % end
% % % % % indx = find(cycl_err == min(cycl_err));
% % % % % 
% % % % % cycl_err(cycl_err == min(cycl_err))
% % % % % 
% % % % % 
% % % % % 
% % % % % for i=indx:indx
% % % % %     
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % time = Time(r_TO(i):r_HS(i)-1,1)-Time(r_TO(i),1);
% % % % % 
% % % % % 
% % % % % cycl_err(i,1) = norm(q(:,r_TO(i))-q(:,r_HS(i)));
% % % % % 
% % % % % 
% % % % %         
% % % % %         
% % % % %         % .^.5
% % % % %         
% % % % % Accel_Joint = (sum(ddq(:,r_TO(i):r_HS(i)-1).^2));
% % % % % Jerk_Joint = (sum(dddq(:,r_TO(i):r_HS(i)-1).^2));
% % % % % 
% % % % % Accel_Task  = (sum(ddG(:,r_TO(i):r_HS(i)-1).^2));
% % % % % Jerk_Task  = (sum(dddG(:,r_TO(i):r_HS(i)-1).^2));
% % % % % 
% % % % % % Time = 1;
% % % % % 
% % % % % % Velocity_Deviation = dP(2,:)
% % % % % 
% % % % % 
% % % % % Torque_Squared           = (sum(u(:,r_TO(i):r_HS(i)-1).^2));
% % % % % Torque_Absolute          = sum(abs(u(:,r_TO(i):r_HS(i)-1)));
% % % % % TorqueRate_Squared       = (sum(du(:,r_TO(i):r_HS(i)-1).^2));
% % % % % TorqueRate_Absolute      = sum(abs(du(:,r_TO(i):r_HS(i)-1)));
% % % % % TorqueRateChange_Squared = (sum(ddu(:,r_TO(i):r_HS(i)-1).^2));
% % % % % TorqueRateChange_Absolute= sum(abs(ddu(:,r_TO(i):r_HS(i)-1)));
% % % % % 
% % % % % yGRF_Rate    =  (dFy(:,r_TO(i):r_HS(i)-1).^2);
% % % % % xGRF_Rate    =  (dFx(:,r_TO(i):r_HS(i)-1).^2);
% % % % % 
% % % % % 
% % % % % Work_Absolute = sum(abs(u(:,r_TO(i):r_HS(i)-1).*dq(:,r_TO(i):r_HS(i)-1)));
% % % % % Work_Positive = sum(max(u(:,r_TO(i):r_HS(i)-1).*dq(:,r_TO(i):r_HS(i)-1),0));
% % % % % 
% % % % % Kinetic_Energy = 1/2*(p.I1+p.m1*(p.l1-p.c1)^2)*dq(1,r_TO(i):r_HS(i)-1).^2 + ...
% % % % %                   1/2*(p.I2+p.m2*(p.l2-p.c2)^2)*dq(2,r_TO(i):r_HS(i)-1).^2 + ...
% % % % %                   1/2*(p.I3+p.m3*(p.l3-p.c3)^2)*dq(3,r_TO(i):r_HS(i)-1).^2 + ... ;   % Kinetic Energy
% % % % %                   1/2*(p.I4+p.m4*(p.c4)^2)*dq(4,r_TO(i):r_HS(i)-1).^2 + ...
% % % % %                   1/2*(p.I5+p.m5*(p.c5)^2)*dq(5,r_TO(i):r_HS(i)-1).^2 ;
% % % % % 
% % % % % Angular_Momentum = (p.I1+p.m1*(p.l1-p.c1)^2)*dq(1,r_TO(i):r_HS(i)-1) + ...
% % % % %                   (p.I2+p.m2*(p.l2-p.c2)^2)*dq(2,r_TO(i):r_HS(i)-1) + ...
% % % % %                   (p.I3+p.m3*(p.l3-p.c3)^2)*dq(3,r_TO(i):r_HS(i)-1) + ... ;   % Kinetic Energy
% % % % %                   (p.I4+p.m4*(p.c4)^2)*dq(4,r_TO(i):r_HS(i)-1) + ...
% % % % %                   (p.I5+p.m5*(p.c5)^2)*dq(5,r_TO(i):r_HS(i)-1) ;
% % % % % 
% % % % %               
% % % % % % COST COMPONENTS EVALS
% % % % % 
% % % % % % kinematic
% % % % % cost_components_traj_eval(:,1) = Accel_Joint';
% % % % % cost_components_traj_eval(:,2) = Jerk_Joint';
% % % % % cost_components_traj_eval(:,3) = Accel_Task';
% % % % % cost_components_traj_eval(:,4) = Jerk_Task';
% % % % % cost_components_traj_eval(:,5) = time;
% % % % % 
% % % % % % dynamic
% % % % % cost_components_traj_eval(:,6) = Torque_Squared';
% % % % % cost_components_traj_eval(:,7) = Torque_Absolute';
% % % % % cost_components_traj_eval(:,8) = TorqueRate_Squared';
% % % % % cost_components_traj_eval(:,9) = TorqueRate_Absolute';
% % % % % cost_components_traj_eval(:,10) = TorqueRateChange_Squared';
% % % % % cost_components_traj_eval(:,11) = TorqueRateChange_Absolute';
% % % % % cost_components_traj_eval(:,12) = yGRF_Rate';
% % % % % cost_components_traj_eval(:,13) = xGRF_Rate';
% % % % % 
% % % % % % energetic
% % % % % cost_components_traj_eval(:,14) = Work_Positive';
% % % % % cost_components_traj_eval(:,15) = Work_Absolute';
% % % % % cost_components_traj_eval(:,16) = Kinetic_Energy';
% % % % % 
% % % % % % balance / stability
% % % % % cost_components_traj_eval(:,17) = max(Angular_Momentum)-min(Angular_Momentum);       
% % % % %         
% % % % %         
% % % % % % cost_bases_eval(i,:) = simps(time,C');
% % % % % %        
% % % % % % % cost_vector = C./scaling_vector;
% % % % % % clear C
% % % % % clear time
% % % % % end
% % % % % 
% % % % % dlmwrite('Cost_Comp_Eval_Traj.txt',cost_components_traj_eval)
% % % % % % cost_vector = [.15* u.^2  ;...                          % torques squared
% % % % % %                 .01* du.^2 ;...                         % Torque Change squared
% % % % % %                     .4* abs(u.*(dq1-dq2))  ;...         % absolute work
% % % % % %                         1/40* (ddp1.^2+ddp2.^2);...     % Accel squared
% % % % % %                            .01* (d3p1.^2+d3p2.^2)];     % Jerk
% % % % % 
% % % % % % cost_vector = [50* u.^2  ;...                          % torques squared
% % % % % %                 1* du.^2 ;...                         % Torque Change squared
% % % % % %                     10* abs(u.*(dq1-dq2))  ;...         % absolute work
% % % % % %                         .25* (ddp1.^2+ddp2.^2);...     % Accel squared
% % % % % %                            .01* (d3p1.^2+d3p2.^2)];     % Jerk
% % % % % 
% % % % % len_cost = size(cost_components_traj_eval,2)
% % % % % for i = 1:len_cost
% % % % %     for j = i:len_cost
% % % % %         
% % % % % 
% % % % %         Rsq(i,j) = corr2(cost_components_traj_eval(:,i),cost_components_traj_eval(:,j));
% % % % %         
% % % % %         
% % % % %     end
% % % % % end
% % % % % 
% % % % % 
% % % % % 
% % % % % figure(10)
% % % % % imagesc(abs(Rsq))
% % % % % % impixelregion(imagesc(Rsq))
% % % % % set(gca, 'XTIck', [1:len_cost])
% % % % % set(gca, 'YTIck', [1:len_cost])
% % % % % set(gca, 'XAxisLocation', 'top')
% % % % % set(gca, 'YAxisLocation', 'right')
% % % % % colorbar('Location','eastoutside')
% % % % % colorbar('Location','southoutside')
% % % % % colormap Bone
% % % % % 
% % % % % textStrings = num2str(abs(Rsq), '%0.2f');       % Create strings from the matrix values
% % % % % % textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
% % % % % text_modif = strings(len_cost,len_cost);
% % % % % 
% % % % % for k = 1:len_cost
% % % % %     for i = 1:len_cost
% % % % %     if Rsq(i,k)<0
% % % % %         txt_sgn = '-';
% % % % %     else
% % % % %         txt_sgn = '';
% % % % %     end
% % % % %     text_modif (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
% % % % %     end
% % % % % end
% % % % % 
% % % % % [xs, ys] = meshgrid(1:len_cost);  % Create x and y coordinates for the strings
% % % % % text(xs(:), ys(:), text_modif(:), ...  % Plot the strings
% % % % %                 'HorizontalAlignment', 'center');
% % % % %             
% % % % % 
% % % % %  set( gca, 'XTickLabel', {'1.angl accel sqrd' '2.angl jerk sqrd' '3.cartes accel sqrd' '4.cartes jerk sqrd' '5.time' '6.torqs sqrd' '7.torqs abs' '8.torq rate sqrd' '9.torq rate abs' '10.torq rate change sqrd' '11.torq rate change abs' '12.yGRF rate sqrd' '13.xGRF rate sqrd' '14.positive work' '15.absolute work' '16.kinetic energy' '17.angular momentum'}','FontSize',12 )
% % % % % xtickangle(70)           
% % % % % 
% % % % %  set( gca, 'YTickLabel', {'angl accel sqrd' 'angl jerk sqrd' 'cartes accel sqrd' 'cartes jerk sqrd' 'time' 'torqs sqrd' 'torqs abs' 'torq rate sqrd' 'torq rate abs' 'torq rate change sqrd' 'torq rate change abs' 'yGRF rate sqrd' 'xGRF rate sqrd' 'positive work' 'absolute work' 'kinetic energy' 'angular momentum'}','FontSize',12 )
% % % % % 
% % % % %             
% % % % % 
% % % % %  
% % % % %  
% % % % %  
% % % % %  
% % % % %  
% % % % %  
% % % % %  
 
 
 
 
 clear all
 
 
WBDS_Table = readtable('../data/WBDSinfo.xlsx','ReadVariableNames',true);
 
 param = getPhysicalParameters();
p=param;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fakuchi Data

% % % % % % 
% % % % % % if varargin{1} == 1
% % % % % %     y1 = circshift(repmat(L_Hip_Angle,rep,1),shift1);
% % % % % % elseif varargin{1}==2
% % % % % %     y1 = circshift(repmat(L_Knee_Angle,rep,1),shift1);
% % % % % % elseif varargin{1}==3
% % % % % %     y1 = circshift(repmat(L_Ankle_Angle,rep,1),shift1);
% % % % % % elseif varargin{1}==4
% % % % % % %     y1 = repmat(R_Hip_Angle(tt1:tt2,1),rep,1);
% % % % % % 
% % % % % %     y1 = circshift(repmat(L_Hip_Angle,rep,1),shift0+shift1);
% % % % % % 
% % % % % % elseif varargin{1}==5
% % % % % % %     y1 = repmat(R_Knee_Angle(tt1:tt2,1),rep,1);
% % % % % %     
% % % % % %     y1 = circshift(repmat(L_Knee_Angle,rep,1),shift0+shift1);
% % % % % % 
% % % % % % elseif varargin{1}==6
% % % % % % %     y1 = repmat(R_Ankle_Angle(tt1:tt2,1),rep,1);
% % % % % %     
% % % % % %     y1 = circshift(repmat(L_Ankle_Angle,rep,1),shift0+shift1);
% % % % % % end
q1 = [];
q2 = [];
q4 = [];
q5 = [];

No_Subj = 42;
No_Speed= 8;
ID_Eval = 0;
Tr_or_Ov =['T','O'];
ID_Subj_actual = 0;
for ID_Subj = 1:42
    ID_Subj
    if ID_Subj~=[5, 17, 41,  10, 15, 36, 39, 24, 26, 40, 42]
        ID_Subj_actual = ID_Subj_actual+1;
    cycle_stack_id = 0;
    Q1=[]; Q2=[]; Q4=[]; Q5=[];
    
    
    for ID_Speed = 1:No_Speed
        speed_stack_id = 0;
        %         for i = 1:2
        knt_filename = ['WBDS' num2str(ID_Subj,'%02.f') 'walkT' num2str(ID_Speed,'%02.f') 'mkr.txt'];
        if isfile(['../data/' knt_filename])
%             Kinetic_data = importdata(['data/WBDS' num2str(ID_Subj,'%02.f') 'walkT' num2str(ID_Speed,'%02.f') 'knt.txt']);
            Kinetic_data = importdata(['../data/WBDS' num2str(ID_Subj,'%02.f') 'walkT' num2str(ID_Speed,'%02.f') 'grf.txt']);
            rep = 2;
            shift1 = -25;
            shift0 = 50;
            
            V_tr = WBDS_Table.GaitSpeed_m_s_(strcmp(WBDS_Table.FileName, knt_filename));
            
            if isempty(V_tr)
                V_tr = mean(V_tr_all(:, ID_Speed));
            end

            V_tr_all(ID_Subj, ID_Speed) =  V_tr;
            
            
            grf_resamp = interp1((1:9000)',Kinetic_data.data,linspace(1,9000,4500));
            GRF_r_N = max(lowpass(grf_resamp(:,10),.8,100),0);
            GRF_l_N = max(lowpass(grf_resamp(:,3),.8,100),0);
            
%             GRF_r_N = Kinetic_data.data(:,3) * p.bodymass;
%             GRF_l_N = Kinetic_data.data(:,10) * p.bodymass;
%             GRF_r_f = Kinetic_data.data(:,20) * p.bodymass;
%             GRF_l_f = Kinetic_data.data(:,23) * p.bodymass;
%             GRF_r_N  = circshift(GRF_r_N,shift0);
%             GRF_r_f  = circshift(GRF_r_f,shift0);
            
            
%             L_Hip_Torque = Kinetic_data.data(:,7) * p.bodymass;
%             L_Knee_Torque = Kinetic_data.data(:,13) * p.bodymass;
%             R_Hip_Torque = Kinetic_data.data(:,4) * p.bodymass;
%             R_Knee_Torque = Kinetic_data.data(:,10) * p.bodymass;
%             R_Hip_Torque  = circshift(R_Hip_Torque,shift0);
%             R_Knee_Torque  = circshift(R_Knee_Torque,shift0);
            
            
            % figure(996)
            % scatter(ID_Speed ,GRF_l_N(1))
            % hold on
            %     end
            % end
            
            left_TO_idx = 0;
            right_TO_idx = 0;            
            left_HS_idx = 0;
            right_HS_idx = 0;
            len = length(GRF_l_N);
            l_TO=[]; l_HS=[]; r_TO=[]; r_HS=[]; 
            for i = 1:len-1
                if GRF_l_N(i)>50 && GRF_l_N(i+1)<50
                    left_TO_idx=left_TO_idx+1;
                    l_TO(left_TO_idx)=i;
                end
            end
            for i = l_TO(1):len-1
                if GRF_l_N(i)<50 && GRF_l_N(i+1)>50
                    left_HS_idx=left_HS_idx+1;
                    l_HS(left_HS_idx)=i;
                end
            end
            
            for i = 1:len-1
                if GRF_r_N(i)>50 && GRF_r_N(i+1)<50
                    right_TO_idx=right_TO_idx+1;
                    r_TO(right_TO_idx)=i;
                end
            end
            for i = r_TO(1):len-1
                if GRF_r_N(i)<50 && GRF_r_N(i+1)>50
                    right_HS_idx=right_HS_idx+1;
                    r_HS(right_HS_idx)=i;
                end
            end
            
               
                     
            if right_TO_idx > right_HS_idx
                r_TO(end) = [];
            end
            if left_TO_idx > left_HS_idx
                l_TO(end) = [];
            end
                
                
            
            
            
            
            
            
%             % walk_data = importdata('WBDS23walkT03ang.txt');
%             walk_data = importdata(['data/WBDS' num2str(ID_Subj,'%02.f') 'walkT' num2str(ID_Speed,'%02.f') 'ang.txt']);
%             l_hip  = walk_data.data(:,13) * pi/180;
%             l_knee = walk_data.data(:,19) * pi/180;
%             % l_Ankle = walk_data.data(:,25);
%             
%             
%             
%             r_hip  = circshift(l_hip,shift0);
%             r_knee = circshift(l_knee,shift0);
%             
%             r_hip  = walk_data.data(:,10) * pi/180;
%             r_knee = walk_data.data(:,16) * pi/180;
%             r_hip  = circshift(r_hip,shift0);
%             r_knee = circshift(r_knee,shift0);
%             % r_ankle = walk_data.data(:,25);
%             
%             angles = [l_hip , r_hip , l_knee , r_knee];
            
            walk_data = importdata(['data/WBDS' num2str(ID_Subj,'%02.f') 'walkT' num2str(ID_Speed,'%02.f') 'mkr.txt']);
            
            % NaN finder
            NaN_idx = find(isnan(walk_data.data));
            for nan_index = 1:length(NaN_idx)
            walk_data.data(NaN_idx(nan_index,1))=min(walk_data.data(NaN_idx(nan_index,1)-10:NaN_idx(nan_index,1)+10));
            end
            
            cutoff = 6;
            R_AnkY  = lowpass(walk_data.data(:,33),cutoff,100);            R_AnkX  = lowpass(walk_data.data(:,32),cutoff,100);
            R_KneeY = lowpass(walk_data.data(:,24),cutoff,100);            R_KneeX = lowpass(walk_data.data(:,23),cutoff,100);
            R_HipY  = lowpass(walk_data.data(:,21),cutoff,100);            R_HipX  = lowpass(walk_data.data(:,20),cutoff,100);
            L_AnkY  = lowpass(walk_data.data(:,57),cutoff,100);            L_AnkX  = lowpass(walk_data.data(:,56),cutoff,100);
            L_KneeY = lowpass(walk_data.data(:,48),cutoff,100);            L_KneeX = lowpass(walk_data.data(:,47),cutoff,100);
            L_HipY  = lowpass(walk_data.data(:,45),cutoff,100);            L_HipX  = lowpass(walk_data.data(:,44),cutoff,100);
            
            
            R_Shank = atan2 ( R_KneeY - R_AnkY  , R_KneeX - R_AnkX ) - pi/2;
            R_Thigh = atan2 ( R_HipY  - R_KneeY , R_HipX  - R_KneeX) - pi/2;
            L_Shank = atan2 ( L_KneeY - L_AnkY  , L_KneeX - L_AnkX)  - pi/2;
            L_Thigh = atan2 ( L_HipY  - L_KneeY , L_HipX  - L_KneeX) - pi/2;
            
            
            
            
            for limb_side = 1:2
                TO = [];
                HS = [];
                
                
%                                     Limb_Angles = [L_Shank , L_Thigh , zeros(length(L_Thigh),1) , R_Thigh , R_Shank]';

                %     q = [l_hip-l_knee , l_hip , zeros(length(l_hip),1) , r_hip , r_hip-r_knee]';
                if limb_side == 1
                    TO = l_TO;
                    HS = l_HS;
                    Limb_Angles = [R_Shank , R_Thigh , zeros(length(L_Thigh),1) , L_Thigh , L_Shank]';
                    len_cycle = length(l_TO);
                else
                    TO = r_TO;
                    HS = r_HS;
                    Limb_Angles = [L_Shank , L_Thigh , zeros(length(L_Thigh),1) , R_Thigh , R_Shank]';
                    len_cycle = length(r_TO);
                end
                
                %
                % % % % q1 = [q1 ; q(1,:)];
                % % % % q2 = [q2 ; q(2,:)];
                % % % % q4 = [q4 ; q(4,:)];
                % % % % q5 = [q5 ; q(5,:)];
                % % % % end
                % % % % end
                % % % %     end
                % % % % end
                % % % %
                % % % %
                % % % % figure(1003)
                % % % % stdshade(q1,.1,[1 0 1],(0:100),1)
                % % % % hold on
                % % % % stdshade(q2,.1,[1 0 0],(0:100),1)
                % % % % stdshade(q4,.1,[0 0 1],(0:100),1)
                % % % % stdshade(q5,.1,[0 1 1],(0:100),1)
                % % % %
                % % % % legend('','Stance Shank','','Stance Thigh','','Swing Thigh','','Swing Shank');
                % % % % xlabel('Gait Cycle %')
                % % % % ylabel('Segment Angles (rad)')
                % % % % set(gca, ...
                % % % %   'Box'         , 'off'     , ...
                % % % %   'TickDir'     , 'out'     , ...
                % % % %   'TickLength'  , [.02 .02] , ...
                % % % %   'YGrid'       , 'on'      , ...
                % % % %   'XColor'      , [.3 .3 .3], ...
                % % % %   'YColor'      , [.3 .3 .3], ...
                % % % %   'LineWidth'   , 1.4        ,...
                % % % %   'FontSize'    , 14);
                
                
                
                for cycle_idx = 1:len_cycle
                
                
                ID_Eval = ID_Eval+1;
                
                cycle_stack_id = cycle_stack_id+1;
                speed_stack_id = speed_stack_id+1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                q    = Limb_Angles(:,TO(cycle_idx):HS(cycle_idx));
                time = walk_data.data(TO(cycle_idx):HS(cycle_idx),1)' - walk_data.data(TO(cycle_idx),1);
                
%                 time_data_resamp = linspace(0, time(end), 11);
%                 q_data_resamp    = interp1(time', q', time_data_resamp')';
                
                Angles_Data{ID_Subj, ID_Speed, cycle_idx, limb_side} = q;
                Time_Data{ID_Subj, ID_Speed, cycle_idx, limb_side}   = time;
                
                GRF_Data_L_N{ID_Subj, ID_Speed, cycle_idx, limb_side} = GRF_l_N(TO(cycle_idx):HS(cycle_idx),1);
                GRF_Data_R_N{ID_Subj, ID_Speed, cycle_idx, limb_side} = GRF_r_N(TO(cycle_idx):HS(cycle_idx),1);
%                 time = (walk_data.data(:,1)-1)*.01;
%                 time = (0:length(q)-1)'*.01;
                

                q1_limb{ID_Eval,:} = q(1,:);
                q2_limb{ID_Eval,:} = q(2,:);
                q4_limb{ID_Eval,:} = q(4,:);
                q5_limb{ID_Eval,:} = q(5,:);
                
                
                q1_limb_Subj_Speed{ID_Subj,cycle_stack_id} = q(1,:);
                q2_limb_Subj_Speed{ID_Subj,cycle_stack_id} = q(2,:);                
                q4_limb_Subj_Speed{ID_Subj,cycle_stack_id} = q(4,:);
                q5_limb_Subj_Speed{ID_Subj,cycle_stack_id} = q(5,:);
                

                Time = time;
                
                t = time;
                
                
                t_resamp = linspace(0,t(end),100);
%                 Angles_Data_resamp {ID_Subj, ID_Speed, cycle_idx, limb_side} = interp1(t',q',t_resamp');
                
                
                 
                 
                 
                 
                 
                t_resamp = linspace(0,t(end),50);
                
                % % % % l_hip   = lowpass(l_hip,6,100);
                % % % % r_hip   = lowpass(r_hip,6,100);
                % % % % l_knee  = lowpass(l_knee,6,100);
                % % % % r_knee  = lowpass(r_knee,6,100);
                % % % %
                % % % % l_knee(l_knee<0)=0;
                % % % % r_knee(r_knee<0)=0;
                
                
                q_resamp = interp1(t,q',t_resamp');
                
%                 Q1 (cycle_stack_id,:) = [Q1; q_resamp(:,1)'];
%                 Q2 (cycle_stack_id,:) = [Q2; q_resamp(:,2)'];
%                 Q4 (cycle_stack_id,:) = [Q4; q_resamp(:,4)'];
%                 Q5 (cycle_stack_id,:) = [Q5; q_resamp(:,5)'];
                
                Q1 (ID_Eval,:) = q_resamp(:,1)';
                Q2 (ID_Eval,:) = q_resamp(:,2)';
                Q4 (ID_Eval,:) = q_resamp(:,4)';
                Q5 (ID_Eval,:) = q_resamp(:,5)';
                
                
                
                
                
                
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
    
% % %     figure(102)
% % % %     clf;
% % %     subplot(2,2,1)
% % %     plot(time,q(2,:),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
% % %     hold on
% % %     subplot(2,2,3)
% % %     plot(time,q(1,:),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
% % %     hold on
% % %     subplot(2,2,2)
% % %     plot(time,q(4,:),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
% % %     hold on
% % %     subplot(2,2,4)
% % %     plot(time,q(5,:),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
% % %     hold on
% % %     
% % %      
% % %    if ID_Speed ~= [3 , 7]
% % %     figure(106)
% % % %     clf;
% % %     subplot(2,2,1)
% % %     plot(time,q(2,:),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
% % %     hold on
% % %     subplot(2,2,3)
% % %     plot(time,q(1,:),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
% % %     hold on
% % %     subplot(2,2,2)
% % %     plot(time,q(4,:),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
% % %     hold on
% % %     subplot(2,2,4)
% % %     plot(time,q(5,:),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
% % %     hold on
% % %     
% % %    end
    
   
% % %    if ID_Speed == 3 || ID_Speed==7
% % %     figure(107)
% % % %     clf;
% % %     subplot(2,2,1)
% % %     plot(time,q(2,:),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
% % %     hold on
% % %     subplot(2,2,3)
% % %     plot(time,q(1,:),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
% % %     hold on
% % %     subplot(2,2,2)
% % %     plot(time,q(4,:),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
% % %     hold on
% % %     subplot(2,2,4)
% % %     plot(time,q(5,:),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
% % %     hold on
% % %     
% % %    end
   
   
   
% % %    if ID_Speed==7
% % %     figure(108)
% % % %     clf;
% % %     subplot(2,2,1)
% % %     plot(time,q(2,:),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
% % %     hold on
% % %     subplot(2,2,3)
% % %     plot(time,q(1,:),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
% % %     hold on
% % %     subplot(2,2,2)
% % %     plot(time,q(4,:),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
% % %     hold on
% % %     subplot(2,2,4)
% % %     plot(time,q(5,:),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
% % %     hold on
% % %     
% % %    end
   

%% All Subjects Trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% % % % % % % %     figure(1033)
% % % % % % % %     subplot(6,7,ID_Subj)
% % % % % % % %     plot(time,q(1,:),'LineWidth',.6, 'Color', [1-ID_Speed/15 1-ID_Speed/15 1-ID_Speed/15])
% % % % % % % %     hold on
% % % % % % % %     title(['Subject  ' num2str(ID_Subj)])
% % % % % % % %     
% % % % % % % %     figure(1044)
% % % % % % % %     subplot(6,7,ID_Subj)
% % % % % % % %     plot(time,q(2,:),'LineWidth',.6, 'Color', [1-ID_Speed/15 1-ID_Speed/15 1-ID_Speed/15])
% % % % % % % %     hold on
% % % % % % % %     title(['Subject  ' num2str(ID_Subj)])
% % % % % % % % 
% % % % % % % %     figure(103)
% % % % % % % %     subplot(6,7,ID_Subj)
% % % % % % % %     plot(time,q(4,:),'LineWidth',.6, 'Color', [1-ID_Speed/15 1-ID_Speed/15 1-ID_Speed/15])
% % % % % % % %     hold on
% % % % % % % %     title(['Subject  ' num2str(ID_Subj)])
% % % % % % % %     
% % % % % % % %     figure(104)
% % % % % % % %     subplot(6,7,ID_Subj)
% % % % % % % %     plot(time,q(5,:),'LineWidth',.6, 'Color', [1-ID_Speed/15 1-ID_Speed/15 1-ID_Speed/15])
% % % % % % % %     hold on
% % % % % % % %     title(['Subject  ' num2str(ID_Subj)])
    
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
                
                % joint kinematics
                
                % % % % [q_c , dq_c , ddq_c dddq_c] = Cubic_Bspline(t , q(4,:)');
                % % % %
                % % % % dq = diff(q')'/(Time(end)/(length(Time)-1));
                % % % % dq = [dq(:,1) dq];
                % % % % ddq = diff(dq')'/(Time(end)/(length(Time)-1));
                % % % % ddq = [ddq(:,1) ddq];
                % % % % dddq = diff(ddq')'/(Time(end)/(length(Time)-1));
                % % % % dddq = [dddq(:,1) dddq];
                % % % %
                % % % %
                % % % % q_f = lowpass(q(4,:)',15,100)';
                % % % % dq_f = diff(q_f')'/(Time(end)/(length(Time)-1));
                % % % % dq_f = [dq_f(:,1) dq_f];
                % % % % ddq_f = diff(dq_f')'/(Time(end)/(length(Time)-1));
                % % % % ddq_f = [ddq_f(:,1) ddq_f];
                % % % % dddq_f = diff(ddq_f')'/(Time(end)/(length(Time)-1));
                % % % % dddq_f = [dddq_f(:,1) dddq_f];
                % % % %
                % % % % figure(202)
                % % % % subplot(1,3,1)
                % % % % plot([dq(4,:)' dq_f' dq_c])
                % % % % subplot(1,3,2)
                % % % % plot([ddq(4,:)' ddq_f' ddq_c])
                % % % % subplot(1,3,3)
                % % % % plot([dddq(4,:)' dddq_f' dddq_c])
                % % % % legend('Direct Derivative' , 'Filtered Derivative', 'Cubic Bspline', 'Filtered Cubic Bspline')
                
                
                
                % % [q , dq , ddq dddq] = Cubic_Bspline(t , q');
                % % q    = q';
                % % dq   = dq';
                % % ddq  = ddq';
                % % dddq = dddq';
                
                
                % filter
                % % % % % q1 = lowpass(q(1,:)',15,100)';
                % % % % % q2 = lowpass(q(2,:)',15,100)';
                % % % % % q4 = lowpass(q(4,:)',15,100)';
                % % % % % q5 = lowpass(q(5,:)',15,100)';
                % % % % % q_f = [q1 ; q2 ; zeros(size(q,2),1)' ; q4 ; q5];
                
                % % figure(14)
                % % plot([q(4,:)' q_f(4,:)'])
                % % legend( 'Thigh', 'Filtered Thigh')
                
                % % % % % q = q_f;
                
                t_t = linspace(0,t(end),1000);
                
                
                dq = diff(q')'/(Time(end)/(length(Time)-1));
                dq = [dq(:,1) dq];
                ddq = diff(dq')'/(Time(end)/(length(Time)-1));
                ddq = [ddq(:,1) ddq];
                dddq = diff(ddq')'/(Time(end)/(length(Time)-1));
                dddq = [dddq(:,1) dddq];
                
                
                t_int = linspace(0, t(end), 99);
                q = interp1(t', q', t_int')';
                dq = interp1(t', dq', t_int')';
                ddq= interp1(t', ddq', t_int')';
                t = t_int;
                
                x = [q ; dq];
                usol = dynamics_forward(t,x,ddq,p);
                
                

                x_t = pwPoly3(t,[q;dq],[dq; ddq],t_t);
                
                q  = x_t(1:5,:);
                dq = x_t(6:10,:);
                ddq = diff(dq')'/(Time(end)/(length(Time)-1));
                ddq = [ddq(:,1) ddq];
                dddq = diff(ddq')'/(Time(end)/(length(Time)-1));
                dddq = [dddq(:,1) dddq];
                
                u = pwPoly2(t,usol,t_t);
                
                % task kinematics
%                 [P, G] = getPoints(q,p);
                [P, G] = getPoints(t_t,q,p,V_tr);
                
                t = t_t;
                
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
                
                
                

                
                %%%%
                
                
                
                
                
                
                [u_c , du_c , ddu_c , dddu_c] = Cubic_Bspline(t , u');
                
                
                % for n_traj = 1:size(u,1)
                % %     u(n_traj,:) = lowpass(u(n_traj,:)',.01,100)';
                %     u_f(n_traj,:) = lowpass(u(n_traj,:)',10,100)';
                %
                % end
                
                % % figure(13)
                % % plot([u(4,:)' u_f(4,:)'])
                % % legend( 'Torque', 'Filtered Torque')
                
                % u = u_f;
                
                du = diff(u')'/(t(end)/(length(t)-1));
                du = [du(:,1) du];
                ddu= diff(du')'/(t(end)/(length(t)-1));
                ddu= [ddu(:,1) ddu];
                
                
                % % figure(202)
                % % subplot(1,3,1)
                % % plot([du(4,:)' du_c(:,4)])
                % % subplot(1,3,2)
                % % plot([ddu(4,:)' ddu_c(:,4)])
                % % legend( 'Filtered Derivative', 'Cubic Bspline')
                
                
                % torques from the Fakuchi data
                % % % u = [L_Knee_Torque(TO:HS,1)' ; L_Hip_Torque(TO:HS,1)' ; zeros(1,length(t)) ; R_Hip_Torque(TO:HS,1)' ; R_Knee_Torque(TO:HS,1)'];
                % % % du = diff(u')'/(t(end)/(length(t)-1));
                % % % du = [du(:,1) du];
                % % % ddu= diff(du')'/(t(end)/(length(t)-1));
                % % % ddu= [ddu(:,1) ddu];
                
                
                [Fx, Fy] = contactForces(q,dq,ddq,p);
                dFx = diff(Fx')'/(t(end)/(length(t)-1));
                dFx = [dFx(:,1) dFx];
                dFy = diff(Fy')'/(t(end)/(length(t)-1));
                dFy = [dFy(:,1) dFy];
                
                %     Fx_f = lowpass(Fx',10,100)';
                %     Fy_f = lowpass(Fy',10,100)';
                % dFx_f = diff(Fx_f')'/(t(end)/(length(t)-1));
                % dFx_f = [dFx_f(:,1) dFx_f];
                
                
                % figure(202)
                % subplot(1,3,1)
                % plot([Fx' Fx_f'])
                % subplot(1,3,2)
                % plot([dFx' dFx_f'])
                % legend( 'GRFx', 'GRFx filtered')
                
                % % % if limb_side == 1
                % % %     Fx = GRF_r_f(TO:HS)';
                % % %     Fy = GRF_r_N(TO:HS)';
                % % % else
                % % %     Fx = GRF_l_f(TO:HS)';
                % % %     Fy = GRF_l_N(TO:HS)';
                % % % end
                % % % dFx = diff(Fx')'/(t(end)/(length(t)-1));
                % % % dFx = [dFx(:,1) dFx];
                % % % dFy = diff(Fy')'/(t(end)/(length(t)-1));
                % % % dFy = [dFy(:,1) dFy];
                
                
                
                
                % cycl_err(i,1) = norm(q(:,r_TO(i))-q(:,r_HS(i)));
                
                
                
                
                
                
                
                
                
                
                
                % .^.5
                
                Accel_Joint = (sum(ddq.^2));
                Jerk_Joint = (sum(dddq.^2));
                
                Accel_Task  = (sum(ddG.^2));
                Jerk_Task  = (sum(dddG.^2));
                
                % Time = 1;
                
                % Velocity_Deviation = dP(2,:)
                
                
                Torque_Squared           = (sum(u.^2));
                Torque_Absolute          = sum(abs(u));
                TorqueRate_Squared       = (sum(du.^2));
                TorqueRate_Absolute      = sum(abs(du));
                TorqueRateChange_Squared = (sum(ddu.^2));
                TorqueRateChange_Absolute= sum(abs(ddu));
                
                yGRF_Rate    =  (dFy.^2);
                xGRF_Rate    =  (dFx.^2);
                
                
                Work_Absolute = sum(abs(u.*dq));
                Work_Positive = sum(max(u.*dq,0));
                
                Kinetic_Energy = 1/2*(p.I1+p.m1*(p.l1-p.c1)^2)*dq(1,:).^2 + ...
                    1/2*(p.I2+p.m2*(p.l2-p.c2)^2)*dq(2,:).^2 + ...
                    1/2*(p.I3+p.m3*(p.l3-p.c3)^2)*dq(3,:).^2 + ... ;   % Kinetic Energy
                    1/2*(p.I4+p.m4*(p.c4)^2)*dq(4,:).^2 + ...
                    1/2*(p.I5+p.m5*(p.c5)^2)*dq(5,:).^2 ;
                
                Angular_Momentum = (p.I1+p.m1*(p.l1-p.c1)^2)*dq(1,:) + ...
                    (p.I2+p.m2*(p.l2-p.c2)^2)*dq(2,:) + ...
                    (p.I3+p.m3*(p.l3-p.c3)^2)*dq(3,:) + ... ;   % Kinetic Energy
                    (p.I4+p.m4*(p.c4)^2)*dq(4,:) + ...
                    (p.I5+p.m5*(p.c5)^2)*dq(5,:) ;
                
                
                % COST COMPONENTS EVALS
                
                % kinematic
                cost_components_traj_eval(:,1) = Accel_Joint';
                cost_components_traj_eval(:,2) = Jerk_Joint';
                cost_components_traj_eval(:,3) = Accel_Task';
                cost_components_traj_eval(:,4) = Jerk_Task';
                % cost_components_traj_eval(:,5) = 1;
                
                % dynamic
                cost_components_traj_eval(:,5) = Torque_Squared';
                cost_components_traj_eval(:,6) = Torque_Absolute';
                cost_components_traj_eval(:,7) = TorqueRate_Squared';
                cost_components_traj_eval(:,8) = TorqueRate_Absolute';
                cost_components_traj_eval(:,9) = TorqueRateChange_Squared';
                cost_components_traj_eval(:,10) = TorqueRateChange_Absolute';
                cost_components_traj_eval(:,11) = yGRF_Rate';
                cost_components_traj_eval(:,12) = xGRF_Rate';
                
                % energetic
                cost_components_traj_eval(:,13) = Work_Positive';
                cost_components_traj_eval(:,14) = Work_Absolute';
                cost_components_traj_eval(:,15) = max(Kinetic_Energy')/t(end);
                
                % balance / stability
                cost_components_traj_eval(:,16) = max(Angular_Momentum)-min(Angular_Momentum)/t(end);
                
                
                
                Cost_Components_Eval_Fakuchi(ID_Eval,:) = trapz(t, cost_components_traj_eval) / t(end);
                
                % cost_bases_eval(i,:) = simps(time,C');
                %
                % % cost_vector = C./scaling_vector;
                % clear C
                clear cost_components_traj_eval
            end
          
        end
        end
        
%     Mean_Q1(ID_Subj,:) = mean(Q1);
%     Mean_Q2(ID_Subj,:) = mean(Q2);
%     Mean_Q4(ID_Subj,:) = mean(Q4);
%     Mean_Q5(ID_Subj,:) = mean(Q5);
    
    for idx_y = ID_Eval-speed_stack_id+1:ID_Eval
        Dev_Q1(idx_y,1) = dtw(q1_limb{idx_y,:},mean(Q1(ID_Eval-speed_stack_id+1:ID_Eval,:)));
        Dev_Q2(idx_y,1) = dtw(q2_limb{idx_y,:},mean(Q2(ID_Eval-speed_stack_id+1:ID_Eval,:)));
        Dev_Q4(idx_y,1) = dtw(q4_limb{idx_y,:},mean(Q4(ID_Eval-speed_stack_id+1:ID_Eval,:)));
        Dev_Q5(idx_y,1) = dtw(q5_limb{idx_y,:},mean(Q5(ID_Eval-speed_stack_id+1:ID_Eval,:)));
        
        
%         Q1_mean_resamp = interp1(1:100,mean(Q1),linspace(1,100,length(q1_limb{idx_y,:})));
%         Q2_mean_resamp = interp1(1:100,mean(Q2),linspace(1,100,length(q2_limb{idx_y,:})));
%         Q4_mean_resamp = interp1(1:100,mean(Q4),linspace(1,100,length(q4_limb{idx_y,:})));
%         Q5_mean_resamp = interp1(1:100,mean(Q5),linspace(1,100,length(q5_limb{idx_y,:})));
%         
%         Dev_Q1(idx_y,1) = rms(q1_limb{idx_y,:}-Q1_mean_resamp);
%         Dev_Q2(idx_y,1) = rms(q2_limb{idx_y,:}-Q2_mean_resamp);
%         Dev_Q4(idx_y,1) = rms(q4_limb{idx_y,:}-Q4_mean_resamp);
%         Dev_Q5(idx_y,1) = rms(q5_limb{idx_y,:}-Q5_mean_resamp);
        
        DTW_Q1(idx_y,1) = Dev_Q1(idx_y,1)/10;
        DTW_Q2(idx_y,1) = Dev_Q2(idx_y,1)/10;
        DTW_Q4(idx_y,1) = Dev_Q4(idx_y,1)/10;
        DTW_Q5(idx_y,1) = Dev_Q5(idx_y,1)/10;

        DTW_Sum(idx_y,1)            =  (Dev_Q1(idx_y,1)+...
                                        Dev_Q2(idx_y,1)+...
                                        Dev_Q4(idx_y,1)+...
                                        Dev_Q5(idx_y,1))/10;
        
                                    
        RMSE_Q1(idx_y,1)             =   rms(Q1(idx_y,:) - mean(Q1(ID_Eval-speed_stack_id+1:ID_Eval,:)));
        RMSE_Q2(idx_y,1)             =   rms(Q2(idx_y,:) - mean(Q2(ID_Eval-speed_stack_id+1:ID_Eval,:)));
        RMSE_Q4(idx_y,1)             =   rms(Q4(idx_y,:) - mean(Q4(ID_Eval-speed_stack_id+1:ID_Eval,:)));
        RMSE_Q5(idx_y,1)             =   rms(Q5(idx_y,:) - mean(Q5(ID_Eval-speed_stack_id+1:ID_Eval,:)));                            
                                    
        RMSE_Sum(idx_y,1)                =   rms(Q1(idx_y,:) - mean(Q1(ID_Eval-speed_stack_id+1:ID_Eval,:)))+...
                                        rms(Q2(idx_y,:) - mean(Q2(ID_Eval-speed_stack_id+1:ID_Eval,:)))+...
                                        rms(Q4(idx_y,:) - mean(Q4(ID_Eval-speed_stack_id+1:ID_Eval,:)))+...
                                        rms(Q5(idx_y,:) - mean(Q5(ID_Eval-speed_stack_id+1:ID_Eval,:)));
                                   
%         Deviation_Eval_mean(idx_y,1) =  Dev_Q1(idx_y,1);
%                                     
%         RMSE(idx_y,1)                =  rms(Q1(idx_y,:) - mean(Q1(ID_Eval-speed_stack_id+1:ID_Eval,:)));                            
                                
        
    end
    
    
    % plot data points on cost features and DTW
%     figure(203)
% %     clf;
%     scatter(DTW_Sum(ID_Eval-speed_stack_id+1:ID_Eval,1)',...
%         ID_Speed/8*ones(1,speed_stack_id)',...
%         'filled','MarkerEdgeColor',[.9-ID_Speed/10 .9-ID_Speed/10 .9-ID_Speed/10],...
%         'MarkerFaceColor',[.9-ID_Speed/10 .9-ID_Speed/10 .9-ID_Speed/10])
%     hold on
    
    
%     
%     figure(202)
% %     clf;
%     scatter3(Cost_Components_Eval_Fakuchi(ID_Eval-speed_stack_id+1:ID_Eval,2)',...
%         Cost_Components_Eval_Fakuchi(ID_Eval-speed_stack_id+1:ID_Eval,5)',...
%         Cost_Components_Eval_Fakuchi(ID_Eval-speed_stack_id+1:ID_Eval,13)',...
%         'filled','MarkerEdgeColor',[1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10],...
%         'MarkerFaceColor',[1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
%     hold on
    
    
    
    
    
    
%             figure(104)
% %     clf;
%     subplot(2,2,1)
%     plot(t_resamp,mean(Q2(ID_Eval-speed_stack_id+1:ID_Eval,:)),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
%     hold on
%     subplot(2,2,3)
%     plot(t_resamp,mean(Q1(ID_Eval-speed_stack_id+1:ID_Eval,:)),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
%     hold on
%     subplot(2,2,2)
%     plot(t_resamp,mean(Q4(ID_Eval-speed_stack_id+1:ID_Eval,:)),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
%     hold on
%     subplot(2,2,4)
%     plot(t_resamp,mean(Q5(ID_Eval-speed_stack_id+1:ID_Eval,:)),'LineWidth',.6, 'Color', [1-ID_Speed/10 1-ID_Speed/10 1-ID_Speed/10])
%     hold on
    
    
                        Cost_Components_Subjects_Speeds{ID_Subj_actual,ID_Speed} = Cost_Components_Eval_Fakuchi(ID_Eval-speed_stack_id+1:ID_Eval,:);
                    DTW_Sum_Subjects_Speeds{ID_Subj_actual,ID_Speed} = DTW_Sum(ID_Eval-speed_stack_id+1:ID_Eval,:);
                    RMSE_Sum_Subjects_Speeds{ID_Subj_actual,ID_Speed} = RMSE_Sum(ID_Eval-speed_stack_id+1:ID_Eval,:);
    
    % add the mean trajectories on top of all traj for subject 1
%     figure(102)
%     subplot(2,2,1)
%     plot(t_resamp,mean(Q2(ID_Eval-speed_stack_id+1:ID_Eval,:)),'LineStyle','-',...
%         'LineWidth',4, 'Color', [.8-ID_Speed/10 .8-ID_Speed/10 .8-ID_Speed/10])
%     
%     subplot(2,2,3)
%     plot(t_resamp,mean(Q1(ID_Eval-speed_stack_id+1:ID_Eval,:)),'LineStyle','-',...
%         'LineWidth',4, 'Color', [.8-ID_Speed/10 .8-ID_Speed/10 .8-ID_Speed/10])
%     
%     subplot(2,2,2)
%     plot(t_resamp,mean(Q4(ID_Eval-speed_stack_id+1:ID_Eval,:)),'LineStyle','-',...
%         'LineWidth',4, 'Color', [.8-ID_Speed/10 .8-ID_Speed/10 .8-ID_Speed/10])
%     
%     subplot(2,2,4)
%         plot(t_resamp,mean(Q5(ID_Eval-speed_stack_id+1:ID_Eval,:)),'LineStyle','-',...
%             'LineWidth',4, 'Color', [.8-ID_Speed/10 .8-ID_Speed/10 .8-ID_Speed/10])
        
        
        
% % %     if ID_Speed==7
% % %         t_resamp = linspace(0,.38,100);
% % %     figure(108)
% % %     subplot(2,2,1)
% % % % %     shadedErrorBar(t_resamp, mean(Q2(ID_Eval-speed_stack_id+1:ID_Eval,:))+ .15*randn(length(t_resamp),100), {@mean,@std}, 'lineprops', 'm-')
% % % %     plot(t_resamp,mean(Q2(ID_Eval-speed_stack_id+1:ID_Eval,:)),'LineStyle','-',...
% % % %         'LineWidth',4, 'Color', [.8-ID_Speed/10 .8-ID_Speed/10 .8-ID_Speed/10])
% % %         shadedErrorBar(t_resamp, mean(Q2(ID_Eval-speed_stack_id+1:ID_Eval,:))+ .05*randn(length(t_resamp),100), {@mean,@std}, 'lineprops', 'c-')
% % % 
% % %     subplot(2,2,3)
% % % %     plot(t_resamp,mean(Q1(ID_Eval-speed_stack_id+1:ID_Eval,:)),'LineStyle','-',...
% % % %         'LineWidth',4, 'Color', [.8-ID_Speed/10 .8-ID_Speed/10 .8-ID_Speed/10])
% % % % %         shadedErrorBar(t_resamp, mean(Q1(ID_Eval-speed_stack_id+1:ID_Eval,:))+ .15*randn(length(t_resamp),100), {@mean,@std}, 'lineprops', 'm-')
% % %         shadedErrorBar(t_resamp, mean(Q1(ID_Eval-speed_stack_id+1:ID_Eval,:))+ .05*randn(length(t_resamp),100), {@mean,@std}, 'lineprops', 'c-')
% % % 
% % %     
% % %     subplot(2,2,2)
% % % %     plot(t_resamp,mean(Q4(ID_Eval-speed_stack_id+1:ID_Eval,:)),'LineStyle','-',...
% % % %         'LineWidth',4, 'Color', [.8-ID_Speed/10 .8-ID_Speed/10 .8-ID_Speed/10])
% % % % %         shadedErrorBar(t_resamp, mean(Q4(ID_Eval-speed_stack_id+1:ID_Eval,:))+ .15*randn(length(t_resamp),100), {@mean,@std}, 'lineprops', 'm-')
% % %         shadedErrorBar(t_resamp, mean(Q4(ID_Eval-speed_stack_id+1:ID_Eval,:))+ .03*randn(length(t_resamp),100), {@mean,@std}, 'lineprops', 'c-')
% % % 
% % %     
% % %     subplot(2,2,4)
% % % %         plot(t_resamp,mean(Q5(ID_Eval-speed_stack_id+1:ID_Eval,:)),'LineStyle','-',...
% % % %             'LineWidth',4, 'Color', [.8-ID_Speed/10 .8-ID_Speed/10 .8-ID_Speed/10])
% % % % %     shadedErrorBar(t_resamp, mean(Q5(ID_Eval-speed_stack_id+1:ID_Eval,:))+ .17*randn(length(t_resamp),100), {@mean,@std}, 'lineprops', 'm-')
% % %     shadedErrorBar(t_resamp, mean(Q5(ID_Eval-speed_stack_id+1:ID_Eval,:))+ .07*randn(length(t_resamp),100), {@mean,@std}, 'lineprops', 'c-')
% % % 
% % %     end
        
        
    end
    





    ID_Subj;
    
    clear HS TO
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    figure(102)
%     colormap gray;
% % c = colorbar( 'Location', 'EastOutside', ...
% %     'Ticks', Z_Bar, 'TickLabels', cellstr( num2str( Z_Bar(:)*Z_Scl, '%.3e' ) ) );
% c = colorbar('direction','reverse', 'Location', 'EastOutside', ...
%     'Ticks', [0, .5, 1], 'TickLabels', {'Fast' , 'Medium', 'Slow'} );
% ylabel( c, 'Walking Speed' );
%     subplot(2,2,1)
%     ylabel('Femur Angle (rad)')
% %     title('Stance Limb')
%     title('\underline{\textbf{Stance Limb}}', 'FontSize', 20, 'Interpreter', 'latex')
%     xlim([0,.7])
%     ylim([-.35, .5])
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 15);
%     subplot(2,2,3)
%     xlabel('Time (s)')
%     ylabel('Tibia Angle (rad)')
%     xlim([0,.7])
%     ylim([-1.1, .5])
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 15);
%     subplot(2,2,2)
%     ylabel('Femur Angle (rad)')
% %     title('Swing Limb')
%     title('\underline{\textbf{Swing Limb}}', 'FontSize', 20, 'Interpreter', 'latex')
%     xlim([0,.7])
%     ylim([-.35, .5])
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 15);
%     subplot(2,2,4)
%     xlabel('Time (s)')
%     ylabel('Tibia Angle (rad)')
%     xlim([0,.7])
%     ylim([-1.1, .5])
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 15);
% %     sgtitle('Mean Limb Angles per Walking Speed : Subject 1', 'FontSize', 20)
%     sgtitle('Segment Angle Trajectories for Subject 1', 'FontSize', 20)
%     set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[15 1 15 12]);
%     colormap gray;
% % c = colorbar( 'Location', 'EastOutside', ...
% %     'Ticks', Z_Bar, 'TickLabels', cellstr( num2str( Z_Bar(:)*Z_Scl, '%.3e' ) ) );
% c = colorbar('direction','reverse', 'Location', 'EastOutside', ...
%     'Ticks', [0, .5, 1], 'TickLabels', {'Fast' , 'Medium', 'Slow'} );
% ylabel( c, 'Walking Speed' );




% figure(106)
%     subplot(2,2,1)
%     ylabel('Femur Angle (rad)')
% %     title('Stance Limb')
%     title('\underline{\textbf{Stance Limb}}', 'FontSize', 20, 'Interpreter', 'latex')
%     xlim([0,.7])
%     ylim([-.35, .5])
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 15);
%     subplot(2,2,3)
%     xlabel('Time (s)')
%     ylabel('Tibia Angle (rad)')
%     xlim([0,.7])
%     ylim([-1.1, .5])
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 15);
%     subplot(2,2,2)
%     ylabel('Femur Angle (rad)')
% %     title('Swing Limb')
%     title('\underline{\textbf{Swing Limb}}', 'FontSize', 20, 'Interpreter', 'latex')
%     xlim([0,.7])
%     ylim([-.35, .5])
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 15);
%     subplot(2,2,4)
%     xlabel('Time (s)')
%     ylabel('Tibia Angle (rad)')
%     xlim([0,.7])
%     ylim([-1.1, .5])
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 15);
% %     sgtitle('Mean Limb Angles per Walking Speed : Subject 1', 'FontSize', 20)
%     sgtitle('Segment Angle Trajectories for Subject 1', 'FontSize', 20)
%     set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[15 1 15 12]);
    
    
    
    
    
    
%     figure(107)
% 
%     subplot(2,2,1)
%     ylabel('Femur Angle (rad)')
% %     title('Stance Limb')
%     title('\underline{\textbf{Stance Limb}}', 'FontSize', 20, 'Interpreter', 'latex')
%     xlim([0,.7])
%     ylim([-.35, .5])
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 15);
%     subplot(2,2,3)
%     xlabel('Time (s)')
%     ylabel('Tibia Angle (rad)')
%     xlim([0,.7])
%     ylim([-1.1, .5])
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 15);
%     subplot(2,2,2)
%     ylabel('Femur Angle (rad)')
% %     title('Swing Limb')
%     title('\underline{\textbf{Swing Limb}}', 'FontSize', 20, 'Interpreter', 'latex')
%     xlim([0,.7])
%     ylim([-.35, .5])
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 15);
%     subplot(2,2,4)
%     xlabel('Time (s)')
%     ylabel('Tibia Angle (rad)')
%     xlim([0,.7])
%     ylim([-1.1, .5])
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 15);
% %     sgtitle('Mean Limb Angles per Walking Speed : Subject 1', 'FontSize', 20)
%     sgtitle('Segment Angle Trajectories for Subject 1', 'FontSize', 20)
%     set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[15 1 15 12]);



    
    
    
%     figure(108)
% 
%     subplot(2,2,1)
%     ylabel('Femur Angle (rad)')
% %     title('Stance Limb')
%     title('\underline{\textbf{Stance Limb}}', 'FontSize', 20, 'Interpreter', 'latex')
%     xlim([0,.7])
%     ylim([-.35, .5])
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 15);
%     subplot(2,2,3)
%     xlabel('Time (s)')
%     ylabel('Tibia Angle (rad)')
%     xlim([0,.7])
%     ylim([-1.1, .5])
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 15);
%     subplot(2,2,2)
%     ylabel('Femur Angle (rad)')
% %     title('Swing Limb')
%     title('\underline{\textbf{Swing Limb}}', 'FontSize', 20, 'Interpreter', 'latex')
%     xlim([0,.7])
%     ylim([-.35, .5])
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 15);
%     subplot(2,2,4)
%     xlabel('Time (s)')
%     ylabel('Tibia Angle (rad)')
%     xlim([0,.7])
%     ylim([-1.1, .5])
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 15);
% %     sgtitle('Mean Limb Angles per Walking Speed : Subject 1', 'FontSize', 20)
%     sgtitle('Segment Angle Trajectories for Subject 1', 'FontSize', 20)
%     set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[15 1 15 12]);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
    
    t_res1 = 0:.01:1;
t_res2 = 0:.01:2;
sin_wave = sin(2*pi*t_res1);
sin_wave_scaled = sin(pi*t_res2);
sin_wave_resamp = interp1(t_res2, sin_wave_scaled, t_res1);

% figure(10)
% plot(sin_wave, 'LineWidth', 1.4)
% hold on
% plot(sin_wave_scaled, 'LineWidth', 1.4)
% sgtitle('Sin wave with 2 different periodicity values', 'FontSize', 17)
% ylabel('Angle (rad)')
% xlabel('Time point')
% legend('sin(2\pi/100 t)','sin(\pi/100 t)')
% set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'out'     , ...
%         'TickLength'  , [.02 .02] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 14);
%     set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[15 1 15 12]);
% 
%     figure(9)
% dtw(sin_wave,sin_wave_scaled)
% ylabel('Angle (rad)')
% xlabel('Time point')
% legend('sin(2\pi/100 t)','sin(\pi/100 t)')
% set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'out'     , ...
%         'TickLength'  , [.02 .02] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 14);
%     set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[15 1 15 12]);
% 
% 
% 
% rms(sin_wave'-sin_wave_resamp')


% sin_wave_resamp = interp1(t_res2, sin_wave_scaled, t_res1*2);
% figure(9)
% plot(sin_wave, 'LineWidth', 1.4)
% hold on
% plot(sin_wave_resamp, 'LineWidth', 1.4)
% rms(sin_wave'-sin_wave_resamp')
% ylabel('Angle (rad)')
% xlabel('Time point')
% title('Aligned Signals for RMS (RMS=0)')
% % legend('sin(2\pi/100 t)','sin(\pi/100 t)')
% set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'out'     , ...
%         'TickLength'  , [.02 .02] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 14);
%     set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[15 1 15 12]);
    
    


    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     figure(101)
%     subplot(2,2,1)
%     plot([Q2]','LineWidth',.6,'Color',[.6 .6 .6])
%     hold on
%     plot(mean(Q2),'k','LineWidth',3)
% %     xlabel('Gait Cycle %')
%     ylabel('Stance Femur (rad)')
%     title('Stance Limb')
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'out'     , ...
%         'TickLength'  , [.02 .02] , ...
%         'YGrid'       , 'on'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 14);
%     subplot(2,2,3)
%     plot([Q1]','LineWidth',.6,'Color',[.6 .6 .6])
%     hold on
%     plot(mean(Q1),'k','LineWidth',3)
%     xlabel('Gait Cycle %')
%     ylabel('Stance Tibia (rad)')
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'out'     , ...
%         'TickLength'  , [.02 .02] , ...
%         'YGrid'       , 'on'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 14);
%     subplot(2,2,2)
%     plot([Q4]','LineWidth',.6,'Color',[.6 .6 .6])
%     hold on
%     plot(mean(Q4),'k','LineWidth',3)
% %     xlabel('Gait Cycle %')
%     ylabel('Swing Femur (rad)')
%     title('Swing Limb')
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'out'     , ...
%         'TickLength'  , [.02 .02] , ...
%         'YGrid'       , 'on'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 14);
%     subplot(2,2,4)
%     plot([Q5]','LineWidth',.6,'Color',[.6 .6 .6])
%     hold on
%     plot(mean(Q5),'k','LineWidth',3)
%     xlabel('Gait Cycle %')
%     ylabel('Swing Tibia (rad)')
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'out'     , ...
%         'TickLength'  , [.02 .02] , ...
%         'YGrid'       , 'on'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 14);
%     
%     colormap gray;

    
%     figure(202)
%     set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[17 1 17 12]);
%     xlabel('Jerk')
%     ylabel('Torque Sq')
%     zlabel('Work^+')
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'out'     , ...
%         'TickLength'  , [.02 .02] , ...
%         'YGrid'       , 'on'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [.3 .3 .3], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 14);
    
%     figure(203)
%     set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[17 1 17 12]);
%     xlabel('DTW')
% %     xlabel('Speed Conditions')
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.01 .01] , ...
%         'YGrid'       , 'off'      , ...
%         'XColor'      , [.3 .3 .3], ...
%         'YColor'      , [1 1 1], ...
%         'LineWidth'   , 1.4        ,...
%         'FontSize'    , 14);
    
    
    
    % c = colorbar( 'Location', 'EastOutside', ...
%     'Ticks', Z_Bar, 'TickLabels', cellstr( num2str( Z_Bar(:)*Z_Scl, '%.3e' ) ) );
% c = colorbar('direction','reverse', 'Location', 'EastOutside', ...
%     'Ticks', [0, .5, 1], 'TickLabels', {'Fast' , 'Medium', 'Slow'} );
% ylabel( c, 'Walking Speed' );
% 
%     sgtitle('Segment Angles from toe-off to heel-strike : Subject 1', 'FontSize', 20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                    Cost_Components_Subjects{ID_Subj_actual,1} = Cost_Components_Eval_Fakuchi(ID_Eval-cycle_stack_id+1:ID_Eval,:);
                    DTW_Sum_Subjects{ID_Subj_actual,1} = DTW_Sum(ID_Eval-cycle_stack_id+1:ID_Eval,:);
                    RMSE_Sum_Subjects{ID_Subj_actual,1} = RMSE_Sum(ID_Eval-cycle_stack_id+1:ID_Eval,:);

    end
end
% % figure(1)
% % set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[17 1 17 12]);
% % % subplot(1,2,1)
% % hist(DTW_Sum, 30)
% % h = findobj(gca,'Type','patch');
% % h.FaceColor = [0 0.5 0.5];
% % h.EdgeColor = 'w';
% % hold on 
% % hist(RMSE_Sum, 30)
% % h = findobj(gca,'Type','patch');
% % h.FaceColor = [0 0.5 0.5];
% % h.EdgeColor = 'w';
% % h.FaceAlpha = .5;
% % xlabel('RMS Measure Value')
% % xlabel('DTW Measure Value')
% % ylabel('Count #')
% % title(['DTW and RMS measures distribution:  DTW_{max}=' num2str(max(DTW_Sum)) '  ;  RMS_{max}=' num2str(max(RMSE_Sum))'])
% % set(gca, ...
% %   'Box'         , 'off'     , ...
% %   'TickDir'     , 'out'     , ...
% %   'TickLength'  , [.01 .01] , ...
% %   'YGrid'       , 'off'      , ...
% %   'XColor'      , [.3 .3 .3], ...
% %   'YColor'      , [.3 .3 .3], ...
% %   'LineWidth'   , 1.3       ,...
% %   'FontSize'    , 12);
% % % subplot(1,2,2)
% % % hist(RMSE_Sum, 30)
% % % h = findobj(gca,'Type','patch');
% % % h.FaceColor = [0 0.5 0.5];
% % % h.EdgeColor = 'w';
% % % xlabel('RMS Measure Value')
% % % ylabel('Count #')
% % % title(['Sum of RMS measure distribution:  Maximum=' num2str(max(RMSE_Sum))])
% % % set(gca, ...
% % %   'Box'         , 'off'     , ...
% % %   'TickDir'     , 'out'     , ...
% % %   'TickLength'  , [.01 .01] , ...
% % %   'YGrid'       , 'off'      , ...
% % %   'XColor'      , [.3 .3 .3], ...
% % %   'YColor'      , [.3 .3 .3], ...
% % %   'LineWidth'   , 1.3       ,...
% % %   'FontSize'    , 12);
% % 
% % figure(2)
% % set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[17 1 17 12]);
% % subplot(4,2,1)
% % hist(DTW_Q1, 20)
% % h = findobj(gca,'Type','patch');
% % h.FaceColor = [0 0.5 0.5];
% % h.EdgeColor = 'w';
% % xlabel('DTW Measure Value')
% % ylabel('Count #')
% % title(['DTW Distribution (Stance Shank):  Maximum=' num2str(max(DTW_Q1))])
% % set(gca, ...
% %   'Box'         , 'off'     , ...
% %   'TickDir'     , 'out'     , ...
% %   'TickLength'  , [.01 .01] , ...
% %   'YGrid'       , 'off'      , ...
% %   'XColor'      , [.3 .3 .3], ...
% %   'YColor'      , [.3 .3 .3], ...
% %   'LineWidth'   , 1.3       ,...
% %   'FontSize'    , 12);
% % subplot(4,2,2)
% % hist(RMSE_Q1, 20)
% % h = findobj(gca,'Type','patch');
% % h.FaceColor = [0 0.5 0.5];
% % h.EdgeColor = 'w';
% % xlabel('RMS Measure Value')
% % ylabel('Count #')
% % title(['RMS Distribution (Stance Shank):  Maximum=' num2str(max(RMSE_Q1))])
% % set(gca, ...
% %   'Box'         , 'off'     , ...
% %   'TickDir'     , 'out'     , ...
% %   'TickLength'  , [.01 .01] , ...
% %   'YGrid'       , 'off'      , ...
% %   'XColor'      , [.3 .3 .3], ...
% %   'YColor'      , [.3 .3 .3], ...
% %   'LineWidth'   , 1.3       ,...
% %   'FontSize'    , 12);
% % subplot(4,2,3)
% % hist(DTW_Q2, 20)
% % h = findobj(gca,'Type','patch');
% % h.FaceColor = [0 0.5 0.5];
% % h.EdgeColor = 'w';
% % xlabel('DTW Measure Value')
% % ylabel('Count #')
% % title(['DTW Distribution (Stance Thigh):  Maximum=' num2str(max(DTW_Q2))])
% % set(gca, ...
% %   'Box'         , 'off'     , ...
% %   'TickDir'     , 'out'     , ...
% %   'TickLength'  , [.01 .01] , ...
% %   'YGrid'       , 'off'      , ...
% %   'XColor'      , [.3 .3 .3], ...
% %   'YColor'      , [.3 .3 .3], ...
% %   'LineWidth'   , 1.3       ,...
% %   'FontSize'    , 12);
% % subplot(4,2,4)
% % hist(RMSE_Q2, 20)
% % h = findobj(gca,'Type','patch');
% % h.FaceColor = [0 0.5 0.5];
% % h.EdgeColor = 'w';
% % xlabel('RMS Measure Value')
% % ylabel('Count #')
% % title(['RMS Distribution (Stance Thigh):  Maximum=' num2str(max(RMSE_Q2))])
% % set(gca, ...
% %   'Box'         , 'off'     , ...
% %   'TickDir'     , 'out'     , ...
% %   'TickLength'  , [.01 .01] , ...
% %   'YGrid'       , 'off'      , ...
% %   'XColor'      , [.3 .3 .3], ...
% %   'YColor'      , [.3 .3 .3], ...
% %   'LineWidth'   , 1.3       ,...
% %   'FontSize'    , 12);
% % subplot(4,2,5)
% % hist(DTW_Q4, 20)
% % h = findobj(gca,'Type','patch');
% % h.FaceColor = [0 0.5 0.5];
% % h.EdgeColor = 'w';
% % xlabel('DTW Measure Value')
% % ylabel('Count #')
% % title(['DTW Distribution (Swing Thigh):  Maximum=' num2str(max(DTW_Q4))])
% % set(gca, ...
% %   'Box'         , 'off'     , ...
% %   'TickDir'     , 'out'     , ...
% %   'TickLength'  , [.01 .01] , ...
% %   'YGrid'       , 'off'      , ...
% %   'XColor'      , [.3 .3 .3], ...
% %   'YColor'      , [.3 .3 .3], ...
% %   'LineWidth'   , 1.3       ,...
% %   'FontSize'    , 12);
% % subplot(4,2,6)
% % hist(RMSE_Q4, 20)
% % h = findobj(gca,'Type','patch');
% % h.FaceColor = [0 0.5 0.5];
% % h.EdgeColor = 'w';
% % xlabel('RMS Measure Value')
% % ylabel('Count #')
% % title(['RMS Distribution (Swing Thigh):  Maximum=' num2str(max(RMSE_Q4))])
% % set(gca, ...
% %   'Box'         , 'off'     , ...
% %   'TickDir'     , 'out'     , ...
% %   'TickLength'  , [.01 .01] , ...
% %   'YGrid'       , 'off'      , ...
% %   'XColor'      , [.3 .3 .3], ...
% %   'YColor'      , [.3 .3 .3], ...
% %   'LineWidth'   , 1.3       ,...
% %   'FontSize'    , 12);
% % subplot(4,2,7)
% % hist(DTW_Q5, 20)
% % h = findobj(gca,'Type','patch');
% % h.FaceColor = [0 0.5 0.5];
% % h.EdgeColor = 'w';
% % xlabel('DTW Measure Value')
% % ylabel('Count #')
% % title(['DTW Distribution (Swing Shank):  Maximum=' num2str(max(DTW_Q5))])
% % set(gca, ...
% %   'Box'         , 'off'     , ...
% %   'TickDir'     , 'out'     , ...
% %   'TickLength'  , [.01 .01] , ...
% %   'YGrid'       , 'off'      , ...
% %   'XColor'      , [.3 .3 .3], ...
% %   'YColor'      , [.3 .3 .3], ...
% %   'LineWidth'   , 1.3       ,...
% %   'FontSize'    , 12);
% % subplot(4,2,8)
% % hist(RMSE_Q5, 20)
% % h = findobj(gca,'Type','patch');
% % h.FaceColor = [0 0.5 0.5];
% % h.EdgeColor = 'w';
% % xlabel('RMS Measure Value')
% % ylabel('Count #')
% % title(['RMS Distribution (Swing Shank):  Maximum=' num2str(max(RMSE_Q5))])
% % set(gca, ...
% %   'Box'         , 'off'     , ...
% %   'TickDir'     , 'out'     , ...
% %   'TickLength'  , [.01 .01] , ...
% %   'YGrid'       , 'off'      , ...
% %   'XColor'      , [.3 .3 .3], ...
% %   'YColor'      , [.3 .3 .3], ...
% %   'LineWidth'   , 1.3       ,...
% %   'FontSize'    , 12);




% [lin_rel, lin_rel_gof] = fit(RMSE_Sum,DTW_Sum,'poly1');
% p1_lin = predint(lin_rel,RMSE_Sum,0.95,'functional', 'off');
% figure(1022)
% set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[17 1 17 12]);
% subplot(1,3,1)
% hold on
% plot(RMSE_Sum,lin_rel.p1*RMSE_Sum+lin_rel.p2,'r','LineWidth',2.5)
% % plot(linspace(.5,3.5,50),p1_cost,'m--')
% 
% fx = [RMSE_Sum', fliplr(RMSE_Sum)'];
% fy = [p1_lin(:,1)', fliplr(p1_lin(:,2)')];
% hold on
% fill(fx,fy,'r', 'FaceAlpha', .2 , 'edgecolor','none');
% 
% scatter(RMSE_Sum , DTW_Sum ,60, 'filled', 'k')
% 
% xlabel('RMS', 'FontSize', 12)
% ylabel('DTW', 'FontSize', 12)
% title('Positive Mechanical Work', 'FontSize', 14)
% text(1.8,110,['DTW_{fit}= ' num2str(lin_rel.p1,'%.2f') ' RMS  +  ' num2str(lin_rel.p2,'%.2f')],'Color','r', 'FontSize', 11)
% text(1.8,90,['r_{adj}= ' num2str(lin_rel_gof.adjrsquare,'%.2f' ) ' , p=' num2str(pvalue_cost(2,1),'%.2f')],'Color','r', 'FontSize', 11)
% 
% set(gca, ...
%   'Box'         , 'off'     , ...
%   'TickDir'     , 'out'     , ...
%   'TickLength'  , [.01 .01] , ...
%   'YGrid'       , 'on'      , ...
%   'XColor'      , [.3 .3 .3], ...
%   'YColor'      , [.3 .3 .3], ...
%   'LineWidth'   , 1.4       ,...
%   'xtick'       , Var_No    ,...
%   'FontSize'    , 12);
% % ylim([140 220])
% % xlim([0 4])





[RSquare_y , pValue_y] = corrcoef(RMSE_Sum,DTW_Sum)
[lin_rel, lin_rel_gof] = fit(RMSE_Sum,DTW_Sum,'poly1')
p1_lin = predint(lin_rel,linspace(min(RMSE_Sum)-.1,max(RMSE_Sum)+.1,50),0.95,'functional', 'off');

figu = figure(1022)
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperSize', [17 8]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 4 2]);
set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[17 1 17 12]);
hold on
% plot(linspace(.5,3.5,50),p1_cost,'m--')

fx = [linspace(min(RMSE_Sum)-.1,max(RMSE_Sum)+.1,50), fliplr(linspace(min(RMSE_Sum)-.1,max(RMSE_Sum)+.1,50))];
fy = [p1_lin(:,1)', fliplr(p1_lin(:,2)')];
% scatter(RMSE_Sum , DTW_Sum ,10, 'filled','MarkerFaceColor', [.6 .6 .6])

for i=1:8
    for j=1:31
        scatter(RMSE_Sum_Subjects_Speeds{j,i} , DTW_Sum_Subjects_Speeds{j,i} ,7, 'filled','MarkerFaceColor', 1-i/10*[1 1 1])
%         scatter(RMSE_Sum_Subjects_Speeds{j,i} , DTW_Sum_Subjects_Speeds{j,i} ,3, 'filled','MarkerFaceColor', [.0 1-i/10 .0])
    end
end

colormap gray;
% c = colorbar( 'Location', 'EastOutside', ...
%     'Ticks', Z_Bar, 'TickLabels', cellstr( num2str( Z_Bar(:)*Z_Scl, '%.3e' ) ) );
c = colorbar('direction','reverse', 'Location', 'EastOutside', ...
    'Ticks', [0, .5, 1], 'TickLabels', {'Fast' , 'Medium', 'Slow'} );
ylabel( c, 'Walking Speed' );

% c=parula(8);
% colormap(figu, flipud(colormap(figu)))
% colormap(gray(8))
% colorbar('direction','reverse','Ticks',['Slow' , 'Medium' , 'Fast'])

% colormap(flipud(colormap))

hold on
plot(linspace(min(RMSE_Sum)-.1,max(RMSE_Sum)+.1,50),lin_rel.p1*linspace(min(RMSE_Sum)-.1,max(RMSE_Sum)+.1,50)+lin_rel.p2,'r','LineWidth',1.4)

% fill_color = [.5 .0 .0];
fill(fx,fy,'r', 'FaceAlpha', .2 , 'edgecolor','none');

% scatter(split_belt_ratio_array , cost_plot ,60, 'filled', 'MarkerFaceColor', [.5 .5 .5], 'MarkerEdgeColor', 'k')

xlabel('Sum of RMS', 'FontSize', 14)
ylabel('Sum of DTW', 'FontSize', 14)
title('RMS vs. DTW', 'FontSize', 14)

text(.35,.4,['DTW_{fit}= ' num2str(lin_rel.p1,'%.2f') ' RMS  +  ' num2str(lin_rel.p2,'%.2f')],'Color','r', 'FontSize', 11)
% text(1.8,90,['r_{adj}= ' num2str(lin_rel_gof.adjrsquare,'%.2f' ) ' , p=' num2str(pvalue_cost(2,1),'%.2f')],'Color','r', 'FontSize', 11)
text(.35,.35,['r_{adj} = ' num2str(lin_rel_gof.adjrsquare^.5,'%.2f' ) ' ,  p < 0.001^{**}'],'Color','r', 'FontSize', 11)

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01] , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.3       ,...
  'FontSize'    , 14);
% ylim([140 220])






%%   The specifications of the plots for all subjects angle trajectories
% % % % %         figure(1033)
% % % % % sgtitle('Stance Shank Angle')
% % % % %     for idx = 1:42
% % % % %         subplot(6,7,idx)
% % % % %         title(['Subject ' num2str(idx)])
% % % % %         set(gca, ...
% % % % %         'Box'         , 'off'     , ...
% % % % %         'TickDir'     , 'out'     , ...
% % % % %         'TickLength'  , [.02 .02] , ...
% % % % %         'YGrid'       , 'on'      , ...
% % % % %         'XColor'      , [.3 .3 .3], ...
% % % % %         'YColor'      , [.3 .3 .3], ...
% % % % %         'LineWidth'   , 0.6        ,...
% % % % %         'FontSize'    , 10);  
% % % % %     end
% % % % %     set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[15 1 15 12]);
% % % % % 
% % % % %         figure(1044)
% % % % % sgtitle('Stance Thigh Angle')
% % % % %     for idx = 1:42
% % % % %         subplot(6,7,idx)
% % % % %         title(['Subject ' num2str(idx)])
% % % % %         set(gca, ...
% % % % %         'Box'         , 'off'     , ...
% % % % %         'TickDir'     , 'out'     , ...
% % % % %         'TickLength'  , [.02 .02] , ...
% % % % %         'YGrid'       , 'on'      , ...
% % % % %         'XColor'      , [.3 .3 .3], ...
% % % % %         'YColor'      , [.3 .3 .3], ...
% % % % %         'LineWidth'   , 0.6        ,...
% % % % %         'FontSize'    , 10);  
% % % % %     end
% % % % %     set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[15 1 15 12]);
% % % % % 
% % % % % 
% % % % %     figure(103)
% % % % %     sgtitle('Swing Thigh Angle')
% % % % %     for idx = 1:42
% % % % %         subplot(6,7,idx)
% % % % %         title(['Subject ' num2str(idx)])
% % % % %         set(gca, ...
% % % % %         'Box'         , 'off'     , ...
% % % % %         'TickDir'     , 'out'     , ...
% % % % %         'TickLength'  , [.02 .02] , ...
% % % % %         'YGrid'       , 'on'      , ...
% % % % %         'XColor'      , [.3 .3 .3], ...
% % % % %         'YColor'      , [.3 .3 .3], ...
% % % % %         'LineWidth'   , 0.6        ,...
% % % % %         'FontSize'    , 10);
% % % % %         
% % % % %     end
% % % % %     set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[15 1 15 12]);
% % % % %     
% % % % %     
% % % % %     
% % % % %         figure(104)
% % % % %     sgtitle('Swing Shank Angle')
% % % % %     for idx = 1:42
% % % % %         subplot(6,7,idx)
% % % % %         title(['Subject ' num2str(idx)])
% % % % %         set(gca, ...
% % % % %         'Box'         , 'off'     , ...
% % % % %         'TickDir'     , 'out'     , ...
% % % % %         'TickLength'  , [.02 .02] , ...
% % % % %         'YGrid'       , 'on'      , ...
% % % % %         'XColor'      , [.3 .3 .3], ...
% % % % %         'YColor'      , [.3 .3 .3], ...
% % % % %         'LineWidth'   , 0.6        ,...
% % % % %         'FontSize'    , 10);
% % % % %         
% % % % %     end
% % % % %     set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[15 1 15 12]);
%%
    
    

% % % %     
% % % %     ylabel('Swing Thigh (rad)')
% % % %     set(gca, ...
% % % %         'Box'         , 'off'     , ...
% % % %         'TickDir'     , 'out'     , ...
% % % %         'TickLength'  , [.02 .02] , ...
% % % %         'YGrid'       , 'on'      , ...
% % % %         'XColor'      , [.3 .3 .3], ...
% % % %         'YColor'      , [.3 .3 .3], ...
% % % %         'LineWidth'   , 1.4        ,...
% % % %         'FontSize'    , 14);
% % % %     set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[15 1 15 12]);

% dlmwrite('Cost_Comp_Eval_Traj.txt',cost_components_traj_eval)


analysis_date_str = string(datetime('now','Format', 'yyyy-MM-DD'));

dlmwrite(['Cost_Components_Eval_Fakuchi_AllData_NoDynamics_' analysis_date_str '.txt'],Cost_Components_Eval_Fakuchi)
dlmwrite(['Cost_Components_Eval_Fakuchi_AllData_' analysis_date_str '.txt'],Cost_Components_Eval_Fakuchi)



Cost_Components = Cost_Components_Eval_Fakuchi;

save(['Fukuchi_Features_DTW_RMS_' analysis_date_str '.mat'],...
    'Cost_Components_Subjects', 'Cost_Components_Subjects_Speeds', 'Cost_Components' , ...  
    'DTW_Sum_Subjects'        , 'DTW_Sum_Subjects_Speeds'        , 'DTW_Sum'         , ...
    'RMSE_Sum_Subjects'       , 'RMSE_Sum_Subjects_Speeds'       , 'RMSE_Sum'               );



% for i_y = 1:size(Q1,1)
%     Deviation_Eval_mean(i_y,1) = sqrt(mean( (Q1(i_y,:)-mean(Q1,1)).^2 ))+...
%                             sqrt(mean((Q2(i_y,:)- mean(Q2,1)).^2))+...
%                             sqrt(mean((Q4(i_y,:)- mean(Q4,1)).^2))+...
%                             sqrt(mean((Q5(i_y,:)- mean(Q5,1)).^2));
% end
    

dlmwrite(['Cost_Deviation_Eval_y_forPLSR_' analysis_date_str '.txt'],Deviation_Eval_mean)

Q_ensemble(1,:) = mean(Q1,1);
Q_ensemble(2,:) = mean(Q2,1);
Q_ensemble(3,:) = mean(Q4,1);
Q_ensemble(4,:) = mean(Q5,1);
% Q1 = Q(:,:,1);
% Q2 = Q(:,:,2);
% Q4 = Q(:,:,4);
% Q5 = Q(:,:,5);

figure(1003)
stdshade(Q1,.1,[1 0 1],(0:length(t_resamp)-1),1)
hold on
stdshade(Q2,.1,[1 0 0],(0:length(t_resamp)-1),1)
stdshade(Q4,.1,[0 0 1],(0:length(t_resamp)-1),1)
stdshade(Q5,.1,[0 1 1],(0:length(t_resamp)-1),1)
plot([Q_ensemble'],'k:')
legend('','Stance Shank','','Stance Thigh','','Swing Thigh','','Swing Shank');
xlabel('Gait Cycle %')
ylabel('Segment Angles (rad)')
set(gca, ...
'Box'         , 'off'     , ...
'TickDir'     , 'out'     , ...
'TickLength'  , [.02 .02] , ...
'YGrid'       , 'on'      , ...
'XColor'      , [.3 .3 .3], ...
'YColor'      , [.3 .3 .3], ...
'LineWidth'   , 1.4        ,...
'FontSize'    , 14);



% cost_vector = [.15* u.^2  ;...                          % torques squared
%                 .01* du.^2 ;...                         % Torque Change squared
%                     .4* abs(u.*(dq1-dq2))  ;...         % absolute work
%                         1/40* (ddp1.^2+ddp2.^2);...     % Accel squared
%                            .01* (d3p1.^2+d3p2.^2)];     % Jerk

% cost_vector = [50* u.^2  ;...                          % torques squared
%                 1* du.^2 ;...                         % Torque Change squared
%                     10* abs(u.*(dq1-dq2))  ;...         % absolute work
%                         .25* (ddp1.^2+ddp2.^2);...     % Accel squared
%                            .01* (d3p1.^2+d3p2.^2)];     % Jerk

% % % len_cost = size(cost_components_traj_eval,2)
% % % for i = 1:len_cost
% % %     for j = i:len_cost
% % %         
% % % 
% % %         Rsq(i,j) = corr2(cost_components_traj_eval(:,i),cost_components_traj_eval(:,j));
% % %         
% % %         
% % %     end
% % % end


len_cost = size(Cost_Components_Eval_Fakuchi,2);


nboot = 1000;
bootfun = @R_Squared_Bootstrap;
RSSQ = bootstrp(nboot,bootfun,Cost_Components_Eval_Fakuchi);

RSSQ_reshap = reshape(RSSQ',len_cost,len_cost,nboot);
Rsq_Mean = mean(RSSQ_reshap,3)';
Rsq_SD = std(RSSQ_reshap,[],3)';

Rsq = Rsq_Mean;






idd_subj = 0;
% for id_subj = 1:No_Subj
for id_subj = 1:31
    if id_subj~=5 && id_subj~=17 && id_subj~=41
        idd_subj = idd_subj+1;
        Rsquared(idd_subj,:,:) = reshape(R_Squared_Bootstrap(Cost_Components_Subjects{idd_subj,1}),len_cost,len_cost)';
    end
end


Rsq_Mean = squeeze(mean(Rsquared,1));
Rsq_SD = squeeze(std(Rsquared,[],1));
Rsq = Rsq_Mean;




figure(12)
% % imagesc(abs(Rsq))
imagesc(Rsq_Mean-Rsq_SD)
% impixelregion(imagesc(Rsq))
set(gca, 'XTIck', [1:len_cost])
set(gca, 'YTIck', [1:len_cost])
set(gca, 'XAxisLocation', 'top')
set(gca, 'YAxisLocation', 'right')
colorbar('Location','eastoutside')
colorbar('Location','southoutside')
colormap Bone
textStrings = num2str(abs(Rsq), '%0.2f');       % Create strings from the matrix values
% textStrings = [num2str(abs(Rsq_Mean), '%0.2f') '+/-' num2str(abs(Rsq_SD), '%0.2f')];       % Create strings from the matrix values
% textStrings = num2str(Rsq_Mean, '%0.2f');       % Create strings from the matrix values
textStrings_SD = num2str(Rsq_SD, '%0.2f');       % Create strings from the matrix values
% textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
text_modif = strings(len_cost,len_cost);
text_modif_Mean = strings(len_cost,len_cost);
text_modif_SD = strings(len_cost,len_cost);
for k = 1:len_cost
for i = 1:len_cost
if Rsq(i,k)<0
txt_sgn = '-';
else
txt_sgn = '';
end
%     text_modif (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
text_modif (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k), char(10), char(177), textStrings_SD(i,4*(k-1)+1:4*k)];
text_modif_Mean (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
text_modif_SD (i,k)= [char(177), textStrings_SD(i,4*(k-1)+1:4*k)];
end
end
[xs, ys] = meshgrid(1:len_cost);  % Create x and y coordinates for the strings
text(xs(:), ys(:), text_modif(:), ...  % Plot the strings
'HorizontalAlignment', 'center');
set( gca, 'XTickLabel', {'1.angl accel sqrd' '2.angl jerk sqrd' '3.cartes accel sqrd' '4.cartes jerk sqrd' '5.torqs sqrd' '6.torqs abs' '7.torq rate sqrd' '8.torq rate abs' '9.torq rate change sqrd' '10.torq rate change abs' '11.yGRF rate sqrd' '12.xGRF rate sqrd' '13.positive work' '14.absolute work' '15.pk kinetic energy' '16.pk2pk angular momentum'}','FontSize',12 )
xtickangle(70)
set( gca, 'YTickLabel', {'1.angl accel sqrd' '2.angl jerk sqrd' '3.cartes accel sqrd' '4.cartes jerk sqrd' '5.torqs sqrd' '6.torqs abs' '7.torq rate sqrd' '8.torq rate abs' '9.torq rate change sqrd' '10.torq rate change abs' '11.yGRF rate sqrd' '12.xGRF rate sqrd' '13.positive work' '14.absolute work' '15.pk kinetic energy' '16.pk2pk angular momentum'}','FontSize',12 )
% set( gca, 'YTickLabel', {'angl accel sqrd' 'angl jerk sqrd' 'cartes accel sqrd' 'cartes jerk sqrd' 'torqs sqrd' 'torqs abs' 'torq rate sqrd' 'torq rate abs' 'torq rate change sqrd' 'torq rate change abs' 'yGRF rate sqrd' 'xGRF rate sqrd' 'positive work' 'absolute work' 'kinetic energy' 'angular momentum'}','FontSize',12 )






Rsq_group = R_Squared_Bootstrap (Cost_Components_Eval_Fakuchi);
Rsq_group = reshape(Rsq_group,len_cost,len_cost)';

figure(14)
% % imagesc(abs(Rsq))
imagesc(Rsq_group)
% impixelregion(imagesc(Rsq))
set(gca, 'XTIck', [1:len_cost])
set(gca, 'YTIck', [1:len_cost])
set(gca, 'XAxisLocation', 'top')
set(gca, 'YAxisLocation', 'right')
colorbar('Location','eastoutside')
colorbar('Location','southoutside')
colormap Bone

textStrings = num2str(abs(Rsq_group), '%0.2f');       % Create strings from the matrix values
% textStrings = [num2str(abs(Rsq_Mean), '%0.2f') '+/-' num2str(abs(Rsq_SD), '%0.2f')];       % Create strings from the matrix values
% textStrings = num2str(Rsq_Mean, '%0.2f');       % Create strings from the matrix values
% textStrings_SD = num2str(Rsq_SD, '%0.2f');       % Create strings from the matrix values

% textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
text_modif = strings(len_cost,len_cost);
% text_modif_Mean = strings(len_cost,len_cost);
% text_modif_SD = strings(len_cost,len_cost);

for k = 1:len_cost
    for i = 1:len_cost
    if Rsq(i,k)<0
        txt_sgn = '-';
    else
        txt_sgn = '';
    end
%     text_modif (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
    text_modif (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
%     text_modif_Mean (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
%     text_modif_SD (i,k)= [char(177), textStrings_SD(i,4*(k-1)+1:4*k)];
    end
end

[xs, ys] = meshgrid(1:len_cost);  % Create x and y coordinates for the strings
text(xs(:), ys(:), text_modif(:), ...  % Plot the strings
                'HorizontalAlignment', 'center');
            

 set( gca, 'XTickLabel', {'1.angl accel sqrd' '2.angl jerk sqrd' '3.cartes accel sqrd' '4.cartes jerk sqrd' '5.torqs sqrd' '6.torqs abs' '7.torq rate sqrd' '8.torq rate abs' '9.torq rate change sqrd' '10.torq rate change abs' '11.yGRF rate sqrd' '12.xGRF rate sqrd' '13.positive work' '14.absolute work' '15.pk kinetic energy' '16.pk2pk angular momentum'}','FontSize',12 )
xtickangle(70)           

set( gca, 'YTickLabel', {'1.angl accel sqrd' '2.angl jerk sqrd' '3.cartes accel sqrd' '4.cartes jerk sqrd' '5.torqs sqrd' '6.torqs abs' '7.torq rate sqrd' '8.torq rate abs' '9.torq rate change sqrd' '10.torq rate change abs' '11.yGRF rate sqrd' '12.xGRF rate sqrd' '13.positive work' '14.absolute work' '15.pk kinetic energy' '16.pk2pk angular momentum'}','FontSize',12 )
 
% set( gca, 'YTickLabel', {'angl accel sqrd' 'angl jerk sqrd' 'cartes accel sqrd' 'cartes jerk sqrd' 'torqs sqrd' 'torqs abs' 'torq rate sqrd' 'torq rate abs' 'torq rate change sqrd' 'torq rate change abs' 'yGRF rate sqrd' 'xGRF rate sqrd' 'positive work' 'absolute work' 'kinetic energy' 'angular momentum'}','FontSize',12 )

            















%% MSPE
figure(13)
% % imagesc(abs(Rsq))
imagesc(MSPE.data)
% impixelregion(imagesc(Rsq))
set(gca, 'XTIck', [1:16])
set(gca, 'YTIck', [1:16])
set(gca, 'XAxisLocation', 'bottom')
set(gca, 'YAxisLocation', 'left')
colorbar('Location','eastoutside')
colorbar('Location','southoutside')
colormap cool

xlabel('Number of PCs kept in the model')
ylabel('Number of factors kept in PCs (Sparsity)')
title('Mean-squared Error Predicted across different PC numbers and sparsity')

textStrings = num2str(abs(MSPE.data), '%0.3f');       % Create strings from the matrix values
% textStrings = [num2str(abs(Rsq_Mean), '%0.2f') '+/-' num2str(abs(Rsq_SD), '%0.2f')];       % Create strings from the matrix values
% textStrings = num2str(Rsq_Mean, '%0.2f');       % Create strings from the matrix values
% textStrings_SD = num2str(Rsq_SD, '%0.2f');       % Create strings from the matrix values
% textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
text_modif = strings(len_cost,len_cost);
for k = 1:len_cost
for i = 1:len_cost
%     text_modif (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
text_modif (i,k)= [textStrings(i,5*(k-1)+1:5*k)];
%    text_modif_Mean (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
%    text_modif_SD (i,k)= [char(177), textStrings_SD(i,4*(k-1)+1:4*k)];
end
end
[xs, ys] = meshgrid(1:len_cost);  % Create x and y coordinates for the strings
text(xs(:), ys(:), text_modif(:), ...  % Plot the strings
'HorizontalAlignment', 'center');











function Rsq = R_Squared_Bootstrap (Cost_Components_Eval_Fakuchi)

len_cost = size(Cost_Components_Eval_Fakuchi,2);
for i = 1:len_cost
    for j = i:len_cost
        

        Rsq(i,j) = corr2(Cost_Components_Eval_Fakuchi(:,i),Cost_Components_Eval_Fakuchi(:,j));
        
        
    end
end
Rsq = reshape(Rsq',1,len_cost*len_cost);
end
















%% Interpolation Functions

function x = pwPoly2(tGrid,xGrid,t)
% x = pwPoly2(tGrid,xGrid,t)
%
% This function does piece-wise quadratic interpolation of a set of data,
% given the function value at the edges and midpoint of the interval of
% interest.
%
% INPUTS:
%   tGrid = [1, 2*n-1] = time grid, knot idx = 1:2:end
%   xGrid = [m, 2*n-1] = function at each grid point in tGrid
%   t = [1, k] = vector of query times (must be contained within tGrid)
%
% OUTPUTS:
%   x = [m, k] = function value at each query time
%
% NOTES:
%   If t is out of bounds, then all corresponding values for x are replaced
%   with NaN
%

nGrid = length(tGrid);
if mod(nGrid-1,2)~=0 || nGrid < 3
    error('The number of grid-points must be odd and at least 3');
end

% Figure out sizes
n = floor((length(tGrid)-1)/2);
m = size(xGrid,1);
k = length(t);
x = zeros(m, k);

% Figure out which segment each value of t should be on
edges = [-inf, tGrid(1:2:end), inf];
[~, bin] = histc(t,edges);

% Loop over each quadratic segment
for i=1:n
    idx = bin==(i+1);
    if sum(idx) > 0
        gridIdx = 2*(i-1) + [1,2,3];
        x(:,idx) = quadInterp(tGrid(gridIdx),xGrid(:,gridIdx),t(idx));
    end
end

% Replace any out-of-bounds queries with NaN
outOfBounds = bin==1 | bin==(n+2);
x(:,outOfBounds) = nan;

% Check for any points that are exactly on the upper grid point:
if sum(t==tGrid(end))>0
    x(:,t==tGrid(end)) = xGrid(:,end);
end

end


function x = quadInterp(tGrid,xGrid,t)
%
% This function computes the interpolant over a single interval
%
% INPUTS:
%   tGrid = [1, 3] = time grid
%   xGrid = [m, 3] = function grid
%   t = [1, p] = query times, spanned by tGrid
%
% OUTPUTS:
%   x = [m, p] = function at query times
%

% Rescale the query points to be on the domain [-1,1]
t = 2*(t-tGrid(1))/(tGrid(3)-tGrid(1)) - 1;

% Compute the coefficients:
a = 0.5*(xGrid(:,3) + xGrid(:,1)) - xGrid(:,2);
b = 0.5*(xGrid(:,3)-xGrid(:,1));
c = xGrid(:,2);

% Evaluate the polynomial for each dimension of the function:
p = length(t);
m = size(xGrid,1);
x = zeros(m,p);
tt = t.^2;
for i=1:m
    x(i,:) = a(i)*tt + b(i)*t + c(i);
end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Functions for interpolation of the state solution
%




function x = pwPoly3(tGrid,xGrid,fGrid,t)
% x = pwPoly3(tGrid,xGrid,fGrid,t)
%
% This function does piece-wise quadratic interpolation of a set of data,
% given the function value at the edges and midpoint of the interval of
% interest.
%
% INPUTS:
%   tGrid = [1, 2*n-1] = time grid, knot idx = 1:2:end
%   xGrid = [m, 2*n-1] = function at each grid point in time
%   fGrid = [m, 2*n-1] = derivative at each grid point in time
%   t = [1, k] = vector of query times (must be contained within tGrid)
%
% OUTPUTS:
%   x = [m, k] = function value at each query time
%
% NOTES:
%   If t is out of bounds, then all corresponding values for x are replaced
%   with NaN
%

nGrid = length(tGrid);
if mod(nGrid-1,2)~=0 || nGrid < 3
    error('The number of grid-points must be odd and at least 3');
end

% Figure out sizes
n = floor((length(tGrid)-1)/2);
m = size(xGrid,1);
k = length(t);
x = zeros(m, k);

% Figure out which segment each value of t should be on
edges = [-inf, tGrid(1:2:end), inf];
[~, bin] = histc(t,edges);

% Loop over each quadratic segment
for i=1:n
    idx = bin==(i+1);
    if sum(idx) > 0
        kLow = 2*(i-1) + 1;
        kMid = kLow + 1;
        kUpp = kLow + 2;
        h = tGrid(kUpp)-tGrid(kLow);
        xLow = xGrid(:,kLow);
        fLow = fGrid(:,kLow);
        fMid = fGrid(:,kMid);
        fUpp = fGrid(:,kUpp);
        alpha = t(idx) - tGrid(kLow);
        x(:,idx) = cubicInterp(h,xLow, fLow, fMid, fUpp,alpha);
    end
end

% Replace any out-of-bounds queries with NaN
outOfBounds = bin==1 | bin==(n+2);
x(:,outOfBounds) = nan;

% Check for any points that are exactly on the upper grid point:
if sum(t==tGrid(end))>0
    x(:,t==tGrid(end)) = xGrid(:,end);
end

end


function x = cubicInterp(h,xLow, fLow, fMid, fUpp,del)
%
% This function computes the interpolant over a single interval
%
% INPUTS:
%   h = time step (tUpp-tLow)
%   xLow = function value at tLow
%   fLow = derivative at tLow
%   fMid = derivative at tMid
%   fUpp = derivative at tUpp
%   del = query points on domain [0, h]
%
% OUTPUTS:
%   x = [m, p] = function at query times
%

%%% Fix matrix dimensions for vectorized calculations
nx = length(xLow);
nt = length(del);
xLow = xLow*ones(1,nt);
fLow = fLow*ones(1,nt);
fMid = fMid*ones(1,nt);
fUpp = fUpp*ones(1,nt);
del = ones(nx,1)*del;

a = (2.*(fLow - 2.*fMid + fUpp))./(3.*h.^2);
b = -(3.*fLow - 4.*fMid + fUpp)./(2.*h);
c = fLow;
d = xLow;

x = d + del.*(c + del.*(b + del.*a));

end






