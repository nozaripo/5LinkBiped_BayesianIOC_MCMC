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

for i=1:length(step_time)
    




time = Time(r_TO(i):r_HS(i)-1,1)-Time(r_TO(i),1);


cycl_err(i,1) = norm(q(:,r_TO(i))-q(:,r_HS(i)));


        
        
        % .^.5
        
Accel_Joint = (sum(ddq(:,r_TO(i):r_HS(i)-1).^2));
Jerk_Joint = (sum(dddq(:,r_TO(i):r_HS(i)-1).^2));

Accel_Task  = (sum(ddG(:,r_TO(i):r_HS(i)-1).^2));
Jerk_Task  = (sum(dddG(:,r_TO(i):r_HS(i)-1).^2));

% Time = 1;

% Velocity_Deviation = dP(2,:)


Torque_Squared           = (sum(u(:,r_TO(i):r_HS(i)-1).^2));
Torque_Absolute          = sum(abs(u(:,r_TO(i):r_HS(i)-1)));
TorqueRate_Squared       = (sum(du(:,r_TO(i):r_HS(i)-1).^2));
TorqueRate_Absolute      = sum(abs(du(:,r_TO(i):r_HS(i)-1)));
TorqueRateChange_Squared = (sum(ddu(:,r_TO(i):r_HS(i)-1).^2));
TorqueRateChange_Absolute= sum(abs(ddu(:,r_TO(i):r_HS(i)-1)));

yGRF_Rate    =  (dFy(:,r_TO(i):r_HS(i)-1).^2);
xGRF_Rate    =  (dFx(:,r_TO(i):r_HS(i)-1).^2);


Work_Absolute = sum(abs(u(:,r_TO(i):r_HS(i)-1).*dq(:,r_TO(i):r_HS(i)-1)));
Work_Positive = sum(max(u(:,r_TO(i):r_HS(i)-1).*dq(:,r_TO(i):r_HS(i)-1),0));

Kinetic_Energy = 1/2*(p.I1+p.m1*(p.l1-p.c1)^2)*dq(1,r_TO(i):r_HS(i)-1).^2 + ...
                  1/2*(p.I2+p.m2*(p.l2-p.c2)^2)*dq(2,r_TO(i):r_HS(i)-1).^2 + ...
                  1/2*(p.I3+p.m3*(p.l3-p.c3)^2)*dq(3,r_TO(i):r_HS(i)-1).^2 + ... ;   % Kinetic Energy
                  1/2*(p.I4+p.m4*(p.c4)^2)*dq(4,r_TO(i):r_HS(i)-1).^2 + ...
                  1/2*(p.I5+p.m5*(p.c5)^2)*dq(5,r_TO(i):r_HS(i)-1).^2 ;

Angular_Momentum = (p.I1+p.m1*(p.l1-p.c1)^2)*dq(1,r_TO(i):r_HS(i)-1) + ...
                  (p.I2+p.m2*(p.l2-p.c2)^2)*dq(2,r_TO(i):r_HS(i)-1) + ...
                  (p.I3+p.m3*(p.l3-p.c3)^2)*dq(3,r_TO(i):r_HS(i)-1) + ... ;   % Kinetic Energy
                  (p.I4+p.m4*(p.c4)^2)*dq(4,r_TO(i):r_HS(i)-1) + ...
                  (p.I5+p.m5*(p.c5)^2)*dq(5,r_TO(i):r_HS(i)-1) ;

              
% COST COMPONENTS EVALS

% kinematic
cost_components_eval(i,1) = simps(time, Accel_Joint');
cost_components_eval(i,2) = simps(time, Jerk_Joint');
cost_components_eval(i,3) = simps(time, Accel_Task');
cost_components_eval(i,4) = simps(time, Jerk_Task');
cost_components_eval(i,5) = time(end);

% dynamic
cost_components_eval(i,6) = simps(time, Torque_Squared');
cost_components_eval(i,7) = simps(time, Torque_Absolute');
cost_components_eval(i,8) = simps(time, TorqueRate_Squared');
cost_components_eval(i,9) = simps(time, TorqueRate_Absolute');
cost_components_eval(i,10) = simps(time, TorqueRateChange_Squared');
cost_components_eval(i,11) = simps(time, TorqueRateChange_Absolute');
cost_components_eval(i,12) = simps(time, yGRF_Rate');
cost_components_eval(i,13) = simps(time, xGRF_Rate');

% energetic
cost_components_eval(i,14) = simps(time, Work_Positive');
cost_components_eval(i,15) = simps(time, Work_Absolute');
cost_components_eval(i,16) = max(Kinetic_Energy)/time(end);

% balance / stability
cost_components_eval(i,17) = (max(Angular_Momentum)-min(Angular_Momentum))/time(end);       
        
        
% cost_bases_eval(i,:) = simps(time,C');
%        
% % cost_vector = C./scaling_vector;
% clear C
clear time
end

dlmwrite('Cost_Comp_Eval.txt',cost_components_eval)
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

len_cost = size(cost_components_eval,2)
for i = 1:len_cost
    for j = i:len_cost
        

        Rsq(i,j) = corr2(cost_components_eval(:,i),cost_components_eval(:,j));
        
        
    end
end



figure(9)
imagesc(abs(Rsq))
% impixelregion(imagesc(Rsq))
set(gca, 'XTIck', [1:len_cost])
set(gca, 'YTIck', [1:len_cost])
set(gca, 'XAxisLocation', 'top')
set(gca, 'YAxisLocation', 'right')
colorbar('Location','eastoutside')
colorbar('Location','southoutside')
colormap Bone

textStrings = num2str(abs(Rsq), '%0.2f');       % Create strings from the matrix values
% textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
text_modif = strings(len_cost,len_cost);

for k = 1:len_cost
    for i = 1:len_cost
    if Rsq(i,k)<0
        txt_sgn = '-';
    else
        txt_sgn = '';
    end
    text_modif (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
    end
end

[xs, ys] = meshgrid(1:len_cost);  % Create x and y coordinates for the strings
text(xs(:), ys(:), text_modif(:), ...  % Plot the strings
                'HorizontalAlignment', 'center');
            

 set( gca, 'XTickLabel', {'1.angl accel sqrd' '2.angl jerk sqrd' '3.cartes accel sqrd' '4.cartes jerk sqrd' '5.time' '6.torqs sqrd' '7.torqs abs' '8.torq rate sqrd' '9.torq rate abs' '10.torq rate change sqrd' '11.torq rate change abs' '12.yGRF rate sqrd' '13.xGRF rate sqrd' '14.positive work' '15.absolute work' '16.kinetic energy' '17.angular momentum'}','FontSize',12 )
xtickangle(70)           

 set( gca, 'YTickLabel', {'angl accel sqrd' 'angl jerk sqrd' 'cartes accel sqrd' 'cartes jerk sqrd' 'time' 'torqs sqrd' 'torqs abs' 'torq rate sqrd' 'torq rate abs' 'torq rate change sqrd' 'torq rate change abs' 'yGRF rate sqrd' 'xGRF rate sqrd' 'positive work' 'absolute work' 'kinetic energy' 'angular momentum'}','FontSize',12 )

            
            
% % %  set( gca, 'XTickLabel', {'torqs sqrd' 'torqs abs' 'torq chng sqrd' 'torq chng abs' 'torq 2der sqrd' 'torq 2der abs' 'abs work' 'positv work' 'cartes accel sqrd' 'cartes jerk sqrd' 'angl accel sqrd' 'angl jerk sqrd' 'time' 'dev from vel_{ref}' 'kinetic energy' 'kinetic energy2'}','FontSize',12 )
% % % xtickangle(60)           
% % % 
% % %  set( gca, 'YTickLabel', {'torqs sqrd' 'torqs abs' 'torq chng sqrd' 'torq chng abs' 'torq 2der sqrd' 'torq 2der abs' 'abs work' 'positv work' 'cartes accel sqrd' 'cartes jerk sqrd' 'angl accel sqrd' 'angl jerk sqrd' 'time' 'dev from vel_{ref}' 'kinetic energy' 'kinetic energy2'}','FontSize',12 )


% figure(10)
% corrplot(cost_components_eval)
%                     
% cost = W*cost_bases_eval'; 

% end




%% correlation based on trajecotry
for i = 1:length(r_TO)
% cycl_err(i,1) = norm(q(:,r_TO(i))-q(5:-1:1,r_HS(i)));
cycl_err(i,1) = norm([l_hip(r_TO(i))-r_hip(r_HS(i)); l_knee(r_TO(i))-r_knee(r_HS(i))]);
end
indx = find(cycl_err == min(cycl_err));

cycl_err(cycl_err == min(cycl_err))



for i=indx:indx
    




time = Time(r_TO(i):r_HS(i)-1,1)-Time(r_TO(i),1);


cycl_err(i,1) = norm(q(:,r_TO(i))-q(:,r_HS(i)));


        
        
        % .^.5
        
Accel_Joint = (sum(ddq(:,r_TO(i):r_HS(i)-1).^2));
Jerk_Joint = (sum(dddq(:,r_TO(i):r_HS(i)-1).^2));

Accel_Task  = (sum(ddG(:,r_TO(i):r_HS(i)-1).^2));
Jerk_Task  = (sum(dddG(:,r_TO(i):r_HS(i)-1).^2));

% Time = 1;

% Velocity_Deviation = dP(2,:)


Torque_Squared           = (sum(u(:,r_TO(i):r_HS(i)-1).^2));
Torque_Absolute          = sum(abs(u(:,r_TO(i):r_HS(i)-1)));
TorqueRate_Squared       = (sum(du(:,r_TO(i):r_HS(i)-1).^2));
TorqueRate_Absolute      = sum(abs(du(:,r_TO(i):r_HS(i)-1)));
TorqueRateChange_Squared = (sum(ddu(:,r_TO(i):r_HS(i)-1).^2));
TorqueRateChange_Absolute= sum(abs(ddu(:,r_TO(i):r_HS(i)-1)));

yGRF_Rate    =  (dFy(:,r_TO(i):r_HS(i)-1).^2);
xGRF_Rate    =  (dFx(:,r_TO(i):r_HS(i)-1).^2);


Work_Absolute = sum(abs(u(:,r_TO(i):r_HS(i)-1).*dq(:,r_TO(i):r_HS(i)-1)));
Work_Positive = sum(max(u(:,r_TO(i):r_HS(i)-1).*dq(:,r_TO(i):r_HS(i)-1),0));

Kinetic_Energy = 1/2*(p.I1+p.m1*(p.l1-p.c1)^2)*dq(1,r_TO(i):r_HS(i)-1).^2 + ...
                  1/2*(p.I2+p.m2*(p.l2-p.c2)^2)*dq(2,r_TO(i):r_HS(i)-1).^2 + ...
                  1/2*(p.I3+p.m3*(p.l3-p.c3)^2)*dq(3,r_TO(i):r_HS(i)-1).^2 + ... ;   % Kinetic Energy
                  1/2*(p.I4+p.m4*(p.c4)^2)*dq(4,r_TO(i):r_HS(i)-1).^2 + ...
                  1/2*(p.I5+p.m5*(p.c5)^2)*dq(5,r_TO(i):r_HS(i)-1).^2 ;

Angular_Momentum = (p.I1+p.m1*(p.l1-p.c1)^2)*dq(1,r_TO(i):r_HS(i)-1) + ...
                  (p.I2+p.m2*(p.l2-p.c2)^2)*dq(2,r_TO(i):r_HS(i)-1) + ...
                  (p.I3+p.m3*(p.l3-p.c3)^2)*dq(3,r_TO(i):r_HS(i)-1) + ... ;   % Kinetic Energy
                  (p.I4+p.m4*(p.c4)^2)*dq(4,r_TO(i):r_HS(i)-1) + ...
                  (p.I5+p.m5*(p.c5)^2)*dq(5,r_TO(i):r_HS(i)-1) ;

              
% COST COMPONENTS EVALS

% kinematic
cost_components_traj_eval(:,1) = Accel_Joint';
cost_components_traj_eval(:,2) = Jerk_Joint';
cost_components_traj_eval(:,3) = Accel_Task';
cost_components_traj_eval(:,4) = Jerk_Task';
cost_components_traj_eval(:,5) = time;

% dynamic
cost_components_traj_eval(:,6) = Torque_Squared';
cost_components_traj_eval(:,7) = Torque_Absolute';
cost_components_traj_eval(:,8) = TorqueRate_Squared';
cost_components_traj_eval(:,9) = TorqueRate_Absolute';
cost_components_traj_eval(:,10) = TorqueRateChange_Squared';
cost_components_traj_eval(:,11) = TorqueRateChange_Absolute';
cost_components_traj_eval(:,12) = yGRF_Rate';
cost_components_traj_eval(:,13) = xGRF_Rate';

% energetic
cost_components_traj_eval(:,14) = Work_Positive';
cost_components_traj_eval(:,15) = Work_Absolute';
cost_components_traj_eval(:,16) = Kinetic_Energy';

% balance / stability
cost_components_traj_eval(:,17) = max(Angular_Momentum)-min(Angular_Momentum);       
        
        
% cost_bases_eval(i,:) = simps(time,C');
%        
% % cost_vector = C./scaling_vector;
% clear C
clear time
end

dlmwrite('Cost_Comp_Eval_Traj.txt',cost_components_traj_eval)
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

len_cost = size(cost_components_traj_eval,2)
for i = 1:len_cost
    for j = i:len_cost
        

        Rsq(i,j) = corr2(cost_components_traj_eval(:,i),cost_components_traj_eval(:,j));
        
        
    end
end



figure(10)
imagesc(abs(Rsq))
% impixelregion(imagesc(Rsq))
set(gca, 'XTIck', [1:len_cost])
set(gca, 'YTIck', [1:len_cost])
set(gca, 'XAxisLocation', 'top')
set(gca, 'YAxisLocation', 'right')
colorbar('Location','eastoutside')
colorbar('Location','southoutside')
colormap Bone

textStrings = num2str(abs(Rsq), '%0.2f');       % Create strings from the matrix values
% textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
text_modif = strings(len_cost,len_cost);

for k = 1:len_cost
    for i = 1:len_cost
    if Rsq(i,k)<0
        txt_sgn = '-';
    else
        txt_sgn = '';
    end
    text_modif (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
    end
end

[xs, ys] = meshgrid(1:len_cost);  % Create x and y coordinates for the strings
text(xs(:), ys(:), text_modif(:), ...  % Plot the strings
                'HorizontalAlignment', 'center');
            

 set( gca, 'XTickLabel', {'1.angl accel sqrd' '2.angl jerk sqrd' '3.cartes accel sqrd' '4.cartes jerk sqrd' '5.time' '6.torqs sqrd' '7.torqs abs' '8.torq rate sqrd' '9.torq rate abs' '10.torq rate change sqrd' '11.torq rate change abs' '12.yGRF rate sqrd' '13.xGRF rate sqrd' '14.positive work' '15.absolute work' '16.kinetic energy' '17.angular momentum'}','FontSize',12 )
xtickangle(70)           

 set( gca, 'YTickLabel', {'angl accel sqrd' 'angl jerk sqrd' 'cartes accel sqrd' 'cartes jerk sqrd' 'time' 'torqs sqrd' 'torqs abs' 'torq rate sqrd' 'torq rate abs' 'torq rate change sqrd' 'torq rate change abs' 'yGRF rate sqrd' 'xGRF rate sqrd' 'positive work' 'absolute work' 'kinetic energy' 'angular momentum'}','FontSize',12 )

            

 
 
 
 
 
 
 
 
 
 
 
 clear all
 
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

param = getPhysicalParameters();
p=param;

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

q = [l_hip-l_knee , l_hip , zeros(length(l_hip),1) , r_hip , r_hip-r_knee]';



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

dq = diff(q')'/(Time(end)/(length(Time)-1));
dq = [dq(:,1) dq];

ddq = diff(dq')'/(Time(end)/(length(Time)-1));
ddq = [ddq(:,1) ddq];

dddq = diff(ddq')'/(Time(end)/(length(Time)-1));
dddq = [dddq(:,1) dddq];


% task kinematics
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
cost_components_traj_eval(:,5) = time;

% dynamic
cost_components_traj_eval(:,6) = Torque_Squared';
cost_components_traj_eval(:,7) = Torque_Absolute';
cost_components_traj_eval(:,8) = TorqueRate_Squared';
cost_components_traj_eval(:,9) = TorqueRate_Absolute';
cost_components_traj_eval(:,10) = TorqueRateChange_Squared';
cost_components_traj_eval(:,11) = TorqueRateChange_Absolute';
cost_components_traj_eval(:,12) = yGRF_Rate';
cost_components_traj_eval(:,13) = xGRF_Rate';

% energetic
cost_components_traj_eval(:,14) = Work_Positive';
cost_components_traj_eval(:,15) = Work_Absolute';
cost_components_traj_eval(:,16) = max(Kinetic_Energy')/t(end);

% balance / stability
cost_components_traj_eval(:,17) = max(Angular_Momentum)-min(Angular_Momentum);       
        
        
% cost_bases_eval(i,:) = simps(time,C');
%        
% % cost_vector = C./scaling_vector;
% clear C


dlmwrite('Cost_Comp_Eval_Traj.txt',cost_components_traj_eval)
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

len_cost = size(cost_components_traj_eval,2)
for i = 1:len_cost
    for j = i:len_cost
        

        Rsq(i,j) = corr2(cost_components_traj_eval(:,i),cost_components_traj_eval(:,j));
        
        
    end
end



figure(11)
imagesc(abs(Rsq))
% impixelregion(imagesc(Rsq))
set(gca, 'XTIck', [1:len_cost])
set(gca, 'YTIck', [1:len_cost])
set(gca, 'XAxisLocation', 'top')
set(gca, 'YAxisLocation', 'right')
colorbar('Location','eastoutside')
colorbar('Location','southoutside')
colormap Bone

textStrings = num2str(abs(Rsq), '%0.2f');       % Create strings from the matrix values
% textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
text_modif = strings(len_cost,len_cost);

for k = 1:len_cost
    for i = 1:len_cost
    if Rsq(i,k)<0
        txt_sgn = '-';
    else
        txt_sgn = '';
    end
    text_modif (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
    end
end

[xs, ys] = meshgrid(1:len_cost);  % Create x and y coordinates for the strings
text(xs(:), ys(:), text_modif(:), ...  % Plot the strings
                'HorizontalAlignment', 'center');
            

 set( gca, 'XTickLabel', {'1.angl accel sqrd' '2.angl jerk sqrd' '3.cartes accel sqrd' '4.cartes jerk sqrd' '5.time' '6.torqs sqrd' '7.torqs abs' '8.torq rate sqrd' '9.torq rate abs' '10.torq rate change sqrd' '11.torq rate change abs' '12.yGRF rate sqrd' '13.xGRF rate sqrd' '14.positive work' '15.absolute work' '16.kinetic energy' '17.angular momentum'}','FontSize',12 )
xtickangle(70)           

 set( gca, 'YTickLabel', {'angl accel sqrd' 'angl jerk sqrd' 'cartes accel sqrd' 'cartes jerk sqrd' 'time' 'torqs sqrd' 'torqs abs' 'torq rate sqrd' 'torq rate abs' 'torq rate change sqrd' 'torq rate change abs' 'yGRF rate sqrd' 'xGRF rate sqrd' 'positive work' 'absolute work' 'kinetic energy' 'angular momentum'}','FontSize',12 )

            

