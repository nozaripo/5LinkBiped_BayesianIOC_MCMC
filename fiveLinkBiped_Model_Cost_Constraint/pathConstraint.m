function [c, ceq, cGrad, ceqGrad] = pathConstraint(x,u,q0,p)
% [c, ceq, cGrad, ceqGrad] = pathConstraint(x)
%
% This function implements a simple path constraint to keep the knee joint
% of the robot from hyer-extending.
%

q1 = x(1,:);
q2 = x(2,:);
q4 = x(4,:);
q5 = x(5,:);

c = [...
    q1'-q2';    %Stance knee joint limit
    q5'-q4'];   %Swing knee joint limit

c_q4_init = [q4(1)+.1 ; -q2(1)+.1];
c_q4_init = [];

% c = [c;
%      q2'-q1'-pi/15;
%      q4'-q5'-pi/3];

% c_u = max(max(abs(diff(u'))))-50;
c_u = [];

ceq_u_limitcycle = u(:,1)-u(5:-1:1,end);
ceq_u_limitcycle = [];

[P, G] = getPoints(x(1:5,:),p);
c_clearance = -P(10,5:end-5)'+.005;
% c_clearance = [];


% ceq_torso_vel = P(3,1)-P(3,end);

% c_foot_foreaft_init = P(9,1)+.1;
c_foot_foreaft_init = [];


c = [c;c_u;c_clearance; c_q4_init; c_foot_foreaft_init];


ceq = [];
% % % ceq = [x(1:5,1)-q0];%
% ceq = q1(1)-q2(1);

% ceq_treadmill = P(3,end)-P(3,1)-.3;
ceq_treadmill = [];

% ceq_foot_height_init = P(10,1);
ceq_foot_height_init = [];

ceq = [ceq ; ceq_u_limitcycle ; ceq_treadmill; ceq_foot_height_init];

if nargout == 4 %Analytic gradients
    % Gradients with respect to:
    % [t,q1,q2,q3,q4,q5,dq1,dq2,dq3,dq4,dq5,u1,u2,u3,u4,u5] = 1+5+5+5
    nCst = 2;   %stance leg ; swing leg
    nGrad = 16;  %time, angles, rates, torques
    nTime = size(x,2);
    cGrad = zeros(nCst,nGrad,nTime);
    cGrad(1,3,:) = -1;  % cst stance wrt q2
    cGrad(1,2,:) = 1; % cst stance wrt q1
    cGrad(2,5,:) = -1;  % cst swing wrt q4
    cGrad(2,6,:) = 1; % cst swing wrt q5
    
    ceqGrad = [];
end

end