function drawRobot(t,q,p, V_tr)
% drawRobot(q,p)
%
% This function draws the robot with configuration q and parameters p
%
% INPUTS:
%   q = [5, 1] = column vector of a single robot configuration
%   p = parameter struct
%


% Compute the points that will be used for plotting
[P, G] = getPoints(t,q,p, V_tr);

offset = - P (3,1);
G (1:2:9,:) = G (1:2:9,:) + offset;
P (1:2:9,:) = P (1:2:9,:) + offset;


x = [-V_tr*t+offset; P(1:2:end,1)];
y = [0; P(2:2:end,1)];

P1 = P(1:2,1);
P2 = P(3:4,1);
P3 = P(5:6,1);
P4 = P(7:8,1);
P5 = P(9:10,1);

G1 = G(1:2,1);
G2 = G(3:4,1);
G3 = G(5:6,1);
G4 = G(7:8,1);
G5 = G(9:10,1);

% Heuristics:
L = (p.l1 + p.l2);  % Maximum extended leg length
xBnd = L*[-1.2,1.2];
yBnd = [-0.2*L, L + p.l3];

% Colors:
colorGround = [100,80,30]/255;
colorGround2 = [150,150,150]/255;
colorStance = [200,60,60]/255;
colorSwing = [60,60,200]/255;
colorTorso = [160, 80, 160]/255;

groundy = [-.06, -.06];

% Set up the figure
hold off;

% Plot the ground:
co = 0;
lim = -3;
while lim <= 3
    if mod(co,2) == 0
        plot([lim lim+.2]-V_tr*t + offset  ,groundy,'LineWidth',10,'Color',colorGround);
    hold on;
    else
        plot([lim lim+.2]-V_tr*t+ offset ,groundy,'LineWidth',10,'Color',colorGround2);
    end
    lim = lim + .2;
    co = co + 1;
end
% plot([-.5 0]-t,groundy,'LineWidth',11,'Color',colorGround2);
% plot([0 .5]-t,groundy,'LineWidth',11,'Color',colorGround);
% plot([.5 1]-t,groundy,'LineWidth',11,'Color',colorGround2);
% plot([1 1.5]-t,groundy,'LineWidth',11,'Color',colorGround);
plot(0, 0,'k|','MarkerSize',30);
% Plot the links:
% plot(x(1:3),y(1:3),'LineWidth',4,'Color',colorStance);
% plot(x(3:4),y(3:4),'LineWidth',4,'Color',colorTorso);
% plot(x([3,5,6]),y([3,5,6]),'LineWidth',4,'Color',colorSwing);

plot(x(1:2),y(1:2),'LineWidth',7,'Color',[1 0 1]);
plot(x(2:3),y(2:3),'LineWidth',7,'Color',[1 0 0]);
plot(x(3:4),y(3:4),'LineWidth',8,'Color',[.3 .3 .3]);
plot(x([3,5]),y([3,5]),'LineWidth',7,'Color',[0 0 1]);
plot(x(5:6),y(5:6),'LineWidth',7,'Color',[0 1 1]);
% ax = gca;
% ax.ColorOrder = [1 0 1; 1 0 0; 0 0 0; 0 0 1; 0 1 1];

% Plot the joints:
plot(-V_tr*t + offset, 0,'k.','MarkerSize',30);
plot(P1(1), P1(2),'k.','MarkerSize',30);
plot(P2(1), P2(2),'k.','MarkerSize',30);
plot(P3(1), P3(2),'k.','MarkerSize',30);
plot(P4(1), P4(2),'k.','MarkerSize',30);
plot(P5(1), P5(2),'k.','MarkerSize',30);

% Plot the CoM:
plot(G1(1), G1(2),'ko','MarkerSize',5,'LineWidth',2);
plot(G2(1), G2(2),'ko','MarkerSize',5,'LineWidth',2);
plot(G3(1), G3(2),'ko','MarkerSize',5,'LineWidth',2);
plot(G4(1), G4(2),'ko','MarkerSize',5,'LineWidth',2);
plot(G5(1), G5(2),'ko','MarkerSize',5,'LineWidth',2);
text(0,2, ['Treadmill Speed = ', num2str(V_tr, '%.1f'), ' m/s'], 'FontSize', 12, 'HorizontalAlignment', 'center')
text(0,1.97, '_______________', 'FontSize', 12, 'HorizontalAlignment', 'center')
text(0,1.8, ['Time = ', num2str(t, '%.3f'), ' s'], 'FontSize', 13, 'HorizontalAlignment', 'center')
% Format the axis:
axis([[-1 1],yBnd]); axis equal; axis off;

ylim([-.06 2])
xlim([-1 1])
end