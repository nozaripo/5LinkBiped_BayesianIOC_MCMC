
ID_Subj = 10;
ID_Speed= 4;

walk_data = importdata(['data/WBDS' num2str(ID_Subj,'%02.f') 'walkT' num2str(ID_Speed,'%02.f') 'ang.txt']);
l_hip  = walk_data.data(:,13) * pi/180;
l_knee = walk_data.data(:,19) * pi/180;
% l_Ankle = walk_data.data(:,25);



r_hip  = circshift(l_hip,shift0);
r_knee = circshift(l_knee,shift0);

r_hip  = walk_data.data(:,10) * pi/180;
r_knee = walk_data.data(:,16) * pi/180;
r_hip  = circshift(r_hip,shift0);
r_knee = circshift(r_knee,shift0);
% r_ankle = walk_data.data(:,25);

angles = [l_hip , r_hip , l_knee , r_knee];




for limb_side = 1:2
    
    ID_Eval = ID_Eval+1;
    
    q = [l_hip-l_knee , l_hip , zeros(length(l_hip),1) , r_hip , r_hip-r_knee]';
if limb_side == 1
    TO = l_TO;
    HS = l_HS;
else
    TO = r_TO;
    HS = r_HS;
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










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = q(:,TO:HS);

time = (walk_data.data(:,1)-1)*.01;
time = (0:length(q)-1)'*.01;

Time = time;

t = time;


l_hip   = lowpass(l_hip,6,100);
r_hip   = lowpass(r_hip,6,100);
l_knee  = lowpass(l_knee,6,100);
r_knee  = lowpass(r_knee,6,100);

l_knee(l_knee<0)=0;
r_knee(r_knee<0)=0;

% joint kinematics

[q_c , dq_c , ddq_c dddq_c] = Cubic_Bspline(t , q(4,:)');

dq = diff(q')'/(Time(end)/(length(Time)-1));
dq = [dq(:,1) dq];
ddq = diff(dq')'/(Time(end)/(length(Time)-1));
ddq = [ddq(:,1) ddq];
dddq = diff(ddq')'/(Time(end)/(length(Time)-1));
dddq = [dddq(:,1) dddq];


q_f = lowpass(q(4,:)',15,100)';
dq_f = diff(q_f')'/(Time(end)/(length(Time)-1));
dq_f = [dq_f(:,1) dq_f];
ddq_f = diff(dq_f')'/(Time(end)/(length(Time)-1));
ddq_f = [ddq_f(:,1) ddq_f];
dddq_f = diff(ddq_f')'/(Time(end)/(length(Time)-1));
dddq_f = [dddq_f(:,1) dddq_f];



% % % q1 = lowpass(q(1,:)',15,100)';
% % % q2 = lowpass(q(2,:)',15,100)';
% % % q3 = lowpass(q(3,:)',15,100)';
% % % q4 = lowpass(q(4,:)',15,100)';
% % % q = [q1 ; q2 ; q3 ; q4];
% % % dq = diff(q')'/(Time(end)/(length(Time)-1));
% % % dq = [dq(:,1) dq];
% % % ddq = diff(dq')'/(Time(end)/(length(Time)-1));
% % % ddq = [ddq(:,1) ddq];
% % % dddq = diff(ddq')'/(Time(end)/(length(Time)-1));
% % % dddq = [dddq(:,1) dddq];

figure(202)
subplot(1,3,1)
plot([dq(4,:)' dq_f' dq_c])
subplot(1,3,2)
plot([ddq(4,:)' ddq_f' ddq_c])
subplot(1,3,3)
plot([dddq(4,:)' dddq_f' dddq_c])
legend('Direct Derivative' , 'Filtered Derivative', 'Cubic Bspline', 'Filtered Cubic Bspline')
