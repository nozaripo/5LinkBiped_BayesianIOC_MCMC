load('Simulated_Speeds_16Costs.mat')
V = [.4, .7, 1, 1.2, 1.4, 1.6];
Time = [];
Torso = [];
Femur_Stance = [];
Femur_Swing  = [];
Tibia_Stance = [];
Tibia_Swing  = [];
Speed     = [];
Cost = [];
for i=1:6
    soln = soln_map(V(i));
    t = soln(end).grid.time;
    tInt   = linspace(t(1),t(end),10*length(t)+1);
    Time = [Time; tInt'];
    xInt      = soln(end).interp.state(tInt);    
    q = xInt(1:5,:);
    Torso = [Torso; q(3,:)'];
    figure(22)
    plot(tInt', Torso')
    Femur_Stance = [Femur_Stance; q(2,:)'];
    Tibia_Stance = [Tibia_Stance; q(1,:)'];
    Femur_Swing  = [Femur_Swing; q(4,:)'];
    Tibia_Swing  = [Tibia_Swing; q(5,:)'];
    Speed        = [Speed; V(i)*ones(length(tInt),1)];
%     Cost = [Cost; repmat({"Individual Costs"},length(tInt),1)];
    Cost = [Cost; 0*ones(length(tInt),1)];

end
% save('Synthetic_Timeseries_Costs_R.mat', 'Time', 'Torso', 'Femur_Stance', 'Tibia_Stance', 'Femur_Swing', 'Tibia_Swing')

load('Simulated_Speeds_5Bases.mat')
for i=1:6
    soln = soln_map(V(i));
    t = soln(end).grid.time;
    tInt   = linspace(t(1),t(end),10*length(t)+1);
    Time = [Time; tInt'];
    xInt      = soln(end).interp.state(tInt);    
    q = xInt(1:5,:);
    Torso = [Torso; q(3,:)'];
    Femur_Stance = [Femur_Stance; q(2,:)'];
    Tibia_Stance = [Tibia_Stance; q(1,:)'];
    Femur_Swing  = [Femur_Swing; q(4,:)'];
    Tibia_Swing  = [Tibia_Swing; q(5,:)'];
    Speed        = [Speed; V(i)*ones(length(tInt),1)];
%     Cost = [Cost; repmat({"Costs Bases"},length(tInt),1)];
    Cost = [Cost; 1*ones(length(tInt),1)];
end
save('Synthetic_Timeseries_Bases_R.mat', 'Time', 'Torso', 'Femur_Stance', 'Tibia_Stance', 'Femur_Swing', 'Tibia_Swing', 'Speed', 'Cost')