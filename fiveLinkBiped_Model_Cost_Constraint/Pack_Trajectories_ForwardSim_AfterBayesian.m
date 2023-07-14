Treadmill_Speed = [.4 .7 1 1.2 1.4 1.6];
Time = [];
Torso = [];
Femur_Stance=[];
Femur_Swing = [];
Tibia_Stance= [];
Tibia_Swing =[];
Simulation = [];
Speed = [];
Cost  = [];
Obs   = [];
error = zeros(2,3);


load('Simulated_Speeds_5Bases.mat')

for i=1:2
    if i==1
        load('LimbAngles_BayesianSolutions.mat')
    else
        load('LimbAngles_BayesianSolutions_Bases.mat')
    end
    for i_sim = 1:3
        for i_sol = 151:201
            X = X_Traj{i_sol, i_sim};
            T = Time_Traj{i_sol, i_sim};
            for i_speed = 1:6
                
                soln = soln_map(Treadmill_Speed(i_speed));
                t = soln(end).grid.time;
                tInt   = linspace(t(1),t(end),10*length(t)+1);
                xInt      = soln(end).interp.state(tInt);
                q = xInt(1:5,:);
                
                summ = 0;
                for kk =1:5
                    summ = summ + norm(q(kk,:)' - X(kk,:,i_speed)');
                end
                error(i, i_sim) = error(i, i_sim) + summ/5;
                
                tmp_Time = linspace(T(1,1,i_speed),T(1,end,i_speed), 10);
                tmp_X    = interp1(T(1,:,i_speed)',X(:,:,i_speed)', tmp_Time')';
                
                tmp_Time_Stride = [tmp_Time,tmp_Time(2:end)+tmp_Time(end)];
                tmp_X_Stride = [tmp_X, tmp_X(:,end-1:-1:1)];
                
                Time            = [Time         ; tmp_Time_Stride'];
                Torso           = [Torso        ; tmp_X_Stride(3,:)'];
                Femur_Stance    = [Femur_Stance ; tmp_X_Stride(2,:)'];
                Femur_Swing     = [Femur_Swing  ; tmp_X_Stride(4,:)'];
                Tibia_Stance    = [Tibia_Stance ; tmp_X_Stride(1,:)'];
                Tibia_Swing     = [Tibia_Swing  ; tmp_X_Stride(5,:)'];
                Simulation      = [Simulation   ; i_sim*ones(19,1)];
                Speed           = [Speed        ; Treadmill_Speed(i_speed)*ones(19,1)];
                Cost            = [Cost         ; (i-1)*ones(19,1)];
                Obs             = [Obs          ; (10*(i_sol-151)+i_speed)*ones(19,1)];
%                 Time            = [Time         ; T(1,:,i_speed)'];
%                 Torso           = [Torso        ; X(3,:,i_speed)'];
%                 Femur_Stance    = [Femur_Stance ; X(2,:,i_speed)'];
%                 Femur_Swing     = [Femur_Swing  ; X(4,:,i_speed)'];
%                 Tibia_Stance    = [Tibia_Stance ; X(1,:,i_speed)'];
%                 Tibia_Swing     = [Tibia_Swing  ; X(5,:,i_speed)'];
%                 Simulation      = [Simulation   ; i_sim*ones(111,1)];
%                 Speed           = [Speed        ; Treadmill_Speed(i_speed)*ones(111,1)];
%                 Cost            = [Cost         ; (i-1)*ones(111,1)];
%                 figure(1)
%                 hold on
%                 plot(tmp_Time, tmp_X(1,:))
            end
        end
    end
end
Error = zeros(length(Time), 1);
for i=1:2
    for i_sim = 1:3
        Error(Cost==(i-1) & Simulation==i_sim) = round(error(i, i_sim)/51/6 , 3);
    end
end
   
% save('Trajectories_Bayesian.mat', 'Time', 'Torso', 'Femur_Stance', 'Femur_Swing', 'Tibia_Stance', 'Tibia_Swing', 'Simulation', 'Speed')
save('Trajectories_Bayesian_Bases.mat', 'Time', 'Torso', 'Femur_Stance', 'Femur_Swing', 'Tibia_Stance', 'Tibia_Swing', 'Simulation', 'Speed', 'Cost', 'Obs', 'Error')


            