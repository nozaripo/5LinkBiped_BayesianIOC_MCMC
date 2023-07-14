Treadmill_Speed = [.4 .7 1 1.2 1.4 1.6];
Treadmill_Speed = [.4 .7 1 1.3 1.6];
Time = [];
Torso = [];
Femur_Stance=[];
Femur_Swing = [];
Tibia_Stance= [];
Tibia_Swing =[];
Time_Demo            = [];
Torso_Demo           = [];
Femur_Stance_Demo    = [];
Femur_Swing_Demo     = [];
Tibia_Stance_Demo    = [];
Tibia_Swing_Demo     = [];

Knee_Right          = [];
Hip_Right           = [];
Knee_Left           = [];
Hip_Left            = [];
Knee_Right_Demo     = [];
Hip_Right_Demo      = [];
Knee_Left_Demo      = [];
Hip_Left_Demo       = [];

Simulation = [];
Subject = [];
Speed = [];
Weight  = [];
Inference=[];
Obs   = [];
RMSE = [];
% error = zeros(2,3);


FolderName = "Results_Forward/2023-04-04";


load('Simulated_Speeds_5Bases.mat')

load('Synthetic_5Weights_5Speeds_10Subjects.mat')


for id_method=1:2
    if id_method==1
        FileExtension = "_Bases";
    else
        FileExtension = "";
    end
    for id_Subject = 1:10
        for id_Weight = 1:5
            fullpath = strcat(FolderName, '/Results_Angles_W', num2str(id_Weight), '_S', num2str(id_Subject), FileExtension, '.mat');
            load(fullpath)
            for i_sol = 1:length(Time_Traj)
                X = X_Traj{1, i_sol};
                T = Time_Traj{1, i_sol};
                
                Mean_Hip_Right = 0;
                Upp_Hip_Right = 0;
                Low_Hip_Right = 0;
                Mean_Hip_Left = 0;
                Upp_Hip_Left = 0;
                Low_Hip_Left = 0;
                Mean_Knee_Right = 0;
                Upp_Knee_Right = 0;
                Low_Knee_Right = 0;
                Mean_Knee_Left = 0;
                Upp_Knee_Left = 0;
                Low_Knee_Left = 0;
                
                for i_speed = 1:5
                    soln_key = strcat('V', num2str(Treadmill_Speed(i_speed)), '_S', num2str(id_Subject), '_W', num2str(id_Weight));
                    Demo = soln_map(soln_key);
                    tInt = Demo(end).grid.time;
                    xInt = Demo(end).grid.state;
%                     tInt   = linspace(t(1),t(end),10*length(t)+1);
%                     xInt      = soln(end).interp.state(tInt);
                    q = xInt(1:5,:);
                    
                    
                    TT = T(1,:,i_speed);
                    XX = X(:,:,i_speed);
                    
                    tmp_Time = linspace(0, TT(end), 11);
                    tmp_X    = interp1(TT', XX', tmp_Time')';
                                      
%                     soln = soln_map(Treadmill_Speed(i_speed));
                    
                    
                    
                    summ = 0;
                    for kk =1:5
                        summ = summ + norm(q(kk,:)' - tmp_X(kk,:)');
                    end
                    error = summ/5;
                    
                    tmp_Time_Stride = [tmp_Time,tmp_Time(2:end)+tmp_Time(end)];
                    tmp_X_Stride = [tmp_X, tmp_X(:,end-1:-1:1)];
                    
                    tInt_Stride = [tInt,tInt(2:end)+tInt(end)];
                    xInt_Stride = [xInt, xInt(:,end-1:-1:1)];
                    
                    Time            = [Time         ; tmp_Time_Stride' ];
                    Torso           = [Torso        ; tmp_X_Stride(3,:)'];
                    Femur_Stance    = [Femur_Stance ; tmp_X_Stride(2,:)'];
                    Femur_Swing     = [Femur_Swing  ; tmp_X_Stride(4,:)'];
                    Tibia_Stance    = [Tibia_Stance ; tmp_X_Stride(1,:)'];
                    Tibia_Swing     = [Tibia_Swing  ; tmp_X_Stride(5,:)'];
                    Time_Demo            = [Time_Demo         ; tInt_Stride'];
                    Torso_Demo           = [Torso_Demo        ; xInt_Stride(3,:)'];                    
                    Femur_Stance_Demo    = [Femur_Stance_Demo ; xInt_Stride(2,:)'];
                    Femur_Swing_Demo     = [Femur_Swing_Demo  ; xInt_Stride(4,:)'];
                    Tibia_Stance_Demo    = [Tibia_Stance_Demo ; xInt_Stride(1,:)'];
                    Tibia_Swing_Demo     = [Tibia_Swing_Demo  ; xInt_Stride(5,:)'];
                    
                    Knee_Right         = [Knee_Right    ; tmp_X_Stride(2,:)' - tmp_X_Stride(1,:)'];
                    Hip_Right          = [Hip_Right     ; tmp_X_Stride(2,:)'];
                    Knee_Left         = [Knee_Left      ; tmp_X_Stride(4,:)' - tmp_X_Stride(5,:)'];
                    Hip_Left          = [Hip_Left       ; tmp_X_Stride(4,:)'];
                    
                    Knee_Right_Demo         = [Knee_Right_Demo  ; xInt_Stride(2,:)' - xInt_Stride(1,:)'];
                    Hip_Right_Demo          = [Hip_Right_Demo   ; xInt_Stride(2,:)'];
                    Knee_Left_Demo         = [Knee_Left_Demo    ; xInt_Stride(4,:)' - xInt_Stride(5,:)'];
                    Hip_Left_Demo          = [Hip_Left_Demo     ; xInt_Stride(4,:)'];
                    
%                     Simulation      = [Simulation   ; i_sim*ones(11,1)];
                    Inference       = [Inference    ; id_method*ones(21,1)];
                    Speed           = [Speed        ; Treadmill_Speed(i_speed)*ones(21,1)];
                    Weight          = [Weight         ; (id_Weight)*ones(21,1)];
                    Subject          = [Subject         ; (id_Subject)*ones(21,1)];
                    RMSE            = [RMSE         ; error*ones(21,1)];

%                     Obs             = [Obs          ; (i_speed)*ones(11,1)];
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
end

save('Trajectories_Bayesian_10Subjects_5Weights.mat', 'Time', 'Torso', 'Femur_Stance', 'Femur_Swing', 'Tibia_Stance', 'Tibia_Swing', 'Time_Demo', 'Torso_Demo', 'Femur_Stance_Demo', 'Femur_Swing_Demo', 'Tibia_Stance_Demo', 'Tibia_Swing_Demo', 'Hip_Right', 'Hip_Right_Demo','Hip_Left','Hip_Left_Demo','Knee_Right','Knee_Right_Demo', 'Knee_Left','Knee_Left_Demo', 'Speed', 'Weight', 'Inference', 'Subject', 'RMSE')



Error = zeros(length(Time), 1);
for i=1:2
    for i_sim = 1:3
        Error(Cost==(i-1) & Simulation==i_sim) = round(error(i, i_sim)/51/6 , 3);
    end
end
   
% save('Trajectories_Bayesian.mat', 'Time', 'Torso', 'Femur_Stance', 'Femur_Swing', 'Tibia_Stance', 'Tibia_Swing', 'Simulation', 'Speed')
save('Trajectories_Bayesian_Bases.mat', 'Time', 'Torso', 'Femur_Stance', 'Femur_Swing', 'Tibia_Stance', 'Tibia_Swing', 'Simulation', 'Speed', 'Cost', 'Obs', 'Error')


            