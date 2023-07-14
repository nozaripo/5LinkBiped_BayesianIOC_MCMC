Study = [];
Ratio = [];
Value = [];
Variable=[];
load('Spatiotemporal_GoodEnough.mat')

sp = 2.75;
% sp = 3;

% Optimal solutions
Ratio = [ Ratio ; 3 ] ;
Value = [Value; .126]; %015
Study = [Study ; 0];
Variable = [Variable; 1];
Ratio = [ Ratio ; 3 ] ;
Value = [Value; .202]; %150
Study = [Study ; 0];
Variable = [Variable; 2];

indices_good_enough = (sp_r_c==1 | sp_r_c==sp);

Ratio = [ Ratio ; sp_r_c(indices_good_enough)];
Ratio = [ Ratio ; sp_r_c(indices_good_enough)];
Value = [ Value ; SLA(indices_good_enough)];
Value = [ Value ; STA(indices_good_enough)];
Study = [ Study; ones(2*sum(indices_good_enough), 1) ];
Variable = [Variable; ones(sum(indices_good_enough), 1)];
Variable = [Variable; 2*ones(sum(indices_good_enough), 1)];


SLA_Sanchez_mean = 0.028;
SLA_Sanchez_std = (0.028 - 0.0018)/1.96;

SLA_Sanchez = SLA_Sanchez_std.*randn(20,1) + SLA_Sanchez_mean;
Ratio = [ Ratio ; 3*ones(20,1) ] ;
Value = [Value; SLA_Sanchez];
Study = [Study ; 2*ones(20,1)];
Variable = [Variable; ones(20, 1)];

STA_Stenum_mean = 0.178;
STA_Stenum_std = .028;
STA_Stenum = STA_Stenum_std.*randn(10,1) + STA_Stenum_mean;
Ratio = [ Ratio ; 3*ones(10,1) ] ;
Value = [Value; STA_Stenum];
Study = [Study ; 3*ones(10,1)];
Variable = [Variable; 2*ones(10, 1)];


SLA_Stenum_mean = -0.034;
SLA_Stenum_std = .025;
SLA_Stenum = SLA_Stenum_std.*randn(10,1) + SLA_Stenum_mean;
Ratio = [ Ratio ; 3*ones(10,1) ] ;
Value = [Value; SLA_Stenum];
Study = [Study ; 3*ones(10,1)];
Variable = [Variable; ones(10, 1)];

% Price Study
Ratio = [ Ratio ; 3 ] ;
Value = [Value; .42];
Study = [Study ; 4];
Variable = [Variable; 1];
Ratio = [ Ratio ; 3 ] ;
Value = [Value; .34];
Study = [Study ; 4];
Variable = [Variable; 2];


% Seethapathi Study
Ratio = [ Ratio ; 3 ] ;
Value = [Value; .4];
Study = [Study ; 5];
Variable = [Variable; 1];
Ratio = [ Ratio ; 3 ] ;
Value = [Value; .25];
Study = [Study ; 5];
Variable = [Variable; 2];

Ratio(Ratio==2.75) = 3;

save('Spatiotemporal_Comparison_ForR.mat', 'Ratio', 'Value', 'Study', 'Variable')




