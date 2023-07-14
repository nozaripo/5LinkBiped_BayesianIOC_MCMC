soln = soln_map('V1_S1_W1');
T = soln.grid.time;

time = linspace(0,T(end),100);
state= soln.interp.state(time);
options = optimoptions(@fmincon, 'Display', 'iter', ...
                        'MaxIterations',5000, 'MaxFunctionEvaluations', 1e8)

f_obj = @(x)RMS_calc(x,state);
nonlcon = @(x)Time_Con(x);

x0 = rand(size(state)+[1 0]);

exitflag = -2;
while exitflag == -2 || exitflag == 0

[x,fval,exitflag,output] = fmincon(f_obj,x0,[],[],[],[],-15,15,nonlcon, options);

fval
exitflag

end

x


function error = RMS_calc(x, state)

dx = diff(x(1:5,:)')'./ diff(x(11,:)')';
% (x(11,end)/100);
dx = [dx dx(:,end)];

xs = [x(1:5,:) ; dx];


error = sum(rms(xs-state, 2));

end

function [c, ceq] = Time_Con(x)
c = [];
ceq = [x(11, :) - linspace(0,x(11, end), 100)]';
end