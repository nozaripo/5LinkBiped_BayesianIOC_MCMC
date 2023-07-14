function [c, ceq, cGrad, ceqGrad] = stepConstraint_behavioral(x0,xF,p,q0)
% [c, ceq, cGrad, ceqGrad] = stepConstraint(x0,xF,p)
%
% This function applies the non-linear boundary constraint to ensure that
% there gait is periodic, with the correct heel-strike map and left-right
% mapping.
%

if nargout == 2 %Numerical gradients
    
    % Gait must be periodic
    ceq1 = cst_heelStrike_behavioral(x0,xF,p,q0);
    ceq1 = [];
    
    % Initial and Final Pose in the vicinity of the start and end from
    % behavioral
%     c_IF = cst_heelStrike_ineq_behavioral(x0,xF,p);
    c_IF = [];
    % Foot collision (toe-off and heel-strike) velocity
    c = cst_footVel(x0,xF,p);
    c = [];
    
    % Step length and height
    ceq2 = cst_stepLength(xF,p); 
    ceq2 = [];
    
%     c2 = -cst_stepLength(xF,p); 
    c2=[];
    
    % Foot Clearance
%     c3=

    % initial configuration realization; find out x(4:5,1) from x(1:2,1)
    ceq3 = init_swing_realization(x0,p);
    ceq3 = [];
    
    
    % Pack up equality constraints:
    ceq = [ceq1;ceq2;ceq3];
    
    % Pack up inequality constraints:
    c = [c ; c2; c_IF];  
    
else %Analytic gradients
    
    % Gait must be periodic
    [ceq1, ceqGrad1] = cst_heelStrike(x0,xF,p);
    
    % Foot collision (toe-off and heel-strike) velocity
    [c, cGrad] = cst_footVel(x0,xF,p);
    
    % Step length and height
    [ceq2, ceqGrad2] = cst_stepLength(xF,p);
    
        % Pack up equality constraints:
    ceq = [ceq1;ceq2];
    ceqGrad = [ceqGrad1;ceqGrad2];
    
end

end