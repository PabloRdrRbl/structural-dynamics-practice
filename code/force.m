function [ F ] = force( t, Ndof, locnod, dof_rem )
% Force applied on the hangar.
% Its time variation is a triangular function af total 
% length equal to half the period of the fourth eigenfrequency
% and has an amplitude of 100 N.

% Half the period of the fourth eigenfrequency
T = 0.5 * (1 / 0.4651);

if (t < 0) || (t > T)
    f = 0;
    
elseif t < T / 2
    f = 100 / (T / 2) * t;

else
    f = 200 - 100 / (T / 2) * t;
    
end
    
% Ndof, number of degrees of freedom after applying boundary conditions
F = zeros(Ndof, 1);

% Force is applied vertically and pointing downwards on Node 9
loc = locnod(9, 3);

F(loc) = - f;

F = F(dof_rem);

end

