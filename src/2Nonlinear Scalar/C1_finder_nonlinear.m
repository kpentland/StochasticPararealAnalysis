function [FU] = C1_finder_nonlinear(DT,u)
% This function estimates the a minimum value for the constant C1 by
% minimising the function given by the 'C1' function below. 

% function handles for the fine and coarse sovlers given an initial value u
F = @(u)(sqrt(2)*sinh(DT + asinh(u/sqrt(2))));
G = @(u)(u + DT*sqrt(u^2 + 2));

[U,FU] = fminsearch(@(u)C1(DT,F,G,u),u,optimset('TolFun',1e-6,'TolX',1e-6,'MaxIter',200,'Display','off'));


end




function value = C1(DT,F,G,u)
% we wish to minimise the following function with respect to two variables
% u(1) and u(2), this corresponds to u and v in the manuscript (see
% equation (4.5) in the manuscript.

value = (abs( (F(u(1)) - G(u(1))) - (F(u(2)) - G(u(2))) )/( (DT*DT*abs(u(1)-u(2)))));

end

