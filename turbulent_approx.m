% This is the function that will be passed into ODE45. It is derived from
% Moran, Chapter 7.  The unknowns in question are Y = [Theta, H], and the 
% function returns F = [dTheta, dH].
function [F] = turbulent_approx(x, Y, xvals, Ve, dVe, Re)
    F = zeros(2,1);

    % Figure out which Ve and dVe values to use, based on the iteration
    i = find(xvals==x);
    if isempty(i)
        % Interpolate the value if x is not in xvals
        Vec = interp1(xvals, Ve, x);
        dVec = interp1(xvals, dVe, x);
    else
        Vec = Ve(i);
        dVec = dVe(i);
    end
    
    % Compute Re theta and the skin friction coefficient using Moran 7-38
    % and Moran 7-65
    ReT = Vec * Y(1) * Re;
    Cf = (0.246 * 10 ^ (-0.678 * Y(2))) * (ReT ^ -0.268);
    
    % Compute and return dTheta and dH using Moran 7-34 and a reworking of
    % Moran 7-63 and Moran 7-64
    F(1) = (0.5 * Cf) - (Y(1) / Vec) * (2 + Y(2)) * dVec;
    if Y(2) <= 1.6
        H1 = 3.3 + 0.8234 * ((Y(2) - 1.1) ^ -1.287);
        C = -1.05972 / ((Y(2) - 1.1) ^ 2.287);
    else
        H1 = 3.3 + 1.5501 * ((Y(2) - 0.6778) ^ -3.064);
        C = -4.74951 / ((Y(2) - 0.6778) ^ 4.064);
    end
    F(2) = (1/C) * (((0.0306 * ((H1-3) ^ -0.6169)) / Y(1))...
        - ((dVec * H1) / Vec) - ((F(1) * H1) / Y(1)));
end