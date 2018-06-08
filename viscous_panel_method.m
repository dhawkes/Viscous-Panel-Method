close all; clear all;

%% Creating the Airfoil %%
% Begin by defining the airfoil and getting its geometry

% Model input parameters
alpha = deg2rad(1);
Re = 100000;
N = 100;

% Load the Lissaman 7769 airfoil as an example
airfoil = load('LISS7769.mat');
airfoil = airfoil.LISS7769;

% Interpolate the desired number of points and plot
x = zeros(N+1,1); y = zeros(N+1,1);

space = linspace(1,N+1,size(airfoil,1));
for i = 1:N+1
    x(i) = interp1(space, airfoil(:,1), i, 'pchip');
    y(i) = interp1(space, airfoil(:,2), i, 'pchip');
end

plot(x,y,'-o'); axis equal;

%% Defining Parameters %%
% Now let's build the parameter vectors and matricies needed to calculate
% the panel coefficients
x1 = x(1:N); x2 = x(2:N+1);
y1 = y(1:N); y2 = y(2:N+1);

% calculate the angles of each panel, relative to the horizontal
theta = atan2(diff(y), diff(x));
theta_diff = bsxfun(@minus, theta, theta');

% calculate the distance between the center of each panel, and each node
mid = 0.5 * [x1 + x2, y1 + y2];
dist = @(a, b) (b - a) .^ 2;
r = sqrt(bsxfun(dist, mid(:,1), x') + bsxfun(dist, mid(:,2), y'));
r1 = r(:,1:N); r2 = r(:,2:N+1);

% calculate the angle subtended by the center of each panel, by the nodes
% composing every other panel
% this fomula is Moran 4-87
vert = bsxfun(@minus,mid(:,2),y2') .* bsxfun(@minus,mid(:,1),x1')...
    - bsxfun(@minus,mid(:,1),x2') .* bsxfun(@minus,mid(:,2),y1');
horiz = bsxfun(@minus,mid(:,1),x2') .* bsxfun(@minus,mid(:,1),x1')...
    + bsxfun(@minus,mid(:,2),y2') .* bsxfun(@minus,mid(:,2),y1');
B = atan2(vert, horiz);
B(logical(eye(size(B)))) = pi;

%% Calculating Coefficients %%
% Now we want to calculate the coefficients defined by Moran 4-89 and 4-93
% this should correspond to an N+1 by N+1 matrix which we can use to
% calculate the N+1 unknown values

% First we calculate the values from 1 to N using Moran 4-90
A = sin(theta_diff) .* log(r2./r1) + cos(theta_diff) .* B;
A = A / (2*pi);

% Now calculate the N+1 column (except A(N+1,N+1)) using Moran 4-91
Acol = sum(cos(theta_diff) .* log(r2./r1) - sin(theta_diff) .* B, 2);
Acol = Acol / (2*pi);

% Next we get the N+1 row (except A(N+1,N+1)) using Moran 4-94
k = [1 N];
Arow = sum(sin(theta_diff(k,:)) .* B(k,:) - cos(theta_diff(k,:))...
    .* log(r2(k,:)./r1(k,:)), 1);
Arow = Arow / (2*pi);

% Finally, we calculate A(N+1,N+1), completing the matrix
Acorner = sum(sum(sin(theta_diff(k,:)) .* log(r2(k,:)./r1(k,:))...
    + cos(theta_diff(k,:)) .* B(k,:), 1), 2);
Acorner = Acorner / (2*pi);

% We can now combine everything into one matrix
A = [A Acol; [Arow Acorner]];

% Now let's calculate the b vector, which contains the solutions of
% equation 4-89 and 4-93
b = sin(theta - alpha);
b = [b; -cos(theta(1) - alpha) - cos(theta(N) - alpha)];

%% Finding Inviscid Velocity and Pressure %%
% Now we can solve the matrix equation in order to find the sink strength
% and vortex strength, and from that, the velocities and pressure
% coefficients

q = (A\b)';

% Now we find velocities using Moran 4-95
V = cos(theta - alpha) +...
    sum(bsxfun(@times,sin(theta_diff).*B - cos(theta_diff) .*...
    log(r2./r1),q(1:N)/(2*pi)),2) +...
    (q(N+1)/(2*pi)) * sum(sin(theta_diff) .* log(r2./r1) +...
    cos(theta_diff) .* B, 2);

% Finally, we get the pressure coefficients from Moran 4-96
Cp_panel = 1 - V.^2;

% Plot the inviscid pressure distribution
figure; plot(mid(:,1), -mid(:,2), mid(:,1), Cp_panel);
axis equal; set(gca,'Ydir','reverse'); xlim([0 1]); ylim([-1 1]);

%% Find Points of Interest %%
% From the inviscid analysis, we can find points of interest along the
% airfoil, and divide the relevant vectors in to upper and lower pieces

% Find the point of stagnation
[a, stag] = min(abs(V));

% Define relevant vectors for clarity
lCp = zeros(stag,1); lVe = zeros(stag,1); lx = zeros(stag,1);
uCp = zeros(N-stag,1); uVe = zeros(N-stag,1); ux = zeros(N-stag,1);
lN = length(lCp); uN = length(uCp);

% Break up the relevant vectors into lower and upper pieces, and find the
% points of minimum pressure
lCp = Cp_panel(1:stag);
lVe = abs(V(1:stag));
lx = mid(1:stag,1);
[a, lminp] = min(lCp);

uCp = Cp_panel(stag+1:end);
uVe = V(stag+1:end);
ux = mid(stag+1:end,1);
[a, uminp] = min(uCp); % Note that uminp is defined relative to the upper
                       % surface

%% Compute Laminar Boundary Layer %%
% Now we can find the displacement thickness for the laminar portions of
% the airfoil. We also want the theta value at the transition points in
% order to find the initial conditions at the transition points
uDelta = zeros(uN,1); lDelta = zeros(lN,1);

% Define relevant equations. These are based on results of the blassius
% solution, found online.
ReX = @(Ve, x_c) Re .* x_c .* Ve;
delta_lam = @(Ve, x_c) (1.72.*x_c) ./ sqrt(ReX(Ve, x_c));
theta_lam = @(Ve, x_c) (0.665.*x_c) ./ sqrt(ReX(Ve, x_c));

% Starting with the lower surface, where the laminar portion is from the
% minimum pressure point to the stagnation point, we can find delta*
span = lminp:lN;
lDelta(span) = delta_lam(lVe(span), lx(span));

% We also want the theta and H values at the transition point for use as
% initial conditions in the turbulent calculations
lThetai = theta_lam(lVe(lminp), lx(lminp));
lHi = lDelta(lminp) / lThetai;

% Now we repeat the process for the upper surface of the airfoil. Here, the
% laminar portion goes from the stagnation point to the point of minimum
% pressure
span = 1:uminp;
uDelta(span) = delta_lam(uVe(span), ux(span));

uThetai = theta_lam(uVe(uminp), ux(uminp));
uHi = uDelta(uminp) / uThetai;

%% Compute Turbulent Boundary Layer %%
% Now we can compute the turbulent boundary layer displacement thickness
% using ODE45. We have to perform ODE45 backwards on the bottom edge, based
% on where we have the initial values

if lminp ~= 1
    % Let's start with the lower by defining relevant variables. Note that the
    % span is in reverse
    span = lminp:-1:1;
    lx_turb = lx(span);
    lVe_turb = lVe(span);
    ldVe_turb = gradient(lVe_turb, lx_turb);
    lYi = [lThetai, lHi];

    % Now let's put the turbulence function into a form that ODE45 can
    % understand. ODE45 wants to specify its own input and output values
    lturb_fun = @(x, Y) turbulent_approx(x, Y, lx_turb, lVe_turb,...
        ldVe_turb, Re);

    % Then we solve
    [lx_turb, lY] = ode45(lturb_fun, lx_turb, lYi);

    % And then get the values that we are looking for through interpolation
    lDelta(span) = lY(:,2) .* lY(:,1);
end

if uminp ~= uN
    % Now we can do the upper layer
    span = uminp:uN;
    ux_turb = ux(span);
    uVe_turb = uVe(span);
    udVe_turb = gradient(uVe_turb, ux_turb);
    uYi = [uThetai, uHi];

    uturb_fun = @(x, Y) turbulent_approx(x, Y, ux_turb, uVe_turb,...
        udVe_turb, Re);

    [ux_turb, uY] = ode45(uturb_fun, ux_turb, uYi);

    % And then get the delta values
    uDelta(span) = uY(:,2) .* uY(:,1);
end

% Combine to get full vectors
delta = [lDelta; uDelta];
Ve = [lVe; uVe];

% Plot the airfoil with the computed boundary layer
yadjusted = mid(:,2) + delta .* sign(mid(:,2));
figure; hold on;
plot(mid(:,1), mid(:,2), mid(:,1), yadjusted, '--');
POIx = [mid(lminp,1) mid(stag,1) mid(uminp+stag,1)];
POIy = [mid(lminp,2) mid(stag,2) mid(uminp+stag,2)];
scatter(POIx, POIy);
axis equal; xlabel('x/c'); ylabel('y/c');
title('Boundary Layer and Minimum Pressure Points');
legend('Airfoil', 'Boundary Layer', 'Points of Interest');

%% Reapply Panel Method %%
% Now we can alter b using the displacement values that we have found, and
% reapply the panel method to get viscous results

b(1:N) = b(1:N) + Ve(1:N) .* gradient(delta(1:N), mid(:,1))...
    + gradient(Ve(1:N), mid(:,1)) .* delta(1:N);

q = (A\b)';

V = cos(theta - alpha) +...
    sum(bsxfun(@times,sin(theta_diff).*B - cos(theta_diff) .*...
    log(r2./r1),q(1:N)/(2*pi)),2) +...
    (q(N+1)/(2*pi)) * sum(sin(theta_diff) .* log(r2./r1) +...
    cos(theta_diff) .* B, 2);

Cp = 1 - V.^2;

%% Plot Results %%
figure; hold on;
plot(mid(:,1), -mid(:,2), mid(:,1), -yadjusted, '--');
plot(mid(:,1), Cp, mid(:,1), Cp_panel);
axis equal; set(gca,'Ydir','reverse'); xlim([0 1]); ylim([-1 1]);
title('Airfoil Pressure Comparison'); xlabel('x/c');
legend('Body', 'Boundary Layer', 'Viscous Pressure',...
    'Inviscid Pressure', 'location', 'SouthEast');
