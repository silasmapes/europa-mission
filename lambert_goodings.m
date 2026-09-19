%% ======= Lambert's Problem - Gooding's Method ======= %%%
% -------------------------------------------------------------------------
% Reference: Gooding, R.H. (1990) "A procedure for the solution of 
% Lambert's orbital boundary-value problem", Celestial Mechanics 48, 145-146.
% -------------------------------------------------------------------------
% Given two position vectors r1, r2 with a flight time dt, this script 
% finds the orbit connecting them by returning the velocity vectors v1, v2
% at each position. Both positions must be in an intertial frame with
% graviational center C.
% -------------------------------------------------------------------------



%% ----- Stage 0: Initial Parameters ----- %%
% -------------------------------------------------------------------------
% The position vectors r1, r2 and the time dt are given as the initial
% values for the script. The user can also change the number of
% iterations used for Halley's Method. This script can be used on it's own,
% or can be called as a function in other scripts. Comment out necessary
% lines in Stage 0 for desired functionality.
% -------------------------------------------------------------------------

%clear; clc; close all

% Gravitational Parameter [km3/s3]
%mu = 398600.4418;

% Position 1,2 [km; row vectors] and flight time [s]
%r1 = [5000 10000 2100];
%r2 = [-14600 2500 7000];
%dt = 20000;

% Direction of motion [bool.]
% True if theta < pi
% False if theta > pi
%smallAngle = true;

% Halley iteration tuning [int]
%nIter = 5;

% Completed revolutions [int.]
%M = 1;

% Function used to call solving method from another script
function [v1, v2] = lambert_goodings(r1,r2,dt,mu,M,smallAngle,nIter)



%% ----- Stage 1: Geometric Calculations ----- %%
% -------------------------------------------------------------------------
% The triangle connecting points C, P1, and P2 is used to determine the 
% three geometric quantities:
%   1. Reduced angle (angle between r1 and r2 without extra revolutions)
%   2. Semi-Perimeter (1/2 the triangle's perimeter)
%   3. Lambert-invariant similarity class (a number that describes 
%      the triangle's shape).
% -------------------------------------------------------------------------

% Scalar distance of r1 and r2
r1m = norm(r1);
r2m = norm(r2);

% Transfer Angle between r1 and r2
if smallAngle == true
    theta = acos(dot(r1,r2)/(r1m*r2m)) + M*2*pi;
elseif smallAngle == false
    theta = 2*pi - acos(dot(r1,r2)/(r1m*r2m)) + M*2*pi;
end
    
% Reduced Angle (Transfer Angle without revolutions)
phi = theta - 2*pi*M;

% Chord Length (scalar distance) between r1 and r2
c = norm(r2-r1);

% Semi-Perimeter (1/2 of triangle perimeter)
s = (r1m + r2m + c)/2;

% L-similarity (represents type of triangle)
q = (sqrt(r1m*r2m)/s)*cos(phi/2);

% Equal to 1 - q^2; solved differently for accuracy
Q = c/s;



%% ----- Stage 2: Non-Dimensional Time ----- %%
% -------------------------------------------------------------------------
% Time is non-dimensionalized so that it's invariant under uniform scaling 
% of mu and legnths. Therefore, the solving process is identical for any
% orbit. The non-dimensional time is solved at x=0 as a anchor value for
% Halley's method of iterative solving.
% -------------------------------------------------------------------------

% Target Non-Dimensional Time (ND Time)
% Can be solved from triangle geometry
Tf = sqrt(8*mu/s^3)*dt;

% Iteration Anchor Value; ND Time at x=0
% Solved with simplifed formula at x=0
T0 = 2*(atan2(sqrt(Q),q) + q*sqrt(Q));

% Difference between Target ND Time and Anchor ND Time
Td = Tf - T0;



%% ----- Stage 3 (M = 0): Initial x Estimate ----- %%
% -------------------------------------------------------------------------
% Finding an initial x estimate is the first step for the no-revolution
% case. Two paths are used to determine the initial x estimate based on the
% Tf and T0. Some patching is needed when Td > 0 to avoid the x = -1 
% asymptote.
% -------------------------------------------------------------------------

if M == 0

    % Creates matrix for the universal variable
    x = zeros(1,(nIter+1));
    
    % Finds the universal variable when Tf > T0 and x < 0 (High Energy case)
    if Td > 0
    
        x01 = -Td/(Td + 4);
        W = x01 + 1.7*sqrt(2 - phi/pi);
        
        if W > 0
            x03 = x01;
        elseif W < 0
            x02 = -sqrt(Td/(Tf + 0.5*T0));
            w = (-W)^(1/16);
            x03 = (1 - w)*x01 + w*x02;
        end
        
        W2 = 4/(4 + Td);
        lambda = 1 + x03*(0.5*W2 - 0.03*x03*sqrt(W2));
        x(1) = lambda*x03;
    
    % Finds the universal variable when Tf < t0 and x > 0 (Low Energy case)
    elseif Td < 0
        
        x(1) = T0*(T0 - Tf)/(4*Tf);
    
    end
end



%% ----- Stage 4: Find Minimum Values (M > 0) ----- %%
% -------------------------------------------------------------------------
% Finding the minimum x, T, and d2T value is the first step for the 
% multi-revolution case. This is accomplished by finding where dT = 0. 
% Halley's method is used to find the zero. This case can be represented as
% a parabola, so it's necessary to first solve for the lowest point on the 
% parabola to find the real x solutions.
% -------------------------------------------------------------------------

if M > 0

    xM = zeros(1,(nIter+1));

    xMpi = 4/(3*pi*(2*M + 1));
    phir = phi/(2*pi);

    if phir < 1/2
        xM(1) = xMpi*(2*phir)^(1/8);
    elseif phir >= 1/2
        xM(1) = xMpi*(2 - (2 - 2*phir)^(1/8));
    end

    for j = 1:nIter
        [~,dT,d2T,d3T] = highEnergy(xM(j),q,Q,M);
        dxM = -dT*d2T/(d2T^2 - dT*d3T/2);
        xM(j+1) = xM(j) + dxM;
    end

    xMf = xM(end);
    [TM,~,d2TM] = highEnergy(xMf,q,Q,M);

end



%% ----- Stage 5: Halley Iteration ----- %%
% -------------------------------------------------------------------------
% Halley's Method is an improved version of Newton's Method that uses cubic
% convergence rather than square convergence. Unlike Newton's Method, it 
% requires the second derivative of the function, but 13 digits of accuracy
% can be achieved in just three iterations. This stage contains both the
% no-revolution and multi-revolution cases.
% -------------------------------------------------------------------------

% Solves for x when there has been no complete revolutions
if M == 0
    for l=1:nIter
            
            if Td > 0
                [T,dT,d2T,~] = highEnergy(x(l),q,Q,0);
            elseif Td < 0
                [T,dT,d2T] = lowEnergy(x(l),q,Q);
            end
    
            e = Tf - T;
            dx = e*dT/(dT^2 + e*d2T/2);
            x(l+1) = x(l) + dx;
    end
    
    xf = x(end);
    nSol = 1;

% Solves for x when there has been a complete revolution
elseif M > 0
    dTM = Tf - TM;

    % One Solution when Tf is TM
    if abs(dTM) <= 1e-12 * abs(TM)
        nSol = 1;
        xf = xMf;

    % No Solutions
    elseif dTM < 0
        nSol = 0;
        xf = [];

    % Two Solutions in all other cases
    else
        x(1,:) = [multiStarter1(xMf,dTM,d2TM,M,phir), multiStarter2(xMf,TM,d2TM,T0,Tf,M,phir)];

        T   = zeros(nIter,2);
        dT  = zeros(nIter,2);
        d2T = zeros(nIter,2);

        for j=1:nIter
            [T(j,:),dT(j,:),d2T(j,:),~] = highEnergy(x(j,:),q,Q,M);

            e = Tf - T(j,:);
            dx = e.*dT(j,:)./(dT(j,:).^2 + e.*d2T(j,:)/2);
            x((j+1),:) = x(j,:) + dx;
        end

        nSol = 2;
        xf = x(end,:);
    end
end


%% ----- Stage 6: Velocity Calculation ----- %%
% -------------------------------------------------------------------------
% Using the x value determined from the Halley Iteration, the radial and
% tangential components of velocity can be found using three equations
% derivived by Gooding. (The tangential equations are the same, but the 
% radial equations have different signs). The velocity vectors are 
% determined by creating a perifocal frame and using the unit position 
% vectors to convert to the inertial frame.
% -------------------------------------------------------------------------

% Gooding's auxiliary quantities (geometry only, independent of x)
gamma = sqrt(mu*s/2);
rho = (r1m - r2m)/c;
sigma = 2*sqrt(r1m*r2m/c^2)*sin(0.5*phi);

% Defines radial direction as the unit position vector (ur)
ur1 = r1/r1m;
ur2 = r2/r2m;

% Defines tangential direction as the unit angular momentum vector (uh)
% cross the unit position vector (ur)
if smallAngle == true
    uh = cross(r1,r2)/norm(cross(r1,r2));
elseif smallAngle == false
    uh = -cross(r1,r2)/norm(cross(r1,r2));
end
ut1 = cross(uh,ur1);
ut2 = cross(uh,ur2);

% Compute velocity vectors for each solution (stored as columns)
v1 = zeros(3,nSol);
v2 = zeros(3,nSol);

for k = 1:nSol
    zk   = sqrt(Q + q^2*xf(k)^2);
    vr1k = gamma*((q*zk - xf(k)) - rho*(q*zk + xf(k)))/r1m;
    vt1k = gamma*sigma*(zk + q*xf(k))/r1m;
    vr2k = -gamma*((q*zk - xf(k)) + rho*(q*zk + xf(k)))/r2m;
    vt2k = gamma*sigma*(zk + q*xf(k))/r2m;
    v1(:,k) = vr1k*ur1' + vt1k*ut1';
    v2(:,k) = vr2k*ur2' + vt2k*ut2';
end



%% ----- Results Display ----- %%
% -------------------------------------------------------------------------
% Displays final velocity vectors along with important variables used
% throughout the process. Huge thanks to Claude for this section :)
% -------------------------------------------------------------------------

%fprintf('\nLambert solution (Gooding''s method, M = %d):\n', M);
%fprintf('  q  = %+.12f       1 - q^2 = %.12f\n', q, Q);
%fprintf('  Tf = %.10f       T0 = %.10f\n', Tf, T0);

if M > 0
    %fprintf('  xM = %+.12f       TM = %.10f\n', xMf, TM);
end

if nSol == 0
    %fprintf('\n  No solution exists for these inputs.\n');
    if M > 0
        %fprintf('  (T_target = %.6f is below T_M = %.6f)\n', Tf, TM);
    end
else
    %fprintf('\n  %d solution(s) found:\n', nSol);
    for k = 1:nSol
        xk  = xf(k);
        v1k = v1(:,k);
        v2k = v2(:,k);

        if M == 0
            if xk > 0
                label = 'low energy  (x > 0, short-time branch)';
            else
                label = 'high energy (x < 0, long-time branch)';
            end
        else
            if abs(xk - xMf) < 1e-9
                label = 'minimum-time orbit (x = xM)';
            elseif xk > xMf
                label = 'fast / smaller a (x > xM)';
            else
                label = 'slow / larger a  (x < xM)';
            end
        end

        %fprintf('\n  Solution %d - %s\n', k, label);
        %fprintf('    x   = %+.12f\n', xk);
        %fprintf('    v1  = [%+10.6f, %+10.6f, %+10.6f] km/s    |v1| = %8.6f\n', v1k, norm(v1k));
        %fprintf('    v2  = [%+10.6f, %+10.6f, %+10.6f] km/s    |v2| = %8.6f\n', v2k, norm(v2k));
    end
end



%% ----- Function: High Energy Case ----- %%
% -------------------------------------------------------------------------
% Solves T and it's derivatives for the High Energy Scenario (Tf >  T0).
% All solutions result in an elliptic orbit. This function is used by both
% the no-revolution and multi-revolution cases.
% -------------------------------------------------------------------------

function [T,dT,d2T,d3T] = highEnergy(x,q,Q,M)

    % Quantities used in the T equations
    % Identical for High and Low Energy
    U = 1 - x.^2;
    Y = sqrt(abs(U));
    Z = sqrt(Q + q^2*x.^2);

    % Compute scalar coefficients for each element of x independently
    for m = 1:length(x)
        if q*x(m) <= 0
            A(m) = Z(m) - q*x(m);
            B(m) = q*Z(m) - x(m);
        else
            A(m) = Q/(Z(m) + q*x(m));
            B(m) = Q*(q^2*U(m) - x(m)^2)/(q*Z(m) + x(m));
        end

        if q*x(m)*U(m) >= 0
            G(m) = x(m)*Z(m) + q*U(m);
        else
            G(m) = (x(m)^2 - q^2*U(m))/(x(m)*Z(m) - q*U(m));
        end

        F(m) = A(m)*Y(m);
    end

    % Compute T; M*pi term vanishes for M = 0
    T = 2.*((M*pi + atan2(F,G))./Y + B)./U;

    % Derivatives (element-wise; same formula for all M)
    dT  = (3.*x.*T  - 4.*(A + q.*x.*Q)./Z)./U;
    d2T = (3.*T  + 5.*x.*dT  + 4.*q^5.*Q./Z.^5)./U;
    d3T = (8.*dT + 7.*x.*d2T - 12.*q^5.*x.*Q./Z.^5)./U;
end



%% ----- Function: Low Energy Case ----- %%
% -------------------------------------------------------------------------
% Solves T and it's derivatives for the Low Energy Scenario (T0 >  Tf).
% The solution can be an elliptic, parabolic, or hyperbolic orbit, and a
% positive x value. Due to poor accuracy near x = 1 a power series is used
% as a more accuracte esitmation when ~0.84 < x < ~1.18. The Power Series
% Function is used to create the power series form of the T equations. This
% function is only used by the no-revolution case.
% -------------------------------------------------------------------------

function [T,dT,d2T] = lowEnergy(x,q,Q)
    
    % Coefficients used in the T equations
    % Identical for High and Low Energy
    U = 1 - x^2;
    Y = sqrt(abs(U));
    Z = sqrt(Q + q^2*x^2);

    if q*x <= 0
        A = Z - q*x;
        B = q*Z - x;
    elseif q*x > 0
        A = Q/(Z + q*x);
        B = Q*(q^2*U - x^2)/(q*Z + x);
    end

    if q*x*U >= 0
        G = x*Z + q*U;
    elseif q*x*U < 0
        G = (x^2 - q^2*U)/(x*Z - q*U);
    end
    
    F = A*Y;

    % Orbit is an ellipse and x << 1
    if x < sqrt(0.7)
        T = 2*(atan2(F,G)/Y + B)/U;
        dT = (3*x*T - 4*(A + q*x*Q)/Z)/U;
        d2T = (3*T + 5*x*dT + 4*q^5*Q/Z^5)/U;

    % Orbit type is transitioning, so a power series is used for accuracy
    elseif (x > sqrt(0.7)) && (x < sqrt(1.4))
        [T,dT,d2T] = powerSeries(x,q,Q);

    % Orbit is a hyperbolia and x >> 1
    elseif x > sqrt(1.4)
        if (F + G) > pi
            T = 2*(log(F+G)/Y + B)/U;
        elseif (F + G) <= pi
            T = 2*(atanh(F/G)/Y + B)/U;
        end
        dT = (3*x*T - 4*(A + q*x*Q)/Z)/U;
        d2T = (3*T + 5*x*dT + 4*q^5*Q/Z^5)/U;

    end
end



%% ----- Function: Power Series Setup ----- %%
% -------------------------------------------------------------------------
% Creates a power series representation of T for improved accuracy near the
% transiption period between elliptical, parabolic, and hyperbolic orbits.
% The symbolic derivative function is used to find the derivatives of T.
% This function is only utilized for the no-revolution case.
% -------------------------------------------------------------------------

function [T,dT,d2T] = powerSeries(x,q,Q)
    U = 1 - x^2;
    N = 20;

    alpha = 4;
    tau = q*Q;
    sig = 1 - q^3;
    aprev = 4/3;
    
    S = 4/3*sig;
    S_U = 0;
    S_UU = 0;

    pU0 = U;
    pUm1 = 1;
    pUm2 = 2;

    for n = 1:N
        alpha = alpha*(2*n-1)/(2*n);
        tau = tau*q^2;
        sig = sig + tau;
        aan = alpha/(2*n + 3);
        termVal = (6*n + 1)*aan*sig/((2*n - 1)*(2*n + 1)) - aprev*tau;
        
        S = S - pU0*termVal;
        S_U = S_U - n*pUm1*termVal;
        S_UU = S_UU - n*(n-1)*pUm2*termVal;

        pUm2 = pUm1;
        pUm1 = pU0;
        pU0 = pU0*U;
        aprev = aan;
    end

    T = S/x^2;
    dT = -2*(S_U*x^2 + S)/x^3;
    d2T = 4*S_UU + 6*S_U/x^2 + 6*S/x^4;
end



%% ----- Function: Multi-Revolution x0's ----- %%
% -------------------------------------------------------------------------
% Finds the two starting x values for the multi-revolution cases. This can
% be thought of as finding an x0 on each side of the parabola where xM is
% defined as the lowest point.
% -------------------------------------------------------------------------

function x0 = multiStarter1(xM,dTM,d2TM,M,phir)
    delta = sqrt(dTM/(d2TM/2 + dTM/(1 - xM)^2));
    P = xM + delta;
    Wp = 4*P/(4 + dTM) + (1 - P)^2;
    K = (1 + M + 1*(phir - 0.5))/(1 + 0.15*M);
    x0 = xM + delta*(1 - K*delta*(0.5*Wp + 0.03*delta*sqrt(Wp)));
end

function x0 = multiStarter2(xM,TM,d2TM,T0,Tf,M,phir)
    dTM = Tf - TM;
    dTfT0 = Tf - T0;

    if dTfT0 <= 0
        dT0M = T0 - TM;
        x0 = xM - sqrt(dTM/(d2TM/2 - dTM*(d2TM/(2*dT0M) - 1/xM^2)));
    else
        x01 = - dTfT0/(dTfT0 + 4);
        P = x01 + 1.7*sqrt(2*(1 - phir));
        if P < 0
            x02 = -sqrt(dTfT0/(Tf + 0.5*T0));
            p = (-P)^(1/16);
            x01 = (1 - p)*x01 + p*x02;
        end
        Wp = 4/(4 + dTfT0);
        K1 = (1 + M + 0.24*(phir - 0.5))/(1 + 0.15*M);
        x0 = x01*(1 + K1*x01*(0.5*Wp - 0.03*x01*sqrt(Wp)));
    end
end



% Gooding Function 'end'
end