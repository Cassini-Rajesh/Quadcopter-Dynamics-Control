%% Problem 2

clear; clc;

g = 9.81;
m = 0.450;
k = 2.98E-7; 
b = 1.14E-6;
l = 0.225;
Ixx = 4.85E-3; Iyy = Ixx; Izz = 8.8E-3;

% Reference variables
zDotRef = 0; zRef = 10;
phiDotRef = 0; phiRef = 0;
thetaDotRef = 0; thetaRef = 0;
psiDotRef = 0; psiRef = 0;

% Autopilot equations
T = @(zDot, z, phi, theta) (g+(zDotRef - zDot) + (zRef - z)) * (m/(cos(phi)*cos(theta)));
L = @(phiDot, phi) Ixx*((phiDotRef - phiDot) + (phiRef - phi));
M = @(thetaDot, theta) Iyy*((thetaDotRef - thetaDot) + (thetaRef - theta));
N = @(psiDot, psi) Izz * ((psiDotRef - psiDot) + (psiRef - psi));

% Rotor Rate calculation
Rotor_RateSqrd = @(T, L, M, N) [T/(4*k) - M/(2*k*l) + N/(4*b);
                                T/(4*k) - L/(2*k*l) - N/(4*b);
                                T/(4*k) + M/(2*k*l) + N/(4*b);
                                T/(4*k) + L/(2*k*l) - N/(4*b)];

% Initial conditions
tspan = [0 120];
y0 = zeros(12,1);
y0(3) = 1;                      
y0(6) = 0;                      
y0(7:9) = deg2rad([10 10 10]); 

[t, Y] = ode45(@(t, y) odefun(t, y), tspan, y0);


% Function for odefun
function dydt = odefun(t, y)
    g = 9.81; m = 0.450; l = 0.225; k = 2.98e-7; b = 1.14e-6;
    Ixx = 4.85e-3; Iyy = Ixx; Izz = 8.8e-3;

    zRef = 10; zDotRef = 0;
    phiRef = 0; phiDotRef = 0;
    thetaRef = 0; thetaDotRef = 0;
    psiRef = 0; psiDotRef = 0;

    z = y(3); zDot = y(6);
    phi = y(7); theta = y(8); psi = y(9);
    p = y(10); q = y(11); r = y(12);

    T_cmd = (g + (zDotRef - zDot) + (zRef - z)) * (m / (cos(phi)*cos(theta)));
    L_cmd = Ixx*((phiDotRef - p) + (phiRef - phi));
    M_cmd = Iyy*((thetaDotRef - q) + (thetaRef - theta));
    N_cmd = Izz*((psiDotRef - r) + (psiRef - psi));

    Omega_sq = [
        T_cmd/(4*k) - M_cmd/(2*k*l) + N_cmd/(4*b);
        T_cmd/(4*k) - L_cmd/(2*k*l) - N_cmd/(4*b);
        T_cmd/(4*k) + M_cmd/(2*k*l) + N_cmd/(4*b);
        T_cmd/(4*k) + L_cmd/(2*k*l) - N_cmd/(4*b)
    ];
    Omega = sqrt(max(Omega_sq, 0));

    tau = [k*l*(-Omega(2)^2 + Omega(4)^2);
           k*l*(-Omega(1)^2 + Omega(3)^2);
           b*(Omega(1)^2 - Omega(2)^2 + Omega(3)^2 - Omega(4)^2)];

    R = [cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta);
         sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi), sin(phi)*cos(theta);
         cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi), cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi), cos(phi)*cos(theta)];
    acc = -g * [0; 0; 1] + (T_cmd / m) * R(:,3);

    A = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
         0, cos(phi), -sin(phi);
         0, sin(phi)/cos(theta), cos(phi)/cos(theta)];
    eulerRates = A \ [p; q; r];

    omegaDot = [((Iyy - Izz)*q*r + tau(1))/Ixx;
                ((Izz - Ixx)*p*r + tau(2))/Iyy;
                ((Ixx - Iyy)*p*q + tau(3))/Izz];

    dydt = zeros(12,1);
    dydt(1:3) = y(4:6);
    dydt(4:6) = acc;
    dydt(7:9) = eulerRates;
    dydt(10:12) = omegaDot;
end


% Rotor speeds
rotor_speeds = zeros(length(t), 4);
for i = 1:length(t)
    z = Y(i,3); zDot = Y(i,6);
    phi = Y(i,7); theta = Y(i,8); psi = Y(i,9);
    p = Y(i,10); q = Y(i,11); r = Y(i,12);

    T_bro = (g + (zDotRef - zDot) + (zRef - z)) * (m / (cos(phi)*cos(theta)));
    L_bro = Ixx*((phiDotRef - p) + (phiRef - phi));
    M_bro = Iyy*((thetaDotRef - q) + (thetaRef - theta));
    N_bro = Izz*((psiDotRef - r) + (psiRef - psi));

    Omega_sq = Rotor_RateSqrd(T_bro, L_bro, M_bro, N_bro);
    rotor_speeds(i, :) = sqrt(max(Omega_sq, 0))';
end

% Plot position in inertial frame
figure; plot(t, Y(:,1), 'b', t, Y(:,2), 'r', t, Y(:,3), 'm');
legend('x', 'y', 'z'); 
title('Position in Inertial Frame'); 
xlabel('Time (s)'); 
ylabel('Position (m)');

% Plot velocity in body frame
figure; plot(t, Y(:,4), 'b', t, Y(:,5), 'r', t, Y(:,6), 'm');
legend('xvel', 'yvel', 'zvel'); 
title('Velocity in Body Frame'); 
xlabel('Time (s)'); 
ylabel('Velocity (m/s)');

% Plot Euler angles
figure; plot(t, (Y(:,7)), 'b', t, (Y(:,8)), 'r', t, (Y(:,9)), 'm');
legend('\phi', '\theta', '\psi'); 
title('Euler Angles'); 
xlabel('Time (s)'); 
ylabel('Angle (deg)');

% Plot angular velocities
figure; plot(t, (Y(:,10)), 'b', t, (Y(:,11)), 'r', t, (Y(:,12)), 'm');
legend('p', 'q', 'r'); 
title('Angular Velocities'); 
xlabel('Time (s)'); 
ylabel('Angular Velocity (deg/s)');

% Plot rotor speeds
figure; plot(t, rotor_speeds(:,1), 'b', t, rotor_speeds(:,2), 'r', t, rotor_speeds(:,3), 'm', t, rotor_speeds(:,4), 'k');
legend('\Omega_1', '\Omega_2', '\Omega_3', '\Omega_4'); 
title('Rotor Angular Velocities'); 
xlabel('Time (s)'); 
ylabel('Rad/s');
