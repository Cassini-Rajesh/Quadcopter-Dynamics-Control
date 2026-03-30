clear; clc;

% Constants
m = 0.450;                                              % Mass (kg)
g = 9.81;                                               % Gravity (m/s^2)
l = 0.225;                                              % Arm length (m)
k = 2.98 * 10^-7;                                       % k coefficient
b = 1.14 * 10^-6;                                       % b coefficient
I = diag([4.85 * 10^-3, 4.85 * 10^-3, 8.80 * 10^-3]);   % Inertia matrix (kg*m^2)

% Omega_hover speed (rad/s)
Omega_hover = sqrt( (m * g) / (4 * k) );

% Initial conditions and time span
y0 = zeros(12,1);
tspan = [0 6];

% Use ode45
[t, Y] = ode45(@(t, y) Motion_ABC(t, y, Omega_hover, I, m, g, l, k, b), tspan, y0);

% Plot position
figure('Name', 'Position');
plot(t, Y(:,1), 'b', t, Y(:,2), 'r', t, Y(:,3), 'm');
legend('x', 'y', 'z', 'Location', 'eastoutside');
title('Position vs Time'); 
xlabel('Time (s)'); 
ylabel('Position (m)');

% Plot velocity
figure('Name', 'Velocity');
plot(t, Y(:,4), 'b', t, Y(:,5), 'r', t, Y(:,6), 'm');
legend('xvel', 'yvel', 'zvel', 'Location', 'eastoutside');
title('Velocity vs Time'); 
xlabel('Time (s)'); 
ylabel('Velocity (m/s)');

% Plot Euler angles
figure('Name', 'Euler Angles');
plot(t, (Y(:,7)), 'b', t, (Y(:,8)), 'r', t, (Y(:,9)), 'm');
legend('\phi (roll)', '\theta (pitch)', '\psi (yaw)', 'Location', 'eastoutside');
title('Euler Angles vs Time'); 
xlabel('Time (s)'); 
ylabel('Angle (deg)');

% Plot angular velocities
figure('Name', 'Angular Velocities');
plot(t, Y(:,10), 'b', t, Y(:,11), 'r', t, Y(:,12), 'm');
legend('p', 'q', 'r', 'Location', 'eastoutside');
title('Angular Velocities vs Time'); 
xlabel('Time (s)'); 
ylabel('Angular Velocity (rad/s)');

       