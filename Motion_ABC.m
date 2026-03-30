function dydt = Motion_ABC(t, y, Omega_hover, I, m, g, l, k, b)

    % Get phi, theta, psi, p, q, and r
    phi = y(7); theta = y(8); psi = y(9);
    p = y(10); q = y(11); r = y(12);

    % Create omega for rotors
    Omega = Omega_hover * ones(4,1);

    % If statements to find Omega
    if t < 1
        Omega = Omega_hover + 70 * sin(2*pi*t/4);
        Omega = ones(4,1) * Omega;
    elseif t < 2
        Omega = Omega_hover - 77 * sin(2*pi*t/4);
        Omega = ones(4,1) * Omega;
    elseif t < 3
        Omega(2) = sqrt(Omega_hover^2 - 70^2 * sin(2*pi*(t-2)/4));
        Omega(4) = sqrt(Omega_hover^2 + 70^2 * sin(2*pi*(t-2)/4));
    elseif t < 4
        Omega(2) = sqrt(Omega_hover^2 + 70^2 * sin(2*pi*(t-2)/4));
        Omega(4) = sqrt(Omega_hover^2 - 70^2 * sin(2*pi*(t-2)/4));
    elseif t < 5
        Omega(1) = sqrt(Omega_hover^2 - 70^2 * sin(2*pi*(t-4)/4));
        Omega(3) = sqrt(Omega_hover^2 + 70^2 * sin(2*pi*(t-4)/4));
    else
        Omega(1) = sqrt(Omega_hover^2 + 70^2 * sin(2*pi*(t-4)/4));
        Omega(3) = sqrt(Omega_hover^2 - 70^2 * sin(2*pi*(t-4)/4));
    end

    % Get T and Tau
    T = k * sum(Omega.^2);
    L = k * l * (-Omega(2)^2 + Omega(4)^2);
    M = k * l * (-Omega(1)^2 + Omega(3)^2);
    N = b * (Omega(1)^2 - Omega(2)^2 + Omega(3)^2 - Omega(4)^2);
    Tau = [L; M; N];

    % Rotation matrix (inertial to body)
    C_IB = [                  cos(theta)*cos(psi),                            cos(theta)*sin(psi),              -sin(theta);
        sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi), sin(phi)*cos(theta);
        cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi), cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi), cos(phi)*cos(theta)];

    % Body to inertial rotation
    C_BI = C_IB';
    Thrust = C_BI(:,3);

    % Get r double dot
    r_dd = [0;0;-g] + (T/m) * Thrust;

    % Get w and Omega_dot
    w = [p; q; r];
    Omega_dot = I \ ( - cross(w, I*w) + Tau);

    % Euler angle rates
    A = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
         0,      cos(phi),           -sin(phi);
         0, sin(phi)/cos(theta), cos(phi)/cos(theta)];
    Euler = A \ w;

    if t > 4
        Euler(1) = 0;
    end

    % dydt
    dydt = zeros(12,1);
    dydt(1:3) = y(4:6);
    dydt(4:6) = r_dd;
    dydt(7:9) = Euler;
    dydt(10:12) = Omega_dot;
end


