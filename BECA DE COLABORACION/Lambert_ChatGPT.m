function [V1, V2] = Lambert_ChatGPT(R1, R2, t, string,mu)
    % LAMBERT Solves Lambert's problem.
    % 
    % Inputs:
    % R1, R2 - Initial and final position vectors (km)
    % t - Time of flight (s)
    % string - 'pro' for prograde, 'retro' for retrograde
    %
    % Outputs:
    % V1, V2 - Initial and final velocity vectors (km/s)

    %global mu
    %mu = 398600; % Standard gravitational parameter for Earth (km^3/s^2)

    % Magnitudes of R1 and R2:
    r1 = norm(R1);
    r2 = norm(R2);

    % Angle between R1 and R2:
    c12 = cross(R1, R2);
    theta = acos(dot(R1, R2) / (r1 * r2));

    % Determine prograde or retrograde trajectory:
    if strcmp(string, 'pro')
        if c12(3) < 0
            theta = 2 * pi - theta;
        end
    elseif strcmp(string, 'retro')
        if c12(3) > 0
            theta = 2 * pi - theta;
        end
    else
        %fprintf('\n ** Prograde trajectory assumed.\n');
    end

    % Equation 5.35:
    A = sin(theta) * sqrt(r1 * r2 / (1 - cos(theta)));

    % Initial guess for z:
    z = -100;
    while F(z, t, A, r1, r2,mu) < 0
        z = z + 0.1;
    end

    % Set tolerance and maximum iterations:
    tol = 1e-8;
    nmax = 5000;

    % Newton-Raphson iteration to solve for z:
    ratio = 1;
    n = 0;
    while abs(ratio) > tol && n <= nmax
        n = n + 1;
        ratio = F(z, t, A, r1, r2,mu) / dFdz(z, A, r1, r2);
        z = z - ratio;
    end

    % Check for convergence:
    if n > nmax
        error('Number of iterations exceeds %d', nmax);
    end

    % Lagrange coefficients:
    f = 1 - y(z, r1, r2, A) / r1;
    g = A * sqrt(y(z, r1, r2, A) / mu);
    gdot = 1 - y(z, r1, r2, A) / r2;

    % Velocity vectors:
    V1 = (1/g)*(R2 - f*R1);
    V2 = (1/g)*((real(gdot)*R2)-R1);
end

% Subfunctions:

function val = y(z, r1, r2, A)
    % Equation 5.38:
    val = r1 + r2 + A * (z * stumpS(z) - 1)/sqrt(stumpC(z));
end

function val = F(z, t, A, r1, r2,mu)
    % Equation 5.40:
    val = ((y(z, r1, r2, A) / stumpC(z))^1.5)*stumpS(z) + A*sqrt(y(z, r1, r2, A)) -sqrt(mu)*t;
end

function val = dFdz(z, A, r1, r2)
    % Equation 5.43:
    if z == 0
        y0 = y(0, r1, r2, A);
        val = sqrt(2) / 40 * y0^1.5 + A / 8 * (sqrt(y0) + A * sqrt(1 / (2 * y0)));
    else
        y_z = y(z, r1, r2, A);
        val = (y_z / stumpC(z))^1.5 * (0.5 / z * (stumpC(z) - 3 * stumpS(z) / (2 * stumpC(z))) + 3 * stumpS(z)^2 / (4 * stumpC(z))) ...
            + A / 8 * (3 * stumpS(z) / stumpC(z) * sqrt(y_z) + A * sqrt(stumpC(z) / y_z));
    end
end

function val = stumpC(z)
    % Stumpff function C(z)
    if z > 0
        val = (1 - cos(sqrt(z))) / z;
    elseif z < 0
        val = (cosh(sqrt(-z)) - 1) / -z;
    else
        val = 1 / 2;
    end
end

function val = stumpS(z)
    % Stumpff function S(z)
    if z > 0
        val = (sqrt(z) - sin(sqrt(z))) / z^1.5;
    elseif z < 0
        val = (sinh(sqrt(-z)) - sqrt(-z)) / (-z)^1.5;
    else
        val = 1 / 6;
    end
end