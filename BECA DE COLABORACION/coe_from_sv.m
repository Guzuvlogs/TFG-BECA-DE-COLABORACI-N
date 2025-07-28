function coe = coe_from_sv(R, V,mu)
    % ORBITAL_ELEMENTS Computes classical orbital elements (coe) from the state vector (R, V).
    %
    % Inputs:
    %   R - Position vector in the geocentric equatorial frame (km)
    %   V - Velocity vector in the geocentric equatorial frame (km/s)
    %
    % Outputs:
    % coe - Vector of orbital elements [h e RA incl w TA a]
    % H - the angular momentum vector (kmˆ2/s)
    % h - the magnitude of H (kmˆ2/s)
    % incl - inclination of the orbit (rad)
    % N - the node line vector (kmˆ2/s)
    % n - the magnitude of N
    % cp - cross product of N and R
    % RA - right ascension of the ascending node (rad)
    % E - eccentricity vector
    % e - eccentricity (magnitude of E)
    % eps - a small number below which the eccentricity is
    % considered to be zero
    % w - argument of perigee (rad)
    % TA - true anomaly (rad)
    % a - semimajor axis (km)

    % Dependencies:
    %   None (self-contained)
    
    % Gravitational parameter (mu) for Earth (km^3/s^2)
    %mu = 398600; % Default gravitational parameter
    
    % Tolerance for small values
    eps = 1e-10;
    
    % Magnitudes of position and velocity vectors
    r = norm(R);
    v = norm(V);
    
    % Radial velocity component
    vr = dot(R, V) / r;
    
    % Angular momentum vector and its magnitude
    H = cross(R, V);
    h = norm(H);
    
    % Inclination (Equation 4.7)
    incl = acos(H(3) / h);
    
    % Node line vector and its magnitude (Equation 4.8)
    N = cross([0, 0, 1], H);
    n = norm(N);
    
    % Right ascension of the ascending node (RA) (Equation 4.9)
    if n ~= 0
        RA = acos(N(1) / n);
        if N(2) < 0
            RA = 2 * pi - RA;
        end
    else
        RA = 0; % Node line is undefined
    end
    
    % Eccentricity vector and its magnitude (Equation 4.10)
    E = (1 / mu) * ((v^2 - mu / r) * R - r * vr * V);
    e = norm(E);
    
    % Argument of perigee (w) (Equation 4.12, including case e = 0)
    if n ~= 0
        if e > eps
            w = acos(dot(N, E) / (n * e));
            if E(3) < 0
                w = 2 * pi - w;
            end
        else
            w = 0; % Circular orbit
        end
    else
        w = 0; % Equatorial orbit
    end
    
    % True anomaly (TA) (Equation 4.13a, including case e = 0)
    if e > eps
        TA = acos(dot(E, R) / (e * r));
        if vr < 0
            TA = 2 * pi - TA;
        end
    else
        cp = cross(N, R);
        if cp(3) >= 0
            TA = acos(dot(N, R) / (n * r));
        else
            TA = 2 * pi - acos(dot(N, R) / (n * r));
        end
    end
    
    % Semimajor axis (a) (Equation 2.61, a < 0 for a hyperbola)
    a = h^2 / mu / (1 - e^2);
    
    % Output vector of orbital elements
    coe= [h, e, RA, incl, w, TA, a];
