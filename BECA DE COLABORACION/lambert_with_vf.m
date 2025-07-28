function [v1_opt, v2_opt, tof_opt] = lambert_with_vf(r1, r2, vf_target, tof_guess, mu)
    % Resuelve el problema de Lambert ajustando el tiempo de vuelo
    % para obtener la velocidad final deseada.
    %
    % Parámetros:
    % r1        - Posición inicial [km]
    % r2        - Posición final [km]
    % vf_target - Velocidad final deseada [km/s]
    % tof_guess - Estimación inicial del tiempo de vuelo [s]
    % mu        - Parámetro gravitacional del cuerpo [km^3/s^2]
    %
    % Retorna:
    % v1_opt  - Velocidad inicial óptima [km/s]
    % v2_opt  - Velocidad final obtenida [km/s] (debe ser cercana a vf_target)
    % tof_opt - Tiempo de vuelo óptimo [s]

    % Función de coste: minimiza la diferencia entre v2 y vf_target
    cost_function = @(tof) norm(lambert_solver(r1, r2, tof, mu,0) - vf_target);

    % Optimización para encontrar el tiempo de vuelo óptimo
    options = optimset('Display', 'iter', 'TolFun', 1e-6);
    tof_opt = fminsearch(cost_function, tof_guess, options);

    % Obtener las velocidades óptimas con el nuevo tof_opt
    [v1_opt, v2_opt] = lambert_solver(r1, r2, tof_opt, mu);
end

function [v1, v2] = lambert_solver(r1, r2, tof, mu, v1_override)
    % Implementación del problema de Lambert (usando Universal Variables)
    % Si se proporciona v1_override, se usa como condición inicial.

    % Resolver el problema de Lambert
    [v1, v2] = Lambert_ChatGPT(r1, r2, tof,string,mu); % Asumiendo prograde 

    % Si se proporciona v1_override, recalcular v2 con ese v1
    if nargin == 5
        v1 = v1_override;
        [~, v2] = Lambert_ChatGPT(r1, r2, tof,string,mu); % Recalcular v2 con v1 dado
    end
end
