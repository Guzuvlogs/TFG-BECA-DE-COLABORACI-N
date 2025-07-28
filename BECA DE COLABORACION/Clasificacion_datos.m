clear;clc
load ("datos_trayectoria_mejorada.mat")
% Vamos a hacer una versión preliminar del código para analizar y
% clasificar los datos de las distintas trayectorias.Lo más importante es
% el análisis de delta_V. Buscamos una muy baja.

 function [valores, ubicaciones] = filtrarValoresCelda(celda, umbral)
    % FILTRARVALORESCELDA Extrae valores menores que el umbral en una celda de matrices
    %   [valores, ubicaciones] = filtrarValoresCelda(celda, umbral)
    %   devuelve un vector con los valores menores que el umbral y una matriz
    %   con el índice de la celda, la fila y la columna de cada uno.
    %   Analiza valores positivos. Aplica un valor absoluto dentro de la
    %   propia función
    
    valores = []; % Inicialización del vector de valores
    ubicaciones = []; % Inicialización de la matriz de ubicaciones
    
    for k = 1:length(celda) % Recorre cada término de la celda
        matriz = celda{k}; % Extrae la matriz en la celda
        [filas, columnas] = find(abs(matriz) < umbral); % Encuentra valores por debajo del umbral
        
        if ~isempty(filas) % Si hay valores que cumplen la condición
            nuevos_valores = matriz(sub2ind(size(matriz), filas, columnas))';
            nuevas_ubicaciones = [repmat(k, length(filas), 1), filas', columnas'];
            
            valores = [valores; nuevos_valores];
            ubicaciones = [ubicaciones; nuevas_ubicaciones];
        end
    end
end

 %Una vez tenemos esta función, basta con usarla en primer lugar para
 %Delta_V. Vamos a imponer primeramente un valor de 10km/s, pues con un
 %análisis no muy exigente se han encontrado valores más pequeños de este

 [D_V_requirements, ubicaciones] = filtrarValoresCelda(Delta_V, 10);
 disp(size(ubicaciones))

 %Ahora dependerá del tamaño de esta muestra. Supongamos que es bastante
 %grande. Lo más interesante ahora consiste en tiempos no muy prolongados de vuelo,
 % así como velocidades no excesivamente grandes de cara a escapar de la
 % Tierra. Molaría imponer límites no demasiados restrictivos. Estos
 % podemos tomarlos en 45km/s, 1.5 años y 6 años respectivamente.

excentricidad_=zeros(length(D_V_requirements),1);
v_infinito_=zeros(length(D_V_requirements),1);
delta__=zeros(length(D_V_requirements),1);
periapside_=zeros(length(D_V_requirements),1);
Delta_V_=zeros(length(D_V_requirements),1);
fecha_earth_=cell(length(D_V_requirements),1);
fecha_marte_=cell(length(D_V_requirements),1);
fecha_jupiter_=cell(length(D_V_requirements),1);

for i=1:length(ubicaciones)
        excentricidad_(i)=excentricidad{ubicaciones(i,1)}(ubicaciones(i,3));
        v_infinito_(i)=v_infinito{ubicaciones(i,1)}(ubicaciones(i,3));
        delta__(i)=deltaa_{ubicaciones(i,1)}(ubicaciones(i,3));
        periapside_(i)=periapside{ubicaciones(i,1)}(ubicaciones(i,3));
        v0(i)=v0_(ubicaciones(i,1));
        
        
        fecha_earth_{i}=fecha_earth{ubicaciones(i,1)}; %Fecha Salida de Marte
        fecha_marte_{i}=fecha_marte{ubicaciones(i,1)}(ubicaciones(i,3)); %Fecha llegada a Marte
        fecha_jupiter_{i}=fecha_jupiter{ubicaciones(i,1)}(ubicaciones(i,3)); %Fecha llegada a Jupiter
            
end

%Tabla
T1 = table(fecha_earth_,fecha_marte_, fecha_jupiter_, excentricidad_, v_infinito_, delta__,v0', periapside_, D_V_requirements, ...
    'VariableNames', {'Fecha de salida Tierra','Fecha de llegada a Marte', 'Fecha de llegada a Jupiter', 'e', 'V_infinity (Km/s)', 'Delta (rad)','V0 (km/s)', 'Rp (Km)', 'Delta_V (Km/s)'});
disp(T);
%A lo mejor este while debe de ir fuera del bucle for. Al final lo que
%quiero yo es que mientras el elemento de la primera columna de
%"ubicaciones y ind_tierra coincida, que se escriba así. Cuando no
%coincidan, tengo que bajar en ind_tierra hasta que los elementos de la
%primera columna coincidan. Una vez hecho eso ya estaría

%La fila de la tabla 368 tiene una pinta increible. El viaje dura menos de
%2,5 años. Las condiciones de delta_V son muy pequeños. Todo resulta ideal.
%Las fechas tienen lugar a inicios de la década siguiente. El único
%problema es v0. Esta se ha conseguido reducirla notoriamente respecto a la
%dada en el paper inicial. El único problema es que los lanzadores actuales
%tienen una potencia de inyección de 16km/s. Esto haría necesario y delta_V
%inicial de 9km/s. 

%La continuación del trabajo recae en correr una simulación y comprobar que
%efectivamente, la tangente de la trayectoria de los problemas de Lambert
%es prácticamente idéntica. Si se corrobora aquello, necesitaríamos
%implementar un simulador que permita, efectivamente, aplicar aquel delta_V

%Si esto se da de una forma correcta, la siguiente fase residiría en
%aplicar control, introducir perturbaciones, y hacer realista el problema.
%No obstante, el problema más complicado, la trayectoria analítica, está
%resuelto. Podrían buscarse nuevas fechas para intentar reducir esta
%velocidad inicial desde la Tierra.