clear; tic; close all;
%% Analysis settings
global figNumber; figNumber = 1;

bodyName = "216 Kleopatra"; % Name of the body (by the input file)
numberOfDivisions = [1, 3, 5, 7]; % Gauss integration order;
spacePoints = 100; % Points in space;
order_spacePoints = 3; % Order of polynomial spacing points;
order_int = [2, 3, 5, 7, 9]; % Gauss integration order;
colors = ['r', 'b', 'c', 'k', 'm', 'g', 'y'];

% Analysis name
analysisName = bodyName + "Comp_" + "n_" + num2str(order_int(1)) + ...
    "to" + num2str(order_int(end)) + "_d_" + ...
    num2str(numberOfDivisions(1)) + "to" + ...
    num2str(numberOfDivisions(end));

% Creates the directory if it doesn't exist
figDir = pwd + "\fig\COMP_" + analysisName; checkDir(figDir);
varDir = pwd + "\var"; checkDir(varDir);

%% Inputs
fprintf('\nStating comparision analysis for %s body.\n', bodyName);
file = importdata(bodyName + ".txt");

% Get the number of vertices
vertices = file.data(cell2mat(file.textdata)=='v',:);
faces = file.data(cell2mat(file.textdata)=='f',:);

% Pre allocation memory
x = linspace(90, 252, spacePoints);
y = 0; z = 0;
Uint = zeros(length(x), length(order_int));
Umas = zeros(length(x), length(numberOfDivisions));

% Figura do corpo
figure(figNumber);
patch('faces', faces, 'vertices', vertices, 'EdgeColor', 'k', ...
    'FaceColor','none'); axis equal; view(30,10);
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
saveas(figure(figNumber), figDir + "\bodyDescription.png");
figNumber = figNumber + 1;

%% Main calculation
for ii = 1:length(numberOfDivisions)
    % Mascons
    fprintf("\t\t Calculation of MASCONS for d = %d\n", ...
        numberOfDivisions(ii));
    MSC = getMASCONS(faces, vertices, numberOfDivisions(ii), figDir, ...
        colors);
    Umas(:,ii) = getPotencialMASCONS(MSC, x, y, z);
end
for ii = 1:length(order_int)
    % Integration
    fprintf("\t\t Calculation of integration method for n = %d\n", ...
        order_int(ii));
    options = struct('alpha', order_int(ii), 'beta', order_int(ii), ...
        'gamma', order_int(ii));
    Uint(:,ii) = getPotentialIntegration(faces, vertices, x, 0,0, options);
end

%% Figures
% Convergence of integration method
figure(figNumber);
for ii = 1:length(order_int)
    plot(x, Uint(:, ii), colors(ii), 'LineWidth', 1.2, 'DisplayName', ...
        "\alpha = \beta = \gamma = " + num2str(order_int(ii))); hold on;
end
G = 6.67259e-20; % [m³/kg*s] Gravitational constant;
M = 2.97e18; % [kg/m³] Mean body density;
plot(x, -G*M./x, '--k', 'LineWidth', 1.2, ...
    'DisplayName', 'Aprox');
grid minor; xlabel("X [km]"); ylabel("U"); legend();
saveas(figure(figNumber), figDir + "\convIntegration.png");
figNumber = figNumber + 1;

% Convergence of integration method (error)
figure(figNumber);
for ii = 1:length(order_int)
    plot(x, (Uint(:, ii) - Uint(:, end))*100./Uint(:, end), colors(ii), ...
        'LineWidth', 1.2, 'DisplayName', "\alpha = \beta = \gamma = " + ...
        num2str(order_int(ii))); hold on;
end
grid minor; xlabel("X [km]"); ylabel("Erro [%%]"); legend();
saveas(figure(figNumber), figDir + "\convIntegrationError.png");
figNumber = figNumber + 1;

% Convergence of MASCONS method
figure(figNumber);
for ii = 1:length(numberOfDivisions)
    plot(x, Umas(:, ii), colors(ii), 'LineWidth', 1.2, 'DisplayName', ...
        num2str(numberOfDivisions(ii)) + "divisions"); hold on;
end
plot(x, -G*M./x, '--k', 'LineWidth', 1.2, ...
    'DisplayName', 'Aprox');
grid minor; xlabel("X [km]"); ylabel("U"); legend();
saveas(figure(figNumber), figDir + "\convMASCONS.png");
figNumber = figNumber + 1;

% Convergence of MASCONS method (error)
figure(figNumber);
for ii = 1:length(numberOfDivisions)
    plot(x, (Umas(:, ii) - Uint(:, end))*100./Uint(:, end), colors(ii), ...
        'LineWidth', 1.2, 'DisplayName', ...
        num2str(numberOfDivisions(ii)) + "divisions"); hold on;
end
grid minor; xlabel("X [km]"); ylabel("Erro [%%]"); legend();
saveas(figure(figNumber), figDir + "\convMASCONSError.png");
figNumber = figNumber + 1;

%% Outs
fprintf("Results: \n");
fprintf("Biggest error value for the last MASCONS evaluation: %.2f", ...
    max(abs((Umas(:, ii) - Uint(:, end))*100./Uint(:, end))));

%% My functions

function MSC = getMASCONS(faces, vertices, numberOfDivisions, figDir, ...
    colors)
    % Function responsable for getting the centroids of each MASCONS and
    % their respective volume. The only input is the body description, un
    % struct with the filds 'vertices' and 'faces'.
    
    % Inputs
    global figNumber;
    
    if ~exist('numberOfDivisions', 'var')
        numberOfDivisions = 1;
    end
    
    len = length(faces);
    L = numberOfDivisions*len;
    MSC = struct();
        MSC.centers = zeros(3, L);
        MSC.volume = zeros(1, L);
    
    figure(figNumber); hold on;
    
    % General calculations
    A = vertices(faces(:,1), :)'/numberOfDivisions;
    B = vertices(faces(:,2), :)'/numberOfDivisions;
    C = vertices(faces(:,3), :)'/numberOfDivisions;
    
    MSC.centers(:, 1:len) = ...
        (A + B + C + zeros(3, len))/4;
    MSC.volume(1:len) = dot(C, cross(A, B))/6;
    
    plot3(MSC.centers(1,:), MSC.centers(2,:), MSC.centers(3,:), 'o', ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 0.8);
    
    for i = 2:numberOfDivisions
        Anew = vertices(faces(:,1), :)'*i/numberOfDivisions;
        Bnew = vertices(faces(:,2), :)'*i/numberOfDivisions;
        Cnew = vertices(faces(:,3), :)'*i/numberOfDivisions;
        
        MSC.centers(:, (1:len) + (i-1)*L) = ...
            (A + B + C + Anew + Bnew + Cnew)/6;
        
        MSC.volume(:, (1:len) + (i-1)*L) = ...
            dot(Cnew, cross(Anew, Bnew))/6 - dot(C, cross(A, B))/6;
        
        % Plot the MASCONS centers
        plot3(MSC.centers(1, (1:len) + (i-1)*L), ...
              MSC.centers(2, (1:len) + (i-1)*L), ...
              MSC.centers(3, (1:len) + (i-1)*L), 'o', ...
              'MarkerFaceColor', colors(i-1), 'MarkerEdgeColor', ...
              colors(i-1), 'MarkerSize', 0.8);
        
        A = Anew; B = Bnew; C = Cnew;
    end
    
    % Plots
    
    axis equal; view(30,10);
    xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
    saveas(figure(figNumber), figDir + "\MASCONS_d" + ...
        num2str(numberOfDivisions) + ".png"); close(figNumber);
    figNumber = figNumber + 1;
    
end

function U = getPotentialIntegration(faces, vertices, X, Y, Z, options)
    % Function to get the total potential of the body at a poit (x, y, z).
    % If the points ar vectors, the potencial are evaluated at each point.
    
    % Constants
    G = 6.67259e-20; % [m³/kg*s] Gravitational constant;
    
    % Inputs
    if ~exist('options', 'var')
       options = struct('alpha', 4, 'beta', 4, 'gamma', 4); 
    end
    alpha = options.alpha; beta = options.beta; gamma = options.gamma;
    M = alpha*beta*gamma;
    rho = 4270*1e9; % [kg/m³] Mean body density;
    
    % Initializing
    U = zeros(length(X), length(Y), length(Z));
    table = gaussTable(alpha, beta, gamma);
    
    % Pass through all points of interest
    for i = 1:length(X)
        for j = 1:length(Y)
            for k = 1:length(Z)

                % Do the summation
                U0 = 0;
                for n = 1:length(faces)
                    % Create the linear transformation matrix T
                    T = [vertices(faces(n,1), :)', ...
                         vertices(faces(n,2), :)', ...
                         vertices(faces(n,3), :)'];
                    sum = 0;
                    for m = 1:M
                        % Take coefficients for quadrature from table
                        xm = table(m,1); ym = table(m,2); zm = table(m,3);
                        
                        % Sum (constant density)
                        sum = sum + table(m,4)*rho/...
                                     norm(T*[xm;ym;zm] - [X(i);Y(j);Z(k)]);
                    end 
                    U0 = U0 - sum*det(T);
                end
                
                % Take the result
                U(i, j, k) = U0;
                
            end % End in z for
        end % End in y for
    end % End in z for
    
    U = U*G;
end

function [U, Ur] = getPotencialMASCONS(MSC, X, Y, Z)
    % Function to get the total potential of the body at a poit (x, y, z).
    % If the points ar vectors, the potencial are evaluated at each point.
    
    G = 6.67259e-20; % [m³/kg*s] Gravitational constant;
    rho = 4270*1e9; % [kg/m³] Mean body density;
    U = zeros(length(X), length(Y), length(Z));
    Ur = zeros(length(X), length(Y), length(Z));
    
    for i = 1:length(X)
        for j = 1:length(Y)
            for k = 1:length(Z)
                % Potencial for each MASCONS divided by G, saved as vector
                r = [X(i); Y(j); Z(k)];
                Uu = MSC.volume*rho ./ sqrt(sum((MSC.centers - r).^2));
                Uur = MSC.volume*rho ./ (sum((MSC.centers - r).^2));

                % Total potencial (sum of each MASCONS potencial)
                U(i, j, k) = -G*sum(Uu);
                Ur(i, j, k) = G*sum(Uur);
            end
        end
    end    
end

function checkDir(nameDir)
    if ~isfolder(nameDir)
        mkdir(nameDir)
    end
end
