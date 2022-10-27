clear; close all; tic;
%% Analysis settings
global figNumber G rho; 
figNumber = 1;

bodyName = "216 Kleopatra"; % Name of the body (by the input file)
order_int = 6; % Gauss integration order;
spacePoints = 50; % Points in space;
maxDist = 700; % [km] Maximum distance avaluation
order_spacePoints = 2; % Order of polynomial spacing points;
rho = 4.270e12; % [kg/km³] Mean body density;
G = 6.67259e-20; % [km³/kg*s] Gravitational constant;

% Analysis name
analysisName = bodyName + "_n" + num2str(order_int) + "_p" + ...
    num2str(spacePoints) + "_o" + num2str(order_spacePoints) + ...
    "_m" + num2str(maxDist);

% Creates the directory if it doesn't exist
figDir = pwd + "\results\Integration_" + analysisName + "\fig"; 
varDir = pwd + "\results\Integration_" + analysisName; 
checkDir(figDir);

%% Main calculations
fprintf('\nStarting potential integration analysis for %s body.\n', ...
    bodyName);

% Inputs
file = importdata(bodyName + ".txt");
vertices = file.data(cell2mat(file.textdata)=='v',:);
faces = file.data(cell2mat(file.textdata)=='f',:);

disp('Analysis settings:')
fprintf('\tFaces number: %d\n', length(faces));
fprintf('\tVertices number: %d\n', length(vertices));

% Evaluation points
x = polySpaced(-maxDist, maxDist, 0, spacePoints, order_spacePoints);
y = polySpaced(-maxDist, maxDist, 0, spacePoints, order_spacePoints);
z = polySpaced(-maxDist, maxDist, 0, spacePoints, order_spacePoints);

% Get the potential by integration
options = struct('alpha',order_int,'beta',order_int,'gamma',order_int);
U = getPotential(faces, vertices, x, y, z, options);

% Plot potential data
figure(figNumber);  [X, Y] = meshgrid(x, y);
surf(X, Y, U(:, :, floor(length(z)/2)), 'EdgeColor', 'none'); colorbar(); 
xlabel('X [km]'); ylabel('Y [km]'); zlabel('U'); view(60, 35);
saveas(figure(figNumber), figDir + "\potential.png");
figNumber = figNumber + 1;

%% Ending
s = toc;
fprintf('Code ended successfully.\n');
fprintf('Total execution time: %dmin e %.2fs.\n', floor(s/60), ...
    s - floor(s/60)*60);
save(varDir + "\var", 'U', 'x', 'y', 'z', ...
    'options', 'bodyName', 's');

%% My functions

function U = getPotential(faces, vertices, X, Y, Z, options)
    % Function to get the total potential of the body at a poit (x, y, z).
    % If the points ar vectors, the potencial are evaluated at each point.

    global G rho;
    
    % Inputs
    if ~exist('options', 'var')
       options = struct('alpha', 4, 'beta', 4, 'gamma', 4); 
    end
    alpha = options.alpha; 
    beta = options.beta; 
    gamma = options.gamma;
    M = alpha*beta*gamma;
    
    % Get the vertices
    A = vertices(faces(:,1), :)';
    B = vertices(faces(:,2), :)';
    C = vertices(faces(:,3), :)';
    
    % Initializing
    U = zeros(length(X), length(Y), length(Z));
    table = gaussTable(alpha, beta, gamma);
    xm = table(:,1); ym = table(:,2); zm = table(:,3); wm = rho*table(:,4);
    
    total = length(X)*length(Y); count = 0;
    % Pass through all points of interest
    for i = 1:length(X)
        for j = 1:length(Y)
            
            x = X(i);
            y = Y(j);
            parfor k = 1:length(Z)
                r = [x; y; Z(k)];
                
                % Do the summation
                U0 = 0;
                for n = 1:length(faces)
                    % Create the linear transformation matrix T
                    T = [A(:,n), B(:,n), C(:,n)];
                    
                    U0 = U0 -det(T)*sum(wm'./ ...
                                          vecnorm(T*[xm'; ym'; zm'] - r));           
                end
                
                % Take the result
                U(i, j, k) = U0;

            end % End in z for
            
            count = count + 1; % Count the number of points evaluated
            fprintf("\t\tProgress: %4.2f %%\n", count*100/total);
        end % End in y for
    end % End in z for
    
    U = U*G;
end

function out = polySpaced(startP, endP, peakP, pointsNumber, exp)
    % Function that creats a polynomial spaced vector, starting in startP
    % point, until endP point, with a peak in peakP and with pointNumber
    % points. Exp input is the exp value for the polynomial. If its not
    % specified, then the functions assumes exp = 2.
    
    if ~exist('exp', 'var')
        exp = 2;
    end
    
    aux1 = floor(pointsNumber/2) - 1;
    aux2 = ceil(pointsNumber/2) - 1;
    out = [peakP + (startP - peakP)*((aux1-(0:aux1))/aux1).^exp, ...
           peakP + (endP - peakP)*((0:aux2)/aux2).^exp];
end

function checkDir(nameDir)
    if ~isfolder(nameDir)
        mkdir(nameDir)
    end
end