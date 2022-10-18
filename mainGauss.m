clear; close all; tic;
%% Analysis settings
global figNumber; figNumber = 1;

bodyName = "216 Kleopatra"; % Name of the body (by the input file)
order_int = 6; % Gauss integration order;
spacePoints = 10; % Points in space;
order_spacePoints = 2; % Order of polynomial spacing points;

% Analysis name
alpha = order_int; beta = order_int; gamma = order_int;
analysisName = bodyName + "_n" + num2str(order_int) + "_p" + ...
    num2str(spacePoints) + "_o" + num2str(order_spacePoints);

% Creates the directory if it doesn't exist
figDir = pwd + "\fig\Integration_" + analysisName; checkDir(figDir);
varDir = pwd + "\var"; checkDir(varDir);

%% Main calculations
fprintf('\nStarting potential integration analysis for %s body.\n', ...
    bodyName);

% Inputs
file = importdata(bodyName + ".txt");
vertices = file.data(cell2mat(file.textdata)=='v',:);
faces = file.data(cell2mat(file.textdata)=='f',:);

disp('Analysis settings:')
fprintf('\tFaces number: %d\n', length(body.faces));
fprintf('\tVertices number: %d\n', length(body.vertices));

% Evaluation points
x = polySpaced(-252, 252, 0, spacePoints, order_spacePoints);
y = polySpaced(-252, 252, 0, spacePoints, order_spacePoints);
z = polySpaced(-252, 252, 0, spacePoints, order_spacePoints);

% Get the potential by integration
options = struct('alpha', alpha, 'beta', beta, 'gamma', gamma);
U = getPotential(faces, vertices, x, y, z);

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
save(varDir + "\Integration_" + analysisName, 'U', 'x', 'y', 'z', 'options', ...
    'bodyName', 's');

%% My functions

function U = getPotential(faces, vertices, X, Y, Z, options)
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
    
    total = length(X)*length(Y)*length(Z); count = 0;
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
                
                count = count + 1; % Count the number of points evaluated
            end % End in z for
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