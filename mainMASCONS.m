clear; tic; close all;
%% Analysis settings
global figNumber G rho; 
figNumber = 1;

bodyName = "216 Kleopatra"; % Name of the body (by the input file)
maxDist = 700; % [km] Maximum distance avaluation
numberOfDivisions = 5; % Gauss integration order;
spacePoints = 50; % Points in space;
order_spacePoints = 2; % Order of polynomial spacing points;
rho = 4.270e12; % [kg/km³] Mean body density;
G = 6.67259e-20; % [km³/kg*s] Gravitational constant;

% Analysis name
analysisName = bodyName + "_d" + num2str(numberOfDivisions) + "_p" + ...
    num2str(spacePoints) + "_o" + num2str(order_spacePoints) + ...
    "_m" + num2str(maxDist);

% Creates the directory if it doesn't exist
figDir = pwd + "\results\MASCONS_" + analysisName + "\fig"; 
varDir = pwd + "\results\MASCONS_" + analysisName; 
checkDir(figDir);

%% Main calculations
fprintf('\nStating MASCONS analysis for %s body.\n', bodyName);
file = importdata(bodyName + ".txt");

% Get the number of vertices
vertices = file.data(cell2mat(file.textdata)=='v',:);
faces = file.data(cell2mat(file.textdata)=='f',:);
    
% Figura do corpo
figure(figNumber);
patch('faces', faces, 'vertices', vertices, 'EdgeColor', 'k', ...
    'FaceColor','none'); axis equal; view(30,10);
saveas(figure(figNumber), figDir + "\bodyDescription.png");
figNumber = figNumber + 1;

disp('Analysis settings:')
fprintf('\tFaces number: %d\n', length(faces));
fprintf('\tVertices number: %d\n', length(vertices));

%% Potencial calculus
MSC = getMASCONS(faces, vertices, numberOfDivisions, figDir);
fprintf('\tMASCONS number: %d\n', length(MSC.centers));

% Points to be calculated the potencial
x = polySpaced(-maxDist, maxDist, 0, spacePoints, order_spacePoints);
y = polySpaced(-maxDist, maxDist, 0, spacePoints, order_spacePoints);
z = polySpaced(-maxDist, maxDist, 0, spacePoints, order_spacePoints);

% Get the potential
[U, a] = getPotencial(MSC, x, y, z);

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
save(varDir + "\var", 'U', 'a', 'x', 'y', 'z', ...
    'bodyName', 's');

%% My functions
function MSC = getMASCONS(faces, vertices, numberOfDivisions, figDir)
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
    
    colors = ['r', 'b', 'g', 'y', 'm', 'k', 'c'];
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
    saveas(figure(figNumber), figDir + "\MASCONS.png");
    figNumber = figNumber + 1;
    
end

function [U, a] = getPotencial(MSC, X, Y, Z)
    % Function to get the total potential of the body at a poit (x, y, z).
    % If the points ar vectors, the potencial are evaluated at each point.
    
    global G rho;
    
    U = zeros(length(X), length(Y), length(Z));
    ax = zeros(length(X), length(Y), length(Z));    
    ay = zeros(length(X), length(Y), length(Z));    
    az = zeros(length(X), length(Y), length(Z));
    mass = rho*MSC.volume;
    centers = MSC.centers;
    
    fprintf('\n');
    % For loop through the inputs points
    total = length(X)*length(Y); count = 0;
    for i = 1:length(X)
        for j = 1:length(Y)
            
            x = X(i); y = Y(j);
            parfor k = 1:length(Z)         
                % Potencial for each MASCONS divided by G, saved as vector
                r = [x; y; Z(k)] - centers;
                Uu =  mass ./ (sum(r.^2).^(1/2));
                au = (mass ./ (sum(r.^2).^(3/2)) ) .* r;

                % Total potencial (sum of each MASCONS potencial)
                U(i, j, k) = -sum(Uu);
                a = sum(au, 2);
                ax(i, j, k) = a(1);
                ay(i, j, k) = a(2);
                az(i, j, k) = a(3);
            end
            
            count = count + 1;
        end
        fprintf("\t\tProgress: %4.2f %%\n", count*100/total);
    end
    U = G*U;
    a = {G*ax, G*ay, G*az};
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