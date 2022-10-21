clear; tic;
%% Analysis settings

global figNumber; figNumber = 1;

bodyName = "216 Kleopatra"; % Name of the body (by the input file)
numberOfDivisions = 5; % Gauss integration order;
spacePoints = 100; % Points in space;
order_spacePoints = 2; % Order of polynomial spacing points;

% Analysis name
analysisName = bodyName + "_d" + num2str(numberOfDivisions) + "_p" + ...
    num2str(spacePoints) + "_o" + num2str(order_spacePoints);

% Creates the directory if it doesn't exist
figDir = pwd + "\fig\MASCONS_" + analysisName; checkDir(figDir);
varDir = pwd + "\var"; checkDir(varDir);

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
x = polySpaced(-252, 252, 0, spacePoints, order_spacePoints);
y = polySpaced(-252, 252, 0, spacePoints, order_spacePoints);
z = polySpaced(-252, 252, 0, spacePoints, order_spacePoints);

% Get the potential
[U, ~] = getPotencial(MSC, x, y, z);

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
save(varDir + "\MASCONS_" + analysisName, 'U', 'x', 'y', 'z', ...
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
    
    colors = ['r', 'b', 'g', 'y', 'm'];
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

function [U, Ur] = getPotencial(MSC, X, Y, Z)
    % Function to get the total potential of the body at a poit (x, y, z).
    % If the points ar vectors, the potencial are evaluated at each point.
    
    G = 6.67259e-20; % [m³/kg*s] Gravitational constant;
    rho = 4270*1e9; % [kg/m³] Mean body density;
    U = zeros(length(X), length(Y), length(Z));
    Ur = zeros(length(X), length(Y), length(Z));
    fprintf('\n');
    % For loop through the inputs points
    count = 1; total = length(X)*length(Y)*length(Z);
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
                count = count + 1;
            end
        end
        fprintf("\t\tProgress: %4.2f %%\n", count*100/total);
    end    
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