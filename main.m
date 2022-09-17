clear; close all; tic;
%% Analysis settings
global figNumber figDir; 
bodyName = "216 Kleopatra"; % Name of the body (by the input file)
numberOfDiv = 5; % Number of MASCONS par tetrahedral;

figDir = pwd + "\fig\";
if ~isfolder('fig')
    mkdir(pwd + "\fig")
end

figNumber = 1;

%% Main calculations
fprintf('Stating analysis for %s body.\n', bodyName);
file = importdata(bodyName + ".txt");

% Get the number of vertices
body = struct();
    body.vertices = file.data(cell2mat(file.textdata)=='v',:);
    body.faces = file.data(cell2mat(file.textdata)=='f',:);
    body.length = length(body.faces);
    
% Figura do corpo
figure(figNumber);
patch('faces', body.faces, 'vertices', body.vertices, ...
    'EdgeColor','k','FaceColor','none');
axis equal; view(30,10);
saveas(figure(figNumber), figDir + "bodyDescription.png");
figNumber = figNumber + 1;

disp('Analysis settings:')
fprintf('\tFaces number: %d\n', length(body.faces));
fprintf('\tVertices number: %d\n', length(body.vertices));

% Building MASCONS data
MSC = getMASCONS(body, 3);
fprintf('\tMASCONS number: %d\n', length(MSC.centers));


%% Potencial
x = polySpaced(-252, 252, 0, 400, 5);
y = polySpaced(-252, 252, 0, 400, 5);
z = 0;
U = getPotencial(MSC, x, y, z);
[X, Y] = meshgrid(x, y);

% Plot potential data
figure(3); 
surf(X, Y, U, 'EdgeColor', 'none');
colorbar(); xlabel('X [km]'); ylabel('Y [km]'); zlabel('U');
view(60, 35);

%% Ending
s = toc;
fprintf('Code ended successfully');
fprintf('Total execution time: %dmin e %.2fs.\n', floor(s/60), ...
    s - floor(s/60)*60);
%close all;

%% My functions
function MSC = getMASCONS(body, numberOfDivisions)
    % Function responsable for getting the centroids of each MASCONS and
    % their respective volume. The only input is the body description, un
    % struct with the filds 'vertices' and 'faces'.
    
    % Inputs
    global figNumber figDir;
    
    if ~exist('numberOfDivisions', 'var')
        numberOfDivisions = 1;
    end
    
    L = numberOfDivisions*body.length;
    MSC = struct();
        MSC.centers = zeros(3, L);
        MSC.volume = zeros(1, L);
    
    colors = ['r', 'b', 'g', 'y', 'm'];
    figure(figNumber); hold on;
    
    % General calculations
    A = body.vertices(body.faces(:,1), :)'/numberOfDivisions;
    B = body.vertices(body.faces(:,2), :)'/numberOfDivisions;
    C = body.vertices(body.faces(:,3), :)'/numberOfDivisions;
    
    MSC.centers(:, 1:body.length) = ...
        (A + B + C + zeros(3, body.length))/4;
    MSC.volume(1:body.length) = dot(C, cross(A, B))/6;
    
    plot3(MSC.centers(1,:), MSC.centers(2,:), MSC.centers(3,:), 'o', ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 0.8);
    
    for i = 2:numberOfDivisions
        Anew = body.vertices(body.faces(:,1), :)'*i/numberOfDivisions;
        Bnew = body.vertices(body.faces(:,2), :)'*i/numberOfDivisions;
        Cnew = body.vertices(body.faces(:,3), :)'*i/numberOfDivisions;
        
        MSC.centers(:, (1:body.length) + (i-1)*L) = ...
            (A + B + C + Anew + Bnew + Cnew)/6;
        
        MSC.volume(:, (1:body.length) + (i-1)*L) = ...
            dot(Cnew, cross(Anew, Bnew))/6 - dot(C, cross(A, B))/6;
        
        % Plot the MASCONS centers
        plot3(MSC.centers(1, (1:body.length) + (i-1)*L), ...
              MSC.centers(2, (1:body.length) + (i-1)*L), ...
              MSC.centers(3, (1:body.length) + (i-1)*L), 'o', ...
              'MarkerFaceColor', colors(i-1), 'MarkerEdgeColor', ...
              colors(i-1), 'MarkerSize', 0.8);
        
        A = Anew; B = Bnew; C = Cnew;
    end
    
    % Plots
    
    axis equal; view(30,10);
    saveas(figure(figNumber), figDir + "MASCONS.png");
    figNumber = figNumber + 1;
    
end

function U = getPotencial(MSC, X, Y, Z)
    % Function to get the total potential of the body at a poit (x, y, z).
    % If the points ar vectors, the potencial are evaluated at each point.
    
    G = 6.67259e-20; % [m³/kg*s] Gravitational constant;
    rho = 4270*1e9; % [kg/m³] Mean body density;
    U = zeros(length(X), length(Y), length(Z));
    
    % For loop through the inputs points
    for i = 1:length(X)
        for j = 1:length(Y)
            for k = 1:length(Z)
                % Potencial for each MASCONS divided by G, saved as vector
                r = [X(i); Y(j); Z(k)];
                Uu = MSC.volume*rho ./ sqrt(sum((MSC.centers - r).^2));

                % Total potencial (sum of each MASCONS potencial)
                U(i, j, k) = -G*sum(Uu);
            end
        end
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