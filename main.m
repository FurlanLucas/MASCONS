clear; close all;
%%
file = importdata("216 Kleopatra.txt");
%file = importdata("1998 KY26.txt");

% Get the number of vertices
body = struct();
    body.vertices = file.data(cell2mat(file.textdata)=='v',:);
    body.faces = file.data(cell2mat(file.textdata)=='f',:);
    body.length = length(body.faces);
MSC = getMASCONS(body);

disp('Analysis settings:')
fprintf('\tFaces number: %d\n', length(body.faces));
fprintf('\tVertices number: %d\n', length(body.vertices));
fprintf('\tMASCONS number: %d\n', length(MSC.centers));

figure(1);
patch('faces', body.faces, 'vertices', body.vertices, ...
    'EdgeColor','k','FaceColor','none');
axis equal; view(30,10);

figure(2)
plot3(MSC.centers(1,:), MSC.centers(2,:), MSC.centers(3,:), 'or', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 0.8);
axis equal; view(30,10)

%% Potencial
x = polySpaced(-252, 252, 0, 800, 5);
y = polySpaced(-252, 252, 0, 800, 5);
z = 0;
U = getPotencial(MSC, x, y, z);
[X, Y] = meshgrid(x, y);

figure(3); 
surf(X, Y, U, 'EdgeColor', 'none');
colorbar(); xlabel('X [km]'); ylabel('Y [km]'); zlabel('U');
view(100,10);

%% My functions
function MSC = getMASCONS(body)
    % Function responsable for getting the centroids of each MASCONS and
    % their respective volume. The only input is the body description, un
    % struct with the filds 'vertices' and 'faces'.
    
    MSC = struct();
    A = body.vertices(body.faces(:,1), :)';
    B = body.vertices(body.faces(:,2), :)';
    C = body.vertices(body.faces(:,3), :)';
    
    MSC.centers = (A + B + C + zeros(3, body.length))/4;
    MSC.volume = dot(C, cross(A, B))/6;
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