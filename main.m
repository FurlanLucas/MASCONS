clear; close all;
%%
file = importdata("1998 KY26.txt");

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
plot3(MSC.centers(:,1), MSC.centers(:,2), MSC.centers(:,3), 'or', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 0.8);
axis equal; view(30,10)




%% My functions
function out = getMASCONS(body)
    out = struct();
    out.centers = (body.vertices(body.faces(:,1), :) + ...
                   body.vertices(body.faces(:,2), :) + ...
                   body.vertices(body.faces(:,3), :) + ...
                   zeros(body.length, 3))/4;
end
