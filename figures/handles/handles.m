addpath(genpath('../../matlab-include')) % path to functions
[V,F] = read_triangle_mesh('../../data/handle.obj'); % read input
V = V-min(V);
V = V./(max(max(V)));
%Set parameters
h = 0.005;
bd = 1/.08;
is_active = [];
dt = 0.001;
% Any point that makes this function "true" will not move
fixed_function = @(V) find((normrow(V(:,1:2)-[0.04,0.12])<0.02) + ...
    (normrow(V(:,1:2)-[0.22,0.12])<0.02) + ...
    (normrow(V(:,1:2)-[0.12,0.20])<0.02) + ...
    (normrow(V(:,1:2)-[0.12,0.04])<0.02) + ...
    (normrow(V(:,1:2)-[0.80,0.12])<0.02) + ...
    (normrow(V(:,1:2)-[0.96,0.12])<0.02) + ...
    (normrow(V(:,1:2)-[0.88,0.04])<0.02) + ...
    (normrow(V(:,1:2)-[0.88,0.20])<0.02));
writeOBJ('handle_input.obj',V,F); % write input
[U,G] = closing_flow(V,F,'Bound',bd,'EdgeLength',h,'TimeStep',dt,...
    'MaxIter',300,'RemeshIterations',1,'Debug',false,'Plot',true,...
    'Write',false,'FixedFunction',fixed_function); % run method
writeOBJ('handle_output.obj',U,G); % write output

% We've already saved input and output. In order to render them with the
% moving part highlighted, we'll do the following to separate the output
% into an "active" part and an "inactive" one. We then render them as in
% ../../render/render-template.blend

clc; clear all; close all;
[Vgt,Fgt] = readOBJ('handle_input.obj');
[V,F] = read_triangle_mesh('handle_output.obj');
[sqrD,I,C] = point_mesh_squared_distance(V,Vgt,Fgt);
A = adjacency_matrix(F) + speye(size(V,1));
moving = find(double(sqrD>1e-6));
%active = find(sum(A(moving,:))>0);
active = moving;
f_active = F(sum(ismember(F,active),2)>2,:);
[I,J,f_active,v_active] = output_sensitive_remove_unreferenced(f_active,V);
f_inactive = F(~(sum(ismember(F,active),2)>2),:);
[Ii,Ji,f_inactive,v_inactive] = output_sensitive_remove_unreferenced(f_inactive,V);
hold off
tsurf(f_active,v_active,'FaceColor',[189,235,252]./255,'EdgeAlpha',0)
hold on
tsurf(f_inactive,v_inactive,'FaceAlpha',0.5,'FaceColor',[.8 .8 .8],'EdgeAlpha',0)
axis equal
camlight
drawnow
%pause
writeOBJ('handle_active.obj',v_active,f_active);
writeOBJ('handle_inactive.obj',v_inactive,f_inactive);
