addpath(genpath('../../matlab-include')) % path to functions
[V,F] = read_triangle_mesh('../../data/fandisk.obj'); % read input
% Set parameters
a = V(:,3);
V(:,3) = V(:,2);
V(:,2) = a;
F = [F(:,2) F(:,1) F(:,3)];
V = V-min(min(V));
V = V./(max(max(V)));
h = 0.002;
bd = 1/.015;
%bd_fun = @(V) 3./(1-sqrt(V(:,1)));
dt = 0.001;
% Parameters set
writeOBJ('fandisk_input.obj',V,F) % write input
[U,G] = closing_flow(V,F,'Bound',bd,'EdgeLength',h,'TimeStep',dt,...
    'MaxIter',300,'RemeshIterations',1,'Debug',false,'Plot',true,...
    'Write',false,'Opening',false);
writeOBJ('fandisk_output.obj',U,G)% write output

% We've already saved input and output. In order to render them with the
% moving part highlighted, we'll do the following to separate the output
% into an "active" part and an "inactive" one. We then render them as in
% ../../render/render-template.blend

clc; clear all; close all;
[Vgt,Fgt] = readOBJ('fandisk_input.obj');
[V,F] = read_triangle_mesh('fandisk_output.obj');
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
writeOBJ('fandisk_active.obj',v_active,f_active);
writeOBJ('fandisk_inactive.obj',v_inactive,f_inactive);