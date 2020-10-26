addpath(genpath('../../matlab-include')) % path to functions
[V,F] = read_triangle_mesh('../../data/nefertiti.obj'); % read input
V = V-min(min(V));
V = V./(max(max(V)));
h = 0.005;
bd = 1/.015;
dt = 0.1;
writeOBJ('nefertiti_input.obj',V,F);
[U,G] = closing_flow(V,F,'Bound',bd,'EdgeLength',h,'TimeStep',dt,...
    'MaxIter',40,'RemeshIterations',1,'Debug',false,'Plot',true,...
    'Write',false,'Opening',false,'AlwaysRecompute',false);
writeOBJ('nefertiti_output.obj',U,G);

% We've already saved input and output. In order to render them with the
% moving part highlighted, we'll do the following to separate the output
% into an "active" part and an "inactive" one. We then render them as in
% ../../render/render-template.blend

clc; clear all; close all;
[Vgt,Fgt] = readOBJ('nefertiti_input.obj');
[V,F] = read_triangle_mesh('nefertiti_output.obj');
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
writeOBJ('nefertiti_active.obj',v_active,f_active);
writeOBJ('nefertiti_inactive.obj',v_inactive,f_inactive);