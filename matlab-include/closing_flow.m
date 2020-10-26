function [U,G,data] = closing_flow(V,F,varargin)





%%%%%%%%%%%%%%%%%%%%%%%%
% Alec's way of defining optional parameters and stuff. Do not touch
% unless you know what you're doing please!
%

% DEFAULTS:
maxiter = 100;
dt = 1;
bd = 1/.4;
h = .05;
remesh_iterations = 1;
debug = false;
make_plot = true;
write = false;
self_intersect = false;
data = [];
bd_fun = [];
normal_distortion = [];
opening = false;
gt = [];
always_recompute = false;
is_active =[];
quadric_curvatures = false;
fixed_function = @(V) [];
%closing = true;
tol = 1e-7;

% ASSIGNED:
params_to_variables = containers.Map({'EdgeLength','TimeStep','Bound','MaxIter',...
    'RemeshIterations','Debug','Plot','Write','SelfIntersect',...
    'BoundFunction','NormalDistortion','Opening','Groundtruth','AlwaysRecompute','IsActive','QuadricCurvatures',...
    'FixedFunction','Tolerance'},...
    {'h','dt','bd','maxiter','remesh_iterations','debug','make_plot','write',...
    'self_intersect','bd_fun','normal_distortion','opening','gt','always_recompute','is_active','quadric_curvatures',...
    'fixed_function','tol'});
v = 1;
while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
        assert(v+1<=numel(varargin));
        v = v+1;
        % Trick: use feval on anonymous funtion to use assignin to this workspace
        feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
        error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%

Vfull = V;
Ffull = F;

if ~isempty(gt)
    data.sqrd = [];
end
if isempty(is_active)
    is_active = @(K) K<-bd;
end
active_size = [];
recompute = true;
iter = 0;
data.active_num = 0;
while iter<maxiter
    
    
    Vprev = Vfull;
    Fprev = Ffull;
    
    
    
    % UPDATE ACTIVE REGION
    
    %recompute = true;
    % UNRAVELING DISCRETE_CURVATURES WITHOUT RECOMPUTING STUFF
    if iter==1 || recompute
        A = adjacency_matrix(Ffull) + speye(size(Vfull,1));
        M = massmatrix(Vfull,Ffull);
        K = discrete_gaussian_curvature(Vfull,Ffull);
        H = discrete_mean_curvature(Vfull,Ffull);
        k = H + [1 -1].*sqrt(H.^2 - M*K);
        if opening
            K = -real(k(:,1)./diag(M));
        else
            K = real(k(:,2)./diag(M));
        end
        if quadric_curvatures
            [K,N,T] = per_vertex_prin_curvature(Vfull,Ffull,'Type','quadric');
            if opening
                K = max(K,[],1);
            else
                K = min(K,[],2);
            end
        end
        if isempty(bd_fun)
            if isempty(normal_distortion)
                moving = find(is_active(K));
            else
                normals = per_vertex_normals(Vfull,Ffull);
                normals = normals./normrow(normals);
                moving = find(K<(-bd.*abs(dot(repmat(normal_distortion,...
                    size(normals,1),1),normals,2))));
            end
        else
            moving = find(K<-bd_fun(Vfull));
        end
        if ~always_recompute
            recompute = false;
        end
    end
    
    
    if size(moving,1)==0
        disp("Active set is empty");
        break
    end
    active = find(sum(A(moving,:))>0);
    active = find(sum(A(active,:))>0); % Tworing
    %fixed_test = setdiff(active,moving);
    % % For two-ring?
    if always_recompute
        active = 1:size(Vfull,1)';
    end
    % % FOR DEBUGGING
    f_active = Ffull(sum(ismember(Ffull,active),2)>1,:);
    if size(f_active,1)==0
        disp("Active set is empty");
        break
    end
    %[v_active,I,J] = remove_unreferenced(Vfull,f_active);
    [I,J,f_active,v_active] = output_sensitive_remove_unreferenced(f_active,Vfull);
    fixed_test = setdiff(1:size(v_active,1),I(moving));
    %f_active = I(f_active);
    f_inactive = Ffull(~(sum(ismember(Ffull,active),2)>1),:);
    %[v_inactive,Ii,Ji] = remove_unreferenced(Vfull,f_inactive);
    if size(f_inactive,1)==0
        Ii = [];
        Ji = [];
        v_inactive = [];
        f_inactive = [];
    else
        [Ii,Ji,f_inactive,v_inactive] = output_sensitive_remove_unreferenced(f_inactive,Vfull);
    end
    V = v_active;
    active_size = [active_size,size(V,1)];
    %plot(active_size)
    %drawnow
    %save('thing.mat','active_size')
    F = f_active;
    data.active_num = data.active_num + size(V,1);
    % END OF UPDATE ACTIVE REGION
    
    %M = massmatrix(V,F);
    M = M(J,J);
    dblA = diag(sparse(doublearea(V,F))/2);
    v_xyz = [V(:,1);V(:,2);V(:,3)];
    
    %     [k,~,~,~] = discrete_curvatures(V,F);
    %     K = real(k(:,2)./diag(M));
    O = outline(F);
    boundary_verts = unique(O(:));
    
    fixed = unique([fixed_test';boundary_verts;fixed_function(V)]);
    
    %     fixed = unique([fixed_test';boundary_verts;find(V(:,3)<0.1)]);
    %     warning('THING');
    %K(boundary_verts) = Inf;
    %    fixed = find(K>-bd);
    %     fixed_full = ones(size(Vfull,1),1);
    %     fixed_full(J') = K>-bd;
    if debug
        if ~is_edge_manifold(F)
            error("NOT EDGE MANIFOLD");
            break
        end
    end
    
    fixed_xyz = [fixed;fixed+size(V,1);fixed+2*size(V,1)];
    if debug
        disp("Entering curvature")
        save('debug_curvature.mat')
        writeOBJ('debug_curvature.obj',V,F);
    end
    [PD1,PD2,~,~]=per_face_prin_curvature_mex(V,F);
    if opening
        PD1 = PD2;
    end
    proy_matrix = [sparse(1:length(PD1(:,1)),1:length(PD1(:,1)),PD1(:,1)),...
        sparse(1:length(PD1(:,2)),1:length(PD1(:,2)),PD1(:,2)),...
        sparse(1:length(PD1(:,3)),1:length(PD1(:,3)),PD1(:,3))];
    proyected_gradient = proy_matrix*grad(V,F);
    int_proy_grad_sq = -proyected_gradient'*dblA*proyected_gradient;
    Q = blkdiag(M,M,M)-0.01.*dt*blkdiag(int_proy_grad_sq,int_proy_grad_sq,int_proy_grad_sq);
    linear = -blkdiag(M,M,M)*v_xyz;
    uu = min_quad_with_fixed(.5*real(Q),real(linear),fixed_xyz,v_xyz(fixed_xyz));
    U = [uu(1:size(V,1)),uu(size(V,1)+1:2*size(V,1)),...
        uu(2*size(V,1)+1:3*size(V,1))];
    %V = V + dot((U-V),per_vertex_normals(V,F),2).*per_vertex_normals(V,F);
        [U,F] = remesh_botsch_mex(U,F,fixed,h.*ones(size(U,1),1),remesh_iterations);
    
    
    
    % Reassign V
    if debug
        save('debug_remesher_before.mat','U','F','fixed','h');
        disp("Entering remesher")
        writeOBJ('debug_remesher_before.obj',U,F);
        writeOBJ('debug_remesher_fixed_points.obj',[fixed zeros(length(fixed),2)],[1 2 3]);
    end
    %[U,F] = remesh_botsch_mex(U,F,fixed,h.*ones(size(U,1),1),remesh_iterations);
    Udup = U;
    %warning('DEAL WITH THE DUPLICATION THING');
    %U = U + 1e-10.*rand(size(U,1),3);
    if debug
        save('debug_remesher_after.mat','U','F','fixed','h');
        disp("Exiting remesher")
        writeOBJ('debug_remesher_after.obj',U,F);
        s = statistics(U,F);
        save('before_pasting.mat')
        if s.num_nonmanifold_edges>0
            error('NON MANIFOLD AFTER REMESHER');
            break
        end
    end
    O = outline(F);
    boundary_verts = unique(O(:));
    interior_verts = setdiff(1:size(U,1),boundary_verts);
    Udup(interior_verts,:) = U(interior_verts,:)+1e-8.*rand(length(interior_verts),3);
    Vfull = [v_inactive;Udup];
    Ffull = [f_inactive;F+size(v_inactive,1)];
    [Vfull,SVI,SVJ] = remove_duplicate_vertices(Vfull,0);
    Ffull = SVJ(Ffull);
    interior_active_full = SVJ(interior_verts+size(v_inactive,1));
    H_interior_active = discrete_mean_curvature(U,F);
    K_interior_active = discrete_gaussian_curvature(U,F);
    Mnew = sparse(1:size(Vfull,1),1:size(Vfull,1),1);
    M_active = massmatrix(U,F);
    m_active = diag(M_active);
    Mnew = Mnew + ...
        sparse(interior_active_full,interior_active_full,m_active(interior_verts)-1,size(Vfull,1),size(Vfull,1));
    M = Mnew;
    k_interior_active = H_interior_active + ...
        [1 -1].*sqrt(H_interior_active.^2 - M_active*K_interior_active);
    if opening
        K_interior_active = -real(k_interior_active(:,1)./diag(M_active));
    else
        K_interior_active = real(k_interior_active(:,2)./diag(M_active));
    end
    if quadric_curvatures
        [PD1,PD2,PV1,PV2]=principal_curvature(U,F);
        if opening
            K_interior_active = -PV2;
        else
            K_interior_active = PV1;
        end
    end
    moving = zeros(size(Vfull,1),1);
    if isempty(bd_fun)
        if isempty(normal_distortion)
            moving_interior_active = is_active(K_interior_active);
        else
            normals = per_vertex_normals(U,F);
            normals = normals./normrow(normals);
            moving_interior_active = K_interior_active<(-bd.*abs(dot(repmat(normal_distortion,...
                size(normals,1),1),normals,2)));
        end
    else
        moving_interior_active = K_interior_active<-bd_fun(U);
    end
    moving(interior_active_full) = moving_interior_active(interior_verts);
    moving = find(moving);
    
    E = SVJ(edges(F+size(v_inactive,1)));
    A = sparse([E(:,1) E(:,2)],[E(:,2) E(:,1)],1,size(Vfull,1),size(Vfull,1));
    A = A + speye(size(Vfull,1));
    
    
    %     target = h.*ones(size(Vfull,1),1);
    %     feature = unique([setdiff(1:size(Vfull,1),active(setdiff(1:length(active),fixed)))]);
    %     [Vfull,Ffull] = remesh_botsch_mex(Vfull,Ffull,feature,target,10);
    %[V,F] = remesh_botsch_mex(Vfull,Ffull,find(fixed_full),h.*ones(size(Vfull,1),1),1);
    %Vfull = V;
    %Ffull = F;
    %     [V,F] = remesh_morphology(V,F,h,...
    %         1.5*bd,1,false,true,[]);
    if ~isempty(gt)
        data.sqrd = [data.sqrd,max(point_mesh_squared_distance(Vfull,gt.V,gt.F))];
        save('data.mat','data');
    end
    %
    if debug || make_plot
        tsurf(F,U,'FaceColor',[189,235,252]./255,'EdgeAlpha',0)
        hold on
        if size(f_inactive,1)>0
            tsurf(f_inactive,v_inactive,'FaceAlpha',0.5,'FaceColor',[.8 .8 .8],'EdgeAlpha',0)
        end
        tsurf(Ffull,Vfull + [max(max(Vfull(:,1))),0,0],'EdgeAlpha',1,fsoft,fphong,'FaceColor',[189,235,252]./255)
        axis equal
        camlight
        view([0 30])
        %view([0 90])
        drawnow
        %figgif('~/Downloads/release_the_kraken.gif')
        hold off
    end
    if self_intersect && mod(iter,10)==0
        disp("Detecting self-intersections...")
        [SV,SF,~,~,IM] = selfintersect(Vfull,Ffull,'StitchAll',true);
        FF = IM(SF);
        [Vfull,IM] = remove_unreferenced(SV,FF);
        Ffull = IM(FF);
        [Vfull,Ffull] = mesh_boolean(Vfull,Ffull,[],[],'union');
        recompute = true;
        disp("...done!")
    end
    if write && mod(iter,2)==0
        writeOBJ(['closing_written_',num2str(iter),'.obj'],Vfull,Ffull);
        disp(iter)
    end
    %disp(iter)
    iter = iter + 1;
    if mod(iter,10)==0
        dist = max(point_mesh_squared_distance(Vfull,Vprev,Fprev));
        if dist<tol
            disp('converged')
            break
        end
        disp(dist)
    end
end

data.active_num = data.active_num /iter;


U = Vfull;
G = Ffull;



end