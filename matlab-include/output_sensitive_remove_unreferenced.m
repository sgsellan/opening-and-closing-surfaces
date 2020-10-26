function [I,U,IF,VU] = output_sensitive_remove_unreferenced(F,n)
% OUTPUT_SENSITIVE_REMOVE_UNREFERENCED Remove any rows in (1:n or V) that are
% not referenced in F. If (V,F) is the original mesh, then (VU,IF) is the new
% mesh.
%
% [I,U,IF] = output_sensitive_remove_unreferenced(F,n)
%   or
% [I,U,IF,VU] = output_sensitive_remove_unreferenced(F,V)
%
% Inputs:
%   F  #F by ss list of element indices into some set of n vertices
%   n  number of vertices in original set
%     or
%   V  n by dim list of vertex positions
% Outputs:
%   I  #V by 1 sparse matrix where I(j) = i ≠ 0 if vertex j in V becomes vertex i
%     in VU. 
%   U  #VU by 1 list of indices so that U(i) = j indicates that vertex j in V
%     becomes vertex i in VU
%   IF  #F by ss list of reindexed faces so that IF = I(F);
%   VU  #VU by dim list of new vertex positions: VU = V(U,:);
%
% Example:
%   
%   tic;[I,U,IF,VU]=output_sensitive_remove_unreferenced(F(end,:),V);toc;
%   tic;[RV,IM,J,IMF]=remove_unreferenced(V,F(end,:));toc;
%   assert(isequal(RV,VU));
%   assert(isequal(IM,sparse(find(I==0),1,max(I)+(1:sum(I==0)),size(I,1),1)+I));
%   assert(isequal(J,U'));
%   assert(isequal(IMF,IF));
%
  if nargout>3
    assert(size(n,1)>1);
    V = n;
    n = size(n,1);
  end
  I = [];
  J = [];
  IF = [];

  % get list of unique vertex indices that occur in faces
  %U = unique(F(:));
  % Slightly faster unique if we don't have infs or nans
  sF = sort(F(:));
  sI = [true;diff(sF)~=0];
  U = sF(sI);
  % allocate space for an indexmap so that I[i] gives new index of vertex i
  I = sparse(U,1,1:size(U,1),n,1);
  if nargout>2
    IF = full(I(F));
  end
  if nargout>3
    VU = V(U,:);
  end
end
