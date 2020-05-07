function [Wp, Wn, V, rank] = project_W( W)
% call: [Wp, Wn, V, rank] = project_W( W);

  M = (W + W')/2;      % should be symmetric
n = size( M, 1);
  
% now compute eigenvalue decomposition of W  
  [ev, lam] = eig(M); lam = diag(lam);

% compute projection Mp and Mn of M onto PSD and NSD
  I = find(lam> 0); j = length(I);
  if j < n/2;
     evp = zeros( n,j);
     for r=1:j;
         ic = I(r); evp(:,r) = ev(:,ic)*sqrt(lam(ic)); 
     end
     if j==0; evp = zeros(n,1); end
     evpt= evp';
     Mp = evp*evpt;     
     %%%     Computation of V         %%%
     In = find(lam < 0); jn = length(In);
     rank = jn;
     pn = zeros(n,jn);
     for r=1:jn;
        icn = In(r); pn(:,r) = ev(:,icn)*sqrt(-lam(icn)); 
     end
     V = -pn;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Mn = M - Mp;
  else
     I = find(lam < 0); j = length( I); % should be <= n/2
     evn = zeros( n,j);
     for r=1:j;
        ic = I(r); evn(:,r) = ev(:,ic)*sqrt( -lam(ic)); 
     end
     if j==0; evn = zeros(n,1); end
     evnt= evn';
     Mn = -evn*evnt;
     V = -evn;
     rank = j;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Mp = M - Mn; 
  end

  Wn = (Mn+Mn')/2; Wp = (Mp+Mp')/2;  % these should be symm.

