function [X, y, Z, V, primal, iter, secs] = DADAL( G_1, G_2, dAAT, b, C, maps, max_iter, tol, sigma, y, V, X)
% function [X, y, Z, V, primal, iter, secs] = DADAL( A, b, C, max_iter, tol, sigma, y, V, X)
% DADAL: a dual step for improving alternating direction 
%        augmented Lagrangian methods for SDP
%  
%         min <C,X> s.t. A*X(:) = b; X psd
%         max b'y s.t. C - A^t(y) = Z psd
%
% A, b, C: A is m by n*n, b is m-vector, C is n by n
% max_iter: max number of iterations
% sigma: penalty parameter for augmented Lagrangian 
% tol: stopping condition
%      we stop once primal and dual relative infeasibility < tol
% 
% output: primal-dual pairing (X,y,Z)
% simplest call: [ X, y, Z] = auglag( A, b, C);
% call: [ X, y, Z, V] = auglag( A, b, C, max_iter, tol, sigma, y, V, X);

% some formal checking
%[ m,n2] = size( A); 
n = size( C, 1);
m=size(b,1);
%assert( length(b)==m, 'mismatch rows of A and length of b');
%assert( n*n==n2, 'mismatch columns of A and size of C');
normb = norm( b); normC = norm( C(:));
tStart = tic;

% initialize
if nargin == 6; max_iter=1000; tol = 1e-5;  end
if nargin <= 7; tol = 1e-5; end
if nargin <= 9; Y = zeros(n); Z = Y; end

% form A*A' (and cholesky, if not part of input) 
%At = A';
%AAT = A* At;           % form A*A'
%aatd = diag(AAT);
%if nargin <= 8;
%   [R,p] = chol(diag(dAAT));
%   Rt = R';    % only for speed-up
   secs = toc(tStart);
%   fprintf(' secs after chol:   %12.5f \n', secs);
%   if p>0 
%     fprintf(' rows of A lin. dep. p = %4.0d, m = %4.0d\n',p,m);
%     error(' nothing done');
%   end
%end

if nargin <=9
% Starting rank of the dual variable V: 
rank = floor((sqrt(1+8*(n*(n+1)/2-m))-1)/2);

% random number generator initialized
rng('default'); rng(1);

%Starting V 
v0=(rand(n*rank,1)-.5);
V = reshape( v0, n, rank); Vt = V'; VVT = V*Vt; Z = VVT;

%Starting y
y0=zeros(m,1); y = y0;
else
    Z = V*V'; Y = X;
    rank = size(V,2);
end

% starting sigma
if(nargin < 7)
    G1 = VVT - C;
    normG1 = norm( G1,'fro');
    z = compute_Ay(G_1,G_2,Y(:),m,maps);
    relp = norm( b - z )^2/(1+normb); 
    sigma = relp*(1+normC)/(normG1^2);
%    fprintf('starting sigma: %10.7f \n', sigma);
end

done = 0;  % stopping condition not satisfied
iter = 0;  % iteration count
iimax = 2; % number of iterations for improved (y,V)
incsig = 0;% counts how often sigma is increased consecutively
redsig = 0;% counts how often sigma is reduced consecutively

% main loop
while done ==0  % while not done
    
% determine new pair (V, y)
    for ii=1:iimax
    [y, V, Z] = iter_AL(G_1,G_2,b,C,maps,sigma,dAAT,Y,y,V,Z);
    end    
    
% determine projection  
    Aty = reshape(compute_ATy(G_1,G_2,y,maps), n, n);    % A^T(y)
    M = Aty - C + Y/sigma;
% compute projections Mp and Mn of M
   [Mp, Mn, V, rank] = project_W(M);
  
% compute new Z and X and primal and dual errors   
   Z = -Mn; X = sigma * Mp;
   g = b-compute_Ay(G_1,G_2,X(:),m,maps);  G = -Z + C - Aty;
   err_0 = [norm(g)/(1+normb), norm(G(:))/(1+normC)];
% fprintf( '%10.4f %10.4f \n', err_0(1), err_0(2));

% update Y
    Y = X;

% check sigma update
  if (err_0(1) > 100*err_0(2))
     redsig = redsig + 1;
  elseif (err_0(1)< 2*err_0(2))
      incsig = incsig + 1;
  elseif (err_0(2)<tol*3)
      sigma = sigma*.9;
  else
      incsig = 0; redsig = 0;
  end
  if incsig>9
%      fprintf('increase sigma\n');
      sigma = sigma *1.3; incsig = 0;
  end
  if redsig>9
%      fprintf('reduce sigma\n'); 
      sigma = sigma/ 1.3; redsig=0;
  end

% check stopping conditions
iter = iter + 1;
secs = toc(tStart);
primal = C(:)'*Y(:); dual = b'*y;
done = (max(err_0)<tol) | (iter>= max_iter);

% print something
%if (mod(iter, 10)==0)
%  fprintf( '%6.0d %8.2f %13.5e %13.5e %8.3f %8.3f %9.6f \n', iter, ...
%   secs, dual, primal, log10(err_0(2)), log10(err_0(1)), sigma);

%end

end   % end while loop

if max(err_0)<tol
  fprintf(' required accuracy reached after %6.0d iterations,',iter)
else
  fprintf(' max iterations, ');
end
fprintf('    total time: %8.2f \n', secs);
end

