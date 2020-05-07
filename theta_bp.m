function [ X,y,Z,sigma] = theta_bp( G, max_iter, tol, X, Z, sigma)
% compute e from G
n=size(G,1);
m=nnz(G)/2;
e=zeros(m,2);
nxt=0;
for i=1:n
    for j=i+1:n
        if(G(i,j)==1)
            nxt=nxt+1;
            e(nxt,1)=i;
            e(nxt,2)=j;
        end
    end
end
%
% boundary point method, based on augmented lagrangian 
% dual version for max-clique
% solves: max <J,X> s.t. A_G*X(:) = b; X psd
%       (means x(i,j) = 0 for all [ij] in E(G), tr(X)=1;) 
% the edge list of the graph is given in e, n = max(max(e))
%
% call: [ X,y,Z,sigma] = theta_bp( e, max_iter, tol, X, Z, sigma);
%
% simplest call: X = theta_bp( e);   % just input the edge list
 
tstart = cputime; 
m = size( e,1) + 1;            % number of constraints is |E|+1
n = max(max(e));               % largest node appearing defines n
L = ones( n);                  % all ones matrix (cost function)
b = zeros( m,1); b(m) = 1;     % right hand side

% some auxiliary data
DA = 2*ones(m,1); DA(m) = n;   % A*A', which is diagonal 
Id = eye(n);                   % identity matrix
I = e(:,1); J = e(:,2);        % edgelist taken apart
m1 = m-1;                      % m1 = number of edges

% initialize
if nargin <= 3;
  if n<= 250
    sigma = .1/n;              % sigma is tuned for random graphs
  else
    sigma = .05/n; 
  end
    X = zeros(n); Z = X;
    if nargin <= 2; tol = 1e-5; end       % default tolerance
    if nargin == 1; max_iter = 1500; end   % default max-iter
end

iter = 1;                      % iteration count

% print some output
fprintf([' it    secs    dual       primal    lg10(r-d   r-p' ...
	 '    sigma) \n']);

% start iterations
while iter <= max_iter;

% determine right hand side of linear system to get y
% we need to form A(.)
  tmp = Z + L + X/sigma; 
  rhs = zeros(m,1);                     % initialize
  rhs(m) = trace(tmp);                  % last component is trace 
  tmp = tmp(:);                         % make tmp a long vector
  rhs(1:m1) = tmp( (I-1)*n + J) * 2; 
  rhs = rhs - b/sigma;                  % the final right hand side

% now compute y
  y = rhs ./DA;    % solve sigma*AA^Ty = sigma*A(L+Z) + A(X)-b

% now compute A^t(y)
  Aty = zeros(n*n,1);
  Aty( (I-1)*n + J) = y(1:m1);
  Aty = reshape(Aty, n,n);
  Aty = Aty + Aty' + y(m)*Id;

% compute W(y)
  W = Aty - L - X/sigma;  

% now compute projections to get Wp and Wn 
   [ev, lam] = eig(W); 
   lam = diag( lam); 
   Ip = find(lam> 0); j = length(Ip); rkz = j;

   if j>n/2           
      evp = zeros( n,j);
      for r=1:j;
         ic = Ip(r); evp( :,r) = ev(:,ic)*sqrt(lam(ic)); 
      end
      Wp = evp*evp';     % the projection
      Wp = (Wp+Wp')/2;   % should be symmetric
      Wn = W - Wp;
         else
      In = find(lam<0); jj = length( In); evp = zeros(n,jj);
      for r=1:jj;    
         ic = In(r); evp( :,r) = ev(:,ic)*sqrt(-lam(ic)); 
      end
      Wn = -evp*evp';    % the projection
      Wn = (Wn+Wn')/2;   % should be symmetric
      Wp = W - Wn; 
  end

% determine V and Z
   Z = Wp; V = - sigma * Wn;

% compute inner error, even though we iterate only once anyway
  tmp = V; 
  rhs = zeros(m,1);                     % initialize
  rhs(m) = trace(tmp);                  % last component is trace 
  tmp = tmp(:);                         % make tmp a long vector
  rhs(1:m1) = tmp( (I-1)*n + J) * 2; 
  res_p = rhs - b;                      % A(V)-b
  err_p = norm(res_p);                  % inner error

% update X
X = V;

% compute outer error ( relative dual error)
res_d = Z - Aty + L;
err_d = norm( res_d,'fro'); 

% compute primal and dual objective (of infeasible points)
primal = sum(sum(X)); dual = b'*y;
secs = cputime-tstart;

% print some output (every 20 iterations)
if (mod(iter,10)==0);
fprintf( '%3.0d %7.2f %11.7f %11.7f %7.2f %7.2f %8.4f \n', iter, ...
   secs, dual, primal, log10(err_d/(n+1)), log10(err_p/2), ...
	 log10(sigma) );
end

if (err_d>5000*err_p )     
  sigma=sigma*1.005; % increase sigma if dual error is much larger 
end                  % than primal error

iter = iter+ 1;
err = max( err_p, err_d/n );
if err < tol; 
 fprintf( '%3.0d %7.2f %11.7f %11.7f %7.2f %7.2f %8.4f \n', iter, ...
   secs, dual, primal, log10(err_d/(n+1)), log10(err_p/2), ...
	 log10(sigma) );
  fprintf('normal termination \n'); 
  fprintf('total time: %10.3f \n', secs); 
  return; end
end   % the stopping condition is set somewhat arbitrarily

fprintf('max iter reached \n');
curr_rel_err = max(err_d/(n+1), err_p/2);
fprintf('current error:  %11.7f,   target: %12.9f \n',curr_rel_err,tol);
fprintf('total time: %10.3f \n', secs);
