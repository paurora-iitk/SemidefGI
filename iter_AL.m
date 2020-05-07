function [ynew, Vnew, VVtnew] = iter_AL(G_1,G_2,b,C,maps,sigma,dAAT,Y,y,V,VVt)
% do one inner iteration
% data: A,b,C,sigma,R,Rt,Y
% input:  y,V,VVt
% output: ynew, Vnew, VVtnew
% call: [ynew, Vnew, VVtnew] = iter_AL(A,b,C,sigma,R,Rt,Y,y,V,VVt);

% compute gradient dV
n = size(C,1);
m=size(b,1);
%At = A';
Aty = reshape(compute_ATy(G_1,G_2,y,maps), n, n);    % A^T(y)
G = VVt - C + Aty;
dV = -2*(Y + sigma*G)*V;      % this is a matrix

% compute second order diagonal scaling hvv
r = size( V,2);
er = ones(r,1);
en = ones(n,1);
M = Y + sigma*(-C + Aty);
Md = max(diag(M),0);
d1 = diag(VVt);
d2 = sum( V.*V); d2 = d2';
kr1 = d1*er'; kr1 = kr1(:); 
kr2 = en*d2'; kr2 = kr2(:); 
kr3 = Md*er'; kr3 = kr3(:);
% hvv2 =  2*(sigma*(V(:).*V(:)+kron(er,d1)+kron(d2,en))+kron(er,Md));
hvv2 = 2*(sigma*(V(:).*V(:) + kr1 + kr2) + kr3);

% search direction wrt to V is dv = dV(:)./hvv2;
dv = dV(:)./hvv2;    % this is a vector
%keyboard
dv = reshape( dv, n, r);


% compute optimal y for dV (y0, y1, y2)
% first do expensive matrix multiplications
DV = dv * V';
DD = dv * dv';

% prepare rhs to solve for y0,y1
Mtmp = Y/sigma - C + VVt;
rs1 = b/sigma - compute_Ay(G_1,G_2,Mtmp(:),m,maps);
DVVt = DV + DV';
rs2 = -compute_Ay(G_1,G_2,DVVt(:),m,maps); 
% now solve
%tmp  = Rt\[rs1 rs2];
%sol = R\tmp;
for i=1:size(dAAT,1)
	y0(i,1)=rs1(i,1)/dAAT(i,1);
	y1(i,1)=rs2(i,1)/dAAT(i,1);
end
%y0 = sol(:,1); 
%y1 = sol(:,2);
y2 = -compute_Ay(G_1,G_2,DD(:),m,maps);

% evaluate at 5 points
alpha = [0 .25 .5 .75 1];
for i = 1:5
    ai = alpha(i);
    yp = y0 + ai*y1 + (ai*ai)*y2;
    Vp = VVt + ai * DVVt + (ai*ai) * DD;
    Aty = reshape(compute_ATy(G_1,G_2,yp,maps), n, n);
    Mp = Vp  - C + Aty; 
    Lag(i) = b'*yp -Mp(:)'*Y(:) -(sigma/2)*Mp(:)'*Mp(:);
end


% interpolation for (alpha, Lag)
P = polyfit( alpha, Lag, 4);

% now do search over interval
xx = [0: 0.0001 : 10];
yy = polyval( P, xx);

% now take maximum
[ymax, imax] = max(yy);

% update V and y
step = xx(imax);
Vnew = V + step * dv;
VVtnew = VVt + step * DVVt + step*step*DD;
ynew = y0 + step*y1 + step*step*y2;


% print some info
% normgrad = sqrt(dV(:)'*dV(:));
% fprintf(' %12.4f %12.4f %6.3f \n', ymax, normgrad, step);
