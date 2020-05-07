function [X,y,Z,val]=algo1(G_1,G_2)

n=size(G_1,1);
C=eye(n^2+1);
C(n^2+1,n^2+1)=0;
C=-1.*C;

maps=-ones(n,1);
val=0;
%for i=2:n
	maps(1)=1; 
	% maps(2)=2; 
	[dAAT,b]=compute_daat(G_1,G_2,maps); 
	[X,y,Z]=DADAL(G_1,G_2,dAAT,b,C,maps);
	if (trace(X)>val)
		val=trace(X);
	end
%end
%[X,y,Z]=mprw2(dAAT,b,C,G_1,G_2,maps);

