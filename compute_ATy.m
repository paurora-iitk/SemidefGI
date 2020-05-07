function [z]=compute_ATy(G_1,G_2,y,maps)

n=size(G_1,1);
Z=zeros(n^2+1);
nxt=1;
% constraints X(i,i)-X(i,n^2+1)=0
for i=1:n
	for j=1:n
		if (maps(i)==j)
			Z(n*(i-1)+j,n*(i-1)+j)=y(nxt,1);
			nxt=nxt+1;
			Z(n*(i-1)+j,n^2+1)=y(nxt,1);
			Z(n^2+1,n*(i-1)+j)=y(nxt,1);
		else
			Z(n*(i-1)+j,n*(i-1)+j)=2*y(nxt,1); 
			Z(n*(i-1)+j,n^2+1)=-y(nxt,1); 
			Z(n^2+1,n*(i-1)+j)=-y(nxt,1);
		end
		nxt=nxt+1;
	end
end
% constraint X(n^2+1,n^2+1)=1
Z(n^2+1,n^2+1)=y(nxt,1);
nxt=nxt+1;
% add the graph constraints
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                if (((i~=k && j==l) || (i==k && j~=l)) && (n*(i-1)+j < n*(k-1)+l))
                    Z(n*(i-1)+j,n*(k-1)+l)=y(nxt,1);
                    Z(n*(k-1)+l,n*(i-1)+j)=y(nxt,1);
                    nxt=nxt+1;
                elseif (((G_1(i,k)==1 && G_2(j,l)==0) || (G_1(i,k)==0 && G_2(j,l)==1)) && (n*(i-1)+j < n*(k-1)+l))
                    Z(n*(i-1)+j,n*(k-1)+l)=y(nxt,1);
                    Z(n*(k-1)+l,n*(i-1)+j)=y(nxt,1);
                    nxt=nxt+1;
                end
            end
        end
    end
end
z=Z(:);
