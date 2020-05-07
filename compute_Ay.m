function [z]=compute_Ay(G_1,G_2,y,m,maps)

n=size(G_1,1);
N=n^2+1;
Y=reshape(y,N,N);
z=zeros(m,1);
nxt=1;
% constraints X(i,i)-X(i,n^2+1)=0
for i=1:n
	for j=1:n
		if (maps(i)==j)
			z(nxt,1)=Y(n*(i-1)+j,n*(i-1)+j);
			nxt=nxt+1;
			z(nxt,1)=Y(n*(i-1)+j,n^2+1)+Y(n^2+1,n*(i-1)+j);
		else
			z(nxt,1)=2*Y(n*(i-1)+j,n*(i-1)+j)-Y(n*(i-1)+j,n^2+1)-Y(n^2+1,n*(i-1)+j);
		end
		nxt=nxt+1;
	end
end
% constraint X(n^2+1,n^2+1)=1
z(nxt,1)=Y(n^2+1,n^2+1);
nxt=nxt+1;
% add the graph constraints
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                if (((i~=k && j==l) || (i==k && j~=l)) && (n*(i-1)+j < n*(k-1)+l))
                    z(nxt,1)=Y(n*(i-1)+j,n*(k-1)+l)+Y(n*(k-1)+l,n*(i-1)+j);
                    nxt=nxt+1;
                elseif (((G_1(i,k)==1 && G_2(j,l)==0) || (G_1(i,k)==0 && G_2(j,l)==1)) && (n*(i-1)+j < n*(k-1)+l))
                    z(nxt,1)=Y(n*(i-1)+j,n*(k-1)+l)+Y(n*(k-1)+l,n*(i-1)+j);
                    nxt=nxt+1;
                end
            end
        end
    end
end

