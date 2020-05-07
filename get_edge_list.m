function e=get_edge_list(G_1,G_2)

n=size(G_1,1);
H=zeros(n^2);
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                if((i==k && j~=l)||(i~=k && j==l)||(G_1(i,k)==1 && G_2(j,l)==0)||(G_1(i,k)==0 && G_2(j,l)==1))
                    H(n*(i-1)+j,n*(k-1)+l)=1;
                end
            end
        end
    end
end
%--------------------------------------------------------------------------
% adding edges to force the vertex labeled 1 in G_1 to map to the vertex labeled 1
% in G_2
i=1;
k=2;
for j=2:n
    for l=1:n
        H(n*(i-1)+j,n*(k-1)+l)=1;
    end
end
%--------------------------------------------------------------------------
m=nnz(H)/2;
e=zeros(m,2);
nxt=0;
for i=1:n^2
    for j=i+1:n^2
        if(H(i,j)==1)
            nxt=nxt+1;
            e(nxt,1)=i;
            e(nxt,2)=j;
        end
    end
end