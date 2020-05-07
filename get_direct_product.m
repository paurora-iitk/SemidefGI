% G = G_1 * G_2
function G = get_direct_product(G_1,G_2)
n=size(G_1,1);
G=zeros(n^2);
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                if ((i==k && G_2(j,l)==1) || (j==l && G_1(i,k)==1) || (G_1(i,k)==1 && G_2(j,l)==1))
                    G(n*(i-1)+j,n*(k-1)+l)=1;
                end
            end
        end
    end
end