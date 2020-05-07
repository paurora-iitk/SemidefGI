
J=ones(size(G_1,1));
I=eye(size(G_1,1));
H=get_direct_product(J-I-G_1,G_2);
L=get_direct_product(G_1,J-I-G_2);
G=L+H;
B=J-I-G_2;
J=ones(size(H,1));  
Y=J-G;
d=nnz(G_1(1,:));
n=size(G_1,1);
A=(1/d).*G_2;
B=(1/(n-1-d)).*B;

for i=1:n  
for j=i+1:n  
if(G_1(i,j)==1)
Y(n*(i-1)+1:n*(i-1)+n,n*(j-1)+1:n*(j-1)+n)=A;
else
Y(n*(i-1)+1:n*(i-1)+n,n*(j-1)+1:n*(j-1)+n)=B;
end;
end;
end;
for i=1:n^2
for j=1:i-1  
Y(i,j)=Y(j,i);                               
end;
end;

