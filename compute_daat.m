function [dAAT,b]=compute_daat(G_1,G_2,maps)

n=size(G_1,1);

% determine the total number of constraints
cnt=0;
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                if ((i~=k && j==l) || (i==k && j~=l))
					cnt=cnt+1;
                elseif ((G_1(i,k)==1 && G_2(j,l)==0) || (G_1(i,k)==0 && G_2(j,l)==1))
					cnt=cnt+1;
                end
            end
        end
    end
end
cnt=cnt/2; % one constraint per upper triangular entry
for i=1:n
	if (maps(i)~=-1)
		cnt=cnt+1;
	end
end
% initialize dAAT=diag(A*A')
m=n^2+1+cnt;
dAAT=zeros(m,1);
b=zeros(m,1);
% constraints X(i,i)-X(i,n^2+1)=0
% let maps be a vector of dimension n where maps(i)=j implies i->j
% also add constraints of the form X_{ij,ij}=1 when maps(i)=j
% by default maps(i)=-1
nxt=1;
for i=1:n
	for j=1:n
		if (maps(i)==j)
			dAAT(nxt,1)=1;
			b(nxt,1)=1;
			nxt=nxt+1;
			dAAT(nxt,1)=2;
			b(nxt,1)=2;
		else 
			dAAT(nxt,1)=6; 
			b(nxt,1)=0;
		end
		nxt=nxt+1;
	end
end
% constraint X(n^2+1,n^2+1)=1
dAAT(nxt,1)=1;
b(nxt,1)=1;
nxt=nxt+1;

% add the graph constraints
for i=nxt:m
    dAAT(i,1)=2;
    b(i,1)=0;
end

