function [flag,itr] = graph_iso_test(G_1,G_2,map,itr)

n=size(G_1,1);
C=eye(n^2+1);
C(n^2+1,n^2+1)=0;
C=-1.*C;
[dAAT,b]=compute_daat(G_1,G_2,map); 
[X,y,Z]=DADAL(G_1,G_2,dAAT,b,C,map);
itr=itr+1;
fprintf('itrn: %d\n',itr);
fprintf('map=> ');
for i=1:n
	fprintf('%d ',map(i));
end
fprintf('\n');
tol=1e-2;
if (abs(trace(X)-n-1)>0.6)
	flag=0;
	fprintf('trace(X): %f\n',trace(X));
	return;
elseif (rank(X,tol)<=3)
	flag=1;
	return;
else
	max_val=0;
	for i=1:n^2
		if (abs(1-X(i,i))>tol && X(i,i)-max_val>tol)
			max_val=X(i,i);
		end
	end
	fprintf('max_val: %f\n',max_val);
	for i=1:n
		t=0;
		for j=1:n 
			if (abs(X(n*(i-1)+j,n*(i-1)+j)-max_val)<tol) 
				map(i)=j; 
				[flag,itr]=graph_iso_test(G_1,G_2,map,itr);
				if (flag==1) 
					return;
				else
					map(i)=-1;
				end
			elseif (X(n*(i-1)+j,n*(i-1)+j)>tol)
				t=1;
			end
		end
		if (t==0)
			flag=0;
			return;
		end
	end
	flag=0;
end
