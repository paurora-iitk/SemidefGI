Y1=Y;
r=1;
c=1;
while(c<=n^2)
	flag=0;
	pv=Y1(r,c);
	if(pv<1e-8)
		for i=r+1:n^2
			if(abs(Y1(i,c))>1e-8)
				if(Y1(i,c)<0)
					Y1(i,c)
				end;
				break;
			end;
		end;
		if(i<=n^2 && abs(Y1(i,c))>1e-8)
			tmp=Y1(i,1:n^2);
			Y1(i,1:n^2)=Y1(r,1:n^2);
			Y1(r,1:n^2)=tmp(1,1:n^2);
			pv=Y1(r,c);
		else
			c=c+1;
			flag=1;
		end;
	end;
	if(flag==0)
		for j=r+1:n^2
			if(abs(Y1(j,c))>1e-8)
				x=Y1(j,c)/pv;
				Y1(j,:)=Y1(j,:)-x.*Y1(r,:);
			end;
		end;
		r=r+1;
		c=c+1;
	end;
end;
