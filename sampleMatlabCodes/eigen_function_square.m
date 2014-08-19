function [e]=eigen(N,m,vec);
[V,D] = eig(m); x=0;y=0;
e=zeros(N);
for x=1:N
    for y=1:N
        e(x,y)=V(get_coor(x,y,N),vec);
    end
end
surf(e);