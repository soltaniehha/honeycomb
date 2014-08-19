function [m]=lattice_square_walls(N,E0,tx,ty,alpha);
%create the matrix of a square lattice with zero boundary condtions.
%Assuming that A=alpha*x*h/(a*e) in the direction of y
x=0;y=0;
m=zeros(N*N);
for x=1:N
    for y=1:N   %stores the interaction in the matrix of the row in X,Y with it neighbours' coloumns.
                %skips the neighbours of the edges.
        m(get_coor(x,y,N),get_coor(x,y,N))=E0;
        if (x+1<N+1)
            m(get_coor(x,y,N),get_coor(x+1,y,N))=tx;
        end
        if (x-1>0)
            m(get_coor(x,y,N),get_coor(x-1,y,N))=tx;
        end
        if (y-1>0)
            m(get_coor(x,y,N),get_coor(x,y-1,N))=ty*exp(-i*x*alpha);;   %phase as a result of the magnetic field
        end
        if (y+1<N+1)
            m(get_coor(x,y,N),get_coor(x,y+1,N))=ty*exp(i*x*alpha);;
        end
    end
end

        
        


