function [m]=lattice_hex_walls(N,E0,t,alpha);
%create the matrix of a hexagonal lattice with zero boundary condtions.
%Assuming that A=Ax in the direction of y
%the exponent is calculated where alpha=A*a^2*sin(60)*e/h
x=0;y=0;
m=zeros(N*(N+1));
for x=1:N
    for y=1:N+1     %stores the interaction in the matrix of the row in X,Y with it neighbours' coloumns.
                    %skips the neighbours of the edges.
        m(get_coor_hex(x,y,N),get_coor_hex(x,y,N))=E0;
        if floor((x+y)/2)*2==x+y    %the two y+1,y-1 neigbours are in the back of the point X,Y
            if (y<N+1) 
                m(get_coor_hex(x,y,N),get_coor_hex(x,y+1,N))=t*exp(i*(get_x_distance_hex(x,y)*alpha-alpha*1/4));
            end                                         
            if (y>1)
                m(get_coor_hex(x,y,N),get_coor_hex(x,y-1,N))=t*exp(i*(-get_x_distance_hex(x,y)*alpha+alpha*1/4));
            end
            if (x<N)
                m(get_coor_hex(x,y,N),get_coor_hex(x+1,y,N))=t;
            end
        else                         %the two y+1,y-1 neigbours are in front of the point X,Y
            if (y<N+1)
                m(get_coor_hex(x,y,N),get_coor_hex(x,y+1,N))=t*exp(i*(get_x_distance_hex(x,y)*alpha+alpha*1/4));
            end
            if (y>1)
                m(get_coor_hex(x,y,N),get_coor_hex(x,y-1,N))=t*exp(i*(-get_x_distance_hex(x,y)*alpha-alpha*1/4));
            end
            if (x>1)
                m(get_coor_hex(x,y,N),get_coor_hex(x-1,y,N))=t;
            end
        end            
    end
end
       
        


