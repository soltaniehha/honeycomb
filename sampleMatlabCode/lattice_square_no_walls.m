function [m]=lattice_square_no_walls(N,E0,tx,ty,alpha);
%create the matrix of a square lattice with perdical boundary condtions.
%Assuming that A=alpha*x*h/(a*e) in the direction of y
x=0;y=0;
m=zeros(N*N);
for x=1:N
    for y=1:N       %stores the interaction in the matrix of the row in X,Y with it neighbours' coloumns.
        m(get_coor(x,y,N),get_coor(x,y,N))=E0;
        m(get_coor(x,y,N),get_coor(x+1,y,N))=tx;    %no phase on the X axis
        m(get_coor(x,y,N),get_coor(x-1,y,N))=tx;
        m(get_coor(x,y,N),get_coor(x,y-1,N))=ty*exp(-i*x*alpha); %phase as a result of the magnetic field
        m(get_coor(x,y,N),get_coor(x,y+1,N))=ty*exp(i*x*alpha);  %- according to the integral on A
    end
end

        
        


