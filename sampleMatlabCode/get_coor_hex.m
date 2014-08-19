function [c]=get_coor_hex(x,y,N)
%gives the location in the matrix of the hexagonal of an X,Y position - X=row
%Y=coloum
if (x>N) x=x-N;
else if (x<1) x=x+N;
    end
end
if (y>N+1) y=y-N-1;
else if (y<1) y=y+N+1;
    end
end
c=x+(y-1)*N;
