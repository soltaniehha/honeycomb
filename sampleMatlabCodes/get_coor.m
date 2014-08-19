function [c]=get_coor(x,y,N)
%gives the location in the matrix of the square of an X,Y position - X=row
%Y=coloum
if (x>N) x=x-N;
else if (x<1) x=x+N;
    end
end
if (y>N) y=y-N;
else if (y<1) y=y+N;
    end
end
c=x+(y-1)*N;
