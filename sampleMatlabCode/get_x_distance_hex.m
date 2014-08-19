function [c]=get_x_distance_hex(x,y)
%the distance in the X axis of a point in the hexagonal lattice from the
%left border of the lattice.
if floor((x+y)/2)*2==x+y 
    c=1.5*(x-1)+0.5;
else
    c=1.5*(x-1);
end
