function [c]=draw_brilouin_hex()
c=zeros(100);
for x=1:100
    for y=1:100
        c(x,y)=10*(cos(2*pi*2/3*y/10)+2*cos(pi*2/3*x/10*sqrt(3))*cos(pi*2/3*y/10));
    end 
end
surf(c);
