function [b]=butterfly_square_walls(N,E0,tx,ty,alpha,Times)
%runs through all magnetic fields from 0 to Times*Alpha, stores and
%displays the energies - with zero boundary condition.
b=zeros(N*N,Times);
for i=1:Times
    m=lattice_square_walls(N,E0,tx,ty,alpha*(i-1));
    b(:,i)=eig(m);
end
plot(b','k.');