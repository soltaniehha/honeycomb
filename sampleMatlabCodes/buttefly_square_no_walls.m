function [b]=buttefly_square_no_walls(N,E0,tx,ty,alpha,Times)
%runs through all magnetic fields from 0 to Times*Alpha, stores and
%displays the energies.
b=zeros(N*N,Times);
for i=1:Times
    m=lattice_square_no_walls(N,E0,tx,ty,alpha*(i-1));
    b(:,i)=eig(m);
end
plot(b','k.');