function [b]=butterfly_hex_walls(N,E0,t,alpha,Times)
%runs through all magnetic fields from 0 to Times*Alpha, stores and
%displays the energies - with zero boundary condition.
b=zeros(N*(N+1),Times);
for i=1:Times
    m=lattice_hex_walls(N,E0,t,alpha*(i-1));
    b(:,i)=real(eig(m));
end
display_contour(b');
%plot(b','k.');