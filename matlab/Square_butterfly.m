% Square Butterfly mapping with magnetic field and potential
%
%   phi = integral(A.dl)
%   A = Bx e_y

clear all

% for B=0.4:0.02:0.6

l=20;       % linear dimension of lattice=2a            
n_states=l; % number of states

% % This is one way to find seed and generate lattice ******
x=-l:1:l;
y=-l:1:l;
[X,Y] = meshgrid(x,y);
Lattice=(abs(X)+abs(Y)<=l);     %Generate Lattice
P=[X(Lattice),Y(Lattice)];
n=size(P,1);
seed = find(sum(P==0,2)==2);
idx=knnsearch(P,P,'k',5); % takes about 5 seconds for l=600
clear x y X Y Lattice
% %*************************************************

V=0;        % Potential
B=0.1;        % Magnetic Field
a0=1;       % length of unit cell

newstates=zeros(n,n_states);    % Initialize states and hopping elements
newstates2=zeros(n,3); %***test*****
b=zeros(n_states-1,1);
a=zeros(n_states-1,1);

tx=1;
ty=1;

% Calculate newstate 1
newstates(seed,1)=1;

% Calculate newstate 2
newstates(idx(seed,2),2)=tx*newstates(seed,1);  % hop in -x
newstates(idx(seed,5),2)=tx*newstates(seed,1);  % hop in +x

x_pos = 0;
t_up = ty*exp(1j*B*2*pi*x_pos*a0);
t_down = ty*exp(-1j*B*2*pi*x_pos*a0);

newstates(idx(seed,3),2) = t_down*newstates(seed,1);    % hop in -y
newstates(idx(seed,4),2) = t_up*newstates(seed,1);      % hop in +y

% Here I normalize as I iterate
newstates(:,2)=newstates(:,2)/norm(newstates(:,2));

% Calculate newstate 3..and so on..
% ************ This method keeps all states ***************
for state = 3:n_states

sites = newstates(:,state-1)~=0;    % finds occupied sites
needs=idx(sites,:);                 % index of nearest neighbors   

    for i=1:size(needs,1)
        active=needs(i,1);
%         hop=[active+1, active-1, active+move, active-move]; % hop to nearest sites
%         newstates(hop,state)=newstates(active,state-1) + newstates(hop,state);
        x_pos = P(active,1);
        y_pos = P(active,2);
        dist = sqrt(x_pos^2+y_pos^2);
        if (active==seed)
            on_ener=0;
        else
            on_ener = -V/dist;  % On site energy
        end
        t_up = ty*exp(1j*B*2*pi*x_pos*a0);
        t_down = conj(t_up);
        newstates(needs(i,2),state) = tx*newstates(active,state-1) + newstates(needs(i,2),state);       % -x
        newstates(needs(i,5),state) = tx*newstates(active,state-1) + newstates(needs(i,5),state);       % +x
        newstates(needs(i,3),state) = t_down*newstates(active,state-1) + newstates(needs(i,3),state);   % -y
        newstates(needs(i,4),state) = t_up*newstates(active,state-1) + newstates(needs(i,4),state);     % +y
        
        newstates(active,state) = on_ener*newstates(active,state-1) + newstates(active,state);  % on site potential;
        
    end
    
%     a(state-1)=newstates(:,state-1)'*newstates(:,state)/(newstates(:,state-1)'*newstates(:,state-1));       % Calculate a, b, and newstate
%     b(state-1)=newstates(:,state-1)'*newstates(:,state-1)/(newstates(:,state-2)'*newstates(:,state-2));

    a(state-1)=newstates(:,state-1)'*newstates(:,state);
    b(state-1)=newstates(:,state-2)'*newstates(:,state);
    newstates(:,state) = newstates(:,state) - a(state-1)*newstates(:,state-1) - b(state-1)*newstates(:,state-2);
    
    % % % Calculate Overlap of previous orbitals
    for n_ov=1:state-1
        over_temp = newstates(:,n_ov)'*newstates(:,state);
        if (abs(over_temp) > 10^-15 )
%             fprintf('%17.15f %d %d\n',over_temp,n_ov,k);
            newstates(:,state) = newstates(:,state) - over_temp*newstates(:,n_ov);
        end
    end
    
    newstates(:,state) = newstates(:,state)/norm(newstates(:,state));


end
% % *******************************************
b=real(b);

% %********check orthogonality**************

for ns=1:n_states-1
    for nt=ns+1:n_states
        overlap = newstates(:,ns)'*newstates(:,nt);
        if (abs(overlap)>10^(-8))
            fprintf('%14.10f %d %d\n', overlap, ns, nt);
        end
    end
end
    



% A=[a'; b'];
% fname = ['Square_Butterfly_Bpi_',num2str(B),'.txt'];
% fid=fopen(fname,'w');
% fprintf(fid,'%15.12f %15.12f\n',A);
% fclose(fid);
% 
% clearvars -except B
% 
% end




% % % Calculate newstate 3..and so on..
% for state = 3:n_states
%     
%     if (mod(state,3)==0)
%         now=3;
%         prev=2;
%         prev_2=1;
%     elseif (mod(state,3)==1)
%         now=1;
%         prev=3;
%         prev_2=2;
%     elseif (mod(state,3)==2)
%         now=2;
%         prev=1;
%         prev_2=3;
%     else
%         error('Something wrong');
%         
%     end
%     
%     newstates(:,now)=0;
% 
% sites = newstates(:,prev)~=0;    % finds occupied sites
% needs=idx(sites,:);                 % index of nearest neighbors   
% 
%     for i=1:size(needs,1)
%         newstates(needs(i,2:5),now) = newstates(needs(i,1),prev) + newstates(needs(i,2:5),now); % hops to nearest neighbors
%     end
%     
%     b(state-1)=newstates(:,prev_2)'*newstates(:,now);    %b=<alpha-2|H|alpha-1>
%     newstates(:,now)=newstates(:,now)-b(state-1)*newstates(:,prev_2); % Calculate next state
%     
%     newstates(:,now)=newstates(:,now)/norm(newstates(:,now)); % Normalize state
% 
% end

