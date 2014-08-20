% Single Site Impurity, Mapping into Lanzcos Basis
clear all

l=601;       % linear dimension of lattice=2a            
n_states=l; % number of states

% % This is one way to find seed and generate lattice ******
x=-l:1:l;
y=-l:1:l;
[X,Y] = meshgrid(x,y);
Lattice=(abs(X)+abs(Y)<=l);     %Generate Lattice
P=[Y(Lattice),X(Lattice)];
n=size(P,1);
seed = find(sum(P==0,2)==2);
[idx,D]=knnsearch(P,P,'k',5); % takes about 5 seconds for l=600
clear x y X Y D P Lattice
% %*************************************************


% % This is other method
% n=(2*l+1)^2;      % number of sites
% seed=ceil(n/2); % find seed position

newstates=zeros(n,n_states);    % Initialize states and hopping elements
newstates2=zeros(n,3); %***test*****
b=zeros(n_states-1,1);
a=zeros(n_states-1,1);

% % For second method
move=2*l+1;     % this hops left and right

% Calculate newstate 1
newstates(seed,1)=1;
newstates2(seed,1)=1; %****test*****

% % Calculte newstate 2 , method 2
% newstates(seed+1,2)=newstates(seed,1);      % hop up
% newstates(seed-1,2)=newstates(seed,1);      % hop down
% newstates(seed+move,2)=newstates(seed,1);   % hop right
% newstates(seed-move,2)=newstates(seed,1);   % hop left

% Calculate newstate 2
newstates(idx(seed,2:5),2)=newstates(seed,1);
newstates2(idx(seed,2:5),2)=newstates2(seed,1);

% Here I normalize as I iterate
newstates(:,2)=newstates(:,2)/norm(newstates(:,2));


% Calculate newstate 3..and so on..
% % ************ This method keeps all states ***************
% for state = 3:n_states
% 
% % sites = find(newstates(:,state-1)); % alternative method
% sites = newstates(:,state-1)~=0;    % finds occupied sites
% needs=idx(sites,:);                 % index of nearest neighbors   
% 
%     for i=1:size(needs,1)
% %         active=sites(i);
% %         hop=[active+1, active-1, active+move, active-move]; % hop to nearest sites
% %         newstates(hop,state)=newstates(active,state-1) + newstates(hop,state);
%         newstates(needs(i,2:5),state) = newstates(needs(i,1),state-1) + newstates(needs(i,2:5),state); % hops to nearest neighbors
%     end
%     
%     b(state-1)=newstates(:,state-2)'*newstates(:,state);    %b=<alpha-2|H|alpha-1>
%     newstates(:,state)=newstates(:,state)-b(state-1)*newstates(:,state-2); % Calculate next state
%     
%     % Note that all a's are zero for single impurity bipartite lattices
%     
%     newstates(:,state)=newstates(:,state)/norm(newstates(:,state)); % Normalize state
% 
% end
% % *******************************************



% % Calculate newstate 3..and so on..
for state = 3:n_states
    
    if (mod(state,3)==0)
        now=3;
        prev=2;
        prev_2=1;
    elseif (mod(state,3)==1)
        now=1;
        prev=3;
        prev_2=2;
    elseif (mod(state,3)==2)
        now=2;
        prev=1;
        prev_2=3;
    else
        error('Something wrong');
        
    end
    
    newstates(:,now)=0;

sites = newstates(:,prev)~=0;    % finds occupied sites
needs=idx(sites,:);                 % index of nearest neighbors   

    for i=1:size(needs,1)
        newstates(needs(i,2:5),now) = newstates(needs(i,1),prev) + newstates(needs(i,2:5),now); % hops to nearest neighbors
    end
    
    b(state-1)=newstates(:,prev_2)'*newstates(:,now);    %b=<alpha-2|H|alpha-1>
    newstates(:,now)=newstates(:,now)-b(state-1)*newstates(:,prev_2); % Calculate next state
    
    newstates(:,now)=newstates(:,now)/norm(newstates(:,now)); % Normalize state

end


