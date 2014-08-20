% Honey-comb (graphene) lattice, one impurity
clear all

n_states=20;    % Number of states to calculate
l=n_states;
tx=1;           % hopping parameter
ty=tx;
n=(2*l+1)^2;    %number of sites in lattice

newstates=zeros(n,n_states);
b=zeros(n_states-1,1);
a=zeros(n_states-1,1);

x=-l:1:l;
y=-l*0.9:0.9:l*0.9;                 % Generate lattice
[X,Y] = meshgrid(x,y);
P=[X(isfinite(X)),Y(isfinite(Y))];

seed = find(sum(P==0,2)==2);    %find the position of the seed
impseed=find(P(:,1)==1 & P(:,2)==0);
newstates(seed,1)=1;            %initialize seed state

[idx,D]=knnsearch(P,P,'k',5);
clear x y P D X Y

% xhop=(roundn(D,-1)==1.0);       % find x and y neighbors 
% Dont need these lines if boundary is never reached
% yhop=(roundn(D,-1)==0.9);

sites=find(newstates(:,1)); %occupied sites in previous state
need_sites=idx(sites,:);

U=0; % on site energy for impurity (ie. nitrogen in graphene)


for i=1:size(need_sites,1)
        
        newstates(need_sites(i,2:3),2)=ty*newstates(need_sites(i,1),1)+newstates(need_sites(i,2:3),2); % hops in y-direction
        newstates(need_sites(i,5),2)=tx*newstates(need_sites(i,1),1)+newstates(need_sites(i,5),2);     % hops in x-direction
end
newstates(seed,2)=U; %seed impurity is substitited 

a(1)=newstates(:,1)'*newstates(:,2);
newstates(:,2)=newstates(:,2)-a(1)*newstates(:,1);

for k=3:n_states
    
    sites=newstates(:,k-1)~=0;  % occupied sites in previous state
    need_sites=idx(sites,:);    % index neighbors
    
    for i=1:size(need_sites,1)
        newstates(need_sites(i,2:3),k)=ty*newstates(need_sites(i,1),k-1)+newstates(need_sites(i,2:3),k);  %hops in y direction
        if(mod(k,2)==0)
            newstates(need_sites(i,5),k)=tx*newstates(need_sites(i,1),k-1)+newstates(need_sites(i,5),k);  %hops in x direction
        end
        if(mod(k,2)==1)
            newstates(need_sites(i,4),k)=tx*newstates(need_sites(i,1),k-1)+newstates(need_sites(i,4),k);  %hops in -x direction
        end
        
    end
    
    a(k-1)=newstates(:,k-1)'*newstates(:,k)/(newstates(:,k-1)'*newstates(:,k-1));       % Calculate a, b, and newstate
    b(k-1)=newstates(:,k-1)'*newstates(:,k-1)/(newstates(:,k-2)'*newstates(:,k-2));
    newstates(:,k)=newstates(:,k)-a(k-1)*newstates(:,k-1)-b(k-1)*newstates(:,k-2);
end

% normalize if wanted
% for ns=1:n_states;
%     newstates(:,ns)=newstates(:,ns)/norm(newstates(:,ns));
% end
    
