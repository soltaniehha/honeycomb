% Butterfly Mapping
clear all


n_states=201;    % Number of states to calculate
l=n_states;

y=-l:1:l;
x=y;        % Generate lattice
[X,Y] = meshgrid(x,y);
P2=[X(isfinite(X)),Y(isfinite(Y))];
n=size(P2,1);
[idx,D] = knnsearch(P2,P2,'k',5);
clear x y X Y

newstates=zeros(n,n_states);
b=zeros(n_states-1,1);
a=zeros(n_states-1,1);

tx=1;
ty=1;

V=0;
B=0.4;
a0=1;

seed = find(sum(P2==0,2)==2);    % find the position of the seed
newstates(seed,1)=1;            % initialize seed state

sites=find(newstates(:,1)); %occupied sites in previous state
need_sites=idx(sites,:);

% ********** NOTES ************
% % Seed is odd 
% % need_sites(i,1) = active site
% % need_sites(i,2) = left site
% % need_sites(i,3) = down site
% % need_sites(i,4) = up site
% % need_sites(i,5) = right site
% % ****************************

for i=1:size(need_sites,1) % loop is unnecessary size=1
    active = need_sites(i,1);
    xpos = P2(active,1); % here xpos is zero becuase seed is odd
    t_up = ty*exp(1j*B*(sqrt(3)*a0*xpos/2 - sqrt(3)*a0*a0/8));
    t_down = ty*exp(-1j*B*(sqrt(3)*a0*xpos/2 - sqrt(3)*a0*a0/8));
        newstates(need_sites(i,3),2) = t_down*newstates(need_sites(i,1),1) + newstates(need_sites(i,3),2); % hops in -y direction
        newstates(need_sites(i,4),2) = t_up*newstates(need_sites(i,1),1) + newstates(need_sites(i,4),2); % hops in +y direction
        newstates(need_sites(i,5),2) = tx*newstates(need_sites(i,1),1) + newstates(need_sites(i,5),2);     % hops in +x direction
end

a(1)=newstates(:,1)'*newstates(:,2);
newstates(:,2)=newstates(:,2)-a(1)*newstates(:,1);

newstates(:,2)=newstates(:,2)/norm(newstates(:,2));

for k=3:n_states
    
    sites=newstates(:,k-1)~=0;  % occupied sites in previous state
    need_sites=idx(sites,:);    % index neighbors
    
    for i=1:size(need_sites,1)
        
        active = need_sites(i,1);
        
        if (active==seed+1)
            active;
        end
        ypos=sqrt(3)*0.5*a0*P2(active,2);
        
        % Check and get x postion in real space ************
        xpos = P2(active,1);
        if (mod(active,2)==0)
            even=true;
        elseif (mod(active,2)==1)
            even=false;
        end
        
        if (even)
            xpos=(xpos*3/2 - 0.5)*a0;
        elseif (~even)
            xpos=xpos*a0*3/2; %*****
        else
            error('Problem with x position');
        end % end of x position finding************
        
        dis=sqrt(xpos^2 + ypos^2);
        
        if (active==seed)
            on_energy = 0;
        else
            on_energy = -V/dis;
        end
        
        if (even) % even site
            
            t_up = ty*exp(1j*B*(sqrt(3)*a0*xpos/2 + sqrt(3)*a0*a0/8));
            t_down = ty*exp(-1j*B*(sqrt(3)*a0*xpos/2 + sqrt(3)*a0*a0/8));
            
            newstates(need_sites(i,2),k) = tx*newstates(active,k-1) + newstates(need_sites(i,2),k);     % hops in -x direction
            newstates(need_sites(i,3),k) = t_down*newstates(active,k-1) + newstates(need_sites(i,3),k); % hops in -y direction
            newstates(need_sites(i,4),k) = t_up*newstates(active,k-1) + newstates(need_sites(i,4),k);   % hops in +y direction
            
        elseif (~even) % odd site
            
            t_up = ty*exp(1j*B*(sqrt(3)*a0*xpos/2 - sqrt(3)*a0*a0/8));
            t_down = ty*exp(-1j*B*(sqrt(3)*a0*xpos/2 - sqrt(3)*a0*a0/8));
            
            newstates(need_sites(i,5),k) = tx*newstates(active,k-1) + newstates(need_sites(i,5),k);     % hops in x direction
            newstates(need_sites(i,3),k) = t_down*newstates(active,k-1) + newstates(need_sites(i,3),k); % hops in -y direction
            newstates(need_sites(i,4),k) = t_up*newstates(active,k-1) + newstates(need_sites(i,4),k);   % hops in +y direction
            
        else
            error('Site not even or odd');
        end
        
        newstates(active,k) = on_energy*newstates(active,k-1) + newstates(active,k);

    end
    
%     a(k-1)=newstates(:,k-1)'*newstates(:,k)/(newstates(:,k-1)'*newstates(:,k-1));       % Calculate a, b, and newstate
%     b(k-1)=newstates(:,k-1)'*newstates(:,k-1)/(newstates(:,k-2)'*newstates(:,k-2));
    a(k-1)=newstates(:,k-1)'*newstates(:,k);
    b(k-1)=newstates(:,k-2)'*newstates(:,k);
    
    newstates(:,k) = newstates(:,k) - a(k-1)*newstates(:,k-1) - b(k-1)*newstates(:,k-2);
    
%     newstates(:,k)=newstates(:,k)/norm(newstates(:,k));
    
    % % ********* Calculate Overlap of previous orbitals
    for n_ov=1:k-1
        over_temp = newstates(:,n_ov)'*newstates(:,k);
        if (abs(over_temp) > 10^-15 )
%             fprintf('%17.15f %d %d\n',over_temp,n_ov,k);
            newstates(:,k) = newstates(:,k) - over_temp*newstates(:,n_ov);
        end
    end
    % % *********************************************
    
    newstates(:,k)=newstates(:,k)/norm(newstates(:,k));
    
end

% for nk=1:n_states
%     newstates(:,nk)=newstates(:,nk)/norm(newstates(:,nk));
% end
