function SPIM_1D_UD(output_location,varargin)
%
% Usage:
%       SPM_1D_UD('/Users/joshw/Documents/MATLAB/1_D Models/Model_Runs/test','m',0.75,'n',1.5,'l',20000,'dx',10,'tmax',1000000,'dt',100,'Ko',3.8*10^-8,'Uo',2.5*10^-4,'contact',true,'c_list',[-1910 10 3.8*10^-8;-2210 10 7.22*10^-7],'uplift',true,'u_list',[15000000 2.5*10^-4],'t_no_write',500000,'dt_write_out',10000)
%
% Description:
%       Function to solve the stream power incision model with the explicit upwind differencing scheme (See supporting information to 
%       Campforts and Govers, 2015). Users establish an initial condition with Hack's Law parameters (Hack, 1957) and define stream
%       power incision model parameters, geometries of planar contacts, and timing and magnitudes of changes in uplift rate.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Written by Josh Wolpert - Updated : 06/3/20 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse Inputs
p = inputParser;
p.FunctionName='SPIM_1D_UD';

addRequired(p,'output_location',@(x) ischar(x)); % File path for model output

addParameter(p,'l',150000,@(x) isnumeric(x)); % Length of model profile [m]
addParameter(p,'dx',1000,@(x) isnumeric(x)); % Node spacing [m]
addParameter(p,'Ai',1000000,@(x) isnumeric(x)); % Critical Drainage Area [m^2]
addParameter(p,'c',6.69,@(x) isnumeric(x)); % Hack coefficient
addParameter(p,'h',1.67,@(x) isnumeric(x)); % Hack exponent
addParameter(p,'m',0.75,@(x) isnumeric(x)); % Area exponent
addParameter(p,'n',1.5,@(x) isnumeric(x)); % Slope exponent
addParameter(p,'Ao',1,@(x) isnumeric(x)); % Reference drainage area [m^2]
addParameter(p,'tmax',10000,@(x) isnumeric(x)); % Model Duration
addParameter(p,'dt',1000,@(x) isnumeric(x)); % Timestep length [yr]
addParameter(p,'Ko',2*10^-5,@(x) isnumeric(x)); % Initial erosional efficiency
addParameter(p,'Uo',0.001,@(x) isnumeric(x)); % Initial uplift rate [m/yr]
addParameter(p,'contacts',false,@(x) islogical(x)); % Logical indicator for contacts
addParameter(p,'c_list',@(x) isnumeric(x)); % mx3 matrix of elevations in column 1, dip angle in column 2, and K value in column 3
addParameter(p,'uplift',false,@(x) islogical(x)); % Logical indicator for imposing uplift rates
addParameter(p,'u_list',[],@(x) isnumeric(x)); % mx2 matrix of time [yr] and uplift rate [m/yr]
addParameter(p,'t_no_write',0,@(x) isnumeric(x)); % Time before beginning writing out results
addParameter(p,'dt_write_out',1000,@(x) isnumeric(x)); % Time [yr] between saving results

parse(p,output_location,varargin{:});
directory=p.Results.output_location;

l=p.Results.l;
dx=p.Results.dx;
Ai=p.Results.Ai;
c=p.Results.c;
h=p.Results.h;
m=p.Results.m;
n=p.Results.n;
Ao=p.Results.Ao;
dt=p.Results.dt;
tmax=p.Results.tmax;
Ko=p.Results.Ko;
Uo=p.Results.Uo;
contacts=p.Results.contacts;
c_list=p.Results.c_list;
uplift=p.Results.uplift;
u_list=p.Results.u_list;
t_no_write=p.Results.t_no_write;
dt_write_out=p.Results.dt_write_out;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make directory for output to write to
mkdir(directory);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate initial condition steady state profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = l/dx; % Number of nodes in model per timestep
x_scale = [dx:dx:l]; % [m] Downstream distances
A = zeros(nn,1);

% Populate drainage area using Hack's Law
for i=1:1:nn
    A(i) = Ai + (c*x_scale(i)^h); % [m^2]
end

% Use chi method to predict elevations at each x_scale interval
% [Perron and Royden, 2013; Mitchell and Yanites, 2019]   
integrand = (Ao./A).^(m/n);
integrand=flipud(integrand); % Flip integrand for upstream integration
chi=cumtrapz(x_scale,integrand); % First index of chi is at the outlet
chi=[flipud(chi)]'; % Convert back to channel head at first index

% Calculate initial condition profile
Zo=chi*0;

Zb=0;

for i=0:numel(chi)-1 
    idx = numel(chi)-i;
    if idx == numel(chi)
        Zo(idx) = Zb;
    else
        Zo(idx) = Zb + Zo(idx+1) + ((Uo/(Ko*(Ao^m)))^(1/n))*(chi(idx)-chi(idx+1));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set contact initial geometries and prepare arrays to store various timestep and output info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

profile_last = zeros(1,nn);

profile = zeros(1,nn);

% Establish a dummy c_list if the user does not include a contact
if ~contacts
    c_list = [0,0,Ko];
else
end

% Initialize matrix to record erosion along profile during model run
E_temp=zeros(1,nn);

% Initial Condition: Steady State Profile
profile = Zo;
profile_last = Zo;

% Initialize vectors of up-to-date uplift rates and K 'felt' by the stream.
Kadj = zeros(1,nn)+Ko;

% Initialize for later use
U=Uo;

% Initialize matrix to record contact elevations
contact = zeros(1,nn,numel(c_list)/3);

c_list = flipud(sortrows(c_list,1)); % sort contact list by elevations at outlet (greatest first)

for ii=1:numel(c_list)/3
    
    % Define initial contact geometries
    contact(1,:,ii) = c_list(ii,1)+(x_scale.*tand(c_list(ii,2)));
    
end

% Initialize array of previous time-step contact elevations
contact_last = contact;

% Sort changes in uplift rate by time of initiation
if uplift
    u_list=sortrows(u_list,1);
else
end

%%%%%%% Populate below for Vertical Stratigraphy %%%%%%%
%Kadj(1:1749)=2.8*10^-7;
%Kadj(1750:2000)=1*10^-7;
%%%%%%% Populate above for Vertical Stratigraphy %%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 2:(tmax/dt)  % This should probably start at 1.
    
    %%%%%%% Disable below if vertical contacts %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find elevation differences between contacts and stream nodes    
    for i=1:numel(c_list)/3
        cadj(i,:)=fliplr(contact(:,:,i))-profile; % Lithology at node is smallest positive value
    end
    
    % Make elevations of contacts stratigraphically beneath stream profile infinity 
    cadj(cadj<0)=inf;
    
    %%%%%%%%%%%% Disable below for Vertical Stratigraphy %%%%%%%%%%%%%%%%%%
    % Assign appropriate K values to stream nodes
    for i=1:numel(Kadj)
        temp_dists=sortrows([cadj(:,i) [1:numel(c_list)/3]'],1);
        if temp_dists(1,1)~=inf
            Kadj(i)=c_list(temp_dists(1,2),3);
        else
        end
    end
    %%%%%%%%%%%% Disable above for vertical contacts %%%%%%%%%%%%%%%%%%%%%%
    
    % Update Uadj from different uplift rates through time
    if uplift
        Us = (k*dt)-(u_list(:,1));
        Us(Us<0) = inf;
        Us = sortrows([Us [1:numel(u_list)/2]'],1); 
        
        if (Us(1,1))~=inf
            U = u_list(Us(1,2),2);
        else
        end
    else
    end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate erosion and check stability with Courant Friedrich-Lax criterion
% Model is stable if 0<=CFL=<1, as per Neumann Stability Analysis.
% [Campforts and Govers, 2015; main and supporting text]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:nn-1
        a = -Kadj(i)*(A(i)^m)*(((profile_last(i)-profile_last(i+1))/dx)^(n-1));
        if (abs(a)*dt)/dx > 0.9
            error('CFL stability criterion failed')
        else
        end
        
        % Set node index
        idx = nn-i;
        
        % Calculate Erosion (explicit forward upwind differencing scheme [Campforts and Govers, 2015])
        E_temp(idx) = (-1*Kadj(idx)*A(idx)^(m)*((-1*((profile_last(idx+1)-profile_last(idx))/dx))^n)*dt);  
    end
    
    % Calculate and apply uplift to timestep k
    profile(1:nn-1) = profile_last(1:nn-1)+(U*dt);
    
    % Apply Erosion to timestep k
    profile(1:nn-1) = profile(1:nn-1)+E_temp(1:nn-1);
    
    % Update previous timestep
    profile_last = profile;
    
    % Update elevations of contacts
    if contacts
        % Apply uplift to contact(s)
        contact = contact_last+(U*dt);
    else
    end
    
    % Update previous contact position
    contact_last=contact;
    
% Save output for each timestep
if ~mod(k*dt,dt_write_out) && (k*dt>t_no_write)
    save(strcat(directory,['\time_',num2str(k*dt)]),'profile');
    save(strcat(directory,['\time_',num2str(k*dt)]),'contact','-append');
    save(strcat(directory,['\time_',num2str(k*dt)]),'chi','-append');
    save(strcat(directory,['\time_',num2str(k*dt)]),'x_scale','-append');
else
end
    
end

% End of main loop


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble results at each timestep into a single matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results=zeros((tmax-t_no_write)/dt_write_out,nn);
c_matrix=zeros((tmax-t_no_write)/dt_write_out,nn,numel(c_list)/3);

j=1;
for i=t_no_write+dt_write_out:dt_write_out:tmax
    load(strcat(directory,['\time_',num2str(i)]),'profile','contact')
    results(j,:)=profile;
    c_matrix(j,:,:)=contact;
    j=j+1;
end

% Save assembled output for all runs
save(strcat(directory,['\full_time']),'results');
save(strcat(directory,['\full_time']),'c_matrix','-append');
save(strcat(directory,['\full_time']),'chi','-append');
save(strcat(directory,['\full_time']),'x_scale','-append');    
