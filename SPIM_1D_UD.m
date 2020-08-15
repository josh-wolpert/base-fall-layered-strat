function SPIM_1D_UD(varargin)

% Examples:
% SPIM_1D_UD('l',20000,'dx',10,'tmax',1000,'dt',100,'Ko',3.80*10^-8,'Uo',5*10^-5,'contact',true,'c_list',[-750,5,7.4*10^-8],'uplift',true,'u_list',[10000 2.5*10^-4],'output_location','\\example_filepath')
% SPIM_1D_UD('l',20000,'dx',10,'tmax',10000000,'dt',100,'Ko',3.80*10^-8,'Uo',5*10^-5,'contact',true,'vert_c_list',[200 7.4*10^-8;10000 2*10^-8],'uplift',true,'u_list',[10000 2.5*10^-4;2000000 1*10^-4],'output_location','\\example_filepath','t_no_write',100000,'dt_write_out',1000)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overview: 
% Function to simulate bedrock stream evolution with block uplift,
% erosion governed by the stream power incision model (solved with an
% explicit upwind differencing scheme), and planar contacts.
% 'extractknickdata' is a companion function that allows users to view the
% model output and gather data on knickpoints formed throughout the
% simulation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% - 'l' (150000): Length of profile (L)
% - 'dx'(1000): Node spacing (L)
% - 'Ai' (1000000): Critical drainage area for channel initiation (L^2)
% - 'c' (6.69): Hack coefficient (Hack's Law)
% - 'h' (1.67): Hack exponent (Hack's Law)
% - 'm' (0.75): Area exponent (SPIM)
% - 'n' (1.5): Slope exponent (SPIM)
% - 'Ao' (1): Reference drainage area (L^2; chi calculations)
% - 'tmax' (10000): Model duration (t)
% - 'dt' (1000): Timestep duration (t)
% - 'Ko' (2*10^-5): Erosional efficiency of initial condition steady state profile (SPIM)
% - 'Uo' (0.001): Initial block rock uplift rate (l/t)
% - 'contacts' (false): Logical indicator for contacts. Must be 'true' to use
%                       'c_list' or 'vert_c_list' parameters
% - 'c_list' ([]): mx3 array of starting elevation in column 1, dip angle in column 2 
%                  (positive = toward outlet, negative = toward divide), and underlying
%                  K in column 3. 'contacts' must be 'true'. Cannot be used with 'vert_c_list'
% - 'vert_c_list' ([]): mx2 array of along stream distance from outlet of a
%                       vertical contact's outcrop in column 1 and K upstream of the contact in
%                       column 2. 'contacts' must be 'true'. Cannot be used with 'c_list'
% - 'uplift_change' (false): Logical indicator for change in uplift rate. Must be 'true' to use 
%                            'u_list'
% - 'u_list' ([]): mx2 array of time (t) of change in uplift rate in column 1 and uplift rate (l/t) 
%                  in column 2
% - 'output_location': File path to model output location
% - 't_no_write' (0): Time before beginning to write output (shortens run time if this is long, but 
%                     you miss time-steps)
% - 'dt_write_out' (1000): Duration of model output timestep (t). Must be a multiple of 'dt'
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
% - 'full_time': ((tmax/dt_write_out)-t_no_write)x(l/dx) matrix of model results with timesteps in 
%                rows and profile elevations in columns
% - 'c_matrix': ((tmax/dt_write_out)-t_no_write)x(l/dx)x(# of contacts) 3-D matrix of timesteps in 
%               rows, contact elevations in columns, and contact # in the dimension 3
% - 'x_scale': 1x(l/dx) array of x-axis of profile in real space
% - 'chi': 1x(l/dx) array of x-axis of profile in chi space
% - 'time_profile': individual time steps of model runs
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.FunctionName='SPIM_1D_UD';

addParameter(p,'l',150000,@(x) isnumeric(x));
addParameter(p,'dx',1000,@(x) isnumeric(x));
addParameter(p,'Ai',1000000,@(x) isnumeric(x));
addParameter(p,'c',6.69,@(x) isnumeric(x));
addParameter(p,'h',1.67,@(x) isnumeric(x));
addParameter(p,'m',0.75,@(x) isnumeric(x));
addParameter(p,'n',1.5,@(x) isnumeric(x));
addParameter(p,'Ao',1,@(x) isnumeric(x));
addParameter(p,'tmax',10000,@(x) isnumeric(x));
addParameter(p,'dt',1000,@(x) isnumeric(x));
addParameter(p,'Ko',2*10^-5,@(x) isnumeric(x));
addParameter(p,'Uo',0.001,@(x) isnumeric(x));
addParameter(p,'contacts',false,@(x) islogical(x));
addParameter(p,'c_list',[],@(x) isnumeric(x));
addParameter(p,'vert_c_list',[],@(x) isnumeric(x));
addParameter(p,'uplift_change',false,@(x) islogical(x));
addParameter(p,'u_list',[],@(x) isnumeric(x));
addParameter(p,'output_location',@(x) ischar(x));
addParameter(p,'t_no_write',0,@(x) isnumeric(x));
addParameter(p,'dt_write_out',1000,@(x) isnumeric(x));

parse(p,varargin{:});
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
U=p.Results.Uo;
contacts=p.Results.contacts;
c_list=p.Results.c_list;
vert_c_list=p.Results.vert_c_list;
uplift_change=p.Results.uplift_change;
u_list=p.Results.u_list;
directory=p.Results.output_location;
t_no_write=p.Results.t_no_write;
dt_write_out=p.Results.dt_write_out;


if numel(c_list)>=1 && numel(vert_c_list)>=1
    error('The stream cannot have vertical and non-vertical stratigraphy.')
else
end

if ~contacts && numel(c_list)>=1 || ~contacts && numel(vert_c_list)>=1
    error('''contact'' parameter must be set to true for the stream to have stratigaphy.')
else
end

% Make directory for output
mkdir(directory);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize variables and Make Initial Condition %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = l/dx; % Number of nodes in model per timestep
x_scale = [dx:dx:l]; % Downstream distances
A = zeros(nn,1);

% Populate drainage area using Hack's Law
for i=1:1:nn
    A(i) = Ai + (c*x_scale(i)^h);
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
        Zo(idx) = Zb + Zo(idx+1) + ((U/(Ko*(Ao^m)))^(1/n))*(chi(idx)-chi(idx+1));
    end
end

% Initial Condition: Steady State Profile
profile = Zo;
profile_last = Zo;

% Establish a dummy c_list if the user does not include a contact
if ~contacts
    c_list = [0,0,Ko];
else
end

% Initialize matrix to record erosion along profile during model run
E_temp=zeros(1,nn);

% Initialize vectors of up-to-date uplift rates and K 'felt' by the stream.
Kadj = zeros(1,nn)+Ko;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare contact geometries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize matrix to record contact elevations
if numel(c_list)>=1
    contact = zeros(1,nn,numel(c_list)/3);
    c_list = flipud(sortrows(c_list,1)); % sort contact list by elevations at outlet (greatest first)    
    for ii=1:numel(c_list)/3
        % Define initial contact geometries
        contact(1,:,ii) = c_list(ii,1)+(x_scale.*tand(c_list(ii,2)));
    end
    % Initialize array of previous time-step contact elevations
    contact_last = contact;

% Set permanent kadj values if there's vertical stratigraphy
elseif numel(vert_c_list)>=1
    vert_c_list(:,1) = vert_c_list(:,1)./dx; % Convert to nodes from distances
    vert_c_list = sortrows(vert_c_list,1);
    for i = 1:numel(vert_c_list(:,1))
        Kadj(vert_c_list(i,1):end) = vert_c_list(i,2);
    end
    Kadj = fliplr(Kadj);
else
end


% Sort changes in uplift by time they are imposed
if uplift_change
    u_list=sortrows(u_list,1);
else
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of model run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:(tmax/dt)
    
    % Assign appropriate K values to each node if non-vertical contacts
    if numel(c_list)>=1
        % Find elevation differences between contacts and stream nodes    
        for i=1:numel(c_list)/3
            cadj(i,:)=fliplr(contact(:,:,i))-profile; % Lithology at node is smallest positive value
        end
    
        % Make elevations of contacts stratigraphically beneath stream profile infinity 
        cadj(cadj<0)=inf;
      
        % Assign appropriate K values to stream nodes
        for i=1:numel(Kadj)
            temp_dists=sortrows([cadj(:,i) [1:numel(c_list)/3]'],1);
            if temp_dists(1,1)~=inf
                Kadj(i)=c_list(temp_dists(1,2),3);
            else
            end
        end
    else
    end
    
    % Update U from different uplift rates through time
    if uplift_change
        Us = (k*dt)-(u_list(:,1));
        Us(Us<0) = inf;
        Us = sortrows([Us [1:numel(u_list)/2]'],1); 
        
        if (Us(1,1))~=inf
            U = u_list(Us(1,2),2);
        else
        end
    else
    end
    
%%% Calculate erosion and check stability with Courant Friedrich-Lax criterion %%%
    % Model is stable if 0<=CFL=<1, as per Neumann Stability Analysis.
    % [Campforts and Govers, 2015; main and supporting text]
    for i = 1:nn-1
        a = -Kadj(i)*(A(i)^m)*(((profile_last(i)-profile_last(i+1))/dx)^(n-1));
        if (abs(a)*dt)/dx > 0.9
            error('CFL stability criterion failed')
        else
        end
        
        % Spatial loop for discretized uplift and erosion calculations
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
    
    % Update elevations of non-vertical contacts
    if numel(c_list)>=1
        % Apply uplift to contact(s)
        contact = contact_last+(U*dt);
        
        % Update previous contact elevation
        contact_last=contact;
    else
    end
    
    
% Save output for each timestep
if ~mod(k*dt,dt_write_out) && (k*dt>t_no_write)
    save(strcat(directory,['/time_',num2str(k*dt)]),'profile');
    if numel(c_list)>=1 
    save(strcat(directory,['/time_',num2str(k*dt)]),'contact','-append');
    else
    end
    save(strcat(directory,['/time_',num2str(k*dt)]),'chi','-append');
    save(strcat(directory,['/time_',num2str(k*dt)]),'x_scale','-append');
else
end
    
end

% End of Run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble results at each timestep into a single matrix

results=zeros((tmax-t_no_write)/dt_write_out,nn);

j=1;
if numel(c_list)>=1
    c_matrix=zeros((tmax-t_no_write)/dt_write_out,nn,numel(c_list)/3);
    for i=t_no_write+dt_write_out:dt_write_out:tmax
        load(strcat(directory,['/time_',num2str(i)]),'profile','contact')
        results(j,:)=profile;
        c_matrix(j,:,:)=contact;
        j=j+1;
    end
else
    for i=t_no_write+dt_write_out:dt_write_out:tmax
        load(strcat(directory,['/time_',num2str(i)]),'profile')
        results(j,:)=profile;
        j=j+1;
    end
end

% Save assembled output for all runs
save(strcat(directory,['/full_time']),'results');
if numel(c_list)>=1
    save(strcat(directory,['/full_time']),'c_matrix','-append');
else
end
save(strcat(directory,['/full_time']),'chi','-append');
save(strcat(directory,['/full_time']),'x_scale','-append');