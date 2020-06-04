function extractknickdata(filepath,output_location,varargin)

% Examples:
% extractknickdata('\\geol-af-nas\jwolpe1\MS Thesis\Additional_Runs\test\full_time','\\geol-af-nas\jwolpe1\MS Thesis\test\test1','show_model',false,'detection_threshold',1000,'model_dt',10000,'viewed_dt',10000,'ck_ksn_node_window',[15 30])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Companion function to view 1-D model runs and extract
% knickpoint data. Takes filepath for output of SPIM_1D_UD function, identifies 
% knickpoints as changes in mean slope of channel steepness arrays, and allows 
% users to pick out knickpoints of interest from time vs. elevation plots. The 
% function then calculates the vertical and horizontal (in chi-space) velocities 
% of the selected knickpoints via linear regression and cross-knickpoint changes
% in ksn at each timestep. Steady state analytical solutions predict knickpoint
% vertical velocity is constant and independent of K in chi and real space,
% while horizontal velocity is predicted to scale with K, be constant in chi space,
% and scale with drainage area in real space (Niemann et al., 2001; Royden and
% Perron, 2013; Mitchell and Yanites, 2019).
%
% Outputs:
%        - Average vertical velocities and chi celerities
%        - Cross knickpoint change in ksn over selected time window
%        - Number of timesteps used in velocity and celerity calculations
%        - Rsq for vertical velocity
%        - Time vs. Elevation plot of knickpoints
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Written by Josh Wolpert - Updated : 05/28/20 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
% Parse Inputs %
%%%%%%%%%%%%%%%%

p = inputParser;
p.FunctionName = 'CalcKnickVertV';
addRequired(p,'filepath',@(x) ischar(x)); % File path to model output
addRequired(p,'output_location',@(x) ischar(x)); % File path of location to store output (includes name of the .mat file that will house the output)

addParameter(p,'show_model',false,@(x) islogical(x)); % Option to show model output while function runs
addParameter(p,'detection_threshold',40000,@(x) isnumeric(x)); % Adjustable sensitivity threshold for detecting knickpoints. Higher is less sensitive. Lower is more sensitive.
addParameter(p,'model_dt',100,@(x) isnumeric(x)); % Timestep in model run output (not necessarily equal to model dt)
addParameter(p,'viewed_dt',500,@(x) isnumeric(x)); % Timestep used when viewing the model output. Must be a factor of tmax/model_dt
addParameter(p,'ck_ksn_node_window',[15 15],@(x) isnumeric(x)) % Node window up and downstream from knickpoint across which ksn is averaged for cross-knickpoint change in ksn calculations

parse(p,filepath,output_location,varargin{:});
filepath=p.Results.filepath;
output_location = p.Results.output_location;

show_model=p.Results.show_model;
detection_threshold=p.Results.detection_threshold;
model_dt=p.Results.model_dt;
viewed_dt=p.Results.viewed_dt;
window=p.Results.ck_ksn_node_window;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input and set up arrays and variables for main loop %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load and prepare model output data
load(filepath);

% Check that viewed_dt is a multiple of model_dt
if mod((viewed_dt),model_dt)
    error('Time step for model run must be a factor of the time step used to view model output. Choose a new ''viewed_dt''.')
else
end

% Check that viewed_dt is a factor of the model's duration
x=size(results);
if mod(model_dt*x(1),viewed_dt)
    error('Time step used to view model output must be a factor of the model''s duration. Choose a new ''viewed_dt''.')
else
end

% Define model output viewing interval
view_step=viewed_dt/model_dt;

% Define number of contacts in model
if numel(size(c_matrix)) < 3
    num_contacts = 1;
else
    num_contacts = size(c_matrix,3);
end

% Initialize output array
knicks_array=[];

% Check if window is compatible
if (numel(window) ~= 1) && (numel(window)) ~= 2
    error('Node window for cross-contact change in ksn requires either a single value or lower and upper bounds.')
else
end

% Duplicate window's value if single position is chosen
if numel(window)==1
    window=[window(1) window(1)];
else
end

% Define length of window (in nodes) for ck_ksn calculations
window_length=(max(window)-min(window))+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% View model output, find knicks, and record knick data at each viewed_dt %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=view_step:view_step:size(results,1)

    if show_model
        figure(1)
        subplot(3,1,1)
        plot(x_scale,fliplr(results(i,:)),'Color','b');
        hold on

        for k=1:num_contacts
            plot(x_scale,c_matrix(i,:,k),'r')
            hold on
        end
        
        ylim([-1000 5000])
        xlabel('Distance from Outlet [m]')
        ylabel('Elevation [m]')
        title([num2str(i*model_dt),' years'])
        hold off

        subplot(3,1,2)
        plot(chi,results(i,:),'Color','b');
        ylim([-1000 5000])
        xlim([0 max(chi)])
        xlabel('\chi')
        ylabel('Elevation [m]')
    else
    end

    dz=[results(i,2:numel(chi))-results(i,1:numel(chi)-1)];
    dchi=[chi(2:numel(chi))-chi(1:numel(chi)-1)]; % Transpose for older runs
    ksn=dz./dchi;

    if show_model
        subplot(3,1,3)
        plot(chi(1:end-1),ksn)
        ylim([0 3000])
        xlim([0 max(chi)])
        xlabel('\chi')
        ylabel('k_{sn}')
    else
    end
    
    % Find knicks and record their elevations, distances (x), and cross-knick ksns
    knicks_ksn=find(ischange(ksn,'linear','Threshold',detection_threshold));
    knicks_elev_temp=results(i,knicks_ksn)';
    knicks_dist_temp=chi(knicks_ksn)';
    
    % Initialize ksn vectors
    ksn_behind_temp=zeros(window_length,numel(knicks_ksn));
    ksn_ahead_temp=zeros(window_length,numel(knicks_ksn));

    % Assemble vectors for ksns ahead and behind knickpoints
    for j=window(1):window(2)
        for ii = 1:numel(knicks_ksn)
            try
                ksn_behind_temp(j-(window_length-1),ii)=ksn(knicks_ksn(ii)+j);
            catch
                ksn_behind_temp(j-(window_length-1),ii)=0;
            end
    
            try
                ksn_ahead_temp(j-(window_length-1),ii)=ksn(knicks_ksn(ii)-j);
            catch
                ksn_ahead_temp(j-(window_length-1),ii)=0;
            end
        end
    end
    
    % Convert 0 values to NaN
    ksn_behind_temp(ksn_behind_temp==0) = NaN;
    ksn_ahead_temp(ksn_ahead_temp==0) = NaN;

    % Calculate mean ksns over defined window for each knickpoint at a given timestep
    ksn_behind_temp = mean(ksn_behind_temp,1,'omitnan')';
    ksn_ahead_temp = mean(ksn_ahead_temp,1,'omitnan')';   
    
    knick_t_temp = zeros(numel(ksn_ahead_temp),1)+(i*model_dt);
    
    % Record: Time | Elevation | Distance | Mean Upstream Ksn | Mean Downstream Ksn
    knicks_array=[knicks_array;knick_t_temp knicks_elev_temp knicks_dist_temp ksn_ahead_temp ksn_behind_temp];
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot, Select, and analyze Knickpoints on t-z plots %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
t = knicks_array(:,1);
z = knicks_array(:,2);
graph=scatter(t,z,'k');

% Save t-z plot of knickpoints
saveas(graph,output_location)

% 1) Draw lines that form data points at every timestep.
% Select locations on the plot.
[x,y] = getpts(figure(2));

% Check that at least two points were selected
if numel(x)<2
    error('Multiple points must be selected')
else
end

% Snap selected points to time-steps
for g = 1:numel(x)
    distance = []; % Initialize and clear distance
    for f = 1:numel(t)
        distance = [distance;t(f) abs(x(g)-t(f))];
    end
    % Sort potential coordinates
    distance=sortrows(distance,2);
    x(g)=distance(1,1);
end    

% Sort selected points by x-coords (time)
xy=sortrows([x y],1);

% Initialize arrays to store x and y coordinates
py=[];
px=[];

% Make arrays of data for each linear segment between two selected points.
for jj=1:(numel(xy)/2)-1
    s = polyfit([xy(jj,1) xy(jj+1,1)],[xy(jj,2) xy(jj+1,2)],1); % First component of s is slope, second component is y-intercept.
    
    % Find y-coordinates in t-z space at each time step along line segment
    for kk=xy(jj,1):viewed_dt:xy(jj+1,1)
        py=[py;(s(1)*kk)+s(2)];
        px=[px;kk];
    end
end

% 2) Find the nearest knickpoint to every point on the line
k = dsearchn([t z],[px py]);

% The output of 'dsearchn' will be the indices of interest in 'knicks_array', so you can then extract info for those knicks.
koi = knicks_array(k,:); % Extract data for knicks of interest

% Remove Repeats
koi = unique(koi,'rows');

% 3) Calculate vertical velocity and Rsq for velocity of knicks of interest
vv = polyfit(koi(:,1),koi(:,2),1);
Vertical_velocity = vv(1);
B = corrcoef(koi(:,1),koi(:,2));
Rsqv = B(1,2)^2;
N = numel(koi(:,1));

% 4) Calculate horizontal velocity in chi-space for knicks
cc = polyfit(koi(:,1),koi(:,3),1);
chi_celerity = cc(1);

% 5) Store: Time | Mean Upstream Ksn | Mean Downstream Ksn | Change in Ksn
ck_ksn = [koi(:,1) koi(:,4) koi(:,5) koi(:,4)-koi(:,5)];


%%%%%%%%%%%%%%%
% Save output %
%%%%%%%%%%%%%%%

save(output_location,'Vertical_velocity','-append');
save(output_location,'chi_celerity','-append');
save(output_location,'Rsqv','-append');
save(output_location,'N','-append');
save(output_location,'ck_ksn','-append');
end