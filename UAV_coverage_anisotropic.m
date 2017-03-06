% Copyright 2016-2017 Sotiris Papatheodorou
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%    http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

clear variables
close all

% TODO
% Partitioning - ok
% Area-objective calculation - ok
% Communication range calculation
% zopt function and Hopt calculation
% Control law, use three parts (planar, altitude, rotational)
%   planar - ok
% Make control law parametric with regards to the jacobian
%   Define sensing pattern through symbolic expression

%%%%%%%%%%%%%%%%%%% Set Simulation Options %%%%%%%%%%%%%%%%%%%
% Network options
% Altitude constraints
zmin = 0.3;
zmax = 2.3;

% Simulation options
% Simulation duration in seconds
Tfinal = 20;
% Time step in seconds
Tstep = 0.1;

% Control law options
% Planar control law gain
axy = 1;
% Altitude control law gain
az = 1;
% Rotational control law gain
ath = 1;

% Network plots to show during simulation
PLOT_STATE_2D = 1;
PLOT_STATE_3D = 0;
PLOT_STATE_QUALITY = 0;
SAVE_PLOTS = 0;

% Save simulation results to file
SAVE_RESULTS = 0;

% Use finite communication range
COMM_RANGE = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%% Sensing Pattern Parametric equation %%%%%%%%%%%%%%%%%%%
% Parametric equation of base sensing patern:
% node at [0 0 zmin] and orientation 0
a = 0.2;
b = 0.1;
syms t gx gy g
assume([t gx gy g],'real');
gx = a * cos(t);
gy = b * sin(t);
g = [gx ; gy];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% Add function path
addpath( genpath('Functions') );

% ---------------- Region ----------------
% Bullo region
Xb=[ 0, 2.125, 2.9325, 2.975, 2.9325, 2.295, 0.85, 0.17 ];
Yb=[ 0, 0, 1.5, 1.6, 1.7, 2.1, 2.3, 1.2 ];
[Xb, Yb] = poly2cw(Xb, Yb);
region = [Xb ; Yb];
region_area = polyarea( Xb, Yb );
axis_scale = [-0.5 3 -0.5 3];

% ---------------- Initial State ----------------
X = [0.40, 0.60, 0.55];
Y = [0.50, 0.60, 0.50];
Z = [0.45, 0.55, 0.50];
TH = [0.1 -0.2 0.5];


% ---------------- Simulation initializations ----------------
% Caclulate optimal altitude
zopt = z_optimal_uniform(zmin, zmax);
% Number of nodes
N = length(X);
% Simulation steps
smax = floor(Tfinal/Tstep);
% Points Per Circle
PPC = 60;
% Radius for points on plots
disk_rad = 0.02;
% Vector for circle parameterization
t = linspace(0, 2*pi, PPC+1);
t = t(1:end-1); % remove duplicate last element
t = fliplr(t); % flip to create CW ordered circles
% Sensing pattern with node at origin, zmin and theta_i = 0
Cb = [a*cos(t) ; b*sin(t)];
Cd_diameter = sensing_diameter(Cb);
% Simulation data storage
Xs = zeros(smax, N);
Ys = zeros(smax, N);
Zs = zeros(smax, N);
THs = zeros(smax, N);
cov_area = zeros(smax,1);
H = zeros(smax,1);
% Initialize (cell) arrays
% Coverage quality
f = zeros(1, N);
% Sensing disks
C = cell([1 N]);
% Sensed space partitioning
W = cell([1 N]);
% Control inputs
uX = zeros(1,N);
uY = zeros(1,N);
uZ = zeros(1,N);
uTH = zeros(1,N);
% Communication range for each node
r_comm = zeros(1,N);
% Adjacency matrix for communication graph
A = zeros(N,N);

% Create sim struct that contains all information, used for plots
sim = struct;
sim.region = region;
sim.axis = axis_scale;
sim.PPC = PPC;
sim.zmin = zmin;
sim.zmax = zmax;
sim.X = X;
sim.Y = Y;
sim.Z = Z;
sim.TH = TH;
sim.N = N;
sim.C = C;
sim.W = W;
sim.f = f;
sim.A = A;
sim.PLOT_COMMS = 0;
sim.PLOT_STATE_3D = PLOT_STATE_3D;
sim.PLOT_STATE_2D = PLOT_STATE_2D;
sim.PLOT_STATE_PHI = 0;
sim.PLOT_STATE_QUALITY = PLOT_STATE_QUALITY;
sim.SAVE_PLOTS = SAVE_PLOTS;






%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%
if PLOT_STATE_3D || PLOT_STATE_2D || PLOT_STATE_QUALITY
	figure
end
tic;
for s=1:smax
	fprintf('%.2f%% complete\n',100*s/smax);

    % ----------------- Partitioning -----------------
    % Coverage quality
    f = fu(Z, zmin, zmax);
    % Sensing patterns
    for i=1:N
        C{i} = ...
            bsxfun(@plus, rot( (Z(i)-zmin+1).*Cb, TH(i) ), [X(i) ; Y(i)]);
    end
    % Communication range %%%%%%%%%%%% FIX THIS %%%%%%%%%%%%
    r_comm = 10*Cd_diameter * ones(size(f));
    
    % Store simulation data
    Xs(s,:) = X;
    Ys(s,:) = Y;
    Zs(s,:) = Z;
	THs(s,:) = TH;
    
    % Sensed space partitioning
    for i=1:N
        if COMM_RANGE
            % Find the nodes in communication range of each node i
            A(i,:) = in_comms_range3( X, Y, Z, i, r_comm(i) );

            % The index of i in the reduced state vector is
            ind = sum(A(i,1:i));

            % Find the cell of each node i based on its neighbors
            W{i} = sensed_partitioning_uniform_cell(region, ...
                C( logical(A(i,:)) ), f( logical(A(i,:)) ), ind);
        else
            % Find the cell of each node i based on all other nodes
            W{i} = sensed_partitioning_uniform_cell(region, C, f, i);
        end
    end
    
    
    % ----------------- Plots -----------------
    sim.X = X;
    sim.Y = Y;
    sim.Z = Z;
    sim.C = C;
    sim.W = W;
    sim.f = f;
    sim.A = A;
    sim.s = s;
    clf
    plot_sim_UAV( sim );
    
    
    % ----------------- Objective -----------------
    % Find covered area and H objective
    for i=1:N
        if ~isempty(W{i})
            cov_area(s) = cov_area(s) + polyarea_nan(W{i}(1,:), W{i}(2,:));
            H(s) = H(s) + f(i) * polyarea_nan(W{i}(1,:), W{i}(2,:));
        end
    end
    
    
    % ----------------- Control law -----------------
    for i=1:N
        % Create anonymous functions for the Jacobians
        % The functions include parameters specific to this node
        Jxy = @(q) J_ellipse_xy(q);
        Jz = @(q) J_ellipse_z(q, X(i), Y(i), Z(i), TH(i), zmin, a, b);
        Jth = @(q) J_ellipse_th(q, X(i), Y(i), Z(i), TH(i), zmin, a, b);
        
        if COMM_RANGE
            % The index of i in the reduced state vector is
            ind = sum(A(i,1:i));

            % Give correct info based on adjacency matrix A
            [uX(i), uY(i), uZ(i)] = ...
                control_uniform(region, zmin, zmax, a, ...
                W(logical(A(i,:))), C(logical(A(i,:))), ...
                f(logical(A(i,:))), ind, X(i), Y(i), Z(i));
        else
            [uX(i), uY(i)] = control_uniform_planar(region, W, C, ...
                f, i, Jxy);
            uZ(i) = control_uniform_altitude(region, W, C, ...
                f, dfu(Z(i), zmin, zmax), i, Jz);
            uTH(i) = control_uniform_rotational(region, W, C, ...
                f, i, Jth);
        end
    end
    
    % Control inputs with gain
    uX = axy * uX;
    uY = axy * uY;
    uZ = az * uZ;
	uTH = ath * uTH;
    
    
    % ----------------- Simulate with ode -----------------
    Tspan = [s*Tstep (s+1)*Tstep];
    IC = [X Y Z TH]';
    u = [uX uY uZ uTH]';
    [T, ode_state] = ode45(@(t,y) DYNAMICS_simple(t, y, u), Tspan, IC);
    
    % Keep the last row of XYZ
    X = ode_state(end, 1:N );
    Y = ode_state(end, N+1:2*N );
    Z = ode_state(end, 2*N+1:3*N );
    TH = ode_state(end, 3*N+1:4*N );
    
end
elapsed_time = toc;
average_iteration = elapsed_time / smax;
fprintf('\nSimulation time: %.4f s\n', elapsed_time)
fprintf('Average iteration time: %.4f s\n', average_iteration)




%%%%%%%%%%%%%%%%%%% Final plots %%%%%%%%%%%%%%%%%%%
% Plot covered area
figure;
plot( Tstep*linspace(1,smax,smax), 100*cov_area/region_area, 'b');
hold on
area_opt = 100 * N * pi * (zopt * tan(a))^2 / region_area;
plot( Tstep*[1 smax], [area_opt area_opt], 'k--');
axis([0 Tstep*smax 0 100]);
% axis([0 Tstep*smax 0 140]);
h = xlabel('$Time ~(s)$');
set(h,'Interpreter','latex')
h = ylabel('$A_{cov}~(\%)$');
set(h,'Interpreter','latex')

% Plot objective
H_opt = N * pi * (zopt * tan(a))^2 * fu(zopt, zmin, zmax);
figure;
plot( Tstep*linspace(1,smax,smax), 100*H / H_opt, 'b');
hold on
plot( Tstep*[1 smax], [100 100], 'k--');
axis([0 Tstep*smax 0 100]);
h = xlabel('$Time ~(s)$');
set(h,'Interpreter','latex')
h = ylabel('$\frac{\mathcal{H}}{\mathcal{H}_{opt}} ~(\%)$');
set(h,'Interpreter','latex')

% Save trajectories
traj = zeros(3,smax,N);
traj(1,:,:) = Xs;
traj(2,:,:) = Ys;
traj(3,:,:) = Zs;

%%%%%%%%%%%%%%%%%%% Save Results %%%%%%%%%%%%%%%%%%%
if SAVE_RESULTS
    filename = ...
        strcat( 'results_uniform_', ...
        datestr(clock,'yyyymmdd_HHMM') , '.mat' );
    save(filename);
end
