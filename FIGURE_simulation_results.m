% Plot results of simulation
clear variables
close all

% Add function path
addpath( genpath('Functions') );

% Load results
% load('uniform_results');
% load('uniform_results_3_nodes');
% load('uniform_results_9_nodes');
load('results_uniform_phi_20170121_2329');


%%%%%%%%%%%%%%% Select plots %%%%%%%%%%%%%%%
AREA = 1;
OBJECTIVE = 1;
TRAJECTORIES = 1;
INITIAL_QUALITY = 1;
FINAL_QUALITY = 1;
INITIAL_STATE = 1;
FINAL_STATE = 1;
USE_PHI = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gridsize = 50;
disk_rad = 0.02;


%%%%%%%%%%%%%%% Plot covered area %%%%%%%%%%%%%%%
if AREA
	figure
	plot( Tstep*linspace(1,smax,smax), 100*cov_area, 'b');
	hold on
	area_opt = 100 * N * pi * (zopt * tan(a))^2 / region_area;
	plot( Tstep*[1 smax], [area_opt area_opt], 'k--');
	axis([0 Tstep*smax 0 100]);
	% axis([0 Tstep*smax 0 140]);
	h = xlabel('$Time ~(s)$');
	set(h,'Interpreter','latex')
	h = ylabel('$A_{cov}~(\%)$');
	set(h,'Interpreter','latex')
end

%%%%%%%%%%%%%%% Plot objective %%%%%%%%%%%%%%%
if OBJECTIVE
    if USE_PHI
        figure
        plot( Tstep*linspace(1,smax,smax), H, 'b');
        axis([0 Tstep*smax 0 ceil(max(H))]);
        h = xlabel('$Time ~(s)$');
        set(h,'Interpreter','latex')
        h = ylabel('$\mathcal{H}$');
        set(h,'Interpreter','latex')
    else
        H_opt = N * pi * (zopt * tan(a))^2 * fi_u(zopt, zmin, zmax);
        figure
        plot( Tstep*linspace(1,smax,smax), 100*H / H_opt, 'b');
        hold on
        plot( Tstep*[1 smax], [100 100], 'k--');
        axis([0 Tstep*smax 0 100]);
        h = xlabel('$Time ~(s)$');
        set(h,'Interpreter','latex')
        h = ylabel('$\frac{\mathcal{H}}{\mathcal{H}_{opt}} ~(\%)$');
        set(h,'Interpreter','latex')
    end
end

%%%%%%%%%%%%%%% Plot trajectories %%%%%%%%%%%%%%%
if TRAJECTORIES
	traj = zeros(3,smax,N);
	traj(1,:,:) = Xs;
	traj(2,:,:) = Ys;
	traj(3,:,:) = Zs;
	figure;
	plot3( [Xb Xb(1)], [Yb Yb(1)], zeros(length(Xb)+1), 'k');
	hold on
	for i=1:N
		% Initial
		plot3( traj(1,1,i), traj(2,1,i), traj(3,1,i), 'rs');
		hold on
		% Trajectory
		plot3( traj(1,:,i), traj(2,:,i), traj(3,:,i), 'r');
		hold on
		% Final
		plot3( traj(1,smax,i), traj(2,smax,i), traj(3,smax,i), 'ro');
		hold on

		% Projections
		% Initial
		plot3( traj(1,1,i), traj(2,1,i), 0, 'ks');
		hold on
		% Trajectory projection
		plot3( traj(1,:,i), traj(2,:,i), zeros(size(traj(3,:,i))), 'k--');
		hold on
		% Final
		plot3( traj(1,smax,i), traj(2,smax,i), 0, 'ko');
		hold on
	end
	axis equal
	view([-16 16])
	% axis([0 3 0 2.5 0 zmax]);
	xlabel('x');
	ylabel('y');
	zlabel('z');
end

%%%%%%%%%%%%%%% Plot initial quality %%%%%%%%%%%%%%%
if INITIAL_QUALITY
	f = zeros(1, N);
	C = cell([1 N]);
	X = squeeze(traj(:,1,:));
	R = squeeze(Rs(1,:));
	% Sensing disks and quality
	for i=1:N
		f(i) = fu(X(3,i), zmin, zmax);
		C{i} = [X(1,i) + R(i) * cos(t) ; X(2,i) + R(i) * sin(t)];
	end
	% Find cells
	[W, ~] = sensed_partitioning_uniform(Xb, Yb, C, f);
	
	figure
	hold on
	% Plot cylinders
	for i=1:N
		plot3_cell_quality(W{i}, f(i), 'r');
	end
	% Plot region
	plot3( [Xb Xb(1)], [Yb Yb(1)], zeros(length(Xb)+1), 'k');
	axis equal
	axis([0 3 0 2.5 0 1]);
	view([-16 26])
	xlabel('x');
	ylabel('y');
	zlabel('f');
end

%%%%%%%%%%%%%%% Plot final quality %%%%%%%%%%%%%%%
if FINAL_QUALITY
	f = zeros(1, N);
	C = cell([1 N]);
	X = squeeze(traj(:,end,:));
	R = squeeze(Rs(end,:));
	% Sensing disks and quality
	for i=1:N
		f(i) = fu(X(3,i), zmin, zmax);
		C{i} = [X(1,i) + R(i) * cos(t) ; X(2,i) + R(i) * sin(t)];
	end
	% Find cells
	[W, ~] = sensed_partitioning_uniform(Xb, Yb, C, f);
	
	figure
	hold on
	% Plot cylinders
	for i=1:N
		plot3_cell_quality(W{i}, f(i), 'r');
	end
	% Plot region
	plot3( [Xb Xb(1)], [Yb Yb(1)], zeros(length(Xb)+1), 'k');
	axis equal
	axis([0 3 0 2.5 0 1]);
	view([-16 26])
	xlabel('x');
	ylabel('y');
	zlabel('f');
end

%%%%%%%%%%%%%%% Plot initial state %%%%%%%%%%%%%%%
if INITIAL_STATE
	f = zeros(1, N);
	C = cell([1 N]);
	X = squeeze(traj(:,1,:));
	R = squeeze(Rs(1,:));
	% Sensing disks and quality
	for i=1:N
		f(i) = fu(X(3,i), zmin, zmax);
		C{i} = [X(1,i) + R(i) * cos(t) ; X(2,i) + R(i) * sin(t)];
	end
	% Find cells
	[W, ~] = sensed_partitioning_uniform(Xb, Yb, C, f);
	
	figure
    % Plot point importance
    if USE_PHI
        plot_phi( phi , region, gridsize );
    end
	hold on
	% Plot sensing
	for i=1:N
		plot_poly(C{i}, 'r--');
		hold on
		tmpc = [X(1,i) + disk_rad * cos(t) ; X(2,i) + disk_rad * sin(t)];
		fill( tmpc(1,:), tmpc(2,:), 'k', 'EdgeColor', 'none' );
	end
	% Plot cells
	for i=1:N
		plot_poly(W{i}, 'k');
		hold on
	end
	% Plot region
	plot( [Xb Xb(1)], [Yb Yb(1)], 'k');
	axis equal
	axis([0 3 0 2.5]);
    view(0, 90);
	axis off
end

%%%%%%%%%%%%%%% Plot final state %%%%%%%%%%%%%%%
if FINAL_STATE
	f = zeros(1, N);
	C = cell([1 N]);
	X = squeeze(traj(:,end,:));
	R = squeeze(Rs(end,:));
	% Sensing disks and quality
	for i=1:N
		f(i) = fu(X(3,i), zmin, zmax);
		C{i} = [X(1,i) + R(i) * cos(t) ; X(2,i) + R(i) * sin(t)];
	end
	% Find cells
	[W, ~] = sensed_partitioning_uniform(Xb, Yb, C, f);
	
	figure
    % Plot point importance
    if USE_PHI
        plot_phi( phi , region, gridsize );
    end
	hold on
	% Plot sensing
	for i=1:N
		plot_poly(C{i}, 'r--');
		hold on
		tmpc = [X(1,i) + disk_rad * cos(t) ; X(2,i) + disk_rad * sin(t)];
		fill( tmpc(1,:), tmpc(2,:), 'k', 'EdgeColor', 'none' );
	end
	% Plot cells
	for i=1:N
		plot_poly(W{i}, 'k');
		hold on
	end
	% Plot region
	plot( [Xb Xb(1)], [Yb Yb(1)], 'k');
	axis equal
	axis([0 3 0 2.5]);
    view(0, 90);
	axis off
end
