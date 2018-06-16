%% Defining thermomechanical processing condition
% This portion will changed on changing the processing condition.
tic;
T = 1373;                                                   % Temperature in Kelvin K;
strain_rate = .001;                                         % in s^-1
k1 = 6.75e8;                                                % Derived from experimental flow curves
k2 = 59.1;                                                  % Derived from experimental flow curves
saturation_stress = 66.28426e6;
steady_state_stress  = 44e6;
start_stress = 55.061e6;

%% Materials Constants

% General Constants
b = 2.5e-10;                                                % burgers vector (m);
T_m = 1609;                                                 % Melting point (K)
Q_act= 450000;                                              % Activation Energy (J/mol)
R = 8.314;                                                  % Universal Gas Constant (J/mol-K)
k = 1.381e-23;                                              % Boltzman Constant (m2 kg s-2 K-1)
G = 8.1e10*(1-(.64*(T-300)/(T_m)));                         % every temp is in Kelvin % Shear modulus (N/m^2)

%Constants for critical dislocation density calculation
gamma_m = .625;                                             % Maximum grain boundary energy (J/m^2)
theta_m = 15;                                               % High angle grain boundary misorientation (degrees)
delta_d_ob = 1.34e-11;                                      % Boundary self diffusion coefficient (m^3/s)
Q_b = 287000;                                               % Boundary seif diffusion activation energy (J/mol)

%Geometrical Constants
Cd = 1;                                                     % In microns Cell dimension

%Derived Constants
M = delta_d_ob*exp(-Q_b/R/T)*b/k/T;                         % Grain Boundary Mobility
tau = 0.5*G*b*b;                                             % Dislocation line energy
rho_not = (start_stress/0.5/b/G)^2;

%% To construct the initial microstructure at t = 0

m = 600;
n = 120;
rng('shuffle')
initial_nuclei = randi([1 m], n, 2);                 %nuclei at random posn needed for crystallization
euler_angles = zeros(m, m, 3);                       %stores RGB values
colors = rand(n, 3);                                 %supplies RGB values
orientation_matrix = zeros(m, m);                    %orientation of each grain
orient = round(rand(n, 1)*90 );
status = zeros(m , m);                               %crystallized or not
label = zeros(m, m);
neighbor_store = zeros(n,15);

%initial pixels which will give birth to grains
for i = 1:n
    euler_angles(initial_nuclei(i, 1), initial_nuclei(i, 2), :) = colors(i, :);
    orientation_matrix(initial_nuclei(i, 1), initial_nuclei(i, 2)) = orient(i,1);
    label(initial_nuclei(i,1), initial_nuclei(i,2)) = i;                           %these are colored and hence will not be again colored in the next block's loop
end
%imshow(euler_angles, 'initialmagnification', 'fit');

%mapping of nuclei with color
mapping_matrix = [initial_nuclei,colors];

%circular growth of the nuclei
w = 1;
while (nnz(euler_angles(:, :, 1)) ~= m*m)
    r = w;
    for x = -r:r
        for y = -r:r
            for i = 1:n
                r_x = initial_nuclei(i, 1);
                r_y = initial_nuclei(i, 2);
                if r_x + x >0 && r_x + x <=m && r_y + y >0 && r_y + y <=m &&...
                        euler_angles(r_x + x, r_y + y, 1) == 0 && x^2 + y^2 <= r^2 ...
                        && status(r_x + x, r_y + y) == 0
                    
                    euler_angles(r_x + x, r_y + y, :) = colors(i, :);
                    orientation_matrix(r_x + x, r_y + y) = orient(i, 1);
                    status(r_x+x, r_y+y) = 1;
                    label(r_x+x, r_y+y) = i;
                end
            end
        end
    end
    w = w+1;
    imshow(euler_angles,'initialmagnification', 'fit');
    pause(0.01);
end

%Storing the neighbor of each grain
for i1 = 2:m-1
    for j1 = 2:m-1
        q = 1;
        dis_similar = [];
        distinct_neighbors = [];
        z = label(i1, j1);
        for a1 = -1:1
            for b1 = -1:1
                if label(i1 + a1, j1 + b1) ~= z
                    dis_similar(1, q) =...
                        label(i1 + a1, j1 + b1);
                    q = q + 1;
                end
            end
        end
        if isempty(dis_similar) == 0                                      %only when u get some dissimilars go in
            
            distinct_neighbors = unique(dis_similar);                     %stores sorted,distinct,dissimilar neighbors 
            
            if nnz(neighbor_store(z, :)) == 0                             %first time filling for that label
                neighbor_store(z, 1:length(distinct_neighbors))...
                    = distinct_neighbors(:, :);
            else                                                          %if not, check out for no repeat filling in store_neighbor
                for i2 = 1:length(distinct_neighbors)
                    if ismember(distinct_neighbors(1, i2), neighbor_store(z, :))
                    else
                        neighbor_store(z, nnz(neighbor_store(z, :)) + 1)...
                            = distinct_neighbors(1, i2);
                    end
                end
            end
        end
        %after here increment j1 up
    end
end
toc;

%% Recovery modelling
max_strain = 0.65;
strain = 0.02;
time_step = (max_strain / strain_rate)/ 100;   
del_strain = strain_rate* time_step;

dislocation_density_matrix(1:m, 1:m) = rho_not;                         %uniform dislocation density for the primary grains
avg_dislocation_density = sum(sum(dislocation_density_matrix))/ (m* m);
critical_dislocation_density = avg_dislocation_density + 1;
t = 1;

while(avg_dislocation_density < critical_dislocation_density)
    
    del_dislocation_density = ((k1* sqrt(avg_dislocation_density))-...
        (k2* avg_dislocation_density))* del_strain;
    
    dislocation_density_matrix(:) =...
        dislocation_density_matrix(:) + del_dislocation_density;
    
    avg_dislocation_density =...
        sum(sum(dislocation_density_matrix))/(m* m);
    
    l = 20/ sqrt(avg_dislocation_density);                              %mean free path for the dislocation movement
    
    flow_stress = 0.5* G* b* sqrt(avg_dislocation_density);
    
    critical_dislocation_density =...
        ((20*gamma_m*strain_rate)/(3*b*l*M*tau^2))^(1/3);
    
    data(t, 1) = strain;
    data(t, 2) = flow_stress;
    strain = strain + del_strain;
    %DD_avg(t,1) = avg_dislocation_density;
    %dislocation_increment(t, 1) = del_dislocation_density ;
    %critical_value(t, 1) = critical_dislocation_density;
    t = t+1;
    
    if strain >= 0.65 && avg_dislocation_density < critical_dislocation_density
        fprintf('criteria not met')
        break
    end
end
figure;
plot(data(:,1), data(:,2), '-o');



