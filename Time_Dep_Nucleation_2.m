clear all;
close all;
clc;
tic;
m=1000;                                           %matrix size
n=20;                                             %nuclei to be picked
label = zeros(m, m);
euler_angles = zeros(m, m, 3);
timestep = 0;                                     %interval for dropping next set of nuclei
R = [];                                           %tracks the radius of each growing nuclei
z = 0;                                            %labeller
growing_nuclei = [];                              %stores address of all the growing nuclei

[new_nuclei, avl_sites] = Uncrystallized_Nuclei_Gen(label, n);         %random nuclei ready for growth
fprintf("\n sites to be nucleated = %d \n", avl_sites);
for i = 1:length(new_nuclei)
    z = z + 1;
    euler_angles(new_nuclei(i, 1), new_nuclei(i, 2), :) = rand(1, 3);
    label(new_nuclei(i, 1), new_nuclei(i, 2)) = z;
end
growing_nuclei = [growing_nuclei; new_nuclei];                         %new ones appended to the list for growing with prior nuclei
R = vertcat(R, ones(length(new_nuclei), 1));                           %appends the growth radius of the new nuclei ie, 1
%after this first gen of nuclei created , z=3

%% time based growth 
while nnz(label(:, :))~= m*m
    timestep = timestep + 1;
    for j = 1:length(growing_nuclei)
        r_x = growing_nuclei(j, 1);
        r_y = growing_nuclei(j, 2);
        for x = -R(j) : R(j)
            for y = -R(j) : R(j)
                if r_x + x > 0 && r_x + x <= m && r_y + y > 0 &&...
                        r_y + y <= m && x^2 + y^2 <= (R(j))^2 && label(r_x + x, r_y + y) == 0
                    label(r_x + x, r_y + y) = label(r_x, r_y);
                    euler_angles(r_x + x, r_y + y, :) = euler_angles(r_x, r_y, :);
                end
            end
        end
        %each iteration grows the nuclei to the target radius
    end
    R(1:z, 1) = R(1:z, 1) + 1;                                          %target radii increment
    f= imshow(euler_angles, 'initialmagnification', 'fit');
    pause(10^-8);
    %fprintf('%d', timestep);
    
    if timestep >= 8
        [new_nuclei,avl_sites] = Uncrystallized_Nuclei_Gen(label, n);   
        fprintf('\n sites to be nucleated == %d \n', avl_sites);
        for i = 1:length(new_nuclei)
            z = z + 1;
            euler_angles(new_nuclei(i, 1), new_nuclei(i, 2), :) = rand(1, 3);
            label(new_nuclei(i, 1), new_nuclei(i, 2)) = z;
        end
        growing_nuclei = [growing_nuclei; new_nuclei];
        R = vertcat(R, ones(length(new_nuclei), 1));          
        timestep = 0;
    end
end
%% plotting
Areas = zeros(z, 1);
dia = zeros(z, 1);
Area_fraction = zeros(z, 1);
for i = 1:z
    [u, v] = find(label == i);
    Areas(i, 1) = length(u);
    Area_fraction(i, 1) = (Areas(i, 1))/ m/m;
    dia(i, 1) = sqrt((4* Areas(i, 1))/ pi);
end
toc;