
tic;
m=500;
n=10;
label = zeros(m, m);
euler_angles = zeros(m, m, 3);
R = [];                               %tracks the radius of each growing nuclei
z = 0;                                %labeller
growing_nuclei = [];                  %stores address of all the growing nuclei

while nnz(label(:, :))~= m*m
    new_nuclei = Uncrystallized_Nuclei_Gen(label, n);     %random nuclei ready for growth
    for i = 1:length(new_nuclei)
        z = z + 1;
        euler_angles(new_nuclei(i, 1), new_nuclei(i, 2), :) = rand(1, 3);
        label(new_nuclei(i, 1), new_nuclei(i, 2)) = z;
    end
    growing_nuclei = [growing_nuclei; new_nuclei];             %new ones appended to the list for growing with prior nuclei
    R = vertcat(R, ones(n, 1));                                %appends the growth radius of the new nuclei ie, 1
    
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
    end
    R(1:z, 1) = R(1:z, 1) + 1;
    f= imshow(euler_angles, 'initialmagnification', 'fit');
    pause(1);
    %after this go for picking the new ones up
end
toc;