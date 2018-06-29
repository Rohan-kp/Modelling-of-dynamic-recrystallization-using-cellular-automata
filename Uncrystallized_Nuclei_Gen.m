function [new_nuclei] = Uncrystallized_Nuclei_Gen(label, n)
new_nuclei = [];                                           %stores randomly selected nuclei to be grown
unlabelled_pixels = [];                                    %stores the unlabelled pixels and then randomly thrown into new_nuclei matrix

[unlabelled_pixels(:, 1), unlabelled_pixels(:, 2)] = find(label==0);
%supply > demand
if  length(unlabelled_pixels) > n
    select_pixels = randperm(length(unlabelled_pixels), n);    %gives random row no.
    new_nuclei(:, 1) = unlabelled_pixels(select_pixels, 1);
    new_nuclei(:, 2) = unlabelled_pixels(select_pixels, 2);
    %supply < demand
else
    new_nuclei(:, 1) = unlabelled_pixels(:, 1);
    new_nuclei(:, 2) = unlabelled_pixels(:, 2);
end
end