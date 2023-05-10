function [data1_new, data2_new] = filter_high_variation_data(data1, data2, eps, step)
% this function takes as inputs 2 data matrices:
% data1 has dimensions m x 1
% data2 has dimensions m x n
%
% a condition is checked on data1 (pressure), to understand if it varies too much
% (more than eps from one data sample before). If this is the case,
% the next data sample is discarded, and also the corresponding elements in
% data2 (capacitance) are removed, to keep the same length of the two arrays.
%
% Author: Ines Sorrentino (i.sorrentino@iit.it), 2020
configurationfile;

list_of_samples = [];
index = 1;
for row = 1 : step : (size(data1, 1)-step)
  if abs(data1(row + step) - data1(row)) < eps
    list_of_samples(index:index+step,1) = row : (row+step);
    index = index + step;
  end
end

data1_new(:,1) = data1(list_of_samples);
data2_new = data2(list_of_samples, :);

end