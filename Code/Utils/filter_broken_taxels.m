function [data_aux, removed_index]=filter_broken_taxels(data)
% This function takes as input a matrix of capacitance data (number datapoints X number of taxels)
% and returns a modified data structure in which the data corresponding to
% mulfunctioning taxels are substituted by zeros. The removed indexes are
% also returned.
%
% author: Italo Belli (i.belli@tudelft.nl), 2021
    
    data_aux = data;
    removed_index = [];
    num_taxels = size(data_aux,2);
    
    % REMOVING LOW VALUES
    % Heuristically, we observed that taxels that were broken presented
    % occasional drops of the recorded capacitance. We have set the minimum
    % threshold to 100 as this allowed to capture and exclude taxels which
    % were behaving incorrectly.
    minimum_val = min(data_aux);
    index = find (minimum_val < 100);
    if(~isempty(index))
        data_aux(:, index) = zeros(size(data_aux(:, index)));
        for i=1:size(index,2)
            if ~any(removed_index==index(i))
                removed_index = [removed_index, index(i)];
            end
        end
    end
    
    % REMOVING SATURATED TAXELS
    % the maximum digital capacitance value must be below 250, 
    % otherwise the taxel is considered saturated
    maximum_val = max(data_aux);
    index = find(maximum_val > 250);
    if(~isempty(index))
        data_aux(:, index)=zeros(size(data_aux(:, index)));
        for i=1:size(index,2)
            if ~any(removed_index==index(i))
                removed_index = [removed_index, index(i)];
            end
        end
    end
    
    % removing inactive taxles
    avg_val = mean(data_aux);
    index = find(abs(avg_val-240)<0.01);
    if(~isempty(index))
        data_aux(:, index)=zeros(size(data_aux(:, index)));
        for i=1:size(index,2)
            if ~any(removed_index==index(i))
                removed_index = [removed_index, index(i)];
            end
        end
    end
    
    % REMOVING TAXELS THAT SUDDENLY STOPPED WORKING
    % We analyze whether the value of the digital capacitance remains
    % constant for too long (establishing a window_span), or changes too quickly. 
    % If that's the case, the taxel stopped working and must be ignored
    window_span = 5;
    for index=1:num_taxels
        index_valid = 1;
        for sample=2:size(data_aux,1)-window_span
            if (index_valid && data_aux(sample-1, index)==data_aux(sample, index))
                problem_detected = 1;
                for i=1:window_span
                    if(data_aux(sample, index)~=data_aux(sample+i, index))
                        problem_detected = 0;
                    end
                end
                if(problem_detected)
                    index_valid = 0;
                    removed_index = [removed_index, index];
                    data_aux(:, index)=zeros(size(data_aux(:, index)));
                end
            end
            if (index_valid && data_aux(sample, index)-data_aux(sample+window_span, index)>8)
                problem_detected = 1;
                index_valid = 0;
                removed_index = [removed_index, index];
                data_aux(:, index)=zeros(size(data_aux(:, index)));
            end
        end
    end

    % REMOVING TAXELS WITH LARGE DISCONTINUITIES
    % Not used for now
    for index=1:num_taxels
        prev_sample= 1;
        index_valid = 1;
        for sample=2:size(data_aux,1)
            if(index_valid && abs(data_aux(sample,index)-data_aux(prev_sample,index))>1000)
                index_valid = 0;
                removed_index = [removed_index, index];
                data_aux(:, index)=zeros(size(data_aux(:, index)));
            end
            prev_sample = sample;
        end
    end
    
    
    removed_index = sort(unique(removed_index));
    
end

