function min_max_N = find_min_of_max_N_in_subfolders(parent_folder)
    % Check if the parent folder exists
    if ~isfolder(parent_folder)
        error('The specified folder does not exist.');
    end
    
    % Get a list of subfolders within the parent folder
    subfolders = dir(parent_folder);
    subfolders = subfolders([subfolders.isdir]);  % Keep only directories
    
    % Initialize a variable to store the minimum of max N values
    min_max_N = inf;
    
    % Loop through each subfolder
    for i = 1:length(subfolders)
        subfolder_name = subfolders(i).name;
        
        % Ignore '.' and '..'
        if strcmp(subfolder_name, '.') | strcmp(subfolder_name, '..')
            continue;
        end
        
        % Create the full path to the subfolder (e.g., C:\...\2_2 or C:\...\i)
        subfolder_path = fullfile(parent_folder, subfolder_name);
        
        % Check if the subfolder name matches the expected pattern (e.g., '2_2' or 'i')
        if ~isempty(regexp(subfolder_name, '^\d+_\d+$', 'once')) || ~isempty(regexp(subfolder_name, '^\d+$', 'once'))
            % Check for the existence of a '0' subfolder inside this subfolder
            zero_subfolder_path = fullfile(subfolder_path, '0');
            if isfolder(zero_subfolder_path)
                % Use the previous function to find the maximum N in this '0' subfolder
                [max_N, ~] = find_max_N_and_job(zero_subfolder_path);
                
                % Update min_max_N if necessary
                if max_N < min_max_N
                    min_max_N = max_N;
                end
            end
        end
    end
    
    % If no valid subfolders were found, set min_max_N to NaN
    if min_max_N == inf
        min_max_N = NaN;
    end
end
