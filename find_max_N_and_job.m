function [max_N, max_job_id] = find_max_N_and_job(folder_path)
    % Check if the input folder path exists
    if ~isfolder(folder_path)
        error('The specified folder does not exist.');
    end
    
    % Get the list of folders in the directory (ignoring '.' and '..')
    files = dir(folder_path);
    folders = files([files.isdir]);  % Keep only directories

    % Initialize variables to hold the maximum N and job_id values
    max_N = -inf;
    max_job_id = -inf;
    
    % Define the pattern to extract N and job_id from the folder names
    pattern = '_(\d+)_job_id_(\d+)';
    
    % Loop over all the directories in the folder
    for i = 1:length(folders)
        folder_name = folders(i).name;
        
        % Ignore '.' and '..'
        if strcmp(folder_name, '.') || strcmp(folder_name, '..')
            continue;
        end
        
        % Apply regular expression to extract N and job_id
        tokens = regexp(folder_name, pattern, 'tokens');
        
        if ~isempty(tokens)
            % Extract N and job_id values
            N = str2double(tokens{1}{1});
            job_id = str2double(tokens{1}{2});
            
            % Update max values if necessary
            if N > max_N
                max_N = N;
            end
            if job_id > max_job_id
                max_job_id = job_id;
            end
        end
    end
    
    % If no valid folders were found, set max_N and max_job_id to NaN
    if max_N == -inf
        max_N = NaN;
    end
    if max_job_id == -inf
        max_job_id = NaN;
    end
    
    % Return the maximum values
end
