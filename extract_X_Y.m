function [X, Y] = extract_X_Y(input_str)
    % Define the pattern to extract X and Y (X_T_Y_P structure)
    pattern = '_(\d+)_T_(\d+)_P';
    
    % Apply regular expression to extract X and Y
    tokens = regexp(input_str, pattern, 'tokens');
    
    % If match is found, extract X and Y
    if ~isempty(tokens)
        X = str2double(tokens{1}{1});
        Y = str2double(tokens{1}{2});
    else
        error('X_T_Y_P structure not found in the input string.');
    end
end
