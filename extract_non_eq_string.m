function non_eq_str = extract_non_eq_string(input_str)
    % Function to extract and form a new string based on the pattern x1_t_x2_A
    
    % Use regular expression to check for the pattern x1_t_x2_A
    % The pattern we're looking for is something like '_x1_t_x2_A', where x1 and x2
    % can be integers or decimals (like 8 or 8_5)
    pattern = '_([\d_]+)_t_([\d_]+)_A';

    % Check if the pattern exists in the string
    tokens = regexp(input_str, pattern, 'tokens');
    
    if ~isempty(tokens)
        % If the pattern is found, raise a flag (implicit since we are returning a result)
        
        % Extract x1 and x2 from the tokens
        x1 = tokens{1}{1};
        x2 = tokens{1}{2};
        
        % Create the new string 'non_eq_str' as required
        non_eq_str = ['non_eq_' x1 '_t_' x2 '_A'];
    else
        % If the pattern is not found, return an empty string or a message
        non_eq_str = '';
        disp('No matching pattern found.');
    end
end
