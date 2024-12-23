function vector_with_interp = interpolate_nans(data_vector)

  % Find indices of non-NaN values
  valid_indices = ~isnan(data_vector);

  % Handle edge cases (all NaNs or single element)
  if all(isnan(data_vector)) || numel(data_vector) == 1
    vector_with_interp = data_vector;  % Return original vector for these cases
    return;
  end

  % Interpolate for each NaN index
  vector_with_interp = data_vector;
  INDEX_not_nan=find(~isnan(data_vector));
  for i = find(isnan(data_vector))
    % Find indices of closest non-NaN values (before and after)
    closest_before = find(INDEX_not_nan<i, 1, 'last');  
    closest_before=INDEX_not_nan(closest_before);
    I2 = find(INDEX_not_nan>i, 1, 'first');
    closest_after=INDEX_not_nan(I2);
  

    % Check for edge cases (NaN at the beginning or end)
    if isempty(closest_before)
      closest_before = 1;
    end
    if isempty(closest_after)
      closest_after = length(data_vector);
    end

    % Perform linear interpolation
    interp_value = interp1([closest_before closest_after], ...
                            [data_vector(closest_before) data_vector(closest_after)], ...
                            i, 'linear');
    vector_with_interp(i) = interp_value;
  end
end

