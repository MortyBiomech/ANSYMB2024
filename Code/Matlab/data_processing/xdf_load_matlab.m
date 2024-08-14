function streams = xdf_load_matlab(subject_id, rawdata_path)
    
    files_path = [rawdata_path, 'sub-', num2str(subject_id), filesep];

    %% Find the numbers of sessions and load streams
    % Get the list of all items in the folder
    items = dir(files_path);
    streams = repmat(cell(1,1), 1, length(items));
   
    % Loop through the items and load xdf files
    for k = 1:length(items)
        
        % Check if the item is a directory and not '.' or '..'
        if items(k).isdir && ~strcmp(items(k).name, '.') ...
                && ~strcmp(items(k).name, '..')
            
            streams_path = [files_path, items(k).name, filesep, ...
                'eeg', filesep];
            xdf_files = dir(streams_path);
            
            % Check if the file name does not contain 'old' in it
            for j = 1:length(xdf_files)
                
                if ~xdf_files(j).isdir ...
                        && ~contains(xdf_files(j).name, 'old')
                    streams{1,k} = struct('dataset', ...
                        load_xdf(fullfile(streams_path, ...
                        xdf_files(j).name)));
                end

            end

        end

    end

    %% Remove empty cell members
    % Create a logical array where true indicates an empty cell
    emptyCells = cellfun(@isempty, streams);
    
    % Remove the empty cells
    streams(emptyCells) = [];

 
end