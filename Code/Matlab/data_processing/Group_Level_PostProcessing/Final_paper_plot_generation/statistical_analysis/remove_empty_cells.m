function data_clean = remove_empty_cells(data)
    non_empty_idx = ~cellfun(@isempty, data);
    data_clean = data(non_empty_idx);
end