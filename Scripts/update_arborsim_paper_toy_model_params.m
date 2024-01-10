function update_arborsim_paper_toy_model_params(model_root_dir, model_name)
    % This function updates the ligament parameters in the model to match the
    % values used in the ArborSim paper. 

    %% Adjustments for the tendon branch count category
    if strcmp(model_name, 'Tendon_Branch_Count_1')
        ligament_names = {'group1_ligament1', 'group1_ligament2'};
        rest_lens = {1.4170, 0.0670};
    end

    if strcmp(model_name, 'Tendon_Branch_Count_2')
        ligament_names = {'group1_ligament1', 'group1_ligament2', 'group1_ligament3', 'group1_ligament4'};
        rest_lens = {0.0670, 0.2250, 1.1920, 0.0670};
    end

    if strcmp(model_name, 'Tendon_Branch_Count_3')
        ligament_names = {'group1_ligament1', 'group1_ligament2', 'group1_ligament3', 'group1_ligament4', 'group1_ligament5', ...
                          'group1_ligament6'};
        rest_lens = {0.0670, 0.2250, 0.2250, 0.0670, 0.0670, ...
                     0.9670};
    end

    if strcmp(model_name, 'Tendon_Branch_Count_4')
        ligament_names = {'group1_ligament1', 'group1_ligament2', 'group1_ligament3', 'group1_ligament4', 'group1_ligament5', ...
                          'group1_ligament6', 'group1_ligament7', 'group1_ligament8'};
        rest_lens = {0.0670, 0.2250, 0.2250, 0.0670, 0.0670, ...
                     0.2250, 0.0670, 0.7420};

    end
    
    if strcmp(model_name, 'Tendon_Branch_Count_5')
        ligament_names = {'group1_ligament1', 'group1_ligament2', 'group1_ligament3', 'group1_ligament4', 'group1_ligament5', ...
                          'group1_ligament6', 'group1_ligament7', 'group1_ligament8', 'group1_ligament9', 'group1_ligament10'};
        rest_lens = {0.0670, 0.2250, 0.2250, 0.0670, 0.0670, ...
                     0.2250, 0.0670, 0.2250, 0.0670, 0.5170};
    end
    
    if strcmp(model_name, 'Tendon_Branch_Count_6')
        ligament_names = {'group1_ligament1', 'group1_ligament2', 'group1_ligament3', 'group1_ligament4', 'group1_ligament5', ...
                          'group1_ligament6', 'group1_ligament7', 'group1_ligament8', 'group1_ligament9', 'group1_ligament10', ...
                          'group1_ligament11', 'group1_ligament12'};
        rest_lens = {0.0670, 0.2250, 0.2250, 0.0670, 0.0670, ...
                     0.2250, 0.0670, 0.2250, 0.0670, 0.2250, ...
                     0.0670, 0.2920};
    end
    
    %% Adjustments for the joint count bewteen insertions category
    if strcmp(model_name, 'Joint_Count_Between_Insertions_1')
        ligament_names = {'group1_ligament1', 'group1_ligament2'};
        rest_lens = {1.1920, 1.4170};
    end
    
    if strcmp(model_name, 'Joint_Count_Between_Insertions_2')
        ligament_names = {'group1_ligament1', 'group1_ligament2'};
        rest_lens = {0.9670, 1.4170};
    end

    if strcmp(model_name, 'Joint_Count_Between_Insertions_3')
        ligament_names = {'group1_ligament1', 'group1_ligament2'};
        rest_lens = {0.7420, 1.4170};
    end

    if strcmp(model_name, 'Joint_Count_Between_Insertions_4')
        ligament_names = {'group1_ligament1', 'group1_ligament2'};
        rest_lens = {1.4170, 0.5170};
    end

    if strcmp(model_name, 'Joint_Count_Between_Insertions_5')
        ligament_names = {'group1_ligament1', 'group1_ligament2'};
        rest_lens = {1.4170, 0.2920};
    end

    if strcmp(model_name, 'Joint_Count_Between_Insertions_6')
        ligament_names = {'group1_ligament1', 'group1_ligament2'};
        rest_lens = {1.4170, 0.0670};
    end
    
    %% Adjustments for the fiber length to MTU length ratio category
    if strcmp(model_name, 'Fiber_Length_MTU_Length_Ratio_1')
        ligament_names = {'group1_ligament1', 'group1_ligament2'};
        rest_lens = {1.1920, 1.4170};
    end
    
    if strcmp(model_name, 'Fiber_Length_MTU_Length_Ratio_2')
        ligament_names = {'group1_ligament1', 'group1_ligament2'};
        rest_lens = {0.9670, 1.1920};
    end

    if strcmp(model_name, 'Fiber_Length_MTU_Length_Ratio_3')
        ligament_names = {'group1_ligament1', 'group1_ligament2'};
        rest_lens = {0.7420, 0.9670};
    end

    if strcmp(model_name, 'Fiber_Length_MTU_Length_Ratio_4')
        ligament_names = {'group1_ligament1', 'group1_ligament2'};
        rest_lens = {0.5170, 0.7420};
    end

    if strcmp(model_name, 'Fiber_Length_MTU_Length_Ratio_5')
        ligament_names = {'group1_ligament1', 'group1_ligament2'};
        rest_lens = {0.2920, 0.5170};
    end
    
    if strcmp(model_name, 'Fiber_Length_MTU_Length_Ratio_6')
        ligament_names = {'group1_ligament1', 'group1_ligament2'};
        rest_lens = {0.0670, 0.2920};
    end
     
    %% Apply the adjustments
    ligament_params_file = [model_root_dir, '/Data/AdjustedLocalData/Ligament_Parameters.csv'];
    ligament_params_tbl = readtable(ligament_params_file, 'Delimiter', ',');
    
    for i = 1 : numel(ligament_names)
        ligament_name = ligament_names{i};
        [ligament_exists, ligament_idx] = ismember(ligament_name, ligament_params_tbl.ligament_name);

        if ~ligament_exists
            error(['Could not find ', ligament_name, ' in Ligament_Parameters.csv']);
        end

        ligament_params_tbl.rest_len(ligament_idx) = rest_lens{i};
    end

    writetable(ligament_params_tbl, [model_root_dir, '/Data/AdjustedLocalData/Ligament_Parameters.csv'], 'WriteMode', 'overwrite');

end

