% This script is used to build the models used for 
% conducting the study in the accompanying paper. The 
% data for model construction is stored in the folder 
% ../Models/ArobrSim_Paper_Toy_Models.

%% Set up the environment
clear; close all;
addpath(genpath('Core'));

%% Set up the model construction parameters
% The master directory where the data for model construction is stored
model_master_dir = '../Models/ArborSim_Paper_Toy_Models/';
categories = {'Tendon_Branch_Count_Category', 'Joint_Count_Between_Insertions_Category', 'Fiber_Length_MTU_Length_Ratio_Category'};
model_name_prefix_per_category = {'Tendon_Branch_Count', 'Joint_Count_Between_Insertions', 'Fiber_Length_MTU_Length_Ratio'};
model_name_suffix_arr_per_category = {1:6, 1:6, 1:6};

% Initialize the individual model root directories and names
individual_model_dirs = {};
individual_model_names = {};

% Generate the individual model root directories and names
for i = 1 : numel(categories)
    model_name_suffix_arr = model_name_suffix_arr_per_category{i};
    for j = 1 : numel(model_name_suffix_arr)
        model_name = [model_name_prefix_per_category{i}, '_', num2str(model_name_suffix_arr(j))];
        
        individual_model_dirs{end+1, 1} = fullfile(model_master_dir, categories{i}, model_name);
        individual_model_names{end+1, 1} = model_name;
    end
    
end

% Define the gravity vector of the model
model_gravity_vec = zeros(1,3);

% Define the length and mass units used in the data. Default units for angle 
% and force are degrees ('deg') and Newtons ('N'), respectively.
units_in_data.len = 'm';
units_in_data.mass = 'kg';

% Set the flag for explicit branching. This cell array contains both
% true and false, meaning that the models will be built with and without
% explicit branching, respectively.
explicit_branching = {false, true};

% Parameters for the floating body, used only if explicit_branching is true.
float_body_props.radius = 0.005;
float_body_props.mass = 0.0005;
float_body_props.geom_radius = 0.005;

% Define the sliding joint degree of freedom, used only if 
% explicit_branching is true.
sliding_joint_dof = {'ty'};

% Set preset values for tendon slack length and ligament rest length, used 
% only if explicit_branching is true.
preset_tendon_slk_len = 0.0001;
preset_ligament_rest_len = 1;

% Loop through the individual model root directories and names
for i = 1 : numel(individual_model_dirs)
    % For each individual model root directory, loop through the explicit branching flag
    for j = 1 : numel(explicit_branching)

        %% Adjust the local data if explicit branching flag is set to true
        if explicit_branching{j}
            % Create a brancher object
            brancher_obj = Brancher(individual_model_dirs{i}, units_in_data, float_body_props);
        
            % Build direct graphs for MTUs in all branch groups
            brancher_obj = brancher_obj.mtu_path_tbl_to_struct();
            brancher_obj = brancher_obj.create_digraphs_for_all_branch_groups();
        
            % Identify the special nodes in the direct graphs
            brancher_obj = brancher_obj.identify_special_nodes();
        
            % Add the floating bodies and the respective sliding joints to the model
            brancher_obj = brancher_obj.augment_body_N_joints(sliding_joint_dof);

            % Adjust the MTU paths to account for the proposed branching muscle-tendon
            % architecture modeling method
            brancher_obj = brancher_obj.adjust_mtu_path_tbl();
            
            % Adjust the MTU parameters based on the adjusted MTU paths
            brancher_obj = brancher_obj.adjust_mtu_params_tbl(preset_tendon_slk_len, preset_ligament_rest_len);
            
            % Add the floating body visualization geometries to the model
            brancher_obj = brancher_obj.adjust_osim_geom_tbl();
            
            % Adjust the wrapping surface pairs
            brancher_obj = brancher_obj.adjust_wrap_surf_pair_tbl();
        
            % Export the adjusted local data to the respective directory
            brancher_obj.export_and_mirror_files();
        end
        
        %% Manually adjust the MTU and Ligament parameters. 
        if explicit_branching{j}
            update_arborsim_paper_toy_model_params(individual_model_dirs{i}, individual_model_names{i});
        end
        
        %% Assembly the model components from either the "local" or "adjusted local" data
        % Create a builder object
        builder_obj = Builder(individual_model_dirs{i}, explicit_branching{j}, units_in_data, model_gravity_vec);
        
        % Add the bodies, joints, MTUs, and ligaments to the model
        builder_obj = builder_obj.add_bodies();
        builder_obj = builder_obj.add_joints();
        builder_obj = builder_obj.add_mtus();
        
        if explicit_branching{j}
            builder_obj = builder_obj.add_ligaments();
        end
        
        % Add the wrap surfaces and link the surfaces to the MTUs and ligaments
        builder_obj = builder_obj.add_wrap_surfs();
        builder_obj = builder_obj.link_surfs_mtus_ligaments();
        
        % Export the model
        builder_obj.finalize_and_export_model(individual_model_names{i});
    end

end

