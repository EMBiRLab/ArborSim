classdef Builder
    % This class implements the Builder, which is used to construct an OpenSim model from data.
    % By specifying a boolean variable, the class can create a model wherein the complex 
    % muscle-tendon architectures are modeled either explicilty or in a conventional manner.
    %
    % For any questions, feedback, issues, please contact author:
    % xunfu@umich.edu (Xun Fu)

    properties
        model_root_dir; 
        data_dir;
        
        explicit_branching;

        osim_model;

        len_unit_in_data;
        mass_unit_in_data;
        force_unit_in_data;
        angle_unit_in_data;

        len_scale;
        mass_scale;
        force_scale;
        angle_scale;
    end
    
    methods
        function obj = Builder(model_root_dir, explicit_branching, unit_struct, gravity_vec)
            % Constructor of the class.
            %
            % Inputs:
            %   - model_root_dir: The base directory for the model (string of 1xN characters, 
            %     where N is the length of the input string)
            %   - explicit_branching: A boolean variable indicating whether the model is
            %     explicitly branched (true) or conventionally branched (false)
            %   - unit_struct: A struct containing the length and mass units used in the data 
            %     (struct with fields 'len' and 'mass')
            %   - gravity_vec: A 3x1 or 1x3 vector containing the gravity vector for the model.
            %     The order of the gravity axes in the vector is x, y, z.
            %
            % Outputs:
            %   - obj: An instance of the Builder class (Builder object)
            %
            % Example usage:
            %   - obj = Builder('../Model', true, struct('len', 'mm', 'mass', 'g'), [0, 0, 9.81])

            % Import OpenSim libraries
            import org.opensim.modeling.*;

            % Set the model root directory
            obj.model_root_dir = model_root_dir;

            % Set the explicit branching flag 
            obj.explicit_branching = explicit_branching;

            % Set the data source directory
            if obj.explicit_branching
                obj.data_dir = [model_root_dir, '/Data/AdjustedLocalData'];
            else
                obj.data_dir = [model_root_dir, '/Data/LocalData'];
            end
            
            % Initialize the OpenSim model
            obj.osim_model = Model();

            % Set the gravity vector
            % If specified, the gravity vector is set to the input value
            if nargin > 3
                obj.osim_model.set_gravity(Vec3(gravity_vec(1), gravity_vec(2), gravity_vec(3)));
            % If not specified, the default value is [0, 0, 9.81]
            else
                obj.osim_model.set_gravity(Vec3(0, 0, 9.81));
            end

            % Set the length, and mass units in the data
            % If specified, the length and mass units are set to the input values
            if nargin > 2
                obj.len_unit_in_data = unit_struct.len;
                obj.mass_unit_in_data = unit_struct.mass;
            % Otherwise, the default values are 'mm' and 'g'
            else
                obj.len_unit_in_data = 'mm';
                obj.mass_unit_in_data = 'g';
            end

            % Set the force and angle units to the default units (Newtons (N) and degrees (deg) respectively)
            % Note: The current tool configuration does not support other options for force and angle units
            obj.force_unit_in_data = 'N';
            obj.angle_unit_in_data = 'deg';
            
            % Get the scale factors for length, mass, force, and angle. This is used to scale the parameters when
            % adding bodies, joints, MTUs, ligaments, and wrapping surfaces to the model.
            obj = obj.get_scale();

        end
        
        function obj = add_bodies(obj)
            % This method adds bodies to the model.
            %
            % Inputs: None
            %
            % Outputs:
            %   - obj: An instance of the Builder class (Builder object)
            %
            % Example usage:
            %   - obj = obj.add_bodies()

            % Import OpenSim libraries
            import org.opensim.modeling.*;

            % Read the body mass properties from the Body_Mass_Properties.csv file
            body_mass_props_file = [obj.data_dir, '/Body_Mass_Properties.csv'];
            body_mass_props_tbl = readtable(body_mass_props_file, 'Delimiter', ',');

            % Check if the body names in the Body_Mass_Properties.csv file are unique
            [is_unique, nonunique_elements] = obj.find_nonunique_elements(body_mass_props_tbl.body_name(:));

            % If the body names in the Body_Mass_Properties.csv file are not unique, throw an error
            if ~is_unique
                error([strjoin(nonunique_elements, ' and '), ' occur(s) multiple times in the body_name column']);
            end

            % Read the artificial OpenSim geometries from the Artificial_OpenSim_Geometries.csv file
            osim_geom_file = [obj.data_dir, '/Artificial_OpenSim_Geometries.csv'];
            osim_geom_imp_opts = detectImportOptions(osim_geom_file);
            osim_geom_imp_opts = setvartype(osim_geom_imp_opts, 'params', 'string');
            osim_geom_imp_opts.Delimiter = ',';
            osim_geom_imp_opts.DataLines = [2, Inf];
            osim_geom_tbl = readtable(osim_geom_file, osim_geom_imp_opts);

            % Check if the body names in the Artificial_OpenSim_Geometries.csv file are unique
            [is_unique, nonunique_elements] = obj.find_nonunique_elements(osim_geom_tbl.body_name(:));
            
            % If the body names in the Artificial_OpenSim_Geometries.csv file are not unique, throw an error
            if ~is_unique
                error([strjoin(nonunique_elements, ' and '), ' occur(s) multiple times in the body_name column']);
            end

            % Iterate through the body mass properties table
            for i = 1 : size(body_mass_props_tbl, 1)
                % Get the body name, mass, center of mass, and inertia from the body mass properties table
                body_name = body_mass_props_tbl.body_name{i};
                mass = body_mass_props_tbl.mass(i) * obj.mass_scale;
                com = Vec3(body_mass_props_tbl.com_x(i) * obj.len_scale, ...
                           body_mass_props_tbl.com_y(i) * obj.len_scale, ...
                           body_mass_props_tbl.com_z(i) * obj.len_scale);
                inertia = Inertia(body_mass_props_tbl.inertia_xx(i) * obj.mass_scale * (obj.len_scale)^2, ...
                                  body_mass_props_tbl.inertia_yy(i) * obj.mass_scale * (obj.len_scale)^2, ...
                                  body_mass_props_tbl.inertia_zz(i) * obj.mass_scale * (obj.len_scale)^2, ...
                                  body_mass_props_tbl.inertia_xy(i) * obj.mass_scale * (obj.len_scale)^2, ...
                                  body_mass_props_tbl.inertia_xz(i) * obj.mass_scale * (obj.len_scale)^2, ...
                                  body_mass_props_tbl.inertia_yz(i) * obj.mass_scale * (obj.len_scale)^2);

                % Define the body to add with the body name, mass, center of mass, and inertia
                body_to_add = Body(body_name, mass, com, inertia);
                
                % Check if the body has an artificial OpenSim geometry
                row_idx_in_osim_geom_tbl = find(strcmp(osim_geom_tbl.body_name, body_name));
                
                % If the body has an artificial OpenSim geometry, associate the artificial geometry with the body
                if ~isempty(row_idx_in_osim_geom_tbl)
                    % Get the artificial OpenSim geometry type and parameters
                    osim_geom_type = osim_geom_tbl.osim_geom_type{row_idx_in_osim_geom_tbl};
                    
                    osim_geom_params = osim_geom_tbl.params{row_idx_in_osim_geom_tbl};
                    osim_geom_params = str2double(strtrim(strsplit(osim_geom_params, ',')));

                    % If the artificial OpenSim geometry is a sphere, add a sphere geometry to the body
                    if strcmp(osim_geom_type, 'Sphere')
                        sphere_radius = osim_geom_params(1) * obj.len_scale;
                        if numel(osim_geom_params) > 1
                            Warning(['For the Sphere geometry of the body ''', body_name, ''', ', ...
                                     'only one parameter, radius, is required.']);
                        end
                        body_to_add.attachGeometry(Sphere(sphere_radius));
                    
                        % If the artificial OpenSim geometry is a cylinder, add a cylinder geometry to the body
                    elseif strcmp(osim_geom_type, 'Cylinder')
                        cylinder_radius = osim_geom_params(1) * obj.len_scale;
                        cylinder_height = osim_geom_params(2) * obj.len_scale;
                        if numel(osim_geom_params) > 2
                            Warning(['For the Sphere geometry of the body ''', body_name, ''', ', ...
                                     'only two parameters, radius and height, are required.']);
                        end
                        body_to_add.attachGeometry(Cylinder(cylinder_radius, cylinder_height));

                    % Otherwise, throw an error
                    % Note the current tool configuration only supports Sphere and Cylinder artificial geometries.
                    % Other types of artifical geometries will be included in the future.
                    else
                        error('An artificial geometry needs to be one of the following: Sphere, Clyinder.')
                    end
                
                % If the body does not have an artificial OpenSim geometry, then it should have a mesh geometry
                else
                    % Get the mesh geometry file name
                    mesh_geom_file = [body_name, '.stl'];
                    
                    % Create a mesh geometry object with the mesh geometry file
                    osim_mesh_geo = Mesh(mesh_geom_file);

                    % Set the mesh geometry properties
                    osim_mesh_geo.set_scale_factors(Vec3(obj.len_scale));
                    osim_mesh_geo.upd_Appearance().set_opacity(1);
                    osim_mesh_geo.upd_Appearance().set_color(Vec3(1));
                    
                    % Associate the mesh geometry with the body
                    body_to_add.attachGeometry(osim_mesh_geo);
                
                end

                % Add the body to the model
                obj.osim_model.addBody(body_to_add);

            end

        end

        function obj = add_joints(obj)
            % This method adds joints to the model.
            %
            % Inputs: None
            %
            % Outputs:
            %   - obj: An instance of the Builder class (Builder object)
            %
            % Example usage:
            %   - obj = obj.add_joints()

            % Import OpenSim libraries
            import org.opensim.modeling.*;
            
            % Read the joint coordinate systems from the JCS.csv file
            jcs_file = [obj.data_dir, '/JCS.csv'];
            jcs_tbl = readtable(jcs_file, 'Delimiter', ',');

            % Check if the joint names in the JCS.csv file are unique
            [is_unique, nonunique_elements] = obj.find_nonunique_elements(jcs_tbl.joint_name(:));

            % If the joint names in the JCS.csv file are not unique, throw an error
            if ~is_unique
                error([strjoin(nonunique_elements, ' and '), ' occur(s) multiple times in the joint_name column']);
            end

            % Read the joint range of motion from the Joint_ROM.csv file
            joint_rom_file = [obj.data_dir, '/Joint_ROM.csv'];
            joint_rom_imp_opts = detectImportOptions(joint_rom_file);
            joint_rom_imp_opts = setvartype(joint_rom_imp_opts, 'dof', 'string');
            joint_rom_imp_opts.Delimiter = ',';
            joint_rom_tbl = readtable(joint_rom_file, joint_rom_imp_opts);

            % Check if the joint names in the Joint_ROM.csv file are unique
            [is_unique, nonunique_elements] = obj.find_nonunique_elements(joint_rom_tbl.joint_name(:));
            
            % If the joint names in the Joint_ROM.csv file are not unique, throw an error
            if ~is_unique
                error([strjoin(nonunique_elements, ' and '), ' occur(s) multiple times in the joint_name column']);
            end

            % Iterate through the joint coordinate systems table
            for i = 1 : size(jcs_tbl, 1)
                % Get the joint name, parent body, child body, location in parent, and orientation in parent
                joint_name = jcs_tbl.joint_name{i};
                par_body = jcs_tbl.par_body{i};
                chd_body = jcs_tbl.chd_body{i};
                
                % Scale the location in parent and orientation in parent using the length scale
                loc_in_par = Vec3(jcs_tbl.loc_in_par_x(i) * obj.len_scale, ...
                                  jcs_tbl.loc_in_par_y(i) * obj.len_scale, ...
                                  jcs_tbl.loc_in_par_z(i) * obj.len_scale);
                
                ori_in_par = Vec3(deg2rad(jcs_tbl.ori_in_par_x(i)), ...
                                  deg2rad(jcs_tbl.ori_in_par_y(i)), ...
                                  deg2rad(jcs_tbl.ori_in_par_z(i)));

                % Scale the location in child and orientation in child using the length scale
                loc_in_chd = Vec3(jcs_tbl.loc_in_chd_x(i) * obj.len_scale, ...
                                  jcs_tbl.loc_in_chd_y(i) * obj.len_scale, ...
                                  jcs_tbl.loc_in_chd_z(i) * obj.len_scale);
                
                ori_in_chd = Vec3(deg2rad(jcs_tbl.ori_in_chd_x(i)), ...
                                  deg2rad(jcs_tbl.ori_in_chd_y(i)), ...
                                  deg2rad(jcs_tbl.ori_in_chd_z(i)));

                % Check if the joint has a range of motion
                row_idx_in_joint_rom_tbl = find(strcmp(joint_rom_tbl.joint_name, joint_name));
                
                % If the joint has a range of motion
                if ~isempty(row_idx_in_joint_rom_tbl)
                    % Get the joint's dof, minimum and maximum values for the dof
                    dof = joint_rom_tbl.dof(row_idx_in_joint_rom_tbl);
                    min_rx = joint_rom_tbl.min_rx(row_idx_in_joint_rom_tbl);
                    max_rx = joint_rom_tbl.max_rx(row_idx_in_joint_rom_tbl);
                    min_ry = joint_rom_tbl.min_ry(row_idx_in_joint_rom_tbl);
                    max_ry = joint_rom_tbl.max_ry(row_idx_in_joint_rom_tbl);
                    min_rz = joint_rom_tbl.min_rz(row_idx_in_joint_rom_tbl);
                    max_rz = joint_rom_tbl.max_rz(row_idx_in_joint_rom_tbl);
                    min_tx = joint_rom_tbl.min_tx(row_idx_in_joint_rom_tbl);
                    max_tx = joint_rom_tbl.max_tx(row_idx_in_joint_rom_tbl);
                    min_ty = joint_rom_tbl.min_ty(row_idx_in_joint_rom_tbl);
                    max_ty = joint_rom_tbl.max_ty(row_idx_in_joint_rom_tbl);
                    min_tz = joint_rom_tbl.min_tz(row_idx_in_joint_rom_tbl);
                    max_tz = joint_rom_tbl.max_tz(row_idx_in_joint_rom_tbl);
                    
                    % Add the joint to the model with the specified dof and range of motion
                    obj = obj.add_single_joint(joint_name, par_body, chd_body, ...
                                               loc_in_par, ori_in_par, loc_in_chd, ori_in_chd, ...
                                               dof, min_rx, max_rx, min_ry, max_ry, min_rz, max_rz, ...
                                               min_tx, max_tx, min_ty, max_ty, min_tz, max_tz);
                
                % If the joint does not have a range of motion, directly add the joint to the model
                % with the default dof and range of motion
                else
                    obj = obj.add_single_joint(joint_name, par_body, chd_body, ...
                                               loc_in_par, ori_in_par, loc_in_chd, ori_in_chd);
                end
                
            end
        
        end
        
        function obj = add_single_joint(obj, joint_name, par_body, chd_body, ...
                                        loc_in_par, ori_in_par, loc_in_chd, ori_in_chd, ...
                                        dof, min_rx, max_rx, min_ry, max_ry, min_rz, max_rz, ...
                                        min_tx, max_tx, min_ty, max_ty, min_tz, max_tz)
            % This method adds a single joint to the model.
            %
            % Inputs:
            %   - joint_name: The name of the joint (string)
            %   - par_body: The parent body of the joint (string)
            %   - chd_body: The child body of the joint (string)
            %   - loc_in_par: The location of the joint in the parent body (1x3 vector)
            %   - ori_in_par: The orientation of the joint in the parent body (1x3 vector)
            %   - loc_in_chd: The location of the joint in the child body (1x3 vector)
            %   - ori_in_chd: The orientation of the joint in the child body (1x3 vector)
            %   - dof: The degrees of freedom of the joint (string)
            %   - min_rx: The minimum value of the rotation about the x-axis (double or nan)
            %   - max_rx: The maximum value of the rotation about the x-axis (double or nan)
            %   - min_ry: The minimum value of the rotation about the y-axis (double or nan)
            %   - max_ry: The maximum value of the rotation about the y-axis (double or nan)
            %   - min_rz: The minimum value of the rotation about the z-axis (double or nan)
            %   - max_rz: The maximum value of the rotation about the z-axis (double or nan)
            %   - min_tx: The minimum value of the translation along the x-axis (double or nan)
            %   - max_tx: The maximum value of the translation along the x-axis (double or nan)
            %   - min_ty: The minimum value of the translation along the y-axis (double or nan)
            %   - max_ty: The maximum value of the translation along the y-axis (double or nan)) 
            %   - min_tz: The minimum value of the translation along the z-axis (double or nan)
            %   - max_tz: The maximum value of the translation along the z-axis (double or nan)
            %
            % Outputs:
            %   - obj: An instance of the Builder class (Builder object)
            %
            % Example usage:
            %   - obj = obj.add_single_joint('joint_name', 'par_body', 'chd_body', ...
            %                                [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], ...
            %                                'rx,ry,rz,tx,ty,tz', -90, 90, -90, 90, -90, 90, ...
            %                                -10, 10, -10, 10, -10, 10)       

            % Import OpenSim libraries
            import org.opensim.modeling.*;
            
            % Get the parent body and child body
            if strcmp(par_body, 'ground')
                par_body = obj.osim_model.getGround();
            else
                par_body = obj.osim_model.getBodySet().get(par_body);
            end
            
            if strcmp(chd_body, 'ground')
                chd_body = obj.osim_model.getGround();
            else
                chd_body = obj.osim_model.getBodySet().get(chd_body);
            end

            % If the joint's dof is specified, add the joint to the model with the specified dof and range of motion
            if nargin > 8
                % If the joint's dof is '0', model the joint as a weld joint
                if strcmp(dof, '0')
                    joint_to_add = WeldJoint(joint_name, par_body, loc_in_par, ori_in_par, ...
                                             chd_body, loc_in_chd, ori_in_chd);
                    obj.osim_model.addJoint(joint_to_add);
                
                % Otherwise, model the joint as a custom joint
                else
                    % Get the joint's dof
                    dof = strtrim(strsplit(dof, ','));

                    % Check if the joint's dof is empty or contains element(s) beyond '0', 'rx', 'ry', 'rz', 'tx', 'ty', 'tz'
                    % If so, throw an error
                    if ~any(ismember(dof, {'rx', 'ry', 'rz', 'tx', 'ty', 'tz'}))
                        error(['Joint ', joint_name, '''s dof cannot be empty or contain element(s) beyond ' ...
                               '''0'', ''rx'', ''ry'', ''rz'', ''tx'', ''ty'', ''tz''']);
                    end

                    % Check if the joint's dof contains duplicate elements
                    [is_unique, nonunique_elements] = obj.find_nonunique_elements(dof);
                    
                    % If the joint's dof contains duplicate elements, throw a warning and remove the duplicate elements
                    if ~is_unique
                        warning([strjoin(nonunique_elements, ' and '), ' occur(s) multiple times in dof for joint ',  joint_name]);
                        dof = unique(dof);
                    end
                    
                    % Create a spatial transform object for the joint
                    spat_trans = SpatialTransform();

                    % Add the joint's dof to the spatial transform object
                    if ismember('rx', dof)
                        spat_trans.upd_rotation1().append_coordinates([joint_name, '_rx']);
                        spat_trans.upd_rotation1().set_axis(Vec3(1, 0, 0));
                        spat_trans.upd_rotation1().set_function(LinearFunction(1, 0));
                    end
                    
                    if ismember('ry', dof)
                        spat_trans.upd_rotation2().append_coordinates([joint_name, '_ry']);
                        spat_trans.upd_rotation2().set_axis(Vec3(0, 1, 0));
                        spat_trans.upd_rotation2().set_function(LinearFunction(1, 0));
                    end
 
                    if ismember('rz', dof)
                        spat_trans.upd_rotation3().append_coordinates([joint_name, '_rz']);
                        spat_trans.upd_rotation3().set_axis(Vec3(0, 0, 1));
                        spat_trans.upd_rotation3().set_function(LinearFunction(1, 0));
                    end

                    if ismember('tx', dof)
                        spat_trans.upd_translation1().append_coordinates([joint_name, '_tx']);
                        spat_trans.upd_translation1().set_axis(Vec3(1, 0, 0));
                        spat_trans.upd_translation1().set_function(LinearFunction(1, 0));
                    end
                    
                    if ismember('ty', dof)
                        spat_trans.upd_translation2().append_coordinates([joint_name, '_ty']);
                        spat_trans.upd_translation2().set_axis(Vec3(0, 1, 0));
                        spat_trans.upd_translation2().set_function(LinearFunction(1, 0));
                    end
                    
                    if ismember('tz', dof)
                        spat_trans.upd_translation3().append_coordinates([joint_name, '_tz']);
                        spat_trans.upd_translation3().set_axis(Vec3(0, 0, 1));
                        spat_trans.upd_translation3().set_function(LinearFunction(1, 0));
                    end

                    % Add the joint to the model
                    joint_to_add = CustomJoint(joint_name, par_body, loc_in_par, ori_in_par, ...
                                               chd_body, loc_in_chd, ori_in_chd, spat_trans);
                    obj.osim_model.addJoint(joint_to_add);

                    % Set the range of motion for the joint
                    custom_order = {'rx', 'ry', 'rz', 'tx', 'ty', 'tz'};
                    [~, order_idx] = ismember(dof, custom_order);
                    [~, sorted_idx] = sort(order_idx);
                    sorted_dof = dof(sorted_idx);

                    for i = 1 : numel(sorted_dof)
                        sel_dof = sorted_dof{i};
                        
                        if strcmp(sel_dof, 'rx')
                            min_coord = min_rx;
                            max_coord = max_rx;
                        elseif strcmp(sel_dof, 'ry')
                            min_coord = min_ry;
                            max_coord = max_ry;
                        elseif strcmp(sel_dof, 'rz')
                            min_coord = min_rz;
                            max_coord = max_rz;
                        elseif strcmp(sel_dof, 'tx')
                            min_coord = min_tx;
                            max_coord = max_tx;
                        elseif strcmp(sel_dof, 'ty')
                            min_coord = min_ty;
                            max_coord = max_ty;
                        else
                            min_coord = min_tz;
                            max_coord = max_tz;
                        end

                        if ~isnan(min_coord)
                            joint_to_add.upd_coordinates(i-1).set_range(0, min_coord * pi / 180);
                        end

                        if ~isnan(max_coord)
                            joint_to_add.upd_coordinates(i-1).set_range(1, max_coord * pi / 180);
                        end

                    end
                 
                end
            
            % If the joint's dof is not specified, add the joint to the model with the default dof and unlimited range of motion
            else
                % Create a spatial transform object for the joint
                spat_trans = SpatialTransform();

                % Add the joint's dof to the spatial transform object
                spat_trans.upd_rotation1().append_coordinates([joint_name, '_rx']);
                spat_trans.upd_rotation1().set_axis(Vec3(1, 0, 0));
                spat_trans.upd_rotation1().set_function(LinearFunction(1, 0));

                spat_trans.upd_rotation2().append_coordinates([joint_name, '_ry']);
                spat_trans.upd_rotation2().set_axis(Vec3(0, 1, 0));
                spat_trans.upd_rotation2().set_function(LinearFunction(1, 0));

                spat_trans.upd_rotation3().append_coordinates([joint_name, '_rz']);
                spat_trans.upd_rotation3().set_axis(Vec3(0, 0, 1));
                spat_trans.upd_rotation3().set_function(LinearFunction(1, 0));

                spat_trans.upd_translation1().append_coordinates([joint_name, '_tx']);
                spat_trans.upd_translation1().set_axis(Vec3(1, 0, 0));
                spat_trans.upd_translation1().set_function(LinearFunction(1, 0));

                spat_trans.upd_translation2().append_coordinates([joint_name, '_ty']);
                spat_trans.upd_translation2().set_axis(Vec3(0, 1, 0));
                spat_trans.upd_translation2().set_function(LinearFunction(1, 0));

                spat_trans.upd_translation3().append_coordinates([joint_name, '_tz']);
                spat_trans.upd_translation3().set_axis(Vec3(0, 0, 1));
                spat_trans.upd_translation3().set_function(LinearFunction(1, 0));
                
                % Add the joint to the model
                joint_to_add = CustomJoint(joint_name, par_body, loc_in_par, ori_in_par, ...
                                           chd_body, loc_in_chd, ori_in_chd, spat_trans);
                obj.osim_model.addJoint(joint_to_add);

            end

        end

        function obj = add_mtus(obj)
            % This method adds MTUs to the model.
            %
            % Inputs: None
            %
            % Outputs:
            %   - obj: An instance of the Builder class (Builder object)
            %
            % Example usage:
            %   - obj = obj.add_mtus()

            % Import OpenSim libraries
            import org.opensim.modeling.*;

            % Read the MTU parameters from the MTU_Parameters.csv file
            mtu_params_file = [obj.data_dir, '/MTU_Parameters.csv'];
            mtu_params_tbl = readtable(mtu_params_file, 'Delimiter', ',');

            % Check if the MTU names in the MTU_Parameters.csv file are unique
            [is_unique, nonunique_elements] = obj.find_nonunique_elements(mtu_params_tbl.mtu_name(:));

            % If the MTU names in the MTU_Parameters.csv file are not unique, throw an error
            if ~is_unique
                error([strjoin(nonunique_elements, ' and '), ' occur(s) multiple times in the mtu_name column']);
            end

            % Read the MTU path from the MTU_Path.csv file
            mtu_path_file = [obj.data_dir, '/MTU_Path.csv'];
            mtu_path_tbl = readtable(mtu_path_file, 'Delimiter', ',');
            
            % If the MTU names in the MTU_Parameters.csv file and MTU_Path.csv file do not match, throw an error
            if ~isequal(sort(mtu_params_tbl.mtu_name(:)), sort(unique(mtu_path_tbl.mtu_name(:))))
                error('MTUs in MTU_Parameters.csv and MTU_Path.csv do not match');
            end

            % Loop over all the MTUs
            for i = 1 : size(mtu_params_tbl, 1)
                % Get the MTU name
                mtu_name = mtu_params_tbl.mtu_name{i};
                
                % Get the MTU path
                mtu_pts_tbl = mtu_path_tbl(strcmp(mtu_path_tbl.mtu_name, mtu_name), :);
                
                % Sort the MTU path table by the point index
                sorted_mtu_pts_tbl = sortrows(mtu_pts_tbl, 'pt_idx', 'ascend');
                
                % If the point indices in the MTU path table do not form a complete, consecutive sequence 
                % from 1 to the number of points, throw an error
                if ~isequal(sorted_mtu_pts_tbl.pt_idx', 1 : size(sorted_mtu_pts_tbl, 1))
                    error(['Incomplete index sequence detected for ', mtu_name, '. ' ...
                           'The indices do not form a complete, consecutive sequence from 1 to the number of points.'])
                end
                
                % Get the MTU's maximum isometric force, optimal fiber length, tendon slack length, pennation angle,
                % Scale the maximum isometric force using the force scale
                max_iso_frc = mtu_params_tbl.max_iso_frc(i) * obj.force_scale;

                % Scale the optimal fiber length, tendon slack length, and pennation angle using the length scale
                opt_fib_len = mtu_params_tbl.opt_fib_len(i) * obj.len_scale;
                tendon_slk_len = mtu_params_tbl.tendon_slk_len(i) * obj.len_scale;

                % Convert the pennation angle from degrees to radians if the angle unit in the data is degrees
                if strcmp(obj.angle_unit_in_data, 'deg')
                    penn_ang = mtu_params_tbl.penn_ang(i) * pi / 180;
                else
                    penn_ang = mtu_params_tbl.penn_ang(i);
                end

                % Scale the default fiber length using the length scale
                def_fib_len = mtu_params_tbl.def_fib_len(i) * obj.len_scale;

                % If any of the MTU parameters are NaN, throw an error
                if any(isnan([max_iso_frc, opt_fib_len, tendon_slk_len, penn_ang, def_fib_len]))
                    error([mtu_name, '''s parameters contain NaN(s).']);
                end

                % Create an MTU object with the MTU name, maximum isometric force, optimal fiber length, tendon slack length,
                mtu_to_add = Millard2012EquilibriumMuscle(mtu_name, max_iso_frc, opt_fib_len, tendon_slk_len, penn_ang);
                mtu_to_add.set_default_fiber_length(def_fib_len);

                % Add MTU path points to the MTU object
                for j = 1 : size(sorted_mtu_pts_tbl, 1)
                    % Get the point name, attached body name, and location in body
                    pt_name = [mtu_name, '_pt', num2str(j)];
                    attached_body = sorted_mtu_pts_tbl.attached_body(j);
                    loc_in_body = Vec3(sorted_mtu_pts_tbl.loc_in_body_x(j) * obj.len_scale, ...
                                       sorted_mtu_pts_tbl.loc_in_body_y(j) * obj.len_scale, ...
                                       sorted_mtu_pts_tbl.loc_in_body_z(j) * obj.len_scale);
                    
                    % Get the attached body
                    attached_body = obj.osim_model.getBodySet().get(attached_body);

                    % Add the MTU path point to the MTU object
                    mtu_to_add.addNewPathPoint(pt_name, attached_body, loc_in_body);

                end

                % Add the MTU to the model
                obj.osim_model.addForce(mtu_to_add);
            
            end

        end

        function obj = add_ligaments(obj)
            % This method adds ligaments to the model.
            %
            % Inputs: None
            %
            % Outputs:
            %   - obj: An instance of the Builder class (Builder object)
            %
            % Example usage:
            %   - obj = obj.add_ligaments()

            % Import OpenSim libraries
            import org.opensim.modeling.*;

            % Read the ligament parameters from the Ligament_Parameters.csv file
            ligament_params_file = [obj.data_dir, '/Ligament_Parameters.csv'];
            ligament_params_tbl = readtable(ligament_params_file, 'Delimiter', ',');

            % Check if the ligament names in the Ligament_Parameters.csv file are unique
            [is_unique, nonunique_elements] = obj.find_nonunique_elements(ligament_params_tbl.ligament_name(:));
            
            % If the ligament names in the Ligament_Parameters.csv file are not unique, throw an error
            if ~is_unique
                error([strjoin(nonunique_elements, ' and '), ' occur(s) multiple times in the ligament_name column']);
            end

            % Check if the ligament names in the Ligament_Parameters.csv file and Ligament_Path.csv file match
            ligament_path_file = [obj.data_dir, '/Ligament_Path.csv'];
            ligament_path_tbl = readtable(ligament_path_file, 'Delimiter', ',');
            
            % If the ligament names in the Ligament_Parameters.csv file and Ligament_Path.csv file do not match, throw an error
            if ~isequal(sort(ligament_params_tbl.ligament_name(:)), sort(unique(ligament_path_tbl.ligament_name(:))))
                error('Ligaments in Ligament_Parameters.csv and Ligament_Path.csv do not match');
            end

            % Loop over all the ligaments
            for i = 1 : size(ligament_params_tbl, 1)
                % Get the ligament name
                ligament_name = ligament_params_tbl.ligament_name{i};
                
                % Get the ligament path
                ligament_pts_tbl = ligament_path_tbl(strcmp(ligament_path_tbl.ligament_name, ligament_name), :);
                
                % Sort the ligament path table by the point index
                sorted_ligament_pts_tbl = sortrows(ligament_pts_tbl, 'pt_idx', 'ascend');

                % If the point indices in the ligament path table do not form a complete, consecutive sequence
                % from 1 to the number of points, throw an error
                if ~isequal(sorted_ligament_pts_tbl.pt_idx', 1 : size(sorted_ligament_pts_tbl, 1))
                    error(['Incomplete index sequence detected for ', ligament_name, '. ' ...
                           'The indices do not form a complete, consecutive sequence from 1 to the number of points.'])
                end
                
                % Scale the ligament's force scale and resting length using the force scale and length scale, respectively
                frc_scale = ligament_params_tbl.frc_scale(i) * obj.force_scale;
                rest_len = ligament_params_tbl.rest_len(i) * obj.len_scale;

                % If any of the ligament parameters are NaN, throw an error
                if any(isnan([frc_scale, rest_len]))
                    error([ligament_name, '''s parameters contain NaN(s).']);
                end

                % Create a ligament object with the ligament name, force scale, and resting length
                ligament_to_add = Ligament();
                ligament_to_add.setName(ligament_name);
                ligament_to_add.set_pcsa_force(frc_scale);
                ligament_to_add.set_resting_length(rest_len);

                % Set the ligament's force-length curve to the tendon force-length curve of the MTUs in the model
                mtu_sample = Millard2012EquilibriumMuscle.safeDownCast(obj.osim_model.getForceSet().getMuscles().get(0));
                tendon_frc_len_curve = mtu_sample.getTendonForceLengthCurve();
                ligament_to_add.set_force_length_curve(tendon_frc_len_curve);

                % Add ligament path points to the ligament object
                for j = 1 : size(sorted_ligament_pts_tbl, 1)
                    % Get the point name, attached body name, and location in body
                    pt_name = [ligament_name, '_pt', num2str(j)];
                    attached_body = sorted_ligament_pts_tbl.attached_body(j);
                    loc_in_body = Vec3(sorted_ligament_pts_tbl.loc_in_body_x(j) * obj.len_scale, ...
                                       sorted_ligament_pts_tbl.loc_in_body_y(j) * obj.len_scale, ...
                                       sorted_ligament_pts_tbl.loc_in_body_z(j) * obj.len_scale);
                    
                    % Get the attached body
                    attached_body = obj.osim_model.getBodySet().get(attached_body);

                    % Add the ligament path point to the ligament object
                    ligament_to_add.updGeometryPath().appendNewPathPoint(pt_name, attached_body, loc_in_body);

                end

                % Add the ligament to the model
                obj.osim_model.addForce(ligament_to_add);
            
            end
        
        end

        function obj = add_wrap_surfs(obj)
            % This method adds wrapping surfaces to the model.
            %
            % Inputs: None
            %
            % Outputs:
            %   - obj: An instance of the Builder class (Builder object)
            %
            % Example usage:
            %   - obj = obj.add_wrap_surfs()

            % Import OpenSim libraries
            import org.opensim.modeling.*;
            
            % Read the wrapping surfaces from the Wrapping_Surfaces.csv file
            wrap_surfs_file = [obj.data_dir, '/Wrapping_Surfaces.csv'];
            wrap_surf_imp_opts = detectImportOptions(wrap_surfs_file);
            wrap_surf_imp_opts = setvartype(wrap_surf_imp_opts, 'geom_params', 'string');
            wrap_surf_imp_opts.Delimiter = ',';
            wrap_surfs = readtable(wrap_surfs_file, wrap_surf_imp_opts);

            % Check if the wrapping surface names in the Wrapping_Surfaces.csv file are unique
            [is_unique, nonunique_elements] = obj.find_nonunique_elements(wrap_surfs.wrap_surf_name(:));
            
            % If the wrapping surface names in the Wrapping_Surfaces.csv file are not unique, throw an error
            if ~is_unique
                error([strjoin(nonunique_elements, ' and '), ' occur(s) multiple times in the mtu_name column']);
            end
            
            % Loop over all the wrapping surfaces
            for i = 1 : size(wrap_surfs, 1)
                % Get the wrapping surface name, attached body, location in body, orientation in body, 
                % type, geometry parameters, and quadrant
                wrap_surf_name = wrap_surfs.wrap_surf_name{i};
                attached_body = wrap_surfs.attached_body{i};

                % Scale the location in body and orientation in body using the length scale
                loc_in_body = Vec3(wrap_surfs.loc_in_body_x(i) * obj.len_scale, ...
                                   wrap_surfs.loc_in_body_y(i) * obj.len_scale, ...
                                   wrap_surfs.loc_in_body_z(i) * obj.len_scale);

                ori_in_body = Vec3(deg2rad(wrap_surfs.ori_in_body_x(i)), ...
                                   deg2rad(wrap_surfs.ori_in_body_y(i)), ...
                                   deg2rad(wrap_surfs.ori_in_body_z(i)));

                type = wrap_surfs.type{i};
                geom_params = wrap_surfs.geom_params{i};
                geom_params = str2double(strtrim(strsplit(geom_params, ',')));
                quadrant = wrap_surfs.quadrant{i};

                % If the wrapping surface type is 'Cylinder'
                if strcmp(type, 'Cylinder')
                    % Get the radius and length of the cylinder
                    radius = geom_params(1) * obj.len_scale;
                    len = geom_params(2) * obj.len_scale;

                    % Create a wrapping cylinder object with the wrapping surface name, radius, length, and quadrant
                    wrap_to_add = WrapCylinder();
                    wrap_to_add.setName(wrap_surf_name);
                    wrap_to_add.set_radius(radius);
                    wrap_to_add.set_length(len);
                    wrap_to_add.set_quadrant(quadrant);
                    
                % If the wrapping surface type is 'Ellipsoid'
                elseif strcmp(type, 'Sphere')
                    % Get the radius of the sphere
                    radius = geom_params * obj.len_scale;

                    % Create a wrapping sphere object with the wrapping surface name, radius, and quadrant
                    wrap_to_add = WrapSphere();
                    wrap_to_add.setName(wrap_surf_name);
                    wrap_to_add.set_radius(radius);
                    wrap_to_add.set_quadrant(quadrant);
                
                % Otherwise, throw an error
                % The current version of the code only supports 'Cylinder' and 'Sphere' wrapping surfaces.
                % Other types of wrapping surfaces can be added in the future.
                else
                    error('Undefined wrapping surface type.')
                
                end

                % Set the wrapping surface's location in body and orientation in body
                wrap_to_add.set_translation(loc_in_body);
                wrap_to_add.set_xyz_body_rotation(ori_in_body);

                % Set the wrapping surface's appearance
                wrap_to_add.upd_Appearance().set_visible(true);
                wrap_to_add.upd_Appearance().set_opacity(0.5);
                wrap_to_add.upd_Appearance().set_color(Vec3(0, 1, 1));
                wrap_to_add.upd_Appearance().upd_SurfaceProperties().set_representation(3);

                % Add the wrapping surface to the model
                body = obj.osim_model.getBodySet().get(attached_body);
                body.addWrapObject(wrap_to_add);

            end
            
        end

        function obj = link_surfs_mtus_ligaments(obj)
            % This method links wrapping surfaces to MTUs and ligaments. Compared to 
            % add_wrap_surfs(), this method activates (or deactivates) the wrapping surfaces
            % the wrapping surfaces existing in the model so that mtus and ligaments can
            % actually wrap around the wrapping surfaces.
            %
            % Inputs: None
            %
            % Outputs:
            %   - obj: An instance of the Builder class (Builder object)
            %
            % Example usage:
            %   - obj = obj.link_surfs_mtus_ligaments()

            % Import OpenSim libraries
            import org.opensim.modeling.*;
            
            % Read the wrapping surfaces from the Wrapping_Surfaces.csv file
            wrap_surfs_file = [obj.data_dir, '/Wrapping_Surfaces.csv'];
            wrap_surf_imp_opts = detectImportOptions(wrap_surfs_file);
            wrap_surf_imp_opts = setvartype(wrap_surf_imp_opts, 'geom_params', 'string');
            wrap_surf_imp_opts.Delimiter = ',';
            wrap_surfs = readtable(wrap_surfs_file, wrap_surf_imp_opts);

            % Check if the wrapping surface names in the Wrapping_Surfaces.csv file are unique
            [is_unique, nonunique_elements] = obj.find_nonunique_elements(wrap_surfs.wrap_surf_name(:));
            
            % If the wrapping surface names in the Wrapping_Surfaces.csv file are not unique, throw an error
            if ~is_unique
                error([strjoin(nonunique_elements, ' and '), ' occur(s) multiple times in the wrap_surf_name column']);
            end
            
            % Read the wrapping surface-MTU/ligament pairings from the Wrapping_Surface_MTU_Ligament_Pairings.csv file
            wrap_pair_file = [obj.data_dir, '/Wrapping_Surface_MTU_Ligament_Pairings.csv'];
            wrap_pair_tbl = readtable(wrap_pair_file, 'Delimiter', ',');
            
            % If the Wrapping_Surface_MTU_Ligament_Pairings is not empty
            if ~isempty(wrap_pair_tbl.wrap_surf_name(:))
                % Check if any wrapping surface in the Wrapping_Surface_MTU_Ligament_Pairings.csv file does not
                % exist in the Wrapping_Surfaces.csv file. If so, throw an error
                if any(~ismember(wrap_pair_tbl.wrap_surf_name(:), wrap_surfs.wrap_surf_name(:)))
                    non_exist_idx = find(~ismember(wrap_pair_tbl.wrap_surf_name(:), wrap_surfs.wrap_surf_name(:)));
    
                    error([strjoin(wrap_pair_tbl.wrap_surf_name(non_exist_idx), ' and '), ...
                           'in Wrapping_Surface_MTU_Ligament_Pairings.csv do(es) not exist in Wrapping_Surfaces.csv.']);
                end
            end

            % Loop over all the wrapping surface-MTU/ligament pairings
            for i = 1 : size(wrap_pair_tbl, 1)
                % Get the wrapping surface name, paired MTU/ligament, and active status
                wrap_surf_name = wrap_pair_tbl.wrap_surf_name{i};
                paired_mtu_ligament = strtrim(strsplit(wrap_pair_tbl.paired_mtu_ligament{i}, ','));
                active = boolean(wrap_pair_tbl.active(i));

                % Get the attached body of the wrapping surface
                [~, idx] = ismember(wrap_surf_name, wrap_surfs.wrap_surf_name);
                attached_body = wrap_surfs.attached_body{idx};
                
                % Get the wrapping surface and set its active status
                attached_body = obj.osim_model.getBodySet().get(attached_body);
                wrap_object = attached_body.getWrapObject(wrap_surf_name);
                wrap_object.set_active(active);

                % Iterate over all the paired MTUs/ligaments and link them to the wrapping surface
                for j = 1 : length(paired_mtu_ligament)
                    % Get the paired MTU/ligament
                    mtu_ligament = obj.osim_model.getForceSet().get(paired_mtu_ligament{j});
                    
                    % Link the MTU/ligament to the wrapping surface
                    concrete_class = char(mtu_ligament.getConcreteClassName);
                    eval(['mtu_ligament = ', concrete_class, '.safeDownCast(mtu_ligament);']);

                    mtu_ligament.updGeometryPath().addPathWrap(wrap_object);

                end

            end
            
        end

        function obj = get_scale(obj)
            % This method gets the scale for length, mass, force, and angle.
            %
            % Inputs: None
            %
            % Outputs:
            %   - obj: An instance of the Builder class (Builder object) with 
            %     updated length, mass, force, and angle scales
            %
            % Example usage:
            %   - obj = obj.get_scale()

            % Set the length scale
            if strcmp(obj.len_unit_in_data , 'm')
                obj.len_scale = 1;
            elseif strcmp(obj.len_unit_in_data, 'cm')
                obj.len_scale = 1/100;
            elseif strcmp(obj.len_unit_in_data, 'mm')
                obj.len_scale = 1/1000;
            else
                error('The length unit in the data needs to be one of the following: m, cm, mm.')
            end

            % Set the mass scale
            if strcmp(obj.mass_unit_in_data, 'kg')
                obj.mass_scale = 1;
            elseif strcmp(obj.mass_unit_in_data, 'g')
                obj.mass_scale = 1/1000;
            else
                error('The mass unit in the data needs to be one of the following: kg, g')
            end

            % Set the force scale and angle scale.
            % Because the current tool configuration does not support other options for force and angle units,
            % the force scale and angle scale are set to 1.
            obj.force_scale = 1;
            obj.angle_scale = 1;
        
        end
        
        function [is_unique, nonunique_elements] = find_nonunique_elements(obj, cell_arr)
            % This method takes a cell array as input and returns a boolean indicating whether 
            % all elements in the array are unique, and a cell array of the non-unique elements. 
            %
            % Inputs:
            %   - cell_arr: A cell array of strings.
            % 
            % Outputs:
            %   - is_unique: A boolean indicating whether all elements in the array are unique.
            %   - nonunique_elements: A cell array of the non-unique elements.
            %
            % Example usage:
            %   [is_unique, nonunique_elements] = obj.find_nonunique_elements({'a', 'b', 'a'});

            % Find the unique elements in the array
            [unique_elements, ~, idx] = unique(cell_arr);

            % Count the number of occurences of each unique element
            occurences = histcounts(idx, 1:numel(unique_elements)+1);
            
            % Identify the non-unique elements
            nonunique_elements = unique_elements(occurences > 1);

            % Determine whether all elements in the array are unique
            is_unique = isempty(nonunique_elements);
        
        end
        
        function obj = finalize_and_export_model(obj, model_name)
            % This method finalizes the model and exports it to the Output folder.
            %
            % Inputs:
            %   - model_name: The name of the model (string)
            %
            % Outputs:
            %   - obj: An instance of the Builder class (Builder object) with the finalized model
            %
            % Example usage:
            %   - obj = obj.finalize_and_export_model('model_name')

            % Set the model name
            obj.osim_model.setName(model_name);
            
            % Set the export directory according to the explicit branching option
            if obj.explicit_branching
                export_dir = [obj.model_root_dir, '/Output/Proposed/', model_name, '.osim'];
            else
                export_dir = [obj.model_root_dir, '/Output/Conventional/', model_name, '.osim'];
            end

            % Finalize the model and export it to the Output folder
            obj.osim_model.finalizeConnections();
            obj.osim_model.print(export_dir);
            
            % Copy the geometry folder to the Output folder
            if obj.explicit_branching
                copyfile([obj.data_dir, '/Geometry'], [obj.model_root_dir, '/Output/Proposed/Geometry']);
            else
                copyfile([obj.data_dir, '/Geometry'], [obj.model_root_dir, '/Output/Conventional/Geometry']);            
            end
        end
    end
end

