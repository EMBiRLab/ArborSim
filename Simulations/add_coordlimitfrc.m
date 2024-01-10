function model = add_coordlimitfrc(model, coordlimitdfrc_transition, coordlimitfrc_stiffness_ub, coordlimitfrc_damping)
    % This function adds a 'CoordinateLimitForce' to each rotational joint
    % in the given model.
    %
    % NOTE: This functionality is intended for extending the toy models in 
    % the accompanying paper for simulation purposes and is not a core
    % component of the current tool's model generation capabilities. 
    % Therefore, there might be some hard-coding involved in this script.
    % Future versions may integrate this feature more comprehensively into 
    % the core scripts, if required.
    
    % Import OpenSim libraries
    import org.opensim.modeling.*;
    
    % Add CoordinateLimitForce to each rotational joint
    for i = 1 : 8
        % Get joint name
        if i == 1 
            joint_name = 'ground_link1';
        else
            joint_name = ['link', num2str(i-1), '_link', num2str(i)];
        end
        
        % Get joint coordinate name and range
        joint = model.updJointSet().get(joint_name);
        coord_name = toString(joint.getCoordinate);
        coord_ub = joint.get_coordinates(0).get_range(1);

        % Create CoordinateLimitForce
        coordlimitfrc_coord_ub = rad2deg(coord_ub) - coordlimitdfrc_transition;
        coordlimitfrc_coord_lb = -coordlimitfrc_coord_ub;
        coordlimitfrc_stiffness_lb = -coordlimitfrc_stiffness_ub;

        % Configure CoordinateLimitForce
        coorlimitfrc_to_add = CoordinateLimitForce(coord_name, ...
                                                   coordlimitfrc_coord_ub, coordlimitfrc_stiffness_ub, ...
                                                   coordlimitfrc_coord_lb, coordlimitfrc_stiffness_lb, ...
                                                   coordlimitfrc_damping, coordlimitdfrc_transition);
        
        % Add CoordinateLimitForce to the model
        joint.addComponent(coorlimitfrc_to_add);
    end

end
