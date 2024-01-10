function model = edit_mtu_properties(model, mtu_fiber_damping, rigid_tendon)
    % This function edits the properties of muscles in the model. It
    % adjusts the fiber damping and tendon compliance properties.
    %
    % NOTE: This functionality is intended for extending the toy models in 
    % the accompanying paper for simulation purposes and is not a core
    % component of the current tool's model generation capabilities. 
    % Therefore, there might be some hard-coding involved in this script.
    % Future versions may integrate this feature more comprehensively into 
    % the core scripts, if required.

    % Import OpenSim libraries
    import org.opensim.modeling.*;
    
    % Get the number of muscles in the model
    mtu_num = model.getForceSet().getMuscles().getSize();

    % Loop through all muscles and set the fiber damping and tendon
    % compliance properties
    for i = 1 : mtu_num
        % Get the muscle
        mtu = Millard2012EquilibriumMuscle.safeDownCast(model.getForceSet().getMuscles().get(i-1));
        
        % Set the fiber damping and tendon compliance properties
        mtu.set_fiber_damping(mtu_fiber_damping);
        mtu.set_ignore_tendon_compliance(rigid_tendon);
        
    end

end
