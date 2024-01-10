function [model, coord_reporter] = add_coordreporter(model, reporter_dt)
    % This function adds a coordinate reporter to the given model.
    %
    % NOTE: This functionality is intended for extending the toy models in 
    % the accompanying paper for simulation purposes and is not a core
    % component of the current tool's model generation capabilities. 
    % Therefore, there might be some hard-coding involved in this script.
    % Future versions may integrate this feature more comprehensively into 
    % the core scripts, if required.

    % Import OpenSim libraries
    import org.opensim.modeling.*;
    
    % Create a coordinate reporter
    coord_reporter = TableReporter();
    coord_reporter.set_report_time_interval(reporter_dt);
    
    % Add coordinates to the reporter
    joint_num = model.getJointSet().getSize();
    for i = 1 : joint_num
        joint_name = model.getJointSet().get(i-1).getName();
        
        joint = model.updJointSet().get(joint_name);
        coord_name = char(joint.getCoordinate);

        coord_reporter.addToReport(joint.getCoordinate().getOutput('value'), coord_name);

    end
    
    % Add the reporter to the model
    model.addComponent(coord_reporter);
    
end
