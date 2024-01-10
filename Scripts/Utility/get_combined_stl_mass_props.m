function [combined_mass, combined_com, combined_inertia] = get_combined_stl_mass_props(stl_files, len_unit, mass_unit, density_vals)
    % This MATLAB function calculates the combined mass, center of mass (com), 
    % and moment of inertia for a composite object represented by multiple STL files.
    % It treats all the input STL files as components of a single entity and computes
    % the overall physical properties of this combined structure. The function outputs
    % the combined mass, center of mass, and moment of inertia based on specified 
    % units of length and mass, and a set of density values for each component.
    %
    % Inputs:
    %   - stl_files: An N by 1 (or 1 by N) cell array of STL files, where each 
    %     element is the path to an STL file. These files collectively represent 
    %     the 3D geometry of the composite object. The function treats these 
    %     files as parts of a single object for the calculations.
    %
    %   - len_unit: A string representing the unit of length used throughout 
    %     this function. This includes:
    %       1. The unit of length in which the STL files were originally created.
    %       2. The unit of length used for the input densities (density_val).
    %       3. The unit of length for the output com and inertia values.
    %
    %   - mass_unit: A string representing the unit of mass used throughout 
    %     this function. This includes:
    %       1. The unit of mass for the input densities (density_val).
    %       2. The unit of mass for the output mass and inertia values.
    %
    %   - density_val: An N by 1 (or 1 by N) array, with each element representing 
    %     the density value of the corresponding part of the object as defined 
    %     by the STL files.
    %
    % Ouputs:
    %   - combined_mass: A double representing the total mass of the combined object.
    %
    %   - combined_com: A 1 by 3 array representing the center of mass of the 
    %     combined object. The elements represent the x, y, and z coordinates.
    %
    %   - combined_inertia: A 3 by 3 matrix representing the moment of inertia 
    %     tensor, measured at the combined_com, of the combined object in the format:
    %     [Ixx, Ixy, Ixz; Iyx, Iyy, Iyz; Izx, Izy, Izz].
    %
    % Examples:
    %   1. If len_unit is 'm' (meters) and mass_unit is 'kg' (kilograms),
    %      then the input density should be in kg/m^3. The output mass, com, 
    %      and inertia will be calculated for the composite object as a whole,
    %      with mass in kg, com in meters, and inertia represented by a 3x3 matrix.
    %
    %   2. If len_unit is 'cm' (centimeters) and mass_unit is 'g' (grams),
    %      then the input density should be in g/cm^3. The output mass, com, 
    %      and inertia will be for the composite object as a whole, with mass
    %      in grams, com in centimeters, and inertia represented by a 3x3 matrix.

    %%
    if length(stl_files) ~= length(density_vals)
        error('The number of provided STL files and the density values do not match.')
    end
    
    %%
    mass_arr = zeros(length(stl_files), 1);
    com_arr = zeros(length(stl_files), 3);
    inertia_arr = cell(length(stl_files), 1);

    for i = 1 : length(stl_files)
        [mass, com, inertia] = get_single_stl_mass_props(stl_files{i}, len_unit, mass_unit, density_vals(i));

        mass_arr(i) = mass;
        com_arr(i, :) = com;
        inertia_arr{i} = inertia;
    
    end
    
    %%
    combined_mass = sum(mass_arr);
    combined_com = mass_arr' * com_arr / combined_mass;
    
    %%
    par_inertia_arr = cell(size(inertia_arr));
    for i = 1 : length(stl_files)
        d = combined_com - com_arr(i, :);

        par_inertia_arr{i} = inertia_arr{i} + mass_arr(i) * (d*d' * eye(3) - kron(d, d'));
    
    end
    
    %%
    combined_inertia = zeros(3);
    for i = 1 : length(stl_files)
        combined_inertia = combined_inertia + par_inertia_arr{i};

    end

end

