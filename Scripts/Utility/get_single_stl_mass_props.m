function [mass, com, inertia] = get_single_stl_mass_props(stl_file, len_unit, mass_unit, density_val)
    % This function calculates the mass, center of mass (com), and moment of inertia
    % of a given STL file, using specified units for length and mass, and density.
    %
    % Inputs:
    %   - stl_file: A string representing the path to the STL file.
    %
    %   - len_unit: A string specifying the unit of length. This affects:
    %       1. The length unit in which the STL file was created.
    %       2. The length unit used for the input density value.
    %       3. The length unit for outputting the center of mass and inertia values.
    %
    %   - mass_unit: A string specifying the unit of mass. This affects:
    %       1. The mass unit for the input density value.
    %       2. The mass unit for outputting the mass and inertia values.
    %
    %   - density_val: A double representing the material's density.
    %
    % Outputs:
    %   - mass: A double indicating the total mass of the object.
    %
    %   - com: A 1x3 array indicating the center of mass position, with 
    %     elements for x, y, z coordinates.
    %
    %   - inertia: A 3x3 matrix representing the moment of inertia tensor, 
    %     measured about the center of mass, formatted as:
    %     [Ixx, Ixy, Ixz; Iyx, Iyy, Iyz; Izx, Izy, Izz].
    %
    % Examples:
    %   1. If len_unit is 'm' (meters) and mass_unit is 'kg' (kilograms),
    %      then the input density should be in kg/m^3. The output mass will
    %      be in kg, com in meters, and inertia in kg*m^2.
    %
    %   2. If len_unit is 'cm' (centimeters) and mass_unit is 'g' (grams),
    %      then the input density should be in g/cm^3. The output mass will
    %      be in grams, com in centimeters, and inertia in g*cm^2.
    %
    % Reference:
    % The algorithm behind this implementation can be found at 
    % http://www.magic-software.com/Documentation/PolyhedralMassProperties.pdf

    allowed_len_unit = {'m', 'cm', 'mm'};
    if ~ismember(len_unit, allowed_len_unit)
        error(['Use one of the following for the length unit: ', strjoin(allowed_len_unit(:), ', ')]);
    end

    allowed_mass_unit_list = {'kg', 'g'};
    if ~ismember(mass_unit, allowed_mass_unit_list)
        error(['Use one of the following for the density unit: ', strjoin(allowed_mass_unit_list(:), ', ')]);
    end

    %%
    TR = stlread(stl_file);
    
    vertices = TR.Points;
    faces = TR.ConnectivityList;
     
    %%
    mult = 1 ./ [6, 24, 24, 24, 60, 60, 60, 120, 120, 120];    
    intg = zeros(1, 10);
    
    %%
    for i = 1 : size(faces, 1)
        i0 = faces(i, 1);
        i1 = faces(i, 2);
        i2 = faces(i, 3);
    
        x0 = vertices(i0, 1);
        y0 = vertices(i0, 2);
        z0 = vertices(i0, 3);
    
        x1 = vertices(i1, 1);
        y1 = vertices(i1, 2);
        z1 = vertices(i1, 3);
    
        x2 = vertices(i2, 1);
        y2 = vertices(i2, 2);
        z2 = vertices(i2, 3);
    
        a1 = x1 - x0;
        b1 = y1 - y0;
        c1 = z1 - z0;
    
        a2 = x2 - x0;
        b2 = y2 - y0;
        c2 = z2 - z0;
        
        d0 = b1 * c2 - b2 * c1;
        d1 = a2 * c1 - a1 * c2;
        d2 = a1 * b2 - a2 * b1;
        
        [f1x, f2x, f3x, g0x, g1x, g2x] = subexpression([x0, x1, x2]);
        [f1y, f2y, f3y, g0y, g1y, g2y] = subexpression([y0, y1, y2]);
        [f1z, f2z, f3z, g0z, g1z, g2z] = subexpression([z0, z1, z2]);
        
        intg(1) = intg(1) + d0 * f1x;
        intg(2) = intg(2) + d0 * f2x;
        intg(3) = intg(3) + d1 * f2y;
        intg(4) = intg(4) + d2 * f2z;
        intg(5) = intg(5) + d0 * f3x;
        intg(6) = intg(6) + d1 * f3y;
        intg(7) = intg(7) + d2 * f3z;
        intg(8) = intg(8) + d0 * (y0 * g0x + y1 * g1x + y2 * g2x);
        intg(9) = intg(9) + d1 * (z0 * g0y + z1 * g1y + z2 * g2y);
        intg(10) = intg(10) + d2 * (x0 * g0z + x1 * g1z + x2 * g2z);
    
    end
    
    for i = 1 : length(mult)
        intg(i) = intg(i) * mult(i);
    
    end

    %%
    volume = intg(1);
    mass = volume * density_val;
    
    com_x = intg(2) / volume;
    com_y = intg(3) / volume;
    com_z = intg(4) / volume;

    inertia_xx = (intg(6) + intg(7) - volume * (com_y * com_y + com_z * com_z)) * density_val;
    inertia_yy = (intg(5) + intg(7) - volume * (com_z * com_z + com_x * com_x)) * density_val;
    inertia_zz = (intg(5) + intg(6) - volume * (com_x * com_x + com_y * com_y)) * density_val;
    
    inertia_xy = - (intg(8) - volume * com_x * com_y) * density_val;
    inertia_yz = - (intg(9) - volume * com_y * com_z) * density_val;
    inertia_xz = - (intg(10) - volume * com_z * com_x) * density_val;
    
    com = [com_x, com_y, com_z];
    
    inertia = [inertia_xx, inertia_xy, inertia_xz;
               inertia_xy, inertia_yy, inertia_yz;
               inertia_xz, inertia_yz, inertia_zz];
end

%%
function [f1, f2, f3, g0, g1, g2] = subexpression(x)
    w0 = x(1);
    w1 = x(2);
    w2 = x(3);

    temp0 = w0 + w1;
    f1 = temp0 + w2;
    temp1 = w0 * w0;
    temp2 = temp1 + w1 * temp0;
    f2 = temp2 + w2 * f1;
    f3 = w0 * temp1 + w1 * temp2 + w2 * f2;
    g0 = f2 + w0 * (f1 + w0);
    g1 = f2 + w1 * (f1 + w1);
    g2 = f2 + w2 * (f1 + w2);

end

