## Data Format Specifications

This document provides guidelines for supplying geometry files and filling out the CSV files within the `Data` directory.

### Outline
   - [CSV Data Format](#csv-data-format)
      - GlobalData
         - [`Artificial_OpenSim_Geometries.csv`](#artificial_opensim_geometriescsv)
         - [`Body_Frame_Landmarks.csv`](#body_frame_landmarkscsv)
         - [`Body_Mass_Properties.csv`](#body_mass_propertiescsv)
         - [`Branch_Groups.csv`](#branch_groupscsv)
         - [`JCS_Landmarks.csv`](#jcs_landmarkscsv)
         - [`Joint_ROM.csv`](#joint_romcsv)
         - [`MTU_Parameters.csv`](#mtu_parameterscsv)
         - [`MTU_Path_Landmarks.csv`](#mtu_path_landmarkscsv)
         - [`Wrapping_Surface_Landmarks`](#wrapping_surface_landmarkscsv)
         - [`Wrapping_Surface_MTU_Ligament_Pairings.csv`](#wrapping_surface_mtu_ligament_pairingscsv)
      - LocalData
         - [`Artificial_OpenSim_Geometries.csv`](#artificial_opensim_geometriescsv)
         - [`Body_Mass_Properties.csv`](#body_mass_propertiescsv)
         - [`Branch_Groups.csv`](#branch_groupscsv)
         - [`JCS.csv`](#jcscsv)
         - [`Joint_ROM.csv`](#joint_romcsv)
         - [`MTU_Parameters.csv`](#mtu_parameterscsv)
         - [`MTU_Path.csv`](#mtu_pathcsv)
         - [`Wrapping_Surface_MTU_Ligament_Pairings.csv`](#wrapping_surface_mtu_ligament_pairingscsv)
         - [`Wrapping_Surfaces.csv`](#wrapping_surfacescsv)
      - AdjustedLocalData
         - [`Artificial_OpenSim_Geometries.csv`](#artificial_opensim_geometriescsv)
         - [`Body_Mass_Properties.csv`](#body_mass_propertiescsv)
         - [`Branch_Groups.csv`](#branch_groupscsv)
         - [`JCS.csv`](#jcscsv)
         - [`Joint_ROM.csv`](#joint_romcsv)
         - [`Ligament_Parameters.csv`](#ligament_parameterscsv)
         - [`Ligament_Path.csv`](#ligament_pathcsv)
         - [`MTU_Parameters.csv`](#mtu_parameterscsv)
         - [`MTU_Path.csv`](#mtu_pathcsv)
         - [`Wrapping_Surface_MTU_Ligament_Pairings.csv`](#wrapping_surface_mtu_ligament_pairingscsv)
         - [`Wrapping_Surfaces.csv`](#wrapping_surfacescsv)
   - [Geometry File Format](#geometry-file-format)
   - [Important Notes](#important-notes)
---

### CSV Data Format

#### `Artificial_OpenSim_Geometries.csv`
   - **Format Specs**:
      | Column Name       | Description                                                                | Data Type                | Unit Type |
      |-------------------|----------------------------------------------------------------------------|--------------------------|-----------|
      | `body_name`       | The body's name.                                                           | Text                     | None      |
      | `osim_geom_type`  | The type of OpenSim artificial geometry used to visualize the body.        | Text                     | None      |
      | `params`          | The parameters for configuring the artificial geometry.                    | Numeric (floating-point) | Length    |
      
   - **Example Entry**:
      ```command-line
      body_name,osim_geom_type,params
      bone1,Sphere,1.0
      bone2,Cylinder,"1.0,2.0"
      ```

   - **Additional Notes**:
      - The current ArborSim setup supports two types of OpenSim artificial geometry: `Sphere` and `Cylinder`.
      - For `Sphere`, users need to supply a single numeric value, indicating the radius. For `Cylinder`, users need to supply two numeric values, with the first indicating the radius and the second indicating the **HALF** length of the cylinder along the axial direction. 


#### `Body_Frame_Landmarks.csv`
   - **Format Specs**:
      | Column Name       | Description                                                                                         | Data Type                | Unit Type |
      |-------------------|-----------------------------------------------------------------------------------------------------|--------------------------|-----------|
      | `body_name`       | The body's name.                                                                                    | Text                     | None      |
      | `org_pt_x`        | X-coordinate of the landmark point, marking the origin of the body's expected coordinate system.    | Numeric (floating-point) | Length    |
      | `org_pt_y`        | Y-coordinate of the landmark point, marking the origin of the body's expected coordinate system.    | Numeric (floating-point) | Length    |
      | `org_pt_z`        | Z-coordinate of the landmark point, marking the origin of the body's expected coordinate system.    | Numeric (floating-point) | Length    |
      | `pos_xaxis_pt_x`  | X-coordinate of the landmark point, indicating the positive X-axis of the body's coordinate system. | Numeric (floating-point) | Length    |
      | `pos_xaxis_pt_y`  | Y-coordinate of the landmark point, indicating the positive X-axis of the body's coordinate system. | Numeric (floating-point) | Length    |
      | `pos_xaxis_pt_z`  | Z-coordinate of the landmark point, indicating the positive X-axis of the body's coordinate system. | Numeric (floating-point) | Length    |
      | `pos_yaxis_pt_x`  | X-coordinate of the landmark point, indicating the positive Y-axis of the body's coordinate system. | Numeric (floating-point) | Length    |
      | `pos_yaxis_pt_y`  | Y-coordinate of the landmark point, indicating the positive Y-axis of the body's coordinate system. | Numeric (floating-point) | Length    |
      | `pos_yaxis_pt_z`  | Z-coordinate of the landmark point, indicating the positive Y-axis of the body's coordinate system. | Numeric (floating-point) | Length    |

   - **Example Entry**:
      ```command-line
      body_name,org_pt_x,org_pt_y,org_pt_z,pos_xaxis_pt_x,pos_xaxis_pt_y,pos_xaxis_pt_z,pos_yaxis_pt_x,pos_yaxis_pt_y,pos_yaxis_pt_z
      ground,20.8287,-51.0093,72.5108,20.7983,-49.9792,74.5576,21.4481,-50.9499,72.4841
      bone1,1.1,-1.2,3.0,2.4,4.0,5.33,22.45,58.2,8.0
      ```

   - **Additional Notes**:
      - The coordinates in this CSV file are with respect to a single "global" coordinate system as this CSV file is part of the "global" data.
      - Please refer to the accompanying paper and the references within for specific instructions on how to obtain the landmarks for each body frame landmark.

#### `Body_Mass_Properties.csv`
   - **Format Specs**:
      | Column Name       | Description                                                | Data Type                | Unit Type     |
      |-------------------|------------------------------------------------------------|--------------------------|---------------|
      | `body_name`       | The body's name.                                           | Text                     | None          |
      | `mass`            | The mass of the body.                                      | Numeric (floating-point) | Mass          |
      | `com_x`           | X-coordinate of the center of mass of the body.            | Numeric (floating-point) | Length        |
      | `com_y`           | Y-coordinate of the center of mass of the body.            | Numeric (floating-point) | Length        |
      | `com_z`           | Z-coordinate of the center of mass of the body.            | Numeric (floating-point) | Length        |
      | `inertia_xx`      | Moment of inertia about the X-axis of the body.            | Numeric (floating-point) | Mass*Length^2 |
      | `inertia_yy`      | Moment of inertia about the Y-axis of the body.            | Numeric (floating-point) | Mass*Length^2 |
      | `inertia_zz`      | Moment of inertia about the Z-axis of the body.            | Numeric (floating-point) | Mass*Length^2 |
      | `inertia_xy`      | Product of inertia about the X and Y axes of the body.     | Numeric (floating-point) | Mass*Length^2 |
      | `inertia_xz`      | Product of inertia about the X and Z axes of the body.     | Numeric (floating-point) | Mass*Length^2 |
      | `inertia_yz`      | Product of inertia about the Y and Z axes of the body.     | Numeric (floating-point) | Mass*Length^2 |

   - **Example Entry**:
      ```command-line
      body_name,mass,com_x,com_y,com_z,inertia_xx,inertia_yy,inertia_zz,inertia_xy,inertia_xz,inertia_yz
      bone1,5.0,0.0,0.0,0.0,1.0,2.0,3.0,0.1,0.2,0.3
      bone2,10.0,1.0,2.0,3.0,4.0,5.0,6.0,0.4,0.5,0.6
      ```

   - **Additional Notes**:
      - If this CSV file is provided as part of the "global" data, the center of mass and inertia provided in the CSV file should be measured with respect to the "global" coordinate system. If this CSV file is provided as part of the "local" data, then the center of mass and inertia provided in the CSV file should be measured relative to each body's "local" coordinate system.
      - The tool incorporates scripts to acquire mass properties, including mass, center of mass, and inertia, from a single STL file or mutliple STL file for a body. Users can calcuate the "global" mass properties accordingly. For more detailed information, please refer [`get_combined_stl_mass_props`](../Scripts/Utility/get_combined_stl_mass_props.m) or [`get_single_stl_mass_props`](../Scripts/Utility/get_combined_stl_mass_props.m).
      - The 'global' mass properties can be transformed into a "local" format based on the transformation matrices acquired from `Body_Frame_Landmarks.csv` by calling function [`apply_trans_mat_to_bmp`](../Scripts/Core/Transformer.m#L602).

#### `Branch_Groups.csv`
   - **Format Specs**:
      | Column Name       | Description                                                | Data Type     | Unit Type     |
      |-------------------|------------------------------------------------------------|---------------|---------------|
      | `group_name`      | The name of the group.                                     | Text          | None          |
      | `mtu_in_group`    | The names of the MTUs that comprise the group.             | Text          | None          |

   - **Example Entry**:
      ```command-line
      group_name,mtu_in_group
      group1,"mtu1,mtu2,mtu10"
      group2,"mtu3,mtu4"
      ```
   - **Additional Notes**:
      Please refer to the accompanying paper for the definition of branching groups.

#### `JCS_Landmarks.csv`
   - **Format Specs**:
      | Column Name                   | Description                                                                                    | Data Type               | Unit Type     |
      |-------------------------------|------------------------------------------------------------------------------------------------|-------------------------|---------------|
      | `joint_name`                  | The name of the joint.                                                                         | Text                    | None          |
      | `par_body`                    | The name of the parent body of the joint.                                                      | Numeric (floating-point)| None          |
      | `chd_body`                    | The name of the child body of the joint.                                                       | Numeric (floating-point)| None          |
      | `par_joint_frm_org_pt_x`      | X-coordinate of the landmark point that represents the origin of the parent joint frame.       | Numeric (floating-point)| Length        |
      | `par_joint_frm_org_pt_y`      | Y-coordinate of the landmark point that represents the origin of the parent joint frame.       | Numeric (floating-point)| Length        |
      | `par_joint_frm_org_pt_z`      | Z-coordinate of the landmark point, representing the origin of the parent joint frame.         | Numeric (floating-point)| Length        |
      | `par_joint_frm_pos_xaxis_pt_x`| X-coordinate of the landmark point that denotes the positive X-axis of the parent joint frame. | Numeric (floating-point)| Length        |
      | `par_joint_frm_pos_xaxis_pt_y`| X-coordinate of the landmark point that denotes the positive X-axis of the parent joint frame. | Numeric (floating-point)| Length        |
      | `par_joint_frm_pos_xaxis_pt_z`| Z-coordinate of the landmark point that denotes the positive X-axis of the parent joint frame. | Numeric (floating-point)| Length        |
      | `par_joint_frm_pos_yaxis_pt_x`| X-coordinate of the landmark point that denotes the positive Y-axis of the parent joint frame. | Numeric (floating-point)| Length        |
      | `par_joint_frm_pos_yaxis_pt_y`| Y-coordinate of the landmark point that denotes the positive Y-axis of the parent joint frame. | Numeric (floating-point)| Length        |
      | `par_joint_frm_pos_yaxis_pt_z`| Z-coordinate of the landmark point that denotes the positive Y-axis of the parent joint frame. | Numeric (floating-point)| Length        |
      | `chd_joint_frm_org_pt_x`      | X-coordinate of the landmark point that represents the origin of the child joint frame.        | Numeric (floating-point)| Length        |
      | `chd_joint_frm_org_pt_y`      | Y-coordinate of the landmark point that represents the origin of the child joint frame.        | Numeric (floating-point)| Length        |
      | `chd_joint_frm_org_pt_z`      | Z-coordinate of the landmark point that represents the origin of the child joint frame.        | Numeric (floating-point)| Length        |
      | `chd_joint_frm_pos_xaxis_pt_x`| X-coordinate of the landmark point that denotes the positive X-axis of the child joint frame.  | Numeric (floating-point)| Length        |
      | `chd_joint_frm_pos_xaxis_pt_y`| Y-coordinate of the landmark point that denotes the positive X-axis of the child joint frame.  | Numeric (floating-point)| Length        |
      | `chd_joint_frm_pos_xaxis_pt_z`| Z-coordinate of the landmark point that denotes the positive X-axis of the child joint frame.  | Numeric (floating-point)| Length        |
      | `chd_joint_frm_pos_yaxis_pt_x`| X-coordinate of the landmark point that denotes the positive Y-axis of the child joint frame.  | Numeric (floating-point)| Length        |
      | `chd_joint_frm_pos_yaxis_pt_y`| Y-coordinate of the landmark point that denotes the positive Y-axis of the child joint frame.  | Numeric (floating-point)| Length        |
      | `chd_joint_frm_pos_yaxis_pt_z`| Z-coordinate of the landmark point that denotes the positive Y-axis of the child joint frame.  | Numeric (floating-point)| Length        |
      
   - **Example Entry**:
      ```command-line
      joint_name,par_body,chd_body,par_joint_frm_org_pt_x,par_joint_frm_org_pt_y,par_joint_frm_org_pt_z,par_joint_frm_pos_xaxis_pt_x,par_joint_frm_pos_xaxis_pt_y,par_joint_frm_pos_xaxis_pt_z,par_joint_frm_pos_yaxis_pt_x,par_joint_frm_pos_yaxis_pt_y,par_joint_frm_pos_yaxis_pt_z,chd_joint_frm_org_pt_x,chd_joint_frm_org_pt_y,chd_joint_frm_org_pt_z,chd_joint_frm_pos_xaxis_pt_x,chd_joint_frm_pos_xaxis_pt_y,chd_joint_frm_pos_xaxis_pt_z,chd_joint_frm_pos_yaxis_pt_x,chd_joint_frm_pos_yaxis_pt_y,chd_joint_frm_pos_yaxis_pt_z
      ground_bone1,ground,bone1,20.8287,-51.0093,72.5108,20.7983,-49.9792,74.5576,21.4481,-50.9499,72.4841,20.8287,-51.0093,72.5108,20.7983,-49.9792,74.5576,21.4481,-50.9499,72.4841
      ```
   - **Additional Notes**:
      - The coordinates in this CSV file are with respect to a single "global" coordinate system as this CSV file is part of the "global" data.
      - Please refer to the accompanying paper and the references within for specific instructions on how to obtain the landmarks for the joints.

#### `Joint_ROM.csv`
   - **Format Specs**:
      | Column Name   | Description                                                                   | Data Type               | Unit Type     |
      |---------------|-------------------------------------------------------------------------------|-------------------------|---------------|
      | `joint_name`  | The name of the joint.                                                        | Text                    | None          |
      | `dof`         | The degrees of freedom of the joint.                                          | Text OR  Integer 0      | None          |
      | `min_rx`      | The minimum value of the rotational degree of freedom around the X-axis.      | Numeric (floating-point)| Angle         |
      | `max_rx`      | The maximum value of the rotational degree of freedom around the X-axis.      | Numeric (floating-point)| Angle         |
      | `min_ry`      | The minimum value of the rotational degree of freedom around the Y-axis.      | Numeric (floating-point)| Angle         |
      | `max_ry`      | The maximum value of the rotational degree of freedom around the Y-axis.      | Numeric (floating-point)| Angle         |
      | `min_rz`      | The minimum value of the rotational degree of freedom around the Z-axis.      | Numeric (floating-point)| Angle         |
      | `max_rz`      | The maximum value of the rotational degree of freedom around the Z-axis.      | Numeric (floating-point)| Angle         |
      | `min_tx`      | The minimum value of the translational degree of freedom along the X-axis.    | Numeric (floating-point)| Length        |
      | `max_tx`      | The maximum value of the translational degree of freedom along the X-axis.    | Numeric (floating-point)| Length        |
      | `min_ty`      | The minimum value of the translational degree of freedom along the Y-axis.    | Numeric (floating-point)| Length        |
      | `max_ty`      | The maximum value of the translational degree of freedom along the Y-axis.    | Numeric (floating-point)| Length        |
      | `min_tz`      | The minimum value of the translational degree of freedom along the Z-axis.    | Numeric (floating-point)| Length        |
      | `max_tz`      | The maximum value of the translational degree of freedom along the Z-axis.    | Numeric (floating-point)| Length        |
   
   - **Example Entry**:
      ```command-line
      joint_name,dof,min_rx,max_rx,min_ry,max_ry,min_rz,max_rz,min_tx,max_tx,min_ty,max_ty,min_tz,max_tz
      ground_bone1,"rx,ry,rz",-10,10,-90,90,-20,20,,,,,,
      bone1_bone2,"rx,tx,",-90,90,,,,,-1,1,,,,
      ```
   - **Additional Notes**:
      - The provided `dof` should be either a non-empty subset of `[rx, ry, rz, tx, ty, tz]`, or `0`, where `0` indicates that this joint is modeled as a WeldJoint in OpenSim.
      - Users do not need to supply the `dof` data for all the joints in the model. For the joints that are not included in this CSV file, they are by default modeled as a 6DOF CustomJoint in OpenSim, with no limited range of motion. If one supplies the `dof` data for a joint, the range of motion data (`min_rx`, `max_rx`, etc.) are not strictly required. One can simply leave the range of motion data empty. However, if one supplies the range of motion data, the user has the responsibility to ensure that the `dof` corresponding to the range of motion is provided in the `dof` data.

#### `JCS.csv`
   - **Format Specs**:
      | Column Name   | Description                                                                                            | Data Type               | Unit Type    |
      |---------------|--------------------------------------------------------------------------------------------------------|-------------------------|--------------|
      | `joint_name`  | The name of the joint.                                                                                 | Text                    | None         |
      | `par_body`    | The name of the parent body of the joint.                                                              | Text                    | None         |
      | `chd_body`    | The name of the child body of the joint.                                                               | Numeric (floating-point)| Length       |
      | `loc_in_par_x`| X-coordinate of the joint in the parent body frame.                                                    | Numeric (floating-point)| Length       |
      | `loc_in_par_y`| Y-coordinate of the joint in the parent body frame.                                                    | Numeric (floating-point)| Length       |
      | `loc_in_par_z`| Z-coordinate of the joint in the parent body frame.                                                    | Numeric (floating-point)| Angle        |
      | `ori_in_par_x`| X-axis rotation in the XYZ body-fixed Euler angles of the joint frame orientation in the parent frame. | Numeric (floating-point)| Angle        |
      | `ori_in_par_y`| Y-axis rotation in the XYZ body-fixed Euler angles of the joint frame orientation in the parent frame. | Numeric (floating-point)| Angle        |
      | `ori_in_par_z`| Z-axis rotation in the XYZ body-fixed Euler angles of the joint frame orientation in the parent frame. | Numeric (floating-point)| Length       |
      | `loc_in_chd_x`| X-coordinate of the joint in the child body frame.                                                     | Numeric (floating-point)| Length       |
      | `loc_in_chd_y`| Y-coordinate of the joint in the child body frame.                                                     | Numeric (floating-point)| Length       |
      | `loc_in_chd_z`| Z-coordinate of the joint in the child body frame.                                                     | Numeric (floating-point)| Length       |
      | `ori_in_chd_x`| X-axis rotation in the XYZ body-fixed Euler angles of the joint frame orientation in the child frame.  | Numeric (floating-point)| Angle        |
      | `ori_in_chd_y`| Y-axis rotation in the XYZ body-fixed Euler angles of the joint frame orientation in the child frame.  | Numeric (floating-point)| Angle        |
      | `ori_in_chd_z`| Z-axis rotation in the XYZ body-fixed Euler angles of the joint frame orientation in the child frame.  | Numeric (floating-point)| Angle        |

   - **Example Entry**:
      ```command-line
      joint_name,par_body,chd_body,loc_in_par_x,loc_in_par_y,loc_in_par_z,ori_in_par_x,ori_in_par_y,ori_in_par_z,loc_in_chd_x,loc_in_chd_y,loc_in_chd_z,ori_in_chd_x,ori_in_chd_y,ori_in_chd_z
      ground_bone1,ground,bone1,0,0,0,-3.9756,-6.361109,7.95138,0,0,0,-3.975693,-6.36110936,7.95138
      bone1_bone2,bone1,bone2,-2.2907,3.55,1.421,-2.5842,3.05333243,-3.0215214,0,0,0,2.484,0,-0
      ```
   - **Additional Notes**:
      As indicated in the data description, the data in this CSV file are measured relative to each body's "local" coordinate system.

#### `MTU_Parameters.csv`
   - **Format Specs**:
      | Column Name     | Description                                                                                 | Data Type                | Unit Type    |
      |-----------------|---------------------------------------------------------------------------------------------|--------------------------|--------------|
      | `mtu_name`      | The name of the MTU.                                                                        | Text                     | None         |
      | `max_iso_frc`   | The maximum isometric force of the muscle.                                                  | Numeric (floating-point) | Force        |
      | `opt_fib_len`   | The optimal fiber length of the muscle                                                      | Numeric (floating-point) | Length       |
      | `tendon_slk_len`| The resting length of the tendon.                                                           | Numeric (floating-point) | Length       |
      | `penn_ang`      | The angle between the tendon and the fiber when the fiber is at its optimal resting length. | Numeric (floating-point) | Angle        |
      | `def_fib_len`   | The default fiber length of the muscle.                                                     | Numeric (floating-point) | Length       |

   - **Example Entry**:
      ```command-line
      mtu_name,max_iso_frc,opt_fib_len,tendon_slk_len,penn_ang,def_fib_len
      MTU1_1,0.666666667,26.08996044,30.3014428,0,32.61245055
      MTU1_2,0.666666667,26.08996044,39.51839526,0,32.61245055
      ```
   
#### `MTU_Path_Landmarks.csv`
   - **Format Specs**:
      | Column Name     | Description                                                    | Data Type                | Unit Type    |
      |-----------------|----------------------------------------------------------------|--------------------------|--------------|
      | `mtu_name`      | The name of the MTU.                                           | Text                     | None         |
      | `pt_idx`        | The index of the path point.                                   | Numeric (integer)        | None         |
      | `pt_type`       | The type of the path point.                                    | Text                     | None         |
      | `attached_body` | The name of the body to which the path point is attached.      | Text                     | None         |
      | `pt_x`          | X-coordinate of the path point.                                | Numeric (floating-point) | Length       |
      | `pt_y`          | Y-coordinate of the path point.                                | Numeric (floating-point) | Length       |
      | `pt_z`          | Z-coordinate of the path point.                                | Numeric (floating-point) | Length       |

   - **Example Entry**:
      ```command-line
      mtu_name,pt_idx,pt_type,attached_body,pt_x,pt_y,pt_z
      MTU1_1,1,origin,bone1,22.7088,-53.5384,70.4372
      MTU1_1,2,viapoint,bone2,20.1314,-38.8696,18.7648
      MTU1_1,3,viapoint,bone2,20.0868,-38.8129,17.705
      MTU1_1,4,insertion,bone3,19.8304,-38.4468,13.0176
      ```
   - **Additional Notes**:
      - For each MTU listed in this CSV file, the associated point indices, i.e., `pt_idx` should form a complete, consecutive integer sequence from 1 to the number of path points associated with the MTU.
      - The `pt_type` value should be one of the following: `origin`, `viapoint`, `insertion`.
      - The location of the path points in the CSV file should be measured with respect to the "global" coordinate system.
      
#### `MTU_Path.csv`
   - **Format Specs**:
      | Column Name     | Description                                                    | Data Type                | Unit Type    |
      |-----------------|----------------------------------------------------------------|--------------------------|--------------|
      | `mtu_name`      | The name of the MTU.                                           | Text                     | None         |
      | `pt_idx`        | The index of the path point.                                   | Numeric (integer)        | None         |
      | `pt_type`       | The type of the path point.                                    | Text                     | None         |
      | `attached_body` | The name of the body to which the path point is attached.      | Text                     | None         |
      | `loc_in_body_x` | X-coordinate of the path point with respect to the body frame. | Numeric (floating-point) | Length       |
      | `loc_in_body_y` | Y-coordinate of the path point with respect to the body frame. | Numeric (floating-point) | Length       |
      | `loc_in_body_z` | Z-coordinate of the path point with respect to the body frame. | Numeric (floating-point) | Length       |

   - **Example Entry**:
      ```command-line
      mtu_name,pt_idx,pt_type,attached_body,loc_in_body_x,loc_in_body_y,loc_in_body_z
      MTU1_1,1,origin,bone1,-3.01385,1.691,-1.9311
      MTU1_1,2,viapoint,bone2,-1.3231,1.6582,-1.2048
      MTU1_1,3,viapoint,bone2,-2.5636,1.7682,-1.1611
      MTU1_1,4,insertion,bone3,-0.8143,1.6356,-1.1952
      ```
   - **Additional Notes**:
      - For each MTU listed in this CSV file, the associated point indices, i.e., `pt_idx` should form a complete, consecutive integer sequence from 1 to the number of path points associated with the MTU.
      - The `pt_type` value should be one of the following: `origin`, `viapoint`, `insertion`.
      - The location of the path points provided in the CSV file should be measured relative to each body's "local" coordinate system.

#### `Wrapping_Surface_Landmarks.csv`
   - **Format Specs**:
      | Column Name       | Description                                                                                                  | Data Type                | Unit Type |
      |-------------------|--------------------------------------------------------------------------------------------------------------|--------------------------|-----------|
      | `wrap_surf_name`  | The name of the wrapping surface.                                                                            | Text                     | None      |
      | `attached_body`   | The name of the body to which the wrapping surface is attached.                                              | Text                     | None      |
      | `org_pt_x`        | X-coordinate of the landmark point, marking the origin of the expected wrapping surface's frame.             | Numeric (floating-point) | Length    |
      | `org_pt_y`        | Y-coordinate of the landmark point, marking the origin of the expected wrapping surface's frame.             | Numeric (floating-point) | Length    |
      | `org_pt_z`        | Z-coordinate of the landmark point, marking the origin of the expected wrapping surface's frame.             | Numeric (floating-point) | Length    |
      | `pos_xaxis_pt_x`  | X-coordinate of the landmark point, indicating the positive X-axis of the expected wrapping surface's frame. | Numeric (floating-point) | Length    |
      | `pos_xaxis_pt_y`  | Y-coordinate of the landmark point, indicating the positive X-axis of the expected wrapping surface's frame. | Numeric (floating-point) | Length    |
      | `pos_xaxis_pt_z`  | Z-coordinate of the landmark point, indicating the positive X-axis of the expected wrapping surface's frame. | Numeric (floating-point) | Length    |
      | `pos_yaxis_pt_x`  | X-coordinate of the landmark point, indicating the positive Y-axis of the expected wrapping surface's frame. | Numeric (floating-point) | Length    |
      | `pos_yaxis_pt_y`  | Y-coordinate of the landmark point, indicating the positive Y-axis of the expected wrapping surface's frame. | Numeric (floating-point) | Length    |
      | `pos_yaxis_pt_z`  | Z-coordinate of the landmark point, indicating the positive Y-axis of the expected wrapping surface's frame. | Numeric (floating-point) | Length    |
      | `type`            | The type of the wrapping surface in OpenSim.                                                                 | Text                     | None      |
      | `geom_params`     | The parameters for configuring the size of the wrapping surface.                                             | Numeric (floating-point) | Length    |
      | `quadrant`        | The name of specific quadrants over which the wrapping surface is effective.                                 | Text                     | None      |

   - **Example Entry**:
      ```command-line
      wrap_surf_name,attached_body,org_pt_x,org_pt_y,org_pt_z,pos_xaxis_pt_x,pos_xaxis_pt_y,pos_xaxis_pt_z,pos_yaxis_pt_x,pos_yaxis_pt_y,pos_yaxis_pt_z,type,geom_params,quadrant
      wrap1,bone1,19.6911,-42.6982,27.2252,19.7525,-44.4717,29.8136,20.3724,-42.6882,27.219,Cylinder,"2,10",+x
      wrap2,bone2,19.6297,-40.9247,24.6368,19.7862,-42.4761,27.4456,20.2993,-40.8902,24.6092,Sphere,2,all
      ```

   - **Additional Notes**:
      - The coordinates in this CSV file are with respect to a single "global" coordinate system as this CSV file is part of the "global" data.
      - The current ArborSim setup supports two types of built-in wrapping surfaces in OpenSim: WrapSphere and WrapCylinder.
      - In OpenSim, the origin of the frame of a WrapSphere or a WrapCylinder is their center.
      - In OpenSim, the Z-axis of a WrapCylinder runs along the axial direction, while the X and Y axes lie in the radial plane of the cylinder's circular cross-section. One shoud keep this information in mind when identifying the landmark points for the wrapping surface frame to ensure the resulting orientation of WrapCylinder matches what is expected.
      - Regarding the `geom_params` in the CSV file, it is required to provide a single numeric value, representing the radius, for WrapSphere, and two numeric values, with the first for the radius, and the second for the length, for WrapCylinder.
      - The `quadrant` value in the CSV file should be one of the following elements: `-x`, `+x`, `x`, `-y`, `+y`, `y`, `-z`, `+z`, `z`, `all`. The value indicate the side over which is wrap object is effective. For example, `+x` indicates the positive X axis side of the wrapping surface is set to be effective to affect the the MTU or Ligament pathways.

#### `Wrapping_Surfaces.csv`
   - **Format Specs**:
      | Column Name      | Description                                                                                                   | Data Type                | Unit Type |
      |------------------|---------------------------------------------------------------------------------------------------------------|--------------------------|-----------|
      | `wrap_surf_name` | The name of the wrapping surface.                                                                             | Text                     | None      |
      | `attached_body`  | The name of the body to which the wrapping surface is attached.                                               | Text                     | None      |
      | `loc_in_body_x`  | X-coordinate of the wrapping surface's frame in the frame of the attached body.                               | Numeric (floating-point) | Length    |
      | `loc_in_body_y`  | Y-coordinate of the wrapping surface's frame in the frame of the attached body.                               | Numeric (floating-point) | Length    |
      | `loc_in_body_z`  | Z-coordinate of the wrapping surface's frame in the frame of the attached body.                               | Numeric (floating-point) | Length    |
      | `ori_in_body_x`  | X-axis rotation in the XYZ body-fixed Euler angles of the wrapping surface orientation in the parent frame.   | Numeric (floating-point) | Length    |
      | `ori_in_body_y`  | Y-axis rotation in the XYZ body-fixed Euler angles of the wrapping surface orientation in the parent frame.   | Numeric (floating-point) | Length    |
      | `ori_in_body_z`  | Z-axis rotation in the XYZ body-fixed Euler angles of the wrapping surface orientation in the parent frame.   | Numeric (floating-point) | Length    |
      | `type`           | The type of the wrapping surface in OpenSim.                                                                  | Text                     | None      |
      | `geom_params`    | The parameters for configuring the size of the wrapping surface.                                              | Numeric (floating-point) | Length    |
      | `quadrant`       | The name of specific quadrants over which the wrapping surface is effective.                                  | Text                     | None      |

   - **Example Entry**:
      ```command-line
      wrap_surf_name,attached_body,loc_in_body_x,loc_in_body_y,loc_in_body_z,ori_in_body_x,ori_in_body_y,ori_in_body_z,type,geom_params,quadrant
      wrap1,bone1,0,0,0,1.987846,-3.180554,9.939233,Cylinder,"2,10",+x
      wrap2,bone2,0,0,0,-4.9696,0,5.963540027,Sphere,"2",all
      ```

   - **Additional Notes**:
      - The coordinates in this CSV file are with respect to a single "global" coordinate system as this CSV file is part of the "global" data.
      - The current ArborSim setup supports two types of built-in wrapping surfaces in OpenSim: WrapSphere and WrapCylinder.
      - In OpenSim, the origin of the frame of a WrapSphere or a WrapCylinder is their center.
      - In OpenSim, the Z-axis of a WrapCylinder runs along the axial direction, while the X and Y axes lie in the radial plane of the cylinder's circular cross-section.
      - Regarding the `geom_params` in the CSV file, it is required to provide a single numeric value, representing the radius, for WrapSphere, and two numeric values, with the first for the radius, and the second for the length, for WrapCylinder.
      - The `quadrant` value in the CSV file should be one of the following elements: `-x`, `+x`, `x`, `-y`, `+y`, `y`, `-z`, `+z`, `z`, `all`. The value indicate the side over which is wrap object is effective. For example, `+x` indicates the positive X axis side of the wrapping surface is set to be effective to affect the the MTU or ligament pathways.

#### `Wrapping_Surface_MTU_Ligament_Pairings.csv`
   - **Format Specs**:
      | Column Name           | Description                                                                       | Data Type                | Unit Type    |
      |-----------------------|-----------------------------------------------------------------------------------|--------------------------|--------------|
      | `wrap_surf_name`      | The name of the wrapping surface.                                                 | Text                     | None         |
      | `paired_mtu_ligament` | The names of the mtus and/or ligaments that are linked with the wrapping surface. | Text                     | None         |
      | `active`              | The flag indicating is the wrapping surface is active or inactive.                | Numeric (integer)        | None         |

   - **Example Entry**:
      ```command-line
      wrap_surf_name,paired_mtu_ligament,active
      wrap1,"MTU1,MTU2,MTU5",0
      wrap2,MTU4,1
      ```
   - **Additional Notes**:
      The `paired_mtu_ligament` value should be either `0` or `1`, with `0` representing inactive, and `1` representing active.

#### `Ligament_Parameters.csv`
   - **Format Specs**:
      | Column Name     | Description                                                                       | Data Type                | Unit Type    |
      |-----------------|-----------------------------------------------------------------------------------|--------------------------|--------------|
      | `ligament_name` | The name of the ligament.                                                         | Text                     | None         |
      | `frc_scale`     | The magnitude for scaling the normalized force-length curve of the ligament.      | Numeric (floating-point) | None         |
      | `rest_len`      | The resting length of the ligament.                                               | Numeric (floating-point) | Length       |

   - **Example Entry**:
      ```command-line
      ligament_name,frc_scale,rest_len
      group1_ligament1,0.666666667,1
      group1_ligament2,1.333333334,1
      ```

#### `Ligament_Path.csv`
   - **Format Specs**:
      | Column Name     | Description                                                    | Data Type                | Unit Type    |
      |-----------------|----------------------------------------------------------------|--------------------------|--------------|
      | `ligament_name` | The name of the ligament.                                      | Text                     | None         |
      | `pt_idx`        | The index of the path point.                                   | Numeric (integer)        | None         |
      | `pt_type`       | The type of the path point.                                    | Text                     | None         |
      | `attached_body` | The name of the body to which the path point is attached.      | Text                     | None         |
      | `loc_in_body_x` | X-coordinate of the path point with respect to the body frame. | Numeric (floating-point) | Length       |
      | `loc_in_body_y` | Y-coordinate of the path point with respect to the body frame. | Numeric (floating-point) | Length       |
      | `loc_in_body_z` | Z-coordinate of the path point with respect to the body frame. | Numeric (floating-point) | Length       |

   - **Example Entry**:
      ```command-line
      ligament_name,pt_idx,pt_type,attached_body,loc_in_body_x,loc_in_body_y,loc_in_body_z
      group1_ligament1,1,origin,float1,0,0,0
      group1_ligament1,2,insertion,caudal5,-0.22847,0.95-1.313
      ```
   - **Additional Notes**:
      - For each ligament listed in this CSV file, the associated point indices, i.e., `pt_idx` should form a complete, consecutive integer sequence from 1 to the number of path points associated with the ligament.
      - The `pt_type` value should be one of the following: `origin`, `viapoint`, `insertion`.
      - The location of the path points provided in the CSV file should be measured relative to each body's "local" coordinate system.

### Geometry File Format
   - The geometry files of the bodies of a model must be supplied under the `Geometry` folders. And the name of each geometry file should match the `body_name` supplied in the CSV files explained above.

   - Geometry files under the `GlobalData` directory are stored in a "global" manner. Meaning, all these geometry files share the same "global" coordinate system. This coordinate system can be whatever default or customized coordinate system in which the geometry files are segmented and exported. In contrast, geometry files under the `LocalData` and `AdjustedLocalData` directories are stored in a "local" manner, with each associated with the respective "local" body frame. (Transformation for the geometry files from "global" to "local" can be realized by function [`apply_trans_mat_to_geom`](../Scripts/Core/Transformer.m#L206).) 


### Important Notes
   - For data that is shared across different CSV files, such as body's name, joint's name, etc., users must ensure consistency across all respective CSV files.

   - For all the CSV files under the `GlobalData` directory, users must ensure that data belonging to the same unit type share the same unit across all CSV files. For example, for data involving the Length unit type in any CSV file, users need to ensure they are all measured in the same unit, such as meters or centimeters. This rule applies to every unit type. Furthermore, the geometry files under the `GlobalData` directory should use the same units as the CSV files. Specifically, the current ArborSim setup supports STL files for body visualizations. Although STL files are by default unitless, users need to ensure the length unit used when generating or exporting these STL files matches the length unit used in the CSV files.

   - The same principles apply to all CSV files and geometry files under the `LocalData` and `AdjustedLocalData` directories. Additionally, the `LocalData` and `AdjustedLocalData` directories should use the same unit for data of the same unit type, as the current ArborSim setup adjusts the "local" data without changing any units.

   - Regarding specific units supported in the current ArborSim setup: for length, the tool supports meters (abbreviated as `m`), centimeters (abbreviated as `cm`), and millimeters (abbreviated as `mm`); for mass, kilograms (abbreviated as `kg`) and grams (abbreviated as `g`); for angle, only degrees (abbreviated as `deg`); and for force, currently only newtons (abbreviated as `N`).

