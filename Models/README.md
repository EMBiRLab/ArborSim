## Directory Structure and Description

This directory is structured to support the construction of OpenSim models using various data formats and methodologies. In the `Models` directory, users can create individual folders for each model or a master folder for groups of models. For individual models, populate `Data` and `Output` subfolders within the model's folder. For multiple models, establish separate subfolders for each model within the master folder. Each model subfolder should also be structured with its own `Data` and `Output` subfolders. 

### Starter Folders in the Directory
The following folders are included in the repository at its initial setup. They are designed to help users understand the structure, get started with model creation, and reproduce the results in the accompanying paper.

#### - `Hello_ArborSim`
  An introductory example for users new to the repository. The folder contains CSV files pertinent to the model construction. It offers a hands-on experience in constructing simple OpenSim models using pre-filled data. Refer to [create_my_first_model.m](../Scripts/create_my_first_model.m) for details in constructing the models from the provided data.

#### - `Empty_Template`
  An empty template folder intended to serve as a model blueprint. It includes CSV files that provide a format guide with column names but no actual data filled in, and placeholders for geometries. These files collectively demonstrate the file structure expected for new models.

#### - `ArborSim_Paper_Toy_Models`
  A collection of models used to conduct the study in the accompanying paper. The CSV files contained include essential data for building the models. One can reproduce the results in the accompanying paper accordingly.

### Model Root Folder Structure and Contents
The term "model root folder" here refers to the parent folder that directly houses the Data and Output subfolders. It doesn't refer to any higher-level directory, if one exists.
The folllowing explains the components of the `Data` and `Output` subfolders within the "model root folder".

#### `Data`
This folder serves as the primary repository of source data for model construction. It is divided into three distinct data formats:

##### - `GlobalData`
   - **Description**: This folder is equivalent in Layer 1 explained in the accompanying paper. It contains data acquired in a single global reference frame, which is the default or customized coordinate systems used in biomedical imaging source data.
   - **Contents**:
     - `Geometry`: Includes segmented body geometries in the global reference frame.
     - CSV Files: 
       1. `Artificial_OpenSim_Geometries.csv`
       2. `Body_Frame_Landmarks.csv`
       3. `Body_Mass_Properties.csv`
       4. `Branch_Groups.csv`
       5. `JCS_Landmarks.csv`
       6. `Joint_ROM.csv`
       7. `MTU_Parameters.csv`
       8. `MTU_Path_Landmarks.csv`
       9. `Wrapping_Surface_Landmarks.csv`
       10. `Wrapping_Surface_MTU_Ligament_Pairings.csv`

##### - `LocalData`
   - **Description**: This folder is equivalent in Layer 2 explained in the accompanying paper. It contains data in local reference systems, tailored for individual bodies within the model.
   - **Contents**:
     - `Geometry`: Segmented body geometries adapted to each body's local coordinate system.
     - CSV Files:
       1. `Artificial_OpenSim_Geometries.csv`
       2. `Body_Mass_Properties.csv`
       3. `Branch_Groups.csv`
       4. `JCS.csv`
       5. `Joint_ROM.csv`
       6. `MTU_Parameters.csv`
       7. `MTU_Path.csv`
       8. `Wrapping_Surface_MTU_Ligament_Pairings.csv`
       9. `Wrapping_Surfaces.csv`

##### - `AdjustedLocalData`
   - **Description**: This folder is equivalent in Layer 3 explained in the accompanying paper. It features data in local format, adapted for use with the proposed branching modeling method in the accompanying paper.
   - **Contents**:
     - `Geometry`: Segmented body geometries adapted to each body's local coordinate system.
     - CSV Files:
       1. `Artificial_OpenSim_Geometries.csv`
       2. `Body_Mass_Properties.csv`
       3. `JCS.csv`
       4. `Joint_ROM.csv`
       5. `Ligament_Parameters.csv`
       6. `Ligament_Path.csv`
       7. `MTU_Parameters.csv`
       8. `MTU_Path.csv`
       9. `Wrapping_Surface_MTU_Ligament_Pairings.csv`
       10. `Wrapping_Surfaces.csv`

#### `Output`
This folder is equivalent to Layer 4 explained in the accompanying paper. It contains the OpenSim models constructed using different branching modeling methods, represented in two subfolders:

##### - `Proposed`
   - **Description**: Includes the model developed using the proposed branching modeling method.
   - **Contents**:
     - `Geometry`: Essential STL files for model visualization.
     - `*.osim`: The constructed OpenSim model using the proposed method.

##### - `Conventional`
   - **Description**: Features the model built using traditional approaches to branching muscle-tendon architectures.
   - **Contents**:
     - `Geometry`: Essential STL files for model visualization.
     - `*.osim`: The constructed OpenSim model using the conventional method.

### Additional Information
For comprehensive details on the proposed branching modeling method, modifications to data files, and model construction, please refer to the accompanying paper and the scripts in the [Scripts](../Scripts/) folder. For instructions on how to supply geometry files, and fill out each CSV file, please refer to the [DATA_FORMAT_SPECS.md](./DATA_FORMAT_SPECS.md) file.

