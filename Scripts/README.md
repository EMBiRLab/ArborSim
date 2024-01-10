## Directory Structure and Description

This directory includes the following default folders and files:

### - `Core`
This folder contains three core classes: `Transformer`, `Brancher`, and `Builder`. The `Transformer` is responsible for transforming "global" data into "local" data. The `Brancher` adjusts the "local" data to align with the branching muscle-tendon architecture method described in the accompanying paper. The `Builder` assembles the components to create respective OpenSim musculoskeletal models.

### - `External`
This folder houses external functions that are not self-written but are used as helper functions. Currently, it includes a function for generating multiple, visibly distinct colors, enhancing the visualization quality in plots.

### - `Utility`
This folder contains self-written helper functions. As of now, it includes scripts for calculating body mass properties from STL files.

### - `create_my_first_model.m`
This script is designed for new users, guiding them in creating introductory example models.

### - `create_arborsim_paper_toy_models.m`
This script generates the toy models used in the accompanying paper, facilitating result reproduction.

### - `update_arborsim_paper_toy_model_params.m`
This script, called within `create_arborsim_paper_toy_models.m`, manually updates some physiological parameters of MTUs and/or ligaments. It's used after modifying the "local" data with the proposed branching modeling method.

### Additional Information
- Each method in the classes and scripts comes with detailed documentation within the respective files. Please consult these, along with the accompanying paper, to familiarize yourself with the tool.
- To create any musculoskeletal model, please place the data in the directory: `../Models`, following the instructions in [README](../Models/README.md), and refer to the default scripts `create_arborsim_paper_toy_models.m` and `create_my_first_model.m` in this directory to develop your script for model construction.
