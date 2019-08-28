#!/bin/bash

# Set your own EuRoC_PATH path to run ice-ba. Use "./bin/ice_ba --help" to get the explanation for all of the flags. Flags [imgs_folder] and [iba_param_path] are necessary.
# Add flag '--save_feature' to save feature message and calibration file for back-end only mode

#005
EuRoC_PATH=~/Dataset/NEAR/005
Dataset_Name=ip7
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

Dataset_Name=p20
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

#014
EuRoC_PATH=~/Dataset/NEAR/014
Dataset_Name=1p5t
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

Dataset_Name=ipxr
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

#018
EuRoC_PATH=~/Dataset/NEAR/018
Dataset_Name=1p5t
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

Dataset_Name=ipxr
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

#022
EuRoC_PATH=~/Dataset/NEAR/022
Dataset_Name=ip7
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

Dataset_Name=p20
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

#025
EuRoC_PATH=~/Dataset/NEAR/025
Dataset_Name=ip7
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

Dataset_Name=p20
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

#010
EuRoC_PATH=~/Dataset/NEAR/010
Dataset_Name=1p5t
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

Dataset_Name=ipxr
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

#037
EuRoC_PATH=~/Dataset/NEAR/037
Dataset_Name=1p5t
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

Dataset_Name=ipxr
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

#039
EuRoC_PATH=~/Dataset/NEAR/039
Dataset_Name=ip7
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

Dataset_Name=p20
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

#045
EuRoC_PATH=~/Dataset/NEAR/045
Dataset_Name=ip7
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

Dataset_Name=p20
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

#052
EuRoC_PATH=~/Dataset/NEAR/052
Dataset_Name=ip7
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd

Dataset_Name=p20
mkdir $EuRoC_PATH/$Dataset_Name/result
cmd="../bin/ice_ba --imgs_folder $EuRoC_PATH/$Dataset_Name --start_idx 0 --end_idx -1 --iba_param_path ../config/config_of_mono.txt  --gba_camera_save_path $EuRoC_PATH/$Dataset_Name/result/$Dataset_Name.txt"
echo $cmd
eval $cmd
