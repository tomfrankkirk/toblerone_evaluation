%%% example

clc
close all
clear all

addpath(genpath('PVC'));
addpath(genpath('../../toblerone/FSmatlab')); 


%%% set paraneters

% fsdir: here refers to the actual fsdir/Subject_name
% the script will assume that the surfaces (xh.white, xh.pial)
% are in fsdir/surf/
fsdir = 'data_test/ctrl_test';

% hemisphere to consider ('lh' or 'rh')
hemi = 'lh';

% Volume that is subject to partial volume
Vol = [fsdir '/mri/T1.mgz'];

% Isotropic voxel resolution in mm
voxel = 3; % 1 mm isotroptic

% ouptut directory
output_path = [fsdir '/output'];
mkdir(output_path);

% number of surfaces that will be extended, adding more 
% surface will give a more precise estimation but will increase the
% computation time
nb_surface = 5; 

% second correction taking into account the normal of the surface
second_correction = false;

%%% Example 2:
% Run step by step (for the left hemisphere):

% Step 1, expending surfaces
% create #nb_surface expended surfaces in the positive and negative
% direction, for pial amd white surfaces; so 2* 2* #nb_surface
% surf created in [fsdir '/surf'], e.g. lh.pial0000, lh.pial0001,
% lh.pialn0001, etc.
% Long step, ~55min per surf, so ~20h on a single i7 processor.

compute_expended_surfaces = 1;
if compute_expended_surfaces
%     pvc_create_expanded_surfaces('lh.pial', [fsdir, '/surf'], nb_surface, voxel);
    pvc_create_expanded_surfaces('lh.white',[fsdir, '/surf'], nb_surface, voxel);
end

% Step 2, compute partial volumes and r-square coefficients of fitting
% ~30 min to compute on a single i7 processor

compute_partial_volumes = 1;
if compute_partial_volumes
    extend = false;
    PVs = pvc_partial_volume_estimation(Vol, [fsdir, '/surf'], hemi, voxel, output_path, nb_surface, extend, second_correction, false, eye(4));
    save_nii(make_nii(PVs, voxel), sprintf('demo%i.nii', voxel)); 
end



