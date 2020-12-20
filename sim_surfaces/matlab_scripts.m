% Produce PV estimates for the simulated surfaces
% This script handles the ground truth, NeuroPVE and RC methods
% Toblerone and resampling are handled in the main python script 
% (run_sim_surfaces.py) that then calls this script in turn, so you
% don't need to run this directly. 

ROOT = char(ROOT); 
addpath(genpath(pwd));

voxSizes = 1 : 0.2 : 3;
bumpysphDs = [ 60 1.05 ]; 

filenameLambda = @(method,v) [ROOT, sprintf('%s_%s.nii',...
  method,sprintf('%1.2f',v))];
surfnameLambda = @(surf) [ROOT '/' surf];

for vS = voxSizes
  
  voxSize = [vS, vS, vS];
  
  Ds = bumpysphDs;
  imgSize = ceil([150 150 150] ./ voxSize);
  
  d1 = Ds(1); d2 = Ds(2); hcpFillDist = mean([d1, d1*d2]);
  bumpyFunc = @(a1,a2) ( @(x,y,z) bumpySphereClassifier(a1,a2,x,y,z) );
  theFunc = bumpyFunc(d1,d2);
  
  inSurfName = surfnameLambda('lh.white');
  outSurfName = surfnameLambda('lh.pial');
  surfExists = exist(inSurfName, 'file');

  %   Check if the simulated surfaces exist already and read in
  if surfExists
    [inPs, inTris] = read_surf(inSurfName);
    [outPs, outTris] = read_surf(outSurfName);
    
    % or generate them now
  else
    
    [inPs, inTris] = spheretribydepth(7);
    [theta, phi, ~] = cart2sph(inPs(:,1), inPs(:,2), inPs(:,3));
    [rIn, rOut] = bumpySphereFunction(d1, d2, inPs(:,1), inPs(:,2), inPs(:,3));
    
    [x1,y1,z1] = sph2cart(theta, phi, rIn);
    [x2,y2,z2] = sph2cart(theta, phi, rOut);
    inPs = [x1, y1, z1];
    outPs = [x2, y2, z2];
    outTris = inTris;
    
    % Write surface files.
    write_surf(inSurfName, inPs, inTris);
    write_surf(outSurfName, outPs, outTris);
    
  end
  
  trufile = filenameLambda('tru',voxSize(1));
  neurofile = filenameLambda('neuro1',voxSize(1));
  reffile = filenameLambda('ref',voxSize(1));
  rcfile = filenameLambda('rc',voxSize(1));
  
  % To ensure unix calls work properly. 
  if ismac
    sourceCmd = 'source ~/.bash_profile;';
  elseif isunix
    sourceCmd = 'source ~/.profile;source ~/.bash_profile;';
  end
  
  % When producing each method's estimates, we load the header from the
  % reference
  refNII = load_nii(reffile);
  refNII.hdr.dime.dim(1) = 4; 
  refNII.hdr.dime.dim(5) = 3;

  vox2world = [ refNII.hdr.hist.srow_x; refNII.hdr.hist.srow_y;
    refNII.hdr.hist.srow_z; 0 0 0 1];
  
  % True solution via numerical method
  if ~exist(trufile,'file')
    disp('Generating true solution');
    super = 2 * round(15*vS/2) + 1;
    voxList = 1:prod(imgSize);
    tru = voxelIntegrator(imgSize, voxSize, vox2world, ...
      super, theFunc, voxList);
    tru = reshape(tru, imgSize(1), imgSize(2), imgSize(3), 3);
    
    % Save as NIFTI
    nii = refNII;
    nii.img = tru;
    save_nii(nii, trufile);
  end
  
  % Neuropoly's estimate
  if ~exist(neurofile, 'file') % && false
    
    Norig = eye(4); Torig = Norig;
    fsdir = fullfile(ROOT, 'surf'); 
    neuro = pvc_partial_volume_estimation(reffile, ROOT, 'lh', voxSize(1), '', 5, false, false, true, vox2world, false, 0);
    nii = refNII;
    nii.img = neuro;
    save_nii(nii, neurofile);
        
  end
    
  % RC's estimate 
  if ~exist(rcfile, 'file')
    
    % Convert binary surface files to GII
    inSurfGII = fullfile(ROOT, 'surf', 'white.surf.gii');
    if ~exist(inSurfGII, 'file')
      convertCmd = [sourceCmd 'mris_convert ' inSurfName ' ' inSurfGII];
      unix(convertCmd);
    end
    
    outSurfGII = fullfile(ROOT, 'surf', 'pial.surf.gii');
    if ~exist(outSurfGII, 'file')
      convertCmd = [sourceCmd 'mris_convert ' outSurfName ' ' outSurfGII];
      unix(convertCmd);
    end
    
    % Generate the midthickness surface
    midSurfGII = fullfile(ROOT, 'surf', 'mid.surf.gii');
    midCmd = [sourceCmd 'wb_command -surface-cortex-layer ' inSurfGII ' ' outSurfGII ' 0.5 ' midSurfGII];
    unix(midCmd); 
    
    RCmethod('LWS', inSurfGII, 'LPS', outSurfGII, 'LMS', midSurfGII, ...
      'reference', reffile, 'outDir', ROOT, 'outName', rcfile);
    
  end
  
end

