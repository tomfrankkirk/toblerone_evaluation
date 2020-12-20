% Ribbon constrained method for estimating cortical PVs 
% With thanks to Tim Coalson, on who's suggestion this method 
% is entirely based 

function [] = RCmethod(varargin)

if ismac
  sourceCmd = 'source ~/.bash_profile;';
elseif isunix
  sourceCmd = 'source ~/.profile;source ~/.bash_profile;';
end

outInd = find(cellfun(@(c) strcmp(c, 'outDir'), varargin));
if ~isempty(outInd)
  outDir = varargin{outInd + 1};
else
  error('RC error: an ouput directory must be provided');
end

% Create output dir if it does not already exist.
if isstring(outDir); outDir = char(outDir); end
if ~exist(outDir, 'dir')
  mkdir(outDir);
end

preInd = find(cellfun(@(c) strcmp(c, 'outName'), varargin));
if ~isempty(preInd)
  outName = varargin{preInd + 1};
else
  outName = '';
end

% Surface paths
LWSind = find(cellfun(@(c) strcmp(c, 'LWS'), varargin));
if ~isempty(LWSind)
  LWS = varargin{LWSind + 1};
else
  LWS = '';
end

LPSind = find(cellfun(@(c) strcmp(c, 'LPS'), varargin));
if ~isempty(LPSind)
  LPS = varargin{LPSind + 1};
else
  LPS = '';
end

LMSind = find(cellfun(@(c) strcmp(c, 'LMS'), varargin));
if ~isempty(LMSind)
  LMS = varargin{LMSind + 1};
else
  LMS = '';
end

RPSind = find(cellfun(@(c) strcmp(c, 'RPS'), varargin));
if ~isempty(RPSind)
  RPS = varargin{RPSind + 1};
else
  RPS = '';
end

RWSind = find(cellfun(@(c) strcmp(c, 'RWS'), varargin));
if ~isempty(RWSind)
  RWS = varargin{RWSind + 1};
else
  RWS = '';
end

RMSind = find(cellfun(@(c) strcmp(c, 'RMS'), varargin));
if ~isempty(RMSind)
  RMS = varargin{RMSind + 1};
else
  RMS = '';
end

if (~strcmp('', RWS) && ~strcmp('', RPS)) ...
    && (strcmp('', LWS) && strcmp('', LPS))
  hemispheres = "R";
  assert(exist(RWS, 'file') == 2, 'Surface path does not exist');
  assert(exist(RPS, 'file') == 2, 'Surface path does not exist');
  % Detect file type of input surfaces
  [~, ~, surfext] = fileparts(RWS);
  
elseif (strcmp('', RWS) && strcmp('', RPS)) ...
    && (~strcmp('', LWS) && ~strcmp('', LPS))
  hemispheres = "L";
  assert(exist(LWS, 'file') == 2, 'Surface path does not exist');
  assert(exist(LPS, 'file') == 2, 'Surface path does not exist');
  % Detect file type of input surfaces
  [~, ~, surfext] = fileparts(LWS);
  
elseif (~strcmp('', RWS) && ~strcmp('', RPS)) ...
    && (~strcmp('', LWS) && ~strcmp('', LPS))
  hemispheres = ["L", "R"];
  assert(exist(RWS, 'file') == 2, 'Surface path does not exist');
  assert(exist(RPS, 'file') == 2, 'Surface path does not exist');
  assert(exist(LWS, 'file') == 2, 'Surface path does not exist');
  assert(exist(LPS, 'file') == 2, 'Surface path does not exist');
  % Detect file type of input surfaces
  [~, ~, surfext] = fileparts(LWS);
end

imgInd = find(cellfun(@(c) strcmp(c, 'reference'), varargin));
if ~isempty(imgInd)
  imagePath = varargin{imgInd + 1};
  assert(exist(imagePath, 'file') == 2, 'Image path does not exist');
else
  error('Toblerone error: a reference image space must be provided');
end

% Test the input image path
[~, fname, ext] = fileparts(imagePath);
if ~strcmp(ext, '.nii')
  error('Input image must be in NIFTI format');
end

% Read the input scale and vox2world matrix
imageNii = load_nii(imagePath);
voxSize = imageNii.hdr.dime.pixdim(2:4);
imgSize = imageNii.hdr.dime.dim(2:4);
vox2world = [
  imageNii.hdr.hist.srow_x;
  imageNii.hdr.hist.srow_y;
  imageNii.hdr.hist.srow_z;
  0 0 0 1;
  ];



for hemi = hemispheres
  h = char(hemi);
  
  v = sprintf('%1.1f', voxSize(1)); 
  
  midCmd = [sourceCmd, 'wb_command -surface-cortex-layer ' eval([h 'WS']) ' ' eval([h 'PS']) ' 0.5 ' eval([h 'MS']) ];
  unix(midCmd); 
  
  % Produce the template metric file that matches topology of the surface.
  templatefile = fullfile(outDir, [h '_metric_template.func.gii']);
  metricCmd = [sourceCmd 'wb_command -surface-vertex-areas ' eval([h 'WS']) ' ' templatefile];
  unix(metricCmd);
  
  % Now produce an equivalent metric file of ones only.
  onesfile = fullfile(outDir, [h '_ones.func.gii']);
  onesCmd = [sourceCmd 'wb_command -metric-math ''1'' ' onesfile ' -var x ' templatefile];
  unix(onesCmd);
  
  % Finally, map metric to volume
  subdiv = max(ceil(max(voxSize) / 0.4), 4);
  gmfile = [outDir '/' h v '_hcpgm.nii']; 
  volCmd = [sourceCmd 'wb_command -metric-to-volume-mapping -ribbon-constrained -voxel-subdiv ' num2str(subdiv)];
  volCmd = [volCmd ' ' eval([h 'WS']) ' ' eval([h 'PS']) ' ' onesfile ' ' eval([h 'WS']) ' ' imagePath ' ' gmfile ];
  unix(volCmd);
  
  
  % And now use volume math to evaluate the PVs in WM and CSF
  % respectively. The voxel-wise expression is (1-GM) * (signdist > ||
  % < 0) depending on whether CSF or GM is desired (ie, for those
  % voxels not wholly GM, label as WM or CSF depending on the sign dist
  % func)
  hcpFillDist = 60;
  sdfile = [outDir, '/' h '_signdist.nii']; 
  sdCmd = [sourceCmd 'wb_command -create-signed-distance-volume ' eval([h 'MS']) ' ' imagePath];
  sdCmd = [sdCmd ' ' sdfile ' -approx-limit ' char(num2str(hcpFillDist)) ' -fill-value ' char(num2str(hcpFillDist)) ];
  unix(sdCmd);
  
  
  % WM
  wmfile = [outDir '/' h v '_hcpwm.nii']; 
  wmCmd = [sourceCmd 'wb_command -volume-math ''(1 - GM) * (sd < 0)'' ' ...
    wmfile ' -var GM ' gmfile ' -var sd ' sdfile ];
  unix(wmCmd);
  
  % CSF
  csffile = [outDir '/' h v '_hcpcsf.nii']; 
  csfCmd = [sourceCmd 'wb_command -volume-math ''(1 - GM) * (sd > 0)'' ' ...
    csffile ' -var GM ' gmfile ' -var sd ' sdfile];
  unix(csfCmd);
  
  unix(['rm ' sdfile ' ' templatefile ' ' onesfile ]);
  
end

for hemi = hemispheres
  h = char(hemi);
  
  wmfile = [outDir '/' h v '_hcpwm.nii'];
  gmfile = [outDir '/' h v '_hcpgm.nii'];
  csffile = [outDir '/' h v '_hcpcsf.nii'];

  imgs = {};
  for f = {gmfile, wmfile, csffile}
    nii = load_nii(char(f));
    imgs{end+1} = nii.img;
  end
  unix(['rm ' gmfile ' ' csffile ' ' wmfile ]); 
  
  eval([h 'PVs = zeros(prod(imgSize),1);']);
  for c = 1:3
    eval([h 'PVs(:,' num2str(c) ') = imgs{' num2str(c) '}(:);']);
  end
  
end



if size(hemispheres,2) > 1
  
  % Initialise final mask as all CSF
  PVs = zeros(size(eval([h 'PVs'])));
  
  % Fix the GM down as the sum of each hemi
  PVs(:,1) = min(1, (LPVs(:,1) + RPVs(:,1))); 
  PVs(:,2) = min(1-PVs(:,1), LPVs(:,2) + RPVs(:,2)); 
  PVs(:,3) = 1 - sum(PVs(:,1:2),2); 
  
  assert(~any(PVs(:) > 1)); 
  assert(~any(PVs(:) < 0)); 
  assert(all(sum(PVs,2) < 1.0001)); 
  
else
  PVs = eval([char(hemispheres(1)) 'PVs']);
end

% Save output
PVs = reshape(PVs, imgSize(1), imgSize(2), imgSize(3), 3);

% Save onto a copy of the original NII
[~,fname,fext] = fileparts(outName); 
outPath = fullfile(outDir, [fname fext]);
fprintf('Saving PVs to %s\n', outPath);
nii = imageNii;
nii.hdr.dime.dim(1) = 4;
nii.hdr.dime.dim(5) = 3;
nii.img = PVs;
save_nii(nii, outPath);

end