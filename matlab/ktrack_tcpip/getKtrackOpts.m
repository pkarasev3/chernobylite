function opts = getKtrackOpts( )

%mode = 'NoCompNoU_LoC';
%mode = 'NoCompNoU_HiC';
%mode = 'YesCompNoU';
mode = 'YesCompNoU_EqC';

saveRate = 0;
showRate = 5;

if strcmp(mode,'NoCompNoU_LoC')
  b_compensateMotion = false(); %true();
  b_computeHorizon   = false();
  C_iters             = 10; %7; % Iffy with 10, 5 fails badly
  i_maxInputFrames    = 200;
elseif strcmp(mode,'NoCompNoU_HiC')
  b_compensateMotion = false(); %true();
  b_computeHorizon   = false();
  C_iters             = 30; % Should work with 30 
  i_maxInputFrames    = 200;  
elseif strcmp(mode,'YesCompNoU')
  b_compensateMotion = true();
  b_computeHorizon   = false();
  i_maxInputFrames   = 200;
  C_iters            = 5; %different count as low, but same *timefactor*
elseif strcmp(mode,'YesCompNoU_EqC')
  b_compensateMotion = true();
  b_computeHorizon   = false();
  i_maxInputFrames   = 1000;
  C_iters            = 10; %same count as "lo"
  showRate           = 1;
end

res_fname           = [mode '_' num2str(C_iters) '_results.mat'];

% TODO: have these options write to the vtk run script too, 
%   so that scenario is setup dynamically... e.g. :
% echo "$PATH/simulator_client -C 0 -W 1  -b galaxy2.jpg 
%          -s 3.0  -F 0  -m 02_Obj_M_Low.obj  -t MiG-35.obj"  > runme.sh 

opts = struct('output_port',5001,'number_of_retries',10,...
    'compensation', b_compensateMotion,...
    'horizon', b_computeHorizon,...
    'max_input_frames',i_maxInputFrames,'contour_iters',C_iters,...
    'getPsiTru',true,...
    'saveImages',saveRate,...
    'showImages',showRate,...
    'result_filename',res_fname);
  
disp(opts);


end
