function opts = ktrackOpts( )

%mode = 'NoCompNoU_HiC';
%mode = 'NoCompNoU_LoC';
mode = 'YesCompNoU';

if strcmp(mode,'NoCompNoU_LoC')
  b_compensateMotion = false(); %true();
  b_computeHorizon   = false();
  C_iters             = 10; % Iffy with 10, 5 fails badly
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
  C_iters            = 5;
end

% TODO: have these options write to the vtk run script too, 
%   so that scenario is setup dynamically... e.g. :
% echo "$PATH/simulator_client -C 0 -W 1  -b galaxy2.jpg 
%          -s 3.0  -F 0  -m 02_Obj_M_Low.obj  -t MiG-35.obj"  > runme.sh 

opts = struct('output_port',5001,'number_of_retries',1000,...
    'compensation', b_compensateMotion,...
    'horizon', b_computeHorizon,...
    'max_input_frames',i_maxInputFrames,'contour_iters',C_iters);
  
disp(opts);


end
