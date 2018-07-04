% @brief examples of using OPA toolbox for different parameters and methods
% Please, run sections one by one to get an impression of how
% different ordinal-patterns-based methods work for different parameters

clear; close all; clc;
addpath( 'Data' ); % https://vis.caltech.edu/~rodri/data.htm (dataset 5, 130-2_c4.asc)

%% compute permutation entropy in sliding windows
load( 'tonicClonic.mat' );
cfg            = [];
cfg.method     = 'PE'; % compute permutation entropy
cfg.order      = 3;    % ordinal pattens of order 3 (4-points ordinal patterns)
cfg.delay      = 2;    % delay 2 between points in ordinal patterns 
                       % (one point between successive points in ordinal patterns)
cfg.windowSize = 512;  % window size = 512 time steps
cfg.time       = 0:1/102.4:179.999; % OPTIONAL time axis for plotting
cfg.units      = 'seconds';         % OPTIONAL units of time for plotting
outdata        = OPanalysis( cfg, indata );

%% compute permutation entropy and ordinal distributions in sliding windows
load( 'tonicClonic.mat' );
cfg            = [];
cfg.method     = 'opdPE'; % compute permutation entropy
cfg.order      = 3;       % ordinal pattens of order 3 (4-points ordinal patterns)
cfg.orderSeq   = 6;       % ordinal pattens of order 6 for plotting their sequence (7-points ordinal patterns)
cfg.delay      = 1;       % delay 1 between points in ordinal patterns (successive points)
cfg.windowSize = 1024;    % window size = 1024 time steps
cfg.time       = 0:1/102.4:179.999; % OPTIONAL time axis for plotting
cfg.units      = 'seconds';         % OPTIONAL units of time for plotting
outdata        = OPanalysis( cfg, indata );

%% compute all the implemented measures simultaneously for comparison
load( 'tonicClonic.mat' );
cfg                = [];
cfg.method         = 'all';  % compute all implemented ordinal-patterns-based measures
cfg.order          = 4;      % ordinal patterns of order 4 (5-points ordinal patterns)
cfg.delay          = 1;      % delay 1 between points in ordinal patterns
cfg.windowSize     = 512;    % window size = 512 time steps
cfg.lowerThreshold = 0.2;    % the distance considered negligible between points
cfg.upperThreshold = 200;    % the distance between points most probably related to artifact
cfg.time           = 0:1/102.4:179.999; % OPTIONAL time axis for plotting
cfg.units          = 'seconds';         % OPTIONAL units of time for plotting
outdata            = OPanalysis( cfg, indata );

%% compute conditional entropy of ordinal patterns in sliding windows
load( 'tonicClonic.mat' );
cfg            = [];
cfg.method     = 'CE'; % we compute conditional entropy of ordinal patterns
cfg.order      = 3;    % ordinal pattens of order 3 (4-points ordinal patterns)
cfg.delay      = 1;    % delay 1 between points in ordinal patterns (successive points)
cfg.windowSize = 512;  % window size = 512 time steps
cfg.time       = 0:1/102.4:179.999; % OPTIONAL time axis for plotting
cfg.units      = 'seconds';         % OPTIONAL units of time for plotting
outdata        = OPanalysis( cfg, indata );

%% compute robust permutation entropy
load( 'tonicClonic.mat' );
cfg                = [];
cfg.method         = 'rePE'; % compute robust permutation entropy
cfg.order          = 6;      % ordinal patterns of order 6 (7-points ordinal patterns)
cfg.delay          = 1;      % delay 1 between points in ordinal patterns
cfg.windowSize     = 2048;   % window size = 2048 time steps
cfg.lowerThreshold = 0.2;    % the distance that is considered negligible between points
cfg.upperThreshold = 100;    % the distance between points most probably related to artifact
cfg.time           = 0:1/102.4:179.999; % OPTIONAL time axis for plotting
cfg.units          = 'seconds';         % OPTIONAL units of time for plotting
outdata            = OPanalysis( cfg, indata );

%% compute permutation entropy for ordinal patterns with tied ranks in sliding windows
load( 'tonicClonic.mat' );
cfg            = [];
cfg.method     = 'PEeq'; % compute permutation entropy for ordinal patterns with tied ranks
cfg.order      = 3;      % ordinal pattens of order 3 (4-points ordinal patterns)
cfg.delay      = 3;      % delay 3 between points in ordinal patterns 
                         % (2 points between successive points in ordinal patterns)
cfg.windowSize = 1024;   % window size = 1024 time steps
cfg.time       = 0:1/102.4:179.999; % OPTIONAL time axis for plotting
cfg.units      = 'seconds';         % OPTIONAL units of time for plotting
outdata        = OPanalysis( cfg, indata );

%% compute permutation entropy for several channels
load( 'tonicClonic.mat' );
indata( 2, : )     = rand( 1, length( indata ) );  
cfg                = [];
cfg.method         = 'PE'; % compute robust permutation entropy
cfg.order          = 3;      % ordinal patterns of order 3 (4-points ordinal patterns)
cfg.delay          = 1;      % delay 1 between points in ordinal patterns
cfg.windowSize     = 1024;   % window size = 1024 time steps
cfg.time           = 0:1/102.4:179.999; % OPTIONAL time axis for plotting
cfg.units          = 'seconds';         % OPTIONAL units of time for plotting
outdata            = OPanalysis( cfg, indata );

%% compute permutation entropy and conditional entropy of ordinal patterns 
% for different parameters of logistic map (we use low-level functions for the example)
orbitLength = 10^4;
% take different r values 
order       = 7;    % for ordinal pattens of order 7 (8-points ordinal patterns)
delay       = 1;    % for delay 1 (successive points in ordinal patterns)
windowSize  = orbitLength - order*delay;
r           = 3.5:5*10^(-4):4; 
peValues    = zeros( 1, length( r ) );
ceValues    = zeros( 1, length( r ) );
leValues    = LEofLogisticMap( 3.5, 4, 5*10^(-4) );
indata      = zeros( 1, orbitLength );
for i = 1:length( r )
  if ( rem( i, 10 ) == 0 )
    disp( [ 'Calculating entropies for r = ' num2str( r( i ) ) ' from 4' ] );
  end
  indata( 1, 1 ) = rand( 1, 1 );
  for j = 2:orbitLength
    indata( j ) = r( i )*indata( j - 1 )*( 1 - indata( j - 1 ) );
  end
  peValues( i ) = PE( indata, delay, order, windowSize );
  ceValues( i ) = CondEn( indata, delay, order, windowSize - delay );
end
figure;
linewidth  = 0.5;
markerSize = 2;
plot( r, leValues, 'k',  'LineWidth',  linewidth ); grid on; hold on;
plot( r, peValues, 'go', 'markerSize', markerSize ); grid on; hold on;
plot( r, ceValues, 'bo', 'markerSize', markerSize ); grid on; hold on;
legend( 'LE', 'PE', 'CE' );
xlabel( 'Values of parameter r for logistic map x(t)=r*x(t-1)*(1-x(t-1))' ); 

%% INEFFICIENT METHOD: compute permutation entropy in sliding windows with an old method
% just for comparison in terms of speed with fast (PE.m) method
load( 'tonicClonic.mat' );
cfg            = [];
cfg.method     = 'oldPE'; % compute permutation entropy
cfg.order      = 6;       % ordinal pattens of order 6 (7-points ordinal patterns)
cfg.delay      = 1;       % delay 1 between points in ordinal patterns (successive points)
cfg.windowSize = 512;     % window size = 512 time steps
cfg.time       = 0:1/102.4:179.999; % OPTIONAL time axis for plotting
cfg.units      = 'seconds';         % OPTIONAL units of time for plotting
outdata        = OPanalysis( cfg, indata );