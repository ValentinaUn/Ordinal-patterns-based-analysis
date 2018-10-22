# Ordinal-patterns-based-analysis
The OPA (ordinal-patterns-analysis) toolbox is intended for nonlinear analysis 
of multivariate time series with becoming more and more popular 
ordinal-patterns-based measures [AKU15,BKP02,KM05,KUU14,ZZR12] which are efficiently  
computed [UK13,U15] and visualized:

 - permutation entropy (cfg.method = 'PE') [BKP02]
 - permutation entropy for ordinal patterns with tied ranks (cfg.method = 'eqPE') [KUU14]
 - permutation entropy and ordinal patterns distributions (cfg.method = 'opdPE') [KM05]
 - conditional entropy of ordinal patterns (cfg.method = cePE') [UK13]
 - robust permutation entropy (cfg.method = 'rePE') [KUU14,U15]
 - new ordinal-patterns-based measures are to be added (suggestions and feedback are welcome!)

The interface of OPA toolbox is provided by a function OPanalysis( cfg, indata ), where 
 - cfg is a configuration structure with method's parameters;
 - indata is data to be analyzed.

Example of use (see more examples for different methods and parameters 
in examples.m and OPanalysis.m help):

cfg            = [];
cfg.method     = 'opdPE'; % try also 'PE', 'CE', 'rePE' and 'all' here
cfg.order      = 3;       % for ordinal pattens of order 3 (4-points ordinal patterns)
cfg.delay      = 1;       % for delay 1 between points in ordinal patterns (successive points)
cfg.windowSize = 512;     % for window size 512 in time points
indata         = rand( 1, 7777 );
for i = 4000:7000         % change of data complexity
  indata( i ) = 4*indata( i - 1 )*( 1 - indata( i - 1 ) );
end 
outdata        = OPanalysis( cfg, indata );

Values of ordinal-patterns-based measures are computed in maximally overlapping 
sliding windows in a fast way [UK13,U15].

examples.m script contains examples of using OPA toolbox for different 
ordinal-patterns-based measures with different parameters
The Code folder contains functions for computing the aforementioned ordinal-patterns-based measures.
The Data folder contains example epileptic EEG datasets from https://vis.caltech.edu/~rodri/data.htm 
(dataset 5, 130-2_c4.asc) used for illustrating in examples.m.
The Tables folder contain precomputed tables as *.mat-files for efficient computing of 
ordinal-patterns-based methods.

NOTE: this is a beta-version of OPA-toolbox, new ordinal-patterns-based measures 
are to be added (your feedback and comments are welcome)!

REFERENCES (alphabetical order):
[AKU15] Amigo, J.M., Keller, K. and Unakafova, V.A., 2015. On entropy, 
entropy-like quantities, and applications. Discrete & Continuous Dynamical 
Systems-Series B, 20(10).
[AZS08] Amigo, J.M., Zambrano, S. and Sanjuan, M.A., 2008. 
Combinatorial detection of determinism in noisy time series. 
EPL (Europhysics Letters), 83(6), p.60005.
[BKP02] Bandt C., Pompe B., Permutation entropy: a natural complexity 
measure for time series. Physical review letters, 2002, APS
[BQM2012] Bian, C., Qin, C., Ma, Q.D. and Shen, Q., 2012. Modified 
permutation-entropy analysis of heartbeat dynamics. Physical Review E, 85(2), p.021906.
[CTG04] Cao, Y., Tung, W.W., Gao, J.B. et al., 2004. Detecting dynamical changes 
in time series using the permutation entropy. Physical Review E, 70(4), p.046217.
[KM05] Keller, K., and M. Sinn. Ordinal analysis of time series. 
Physica A: Statistical Mechanics and its Applications 356.1 (2005): 114--120
[KUU14] Keller, K., Unakafov, A.M. and Unakafova, V.A., 2014. 
Ordinal patterns, entropy, and EEG. Entropy, 16(12), pp.6212-6239.
[RMW13] Riedl, M., Muller, A. and Wessel, N., 2013. Practical considerations 
of permutation entropy. The European Physical Journal Special Topics, 222(2), pp.249-262.
[UK13] Unakafova, V.A., Keller, K., 2013. Efficiently measuring 
complexity on the basis of real-world Data. Entropy, 15(10), 4392-4415.
[U15] Unakafova, V.A., 2015. Investigating measures of complexity 
for dynamical systems and for time series (Doctoral dissertation, 
University of Luebeck).
[ZZR12] Zanin, M., Zunino, L., Rosso, O.A. and Papo, D., 2012. 
Permutation entropy and its main biomedical and econophysics applications: a review. 
Entropy, 14(8), pp.1553-1577.
 
@author Valentina Unakafova
@email UnakafovaValentina(at)gmail.com
@date 18.07.2017
