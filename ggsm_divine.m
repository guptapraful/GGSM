function q = ggsm_divine(im)

% ========================================================================
% 
% -----------COPYRIGHT NOTICE STARTS WITH THIS LINE------------
% Copyright (c) 2017 The University of Texas at Austin
% All rights reserved.
% 
% Permission is hereby granted, without written agreement and without license or royalty fees, to use, copy, 
% modify, and distribute this code (the source files) and its documentation for
% any purpose, provided that the copyright notice in its entirety appear in all copies of this code, and the 
% original source of this code, Laboratory for Image and Video Engineering (LIVE, http://live.ece.utexas.edu)
% and Center for Perceptual Systems (CPS, http://www.cps.utexas.edu) at the University of Texas at Austin (UT Austin, 
% http://www.utexas.edu), is acknowledged in any publication that reports research using this code. The research
% is to be cited in the bibliography as:
% 
% 1. P. Gupta, A. K. Moorthy, R. Soundararajan and A. C. Bovik, "DIIVINE-GGSM Software Release", 
% URL: http://live.ece.utexas.edu/research/quality/diivine-ggsm.zip, 2017
% 
% IN NO EVENT SHALL THE UNIVERSITY OF TEXAS AT AUSTIN BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, 
% OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OF THIS DATABASE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF TEXAS
% AT AUSTIN HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% THE UNIVERSITY OF TEXAS AT AUSTIN SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE DATABASE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS,
% AND THE UNIVERSITY OF TEXAS AT AUSTIN HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
% 
% -----------COPYRIGHT NOTICE ENDS WITH THIS LINE------------%
% 
% Author  : Praful Gupta
% Version : 1.1
% 
% The authors are with the Laboratory for Image and Video Engineering
% (LIVE), Department of Electrical and Computer Engineering, The
% University of Texas at Austin, Austin, TX.
% 
% Kindly report any suggestions or corrections to praful_gupta@utexas.edu
% 
% ========================================================================


addpath(genpath('./matlabPyrTools/'))
addpath(genpath('./utils/'))
% Feature extraction
f = gdivine_feature_extract(im);
% DIIVINE-GGSM quality score
q = gdivine_overall_quality(f);