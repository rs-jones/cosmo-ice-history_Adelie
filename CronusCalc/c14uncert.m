%
% uncert=c14uncert(c)
%
% Given a 14-C concentration in atoms/gram, computes the
% uncertainty based on interpolation of the standard deviations of
% the intercomparison sample A uncertainties.
%
% Inputs:
%       c            concentration (atoms/gram)
%
% Output
%       uncert       uncertainty (1-sigma) atoms/gram)
%
function uncert=c14uncert(c)
% From summary in Phillips et al. (2015) "The CRONUS-Earth Project: A Synthesis"
uncert=( 100 * 0.44e5 / 6.93e5 ) * c;
