% FUNCTION:     CalcPosition.m
% DESCRIPTION:  Calculates the unconstrained optimal position of probes
%               based on the measured frequencies.
% INPUTS:       freqs:   A 3xn matrix where n is the number of probes and
%                        the rows are the initial x-, y-, z-positions.
%               grads:   A 3xn matrix where each row is the coeffs of
%                        the sph harm decomposition of the x-, y- and
%                        z-gradients respectively.
% OUTPUTS:      opt_pos: A 3xn matrix with the optimum positions.
% DEPENDENCIES: spha.mex64w (<project dir>/common)

function opt_pos = CalcPosition( freqs, grads )

% Check input arguments
if nargin < 1
    error( 'Requires ''freqs'' input variable.' );
else
    if size( freqs, 1 ) == 3, freqs = freqs'; end
    if nargin > 1
        if size( grads, 1 ) == 3, grads = grads'; end
    else
        grads = zeros( 3,4 );
        grads(3,2) = 1; grads(1,3) = 1; grads(2,4) = 1;
    end
end

%% Parameters

gam     = 2.67513e8;
deltaTE = 0.025e-3;

options             = optimset( 'Display', 'off' );
options.TolFun      = 1e-16;
options.TolX        = 1e-8;
options.MaxFunEvals = 3200;

%% Calculate

[npos, nfields] = size(freqs);
nterms = size(grads,1);

opt_pos = zeros( 3,npos );
for k = 1 : npos
    if nfields > 3, x = freqs( k,2:4 );
    else x = freqs( k,: ); end
    optfun = @(pos) sum(sum(( spha(1:nterms,pos)*grads(:,1:nfields) - freqs(k,:) ).^2));
    [x, ~] = fminsearch( optfun, x, options );
    opt_pos( :,k ) = x';
end

end