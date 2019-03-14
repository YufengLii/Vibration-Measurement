function G2 = steerGaussFilterG2(theta,sigma)
% Process input arguments (according to usage mode).
if ~exist('arg1','var') || ~isstruct(arg1)
   
    % Assign default filter orientation (if not provided).
    if ~exist('theta','var') || isempty(theta)
      theta = 0;
    end
    theta = -theta*(pi/180);
    % Assign default standard deviation (if not provided).
    if ~exist('sigma','var') || isempty(sigma)
       sigma = 1;
    end        
end % End of input pre-processing.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part II: Evaluate separable filter kernels.

% Determine necessary filter support (for Gaussian).
Wx = floor((8/2)*sigma); 
if Wx < 1
  Wx = 1;
end
x = [-Wx:Wx];
[xx,yy] = meshgrid(x,x);
g0 = exp(-(xx.^2+yy.^2)/(2*sigma^2))/(sigma*sqrt(2*pi));
G2a = -g0/sigma^2+g0.*xx.^2/sigma^4;
G2b =  g0.*xx.*yy/sigma^4;
G2c = -g0/sigma^2+g0.*yy.^2/sigma^4;
G2 = (cos(theta))^2*G2a+sin(theta)^2*G2c-2*cos(theta)*sin(theta)*G2b;

