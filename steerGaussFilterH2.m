function H2 = steerGaussFilterH2(theta,sigma)


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

H0 = 0.9780*exp(-(xx.^2+yy.^2));

H2a = H0.*(-2.254.*xx + xx.^3);
H2d = H0.*(-2.254.*yy + xx.^3);
H2b = H0.*(-0.7515 + xx.^2).*yy;
H2c = H0.*(-0.7515 + yy.^2).*xx;

H2 = (cos(theta))^3*H2a - 3*(cos(theta))^2*(sin(theta))*H2b + 3*(cos(theta))*(sin(theta))^2*H2c -(sin(theta))^3*H2d;



