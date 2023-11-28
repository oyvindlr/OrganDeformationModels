function [x, fv] = quasinewton(f, xstart, options)
%QUASINEWTON Minimize the function f using a Quasi-Newton method with a
%BFGS approximation to the Hessian. See the arguments block for optional
%arguments. 
% Outputs:
%  x: Point at which f is minimized
%  fv: Function value at x
arguments
    %Function to minimize
    f function_handle
    %Starting point
    xstart double
    %Function change tolerance.
    options.ftol (1, 1) double = 1e-6
    %Gradient norm tolerance
    options.gtol (1, 1) double = 1e-4
    %Max number of iterations
    options.maxiter (1,1) double = 200
    %If verbose is true, detailed information is printed every iteration
    options.verbose logical = false
    %Delta value for the finite difference derivative approximation
    options.finitediffdelta double = 1e-6
    %Max number of iterations for line search
    options.linesearchmaxiter double = 30
    %At the first iteration, the Hessian is computed based on a finite
    %difference approximation. If computefullhessian is false, only the
    %diagonal elements are computed. 
    options.computefullhessian logical = false    
end
%Turn singular matrix-warning into an error so we can catch it
w = warning('error', 'MATLAB:singularMatrix'); 
w2 = warning('error', 'MATLAB:nearlySingularMatrix');

fv_prev = inf;
x = xstart;
i = 0;
fv = f(x);
neval = 1;
reset = true;
g = inf;
while abs(fv-fv_prev) > options.ftol && norm(g) > options.gtol && i < options.maxiter
    fv_prev = fv;
    if reset
        %[H, g, ne] = hessian(f, x, fv, options.finitediffdelta, options.computefullhessian);
        H = eye(length(x));
        [g, ne] = gradient(f, x, fv, options.finitediffdelta);
        reset = false;
    else
        [H, g, ne] = quasi_hessian(f, x, g, dx, H, fv, options.finitediffdelta);        
    end
    neval = neval + ne;
    try
        ddx = -H\g;
    catch
        mydisp('Hessian is singular to working precision, adding regularization');
        %H = H + eye(length(x));
        %ddx = -H\g;
        H = eye(length(x));
        ddx = -g;
    end
    if ddx'*g > 0 
        fprintf('At iteration %d: Direction is not a descent direction. Resorting to steepest descent.\n', i + 1);
        ddx = 0.8*-norm(ddx)*g/norm(g) + 0.2*ddx;                
        reset = true;
    end
    [step, fv, ne] = linesearch(f, x, fv, g, ddx, 1, 20, 'maxiter',options.linesearchmaxiter, 'finitediffdelta', options.finitediffdelta, 'verbose', options.verbose);      
    
    neval = neval + ne;
    dx = step*ddx;
    x = x + dx;
    i = i+1;
    disp(['Iteration ' num2str(i) ': function value=' num2str(fv) ', improvement=' num2str(fv_prev-fv) ', x(1): ' num2str(x(1)) ', function evaluations: ' num2str(neval)]);
end
if i >= options.maxiter 
    disp('Optimization terminated: Maximum number of iteration reached');
elseif abs(fv-fv_prev) <= options.ftol
    disp('Optimization terminated: Function value change less than tolerance');
else
    disp('Optimization terminated: Norm of gradient less than tolerance');
end
%Turn error back to warning
warning(w);
warning(w2);


    function mydisp(str)
        if options.verbose
            disp(str)
        end
    end

end


function [H, g, neval] = hessian(f, x, fv, d, fullhessian)
  if nargin < 4
      d = 1e-5;
  end
  n = length(x);
  g = zeros(size(x));
  H = zeros(n);
  neval = 0;
  for i = 1:n
      xi = x;
      xi(i) = x(i) + d;
      fvv = f(xi);      
      g(i) = (fvv-fv)/d;
      
      xi(i) = xi(i) + d;
      fvvv = f(xi);
      H(i, i) = ((fvvv-fvv)/d - g(i))/d;
      neval = neval + 2;
      if fullhessian
          for j = (i+1):n
              xj = x;
              xj(j) = xj(j) + d;
              fvj = f(xj);
              xj(i) = xj(i) + d;
              fvvj = f(xj);
              H(i, j) = ((fvvj-fvj)/d - g(i))/d;
              H(j, i) = H(i, j);
              neval = neval + 2;
          end
      end    
  end
end

function [g, neval] = gradient(f, x, fv, d)
  n = length(x);
  g = zeros(size(x));
  for i = 1:n
      xi = x;
      xi(i) = x(i) + d;
      fvv = f(xi);
      g(i) = (fvv-fv)/d;      
  end  
  neval = n;
end

function [H, g, neval] = quasi_hessian(f, x, g_prev, dx, H, fv, d)
  if nargin < 7
      d = 1e-5;
  end
  n = length(x);
  g = zeros(size(x));
  for i = 1:n
      xi = x;
      xi(i) = x(i) + d;
      fvv = f(xi);
      g(i) = (fvv-fv)/d;      
  end  
  neval = n;
  y = g - g_prev;
  Hb = H;
  H = H + (y*y')./(y'*dx);
  dxh = Hb*dx;
  H = H-(dxh*dxh')./(dx'*Hb*dx);  
end
