function [step,fvnew, neval] = linesearch(f, x, fv, g, dir, startstep, maxstep, options)
%LINESEARCH Perform a line search algorithm to identify
% a point satisfying the strong Wolfe-conditions.
% The algorithm is adapted from Nocedal's "numerical optimization".

% Input args:
%  f: Objective function
%  x: Current point 
%  fv: Function value at current point
%  g: Gradient at current point
%  dir: Direction to search (vector, same size as x)
%  startstep: Initial value of step size
%  maxstep: Maximum value of the step size
% 
% Output args:
%  step: Step size satisfying the wolfe conditions
%  fvnew: Function value at new point
%  neval: Number of function evaluations performed

arguments
    f function_handle
    x double    
    fv (1, 1) double
    g double
    dir double
    startstep double
    maxstep double
    options.maxiter (1,1) double = 30
    options.verbose logical = false
    options.finitediffdelta double = 1e-6
    options.tol_derivative double = 1e-5
end


psi = @(a)f(x+a*dir);

c1 = 1e-4;
c2 = 0.9;
increasefactor = 2; %Golden ratio
fvd0 = g'*dir;
if fvd0 > 0
    mydisp('Linesearch: Chosen direction is not a descent direction. Returning starting step size');
    step = startstep;
    fvnew = psi(step);
    neval = 1;
    return;
end
neval = 0;

i = 1;
step = startstep;
smallerstep = 0;
smallerstepfv = fv;
smallerstepd = fvd0;
smallesttoobigstep = inf;

while true       
    i = i + 1;
    %If step greater than maximum, or too many iterations, exit
    if step > maxstep
        step = maxstep;
        fvnew = psi(step);
        neval = neval + 1;
        mydisp('Linesearch: Maximum step size for line search reached');
        return;
    end
    if i > options.maxiter
        mydisp('Linesearch: Maximum number of iterations reached for line search');
        return;
    end
    
    fvnew = psi(step);
    neval = neval + 1;
    
    %If step too large, reduce step and continue
    if fvnew > fv + c1*step*fvd0
        prevstep = step;
        if step < smallesttoobigstep
            smallesttoobigstep = step;
        end
        step = quadraticmin(smallerstepfv, smallerstepd, fvnew, step-smallerstep) + smallerstep;
        mydisp(append('Linesearch: Step size ', num2str(prevstep),  ' too large. New step size: ', num2str(step)));
        continue;
    end
    %Find derivative
    fvdnew = finitediff(psi, fvnew, step, options.finitediffdelta);
    neval = neval + 1;
    %If step satisfies curvature conditions, we are done
    if abs(fvdnew) <= -c2*fvd0
        mydisp(append('Linesearch: Step size ', num2str(step), ' satisfies the strong Wolfe conditions. Line search complete.'));
        return;
    elseif abs(fvdnew) < options.tol_derivative
        mydisp('Linesearch terminated: Magnitude of derivative less than tolerance.');
        return;
    end
    mydisp(append('Linesearch: Step size ', num2str(step), ' does not satisfy curvature conditions'));
    
    %If derivative is positive, we've moved to far. Interpolate to reduce
    %step and continue
    if fvdnew > 0
        if step < smallesttoobigstep
            smallesttoobigstep = step;
        end   
        step = quadraticmin(fvnew, fvdnew, smallerstepfv, smallerstep-step) + step;
        mydisp(append('Linesearch: Step too large. New step: ', num2str(step)));
        continue;
    end    
    %Otherwise, we have not moved far enough. No guarantee that polynomial
    %interpolation will work in this case, so increase by a constant
    %factor, but not beyond any value that has previously been shown to be
    %too large. 
    smallerstep = step;
    smallerstepfv = fvnew;
    smallerstepd = fvdnew;    
    step =  increasefactor*step;
    if step > smallesttoobigstep
        step = smallerstep + 0.5*(smallesttoobigstep - smallerstep);
    end
    mydisp(append('Linesearch: Step too small. New step: ', num2str(step)));           
end

    function mydisp(str)
        if options.verbose
            disp(str);
        end
    end

end

%Finite difference derivative of univariate function f
function d = finitediff(f,fv, a, delta)
    fvv = f(a+delta);
    d = (fvv-fv)/delta;    
end

function r = quadraticmin(f0, f0d, fa, a)
r = -a^2*f0d;
r = r/(2*(fa-f0-f0d*a));
end
