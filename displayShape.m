function p1 = displayShape(shape, color, alpha, varargin)
% DISPLAYSHAPE Create a 3D figure of a polyhedron
%
%function p1 = displayShape(fracData, color, alpha, varargin)
%
%Displays a shape defined by a polyhedron in 3D with lighting.
% Input args:
%   shape: Either a polyhedron structure with two fields: vertices and
%             faces (see matlab's isosurface function for details) or just a N by 3
%             matrix of points. If the second type, this function will look for an
%             array called "faces" in the matlab workspace for the definition of the
%             polyhedron faces.
%  color (optional):     Color, same type as matlab uses in e.g. plot. For example
%                        'blue', '#ffccbb' or [0.2 0.9 0.1].
%  alpha (optional):      Alpha, or opaqueness value between 0 and 1.
% other parameters (optional): All parameters and name/value pairs that can
% be used in matlab's patch function can be used here. 
%
% Output args:
%  pl:  The patch object, can be used to manipulate the plot later.
%
% See also PATCH, ISOSURFACE, PLOT.
    

if isstruct(shape)
    faces = shape.faces;
    shape = shape.vertices;
else
    faces = evalin('base', 'faces');
end
if size(shape, 2) == 1
    shape = vec(shape);
end

if nargin == 2 && isa(color, 'matlab.graphics.primitive.Patch')
    p1 = color;
    p1.Vertices = shape;
    return;
end

if nargin < 2 || isempty(color)
    color = 'blue';
end
if nargin < 3 || isempty(alpha)
    alpha = 0.6;
end

ax = gca;
if strcmp(ax.NextPlot, 'replace')
    cla(ax);
end
if isempty(varargin)
    p1 = patch('Faces', faces, 'Vertices', shape,'EdgeColor','none', 'FaceColor', color);
else
    p1 = patch('Faces', faces, 'Vertices', shape,'EdgeColor','none', 'FaceColor', 'interp', varargin{:});
end

p1.FaceAlpha = alpha;

view(3)
%daspect([1,1,1])
axis equal
lighting gouraud
%axis tight
if ~any(cellfun(@(x)eq(x, "matlab.graphics.primitive.Light"), arrayfun(@class, gca().Children, 'UniformOutput', false)))
    camlight
    camlight(-100,-10)
    %camlight(80,-10)
end