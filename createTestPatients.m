function [patientData, polyhedronFaces] = createTestPatients()
%%
%define the grid
xx = -10:0.1:10;
yy = xx;
zz = xx;

[X, Y, Z] = meshgrid(xx, yy, zz);

%%
patientData = struct;

%First patient has an organ that looks like a ball of diameter 3
%This generates its binary map
P1_1 = sqrt(X.^2+Y.^2+Z.^2)< 3;

%Generate polyhedron from binary map
P1_1 = smooth3(P1_1, 'box', 5); 
poly = isosurface(X, Y, Z, P1_1, 0.5);

%First shape for this patient is just the ball
patientData(1).contourPoints{1} = poly.vertices;

%Faces are the same for all shapes 
polyhedronFaces = poly.faces;


% Second shape for that patient is a slightly larger ball with a bump
% inwards
vertices_s2 = poly.vertices * 1.1;
mx = -2.5;
ys = vertices_s2((poly.vertices(:, 2) < mx), 2);
ys = ys - 2*(ys-mx*1.1);
vertices_s2((poly.vertices(:, 2) < mx), 2) = ys;
patientData(1).contourPoints{2} = vertices_s2;

% Third shape is a slightly smaller ball same but with the bump outwards
vertices_s3 = poly.vertices *0.9;
ys = vertices_s3((poly.vertices(:, 2) < mx), 2);
ys = ys + 2*(ys-mx*0.9);
vertices_s3((poly.vertices(:, 2) < mx), 2) = ys;
patientData(1).contourPoints{3} = vertices_s3;

% second patient is a bigger version of the first patient
patientData(2).contourPoints{1} = patientData(1).contourPoints{1}*1.5;
patientData(2).contourPoints{2} = patientData(1).contourPoints{2}*1.5;
patientData(2).contourPoints{3} = patientData(1).contourPoints{3}*1.5;

% Third patient has the bump outward always, and rotates a bit
patientData(3).contourPoints{1} = patientData(1).contourPoints{3};
angle = 10/360*2*pi; %10 degrees
rotmat = [1 0 0; 0 cos(angle) -sin(angle); 0 sin(angle) cos(angle)];
patientData(3).contourPoints{2} = patientData(3).contourPoints{1} * rotmat;
rotmat = [1 0 0; 0 cos(-angle) -sin(-angle); 0 sin(-angle) cos(-angle)];
patientData(3).contourPoints{3} = patientData(3).contourPoints{1} * rotmat;


