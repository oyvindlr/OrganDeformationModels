%Vectorization/de-vectorization of a N by 3 matrix
% Returns a vector if the input is a matrix, and a matrix if the input is a
% vector.

function y=vec(x)

if size(x, 2) == 3
    y = reshape(x, 3* size(x, 1), 1);
else
    y = reshape(x, size(x, 1)/3, 3);
end