function polys = matrix2polynomials(matrix, mons)
% polys = matrix2polynomials(matrix, mons) creates a cell array of multipol
% polynomials from the coefficients in matrix using the monomials in the
% matrix mons.

polys = {};

for k = 1 : size(matrix, 1)
    polys{end + 1} = multipol(matrix(k, :), mons);
end
