out_csv = '10 PSI FY 2IA.csv';

% Column names: c_n (highest order) ... c_0 (constant), matching polyval order
n = poly_order;
coeff_names = arrayfun(@(k) sprintf('c%d', n-k+1), 1:(n+1), 'UniformOutput', false);
% coeff_names{1} corresponds to highest-power coefficient, last is constant term


polyTable = table(B_poly(:), C_poly(:), D_poly(:), E_poly(:), ...
    'VariableNames', {'B', 'C', 'D', 'E'});



writetable(polyTable, out_csv);

opts = detectImportOptions(out_csv);
T = readtable(out_csv, opts);
polyval(T.B,12)