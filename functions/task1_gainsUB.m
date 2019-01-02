function G = task1_gainsUB(task1_ub)
G = cell(1, 5);

G{1} = -diag([task1_ub task1_ub]);
G{2} = [-1 0;0 0];
G{3} = [0 0;0 -1];

end

