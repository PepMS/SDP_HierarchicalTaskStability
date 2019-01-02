function G = task1_gainsLB(task1_lb)
G = cell(1, 5);

G{1} = diag([task1_lb task1_lb]);
G{2} = [1 0;0 0];
G{3} = [0 0;0 1];

end

