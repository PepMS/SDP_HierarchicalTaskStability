function [A] = minors(P)
    disp('Correctness 12')
    norm(P(2)) < sqrt(P(1)*P(5))
    
    disp('Correctness 13')
    norm(P(3)) < sqrt(P(1)*P(9))
    
    disp('Correctness 23')
    norm(P(6)) < sqrt(P(5)*P(9))
    

end
