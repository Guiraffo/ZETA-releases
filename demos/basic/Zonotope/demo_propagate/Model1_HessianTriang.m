function HessianTriang = ChapterModel1_HessianTriang(x)


% System triangular Hessian


x1 = x(1);
x2 = x(2);

HessianTriang = cell(size(x,1),1);



HessianTriang{1} = [(1/2)*((32*x2)/(x1 + 4)^3 - 2/7),       -16/(x1 + 4)^2; 
                                                  0,                    0];
                                 
HessianTriang{2} = [(1/2)*(-(24*x2)/(x1 + 4)^3),       12/(x1 + 4)^2; 
                                             0,                   0];      
                                              
                                            
end

