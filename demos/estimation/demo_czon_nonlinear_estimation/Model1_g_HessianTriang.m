function g_HessianTriang = Model1_g_HessianTriang(x,v)

% Output equation (triangular Hessian)
                  
g_HessianTriang = cell(2,1);

g_HessianTriang{1} = [      0,                   0,   0,   0;
                            0,   (1/4)*sin(x(2)/2),   0,   0;
                            0,                   0,   0,   0;
                            0,                   0,   0,   0];
                  
g_HessianTriang{2} = [   0,  -1,   0,   0;
                        -1,   0,   0,   0;
                         0,   0,   0,   0;
                         0,   0,   0,   0];                
   

end

