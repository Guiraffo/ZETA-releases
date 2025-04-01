function g_Jacob = Model1_g_Jacob(x,v)


% Output equation (Jacobian)


g_Jacob = [      1,  -(1/2)*cos(x(2)/2), 1, 0;
                     -x(2),  - x(1) + 1, 0, 1]; 
 

end

