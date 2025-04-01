function y = Model1_g(x,v)

% Output equation

if(isempty(v))
    v = zeros(2,1);
end

y = [     x(1) - sin(x(2)/2) + v(1);
           (- x(1) + 1)*x(2) + v(2)];

          

end

