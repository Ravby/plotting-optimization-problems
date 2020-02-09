function [out] = easom(x)

x1 = x(1);
x2 = x(2);

fact1 = -cos(x1)*cos(x2);
fact2 = exp(-(x1-pi)^2-(x2-pi)^2);

out = fact1*fact2;

end
