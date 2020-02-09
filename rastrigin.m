function [out] = rastrigin(x)

d = length(x);
sum = 0;
for ii = 1:d
	xi = x(ii);
	sum = sum + (xi^2 - 10*cos(2*pi*xi));
end

out = 10*d + sum;

end