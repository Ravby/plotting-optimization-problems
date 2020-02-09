function [out] = sphere(x)

d = length(x);
sum = 0;
for ii = 1:d
	xi = x(ii);
	sum = sum + xi^2;
end

out = sum;

end