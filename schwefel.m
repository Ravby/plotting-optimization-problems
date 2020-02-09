function [out] = schwefel(x)

d = length(x);
sum = 0;
for ii = 1:d
	xi = x(ii);
	sum = sum + xi*sin(sqrt(abs(xi)));
end

out = 418.9829*d - sum;

end