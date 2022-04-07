function y = makepowerof2 (x)
    n1 = length(x);
	base2 = fix(log(n1)/log(2) + 0.4999);   % power of 2 nearest to N
	y = [x,zeros(1,2^(base2+1)-n1)];
end

