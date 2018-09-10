function [x, res] = fitFunctionT1IR_3p(data,TI)
%   [x,res] = fitFunctionT1IR_3p(data,TI)
%
%   Fit a 3-parameter model to inversion recovery data
%
%   Chong Duan, Oct 2017

opts = optimset('lsqcurvefit');
opts = optimset(opts,'Display','off');

sz = size(data);
x = zeros([3 sz(2:end)]);
res = zeros([sz(2:end) 1]);
% errorMessage = cell(0);

fun = @(x,xd) x(1) - x(2)*exp(-xd(:)/x(3)); 

parfor i = 1:prod(sz(2:end))
    [x(:,i), res(i)] = lsqcurvefit(@(x,xd) fun(x,xd),[max(data(:,i)) 2*max(data(:,i)) 1000],TI,data(:,i),[0 -Inf 0],[Inf Inf Inf],opts);
end

end