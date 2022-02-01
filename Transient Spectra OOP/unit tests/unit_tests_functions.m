%Unit tests for functions
%% findGridIntercepts.m
x = [0:0.01:2.99 3.000001];
y = sin(2*pi*x);
xInt = findGridIntercepts(x,y,0);
soln = [0, 0.5, 1, 1.5, 2, 2.5, 3];
assert(all((xInt(:)-soln(:)) < 1e-6), 'All intercepts were not found for sin');

y = [zeros(1,length(x)-1) 1];
xInt = findGridIntercepts(x,y,0);