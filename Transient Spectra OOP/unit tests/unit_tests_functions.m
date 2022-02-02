%Unit tests for functions
%% findGridIntercepts.m
x = [0:0.01:2.99 3.000001];
y = sin(2*pi*x);
xInt = findGridIntercepts(x,y,0);
soln = [0, 0.5, 1, 1.5, 2, 2.5, 3];
assert(all((xInt(:)-soln(:)) < 1e-6), 'All intercepts were not found for sin');

y = [zeros(1,length(x)-1) 1];
xInt = findGridIntercepts(x,y,0);

%% nearestVal.m
x = 1:5;
[vals,ind] = nearestVal(x,1.1);
assert(vals==1 && ind==1,'failed to find nearest val for a vector with 1 target val');

[vals,ind] = nearestVal(x,[1.1, 2.1]);
assert(all(vals(:)==[1;2]) && all((ind)==[1;2]),...
    'failed to find nearest val for a vector with 2 target vals');

[vals,ind] = nearestVal(repmat(x(:),1,2),[1.1, 2.1]);
assert(all(vals==[1,1;2,2],'all') && all(ind==[1,1;2,2],'all'),...
    'failed to find nearest val for a matrix with 2 target vals');

[vals,ind] = nearestVal(repmat(x(:),1,2,2),1.2);
assert(all(size(vals)==[1,2,2]) && all(size(ind)==[1,2,2]),...
    'failed to return correct vals and ind size for a 3d array with 1 target val');
assert(all(vals(:)==ones(4,1)) && all(ind(:)==ones(4,1)),...
    'failed to find nearest val for a 3d array with 1 target val');

[vals,ind] = nearestVal(repmat(x(:),1,2,2),[1.2, 2.2]);
assert(all(size(vals)==[2,2,2]) && all(size(ind)==[2,2,2]),...
    'failed to return correct vals and ind size for a 3d array with 2 target val');
assert(all(vals(:)==[1,2,1,2,1,2,1,2]') && all(ind(:)==[1,2,1,2,1,2,1,2]'),...
    'failed to find nearest val for a 3d array with 2 target val');

[vals,ind] = nearestVal(repmat(x(:),1,2,2),[1.2, 1.3],'unique',false,'trim',false);
assert(all(size(vals)==[2,2,2]) && all(size(ind)==[2,2,2]),...
    'failed to return correct vals and ind size for a 3d array with 2 target val with unique set to false');
assert(all(vals(:)==[1,1,1,1,1,1,1,1]') && all(ind(:)==[1,1,1,1,1,1,1,1]'),...
    'failed to find nearest val for a 3d array with 2 target val with unique set to false');

[vals,ind] = nearestVal(repmat(x(:),1,2,2),[1.2, 1.3]);
assert(all(size(vals)==[1,2,2]) && all(size(ind)==[1,2,2]),...
    'failed to return correct vals and ind size for a 3d array with 2 non-unique target val');
assert(all(squeeze(vals(1,:,:))==[1,1;1,1],'all') && all(squeeze(vals(1,:,:))==[1,1;1,1],'all'),...
    'failed to find nearest val for a 3d array with 2 non-unique target val');