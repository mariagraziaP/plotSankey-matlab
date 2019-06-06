function [] = plotSankeyDiagram(class1, class2, colors, opts)
% plotSankeyDiagram
% 
%
%
%

% make sure class1 and class2 are colums
if isrow(class1); class1 = class1'; end
if isrow(class2); class2 = class2'; end

cl1 = unique(class1);
cl1flip = flipud(unique(class1));

cl2 = unique(class2);
cl2flip = flipud(unique(class2));

if nargin<3
    colors = parula(length(unique([cl1, cl2])));
end

if nargin <4
    % set parameters
    opts.dist = 0.1; % half of the distance between blocks 
    opts.width = 1;  % width of the blocks
    opts.unit = 1;   % scale factor of the nodes
end

figure;
xlim([0 11])
ylim([1 opts.unit*length(class1)+1])
bly(1) = 1 + opts.dist;
for i=1:length(cl1)
    hl(i) = opts.unit*length(find(class1==cl1flip(i)));
    rl{i} = rectangle('Position', [1,bly(i),opts.width,hl(i)-2*opts.dist],...
        'EdgeColor', 'none',...
        'FaceColor', colors(cl1flip(i),:));
    bly(i+1) = bly(i)+hl(i); 
    lenL(i) = bly(i+1)-bly(i)-2*opts.dist;
end
clear i
% old
% lenL = fliplr(lenL);
% hl = fliplr(hl);
% end old
% new
lenLt = fliplr(lenL);
lenL = zeros(1,max(cl1));
lenL(cl1) = lenLt;
hlt = fliplr(hl);
hl = zeros(1,max(cl1));
hl(cl1) = hlt;
% end new

hold on
bry(1) = 1 + opts.dist;
for i=1:length(cl2)
    hr(i) = opts.unit*length(find(class2==cl2flip(i)));
    rr{i} = rectangle('Position', [9,bry(i),opts.width,hr(i)-2*opts.dist],...
        'EdgeColor', 'none',...
        'FaceColor', colors(cl2flip(i),:));
    bry(i+1) = bry(i)+hr(i); 
    lenR(i) = bry(i+1)-bry(i)-2*opts.dist;
end
clear i

lenRt = fliplr(lenR);
lenR = zeros(1,max(cl2));
lenR(cl2) = lenRt;
hrt = fliplr(hr);
hr = zeros(1,max(cl2));
hr(cl2) = hrt;
grid on

% new
plt = bly(1:end-1)';
plt(1:end, 2) = 0;
plt = flipud(plt);
pl = zeros(max(cl1),2);
pl(cl1,:) = plt;
prt = bry(1:end-1)';
prt(1:end, 2) = 0;
prt = flipud(prt);
pr = zeros(max(cl2),2);
pr(cl2,:) = prt;
% end new

x = 2:0.1:9;
xf = fliplr(x);
y = sigmf(x,[2 5]); y(1) = 0; y(end) = 1;

for i=1:length(cl1)
    
    L = cl1flip(i);
    posL = find(class1==cl1flip(i));
    % posR = find(class2==cl1flip(i));
    
    newlab = class2(posL);
    newlabord = sort(unique(newlab), 'descend');
    
    for n=1:length(newlabord)
        R = newlabord(n);
        currlen = length(find(newlab==R));
        pl(L,2) = pl(L,1) + currlen*lenL(L)/hl(L);
        pr(R,2) = pr(R,1) + currlen*lenR(R)/hr(R);
        if pl(L,1)<pr(R,1) % sigmoide crescente
            y1 = (pr(R,1)-pl(L,1)) * y + pl(L,1);
        elseif pl(L,1)>pr(R,1) % sigmoide decrescente
            y1 = -(pl(L,1)-pr(R,1)) * y + pl(L,1);
        elseif pl(L,1)==pr(R,1) % retta
            y1 = pl(L,1)*ones(1, length(x));
        end
        if pl(L,2)<pr(R,2)
            y2 = (pr(R,2)-pl(L,2)) * y + pl(L,2);
            x2 = [2 x 9];
        elseif pl(L,2)>pr(R,2)
            y2 = -(pl(L,2)-pr(R,2)) * y + pl(L,2);
        elseif pl(L,2)==pr(R,2)
            y2 = pl(L,2)*ones(1, length(xf));
        end
        y2 = fliplr(y2);
        ar = fill([x xf], [y1 y2], [0.7 0.7 0.7]);
        ar.EdgeColor = 'none';
        ar.FaceAlpha = 0.5;
        
        pl(L,1) = pl(L,2);
        pl(L,2) = 0;
        pr(R,1) = pr(R,2);
        pr(R,2) = 0;
    end
    clear n
    
end

axis off

end

function y = sigmf(x, params)
%       Roger Jang, 6-29-93, 4-17-93.
%   Copyright (c) 1994-96 by The MathWorks, Inc.
%   $Revision: 1.10 $  $Date: 1996/04/02 23:16:02 $

if nargin ~= 2
    error('Two arguments are required by the sigmoidal MF.');
elseif length(params) < 2
    error('The sigmoidal MF needs at least two parameters.');
end

a = params(1); c = params(2);
y = 1./(1 + exp(-a*(x-c)));

end