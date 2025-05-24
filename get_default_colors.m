function clrs = get_default_colors(k)
if nargin < 1
    k = 8;
end
figure;
plot(rand(10, k));
h = gca;
clrs = flipud(cat(1, h.Children(:).Color));
close;
end