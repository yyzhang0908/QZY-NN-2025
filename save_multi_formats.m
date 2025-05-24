function fig = save_multi_formats(fig, p, types)
if any(ismember(types, 'pdf'))
    print(fig, p, '-dpdf', '-r300');
end
if any(ismember(types, 'jpg'))
    print(fig, p, '-djpeg', '-r300');
end
if any(ismember(types, 'tiff'))
    print(fig, p, '-dtiff', '-r300');
end
if any(ismember(types, 'fig'))
    saveas(fig, p, 'fig');
end
end