function stars = get_stars(p)
if p<0.0001
    stars = '****';
elseif p<0.001
    stars = '***';
elseif p<0.01
    stars = '**';
elseif p<0.05
    stars = '*';
else
    stars = 'ns';
end
end