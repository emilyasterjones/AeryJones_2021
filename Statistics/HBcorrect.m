function stars = HBcorrect(p)

% given a double of NxG p-values, corrects over N comparisons
% returns stars for significance value if meets the threshold, ns otherwise

stars = strings(size(p,1),size(p,2));
for g = 1:size(p,2)
    
    % Holm-Bonferroni correction
    multcorrect = false(length(p),1);
    significant = zeros(length(p),1);
    [p_sort,ind] = sort(p);
    for k = 1:length(p_sort)
        multcorrect(k) = p_sort(k)<0.05/(length(p_sort)+1-k);
    end
    sigind = ind(multcorrect);
    significant(sigind) = 1;
    
    for i = 1:length(p)
        if significant(i)
            stars(i,g) = string(get_stars(p(i,g)));
        else
            stars(i,g) = "ns";
        end
    end
end

end