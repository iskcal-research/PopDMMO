function niching = fake_niching(fit)
    [~, idx] = sort(fit, 'descend');
    niching = {idx};
end