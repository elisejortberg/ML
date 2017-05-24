function score_out = test(model, data)
    %use predict function after model
    [~, score] = predict(model, data);
    score_out = score(:,2); %assigned true class as 1 (2nd column)
end
