function H = hankel(w, t1, t2, r)
    H = zeros(r*t1, t2);
    for i = 1:t2
        H(:, i) = w((i-1)*2+1:(i+t1-1)*2, :);
    end
    
    if size(H, 1) == 1
        H = H';
    end
end