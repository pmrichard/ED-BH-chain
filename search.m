function index = search(value)

global  D
global  address

L = 1;
R = D;
M = floor((L+R)/2);

while (L <= R)
    if (value > address(M))
        L = M+1;
        M = floor((L+R)/2);
    elseif (value < address(M))
        R = M-1;
        M = floor((L+R)/2);
    elseif (value == address(M))
        index = M;
        return
    end
end

end
