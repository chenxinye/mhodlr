function check_point = product(varargin)
    if nargin <= 1
        seed = 0;
    else
        seed = varargin{1};
    end

    check_point1 = product1(seed);
    check_point2 = product2(seed);

    check_point = (check_point1) + (check_point2 == 1);
end


function check_point = product1(seed)
    rng(seed);
    addpath("../mhodlr/")
    rng(0); %fix randomness
    A = rand(30, 20);
    B = rand(20, 50);

    depth = 99;
    min_block_size = 2;
    epsilon = 1e-16;
    hA = hodlr(A, depth, min_block_size, 'svd', epsilon); 
    hB = hodlr(B, depth, min_block_size, 'svd', epsilon); 

    C = A * B;
    hC = dot(hA, hB); 
    disp("checkpoint1 error:")
    disp(norm(C - hC.dense));

    check_point = sum(hC.shape == size(C)) == 2;
end



function check_point = product2(seed)
    rng(seed);
    addpath("../mhodlr/")

    u = precision('h'); % or ther precision customization
    set_prec(u);

    rng(0); %fix randomness
    A = rand(30, 20);
    B = rand(20, 50);

    depth = 99;
    min_block_size = 2;
    epsilon = 1e-16;
    hA = hodlr(A, depth, min_block_size, 'svd', epsilon); 
    hB = hodlr(B, depth, min_block_size, 'svd', epsilon); 

    C = A * B;
    hC = mdot(hA, hB); 
    disp("checkpoint2 error:")
    disp(norm(C - hC.dense));

    check_point = sum(hC.shape == size(C)) == 2;
end



function check_point = product3(seed)
    rng(seed);
    addpath("../mhodlr/")
    
    u1 = precision('d');
    u2 = precision('s');
    u3 = precision('h');
    u4 = precision('b');
    
    u_chain = prec_chain(u1, u2, u3, u4);
    depth = 5;
    eps = 1e-12;
    A = rand(100, 80);
    B = rand(80, 120);

    aphA = amphodlr(u_chain, A, depth, 10, 'svd', eps); 
    aphB = amphodlr(u_chain, B, depth, 10, 'svd', eps); 
    
    C = hdot(aphA, aphB);
    disp("checkpoint3 error:")
    disp(norm(C - hC.dense));
    check_point = sum(hC.shape == size(C)) == 2;
end


function check_point = product4(seed)
    rng(seed);
    addpath("../mhodlr/")
    
    u1 = precision('d');
    u2 = precision('s');
    u3 = precision('h');
    u4 = precision('b');
    
    u_chain = prec_chain(u1, u2, u3, u4);
    depth = 5;
    eps = 1e-12;
    A = rand(100, 80);
    B = rand(80, 120);
    aphA = amphodlr(u_chain, A, depth, 10, 'svd', eps); 
    aphB = amphodlr(u_chain, B, depth, 10, 'svd', eps); 
    
    C = hdot(aphA, aphB);
    disp("checkpoint4 error:")
    disp(norm(C - hC.dense));
    check_point = sum(hC.shape == size(C)) == 2;
end
