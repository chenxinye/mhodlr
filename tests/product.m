function check_point = product(seed)
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

    check_point = sum(hC.shape == size(C)) == 2;
end