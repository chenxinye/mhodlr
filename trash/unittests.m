addpath('tests');

function verbose = unittests()
    if ~test_recover
        error('Do not pass test 1.')
    end

    if ~test_chol
        error('Do not pass test 2.')
    end
    
    if ~test_lu
        error('Do not pass test 3.')
    end
    
    if ~test_dot
        error('Do not pass test 4.')
    end

    verbose = true;
end
