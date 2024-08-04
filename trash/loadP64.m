% read matrix

function data = loadP64()
    path = 'data';
    matrix = 'P64';
    cs = 128;
    file = sprintf('%s/root_%s_cs%d.bin', path, matrix, cs);
    f = fopen(file);
    data = fread(f,'double');
    fclose(f);
    n = size(data,1);
    n = sqrt(n);
    data = reshape(data,n,n);
end