format long

M = 16 * 1024; % cache size in cache lines
p = (M - 1) / M; % probability a particular line is chosen
m = 1; % # writes in stencil

for n = 2:6:8

    n % number of trailing reads (planes) in stencil access pattern

    % k = # cache lines between accesses in access pattern
    bx = [1 2 4 8 16 32 64 96 128 192 384].'; by = 1;
    k = bx * by / 8;
    bx = 384; by = [2 4 8 16 32 48 64 96 128 192 384].';
    k = [k; bx * by / 8];

    % solve alpha = 1 - ((M-1)/M)^(k*(1+m+n*alpha)) for alpha
    result = -1./(n*k*log(p)).*lambertw(n*k*log(p).*p.^(k*(m+1)+n*k)) + 1;

    disp(result);

    % check answer
    fprintf('Error = %g\n', norm(result - (1 - p.^(k.*(1+m+n*result)))));

end
