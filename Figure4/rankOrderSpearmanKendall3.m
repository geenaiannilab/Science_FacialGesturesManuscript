function out = rankOrderSpearmanKendall3(R1,R2,R3,varargin)
% Rank-order similarity between three MxM Pearson-R matrices (upper triangle).
% Computes BOTH Spearman rho and Kendall tau with:
%  - Mantel permutation test (neuron-label shuffles) for p-values
%  - Bootstrap 95% CIs by resampling neurons with replacement
%
% Usage:
%  out = rankOrderSpearmanKendall3(R1,R2,R3,'NPerm',10000,'NBoot',2000,'Tail','both','Seed',123,'UseDiag',false)
%
% Output fields:
%  .pairs         {'R1-R2','R1-R3','R2-R3'}, gesture comparisons 
%  .spearman.r    [3x1], .spearman.p [3x1], .spearman.CI [2x3], .spearman.null {3x1}
%  .kendall.tau   [3x1], .kendall.p  [3x1], .kendall.CI  [2x3], .kendall.null  {3x1}
%  .idxUpper      logical mask used
%  .summary       table with both metrics

    % ---- options ----
    p = inputParser;
    addParameter(p,'NPerm',10000,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(p,'NBoot',2000, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(p,'Tail','both', @(s)ischar(s)||isstring(s));
    addParameter(p,'Seed',0, @(x)isempty(x)||(isnumeric(x)&&isscalar(x)));
    addParameter(p,'UseDiag',false, @(x)islogical(x)&&isscalar(x));
    parse(p,varargin{:});
    NPerm  = p.Results.NPerm;
    NBoot  = p.Results.NBoot;
    Tail   = lower(string(p.Results.Tail));
    Seed   = p.Results.Seed;
    UseDiag= p.Results.UseDiag;

    if ~isequal(size(R1),size(R2),size(R3))
        error('R1, R2, R3 must be the same size MxM.');
    end
    M = size(R1,1);

    if ~isempty(Seed), rng(Seed); end

    % Symmetrize defensively
    R1=(R1+R1')/2; R2=(R2+R2')/2; R3=(R3+R3')/2;

    % Vectorization mask
    if UseDiag
        idx = triu(true(M),0);
    else
        idx = triu(true(M),1);
    end

    Rs = {R1,R2,R3};
    pairIdx  = [1 2; 1 3; 2 3];
    pairName = {'R1-R2','R1-R3','R2-R3'};

    % Prealloc
    rS  = nan(3,1); pS = nan(3,1); CI_S = nan(2,3); nullS = cell(3,1);
    rK  = nan(3,1); pK = nan(3,1); CI_K = nan(2,3); nullK = cell(3,1);

    for k = 1:3
        A = Rs{pairIdx(k,1)}; B = Rs{pairIdx(k,2)};
        a = A(idx); b = B(idx);
        good = isfinite(a) & isfinite(b);
        a = a(good); b = b(good);

        % ---- observed metrics ----
        rS(k) = corr(a, b, 'type','Spearman', 'rows','pairwise');
        rK(k) = corr(a, b, 'type','Kendall',  'rows','pairwise'); % tau-b

        % ---- Mantel permutations (shared for S & K) ----
        if NPerm>0
            ns = nan(NPerm,1); nk = nan(NPerm,1);
            if ~isempty(Seed), rng(Seed+101+k); end
            for i = 1:NPerm
                pperm = randperm(M);
                Bp = B(pperm,pperm);
                bp = Bp(idx);
                goodp = isfinite(a) & isfinite(bp);
                ap = a(goodp); bp = bp(goodp);
                % compute both with same permutation
                ns(i) = corr(ap, bp, 'type','Spearman','rows','pairwise');
                nk(i) = corr(ap, bp, 'type','Kendall', 'rows','pairwise');
            end
            pS(k) = tail_p(ns, rS(k), Tail);
            pK(k) = tail_p(nk, rK(k), Tail);
            nullS{k} = ns; nullK{k} = nk;
        end

        % ---- bootstrap CIs (shared resamples; compute both) ----
        if NBoot>0
            bs = nan(NBoot,1); bk = nan(NBoot,1);
            if ~isempty(Seed), rng(Seed+202+k); end
            for bIter = 1:NBoot
                ii = randsample(M, M, true);
                Aii = (A(ii,ii)+A(ii,ii)')/2;
                Bii = (B(ii,ii)+B(ii,ii)')/2;
                a2 = Aii(idx); b2 = Bii(idx);
                g2 = isfinite(a2) & isfinite(b2);
                a2 = a2(g2); b2 = b2(g2);
                bs(bIter) = corr(a2, b2, 'type','Spearman','rows','pairwise');
                bk(bIter) = corr(a2, b2, 'type','Kendall', 'rows','pairwise');
            end
            CI_S(:,k) = prctile(bs,[2.5 97.5])';
            CI_K(:,k) = prctile(bk,[2.5 97.5])';
        end
    end

    % Summary table
    T = table(pairName', ...
        rS, pS, CI_S(1,:)', CI_S(2,:)', ...
        rK, pK, CI_K(1,:)', CI_K(2,:)', ...
        'VariableNames', {'Pair', ...
            'Spearman_r','Spearman_p','Spearman_CI_low','Spearman_CI_high', ...
            'Kendall_tau','Kendall_p','Kendall_CI_low','Kendall_CI_high'});

    % Pack output
    out.pairs        = pairName;
    out.idxUpper     = idx;
    out.spearman.r   = rS;   out.spearman.p = pS;   out.spearman.CI = CI_S; out.spearman.null = nullS;
    out.kendall.tau  = rK;   out.kendall.p  = pK;   out.kendall.CI  = CI_K; out.kendall.null  = nullK;
    out.summary      = T;
end

function p = tail_p(nullDist, obs, Tail)
    switch Tail
        case "both",  p = mean(abs(nullDist) >= abs(obs));
        case "right", p = mean(nullDist >= obs);
        case "left",  p = mean(nullDist <= obs);
        otherwise,    error('Tail must be ''both'',''right'',''left''.');
    end
end
