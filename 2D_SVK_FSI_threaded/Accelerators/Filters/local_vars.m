function [Vprev, Wprev, Tprev, T2prev, Kprev, nRprev, nXprev, nVprev, ...
    nWprev, V, W, T, T2, K, nR, nX, nV, nW] = local_vars(M)

    Vprev = M.Vprev;
    Wprev = M.Wprev;
    Tprev = M.Tprev;
    
    Kprev = M.Kprev;
    T2prev = M.T2prev;
    nRprev = M.nRprev;
    nXprev = M.nXprev;
    nVprev = M.nVprev;
    nWprev = M.nWprev;

    V = M.V;
    W = M.W;
    T = M.T;
    
    T2 = M.T2;
    K = M.K;
    nR = M.nR;
    nX = M.nX;
    nV = M.nV;
    nW = M.nW;
    
end
