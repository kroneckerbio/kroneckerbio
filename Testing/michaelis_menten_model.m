function [m_kron, con, obj, opts] = michaelis_menten_model()
syms Km kcat S0 E S P

m.kNames = {'Km';'kcat'};
m.kSyms = [kcat;Km];
m.k = [10;2];

m.sNames = {'S0'};
m.sSyms = S0;
m.s = 30;

m.uNames = {'E'};
m.uSyms = E;
m.u = 1;

m.xNames = {'S';'P'};
m.xSyms = [S;P];
m.x0 = [S0;0];

r = kcat*E*S/(Km+S);
m.f = [-r;r];

m.yNames = {'S','P'};
m.y = [S;P];

m_kron = symbolic2PseudoKronecker(m);
con = InitialValueExperiment(m_kron);
obj = [];
opts.Verbose = false;
end