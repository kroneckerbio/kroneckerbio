function m = initializeModelBase(m)
% Create common required fields in all model types

if nargin < 1
    m = [];
end

m.nv = 0;
m.nk = 0;
m.ns = 0;
m.nu = 0;
m.nx = 0;
m.nr = 0;
m.ny = 0;
m.nz = 0;

m.k    = zeros(0,1);
m.s    = zeros(0,1);
m.u    = zeros(0,1);

m.dv   = zeros(0,1);
m.vxInd  = zeros(0,1);
m.vuInd  = zeros(0,1);
m.rOrder = zeros(0,1);
m.krInd  = zeros(0,1);

m.d2x0ds2   = @(s)zeros(0,0);
m.dx0ds     = @(s)zeros(0,0);
m.x0        = @(s)zeros(0,1);

m.A1 = zeros(0,0);
m.A2 = zeros(0,0);
m.A3 = zeros(0,0);
m.A4 = zeros(0,0);
m.A5 = zeros(0,0);
m.A6 = zeros(0,0);
m.a  = zeros(0,0);

m.B1 = zeros(0,0);
m.B2 = zeros(0,0);
m.b  = zeros(0,1);

m.C1      = zeros(0,0);
m.C2      = zeros(0,0);
m.c       = zeros(0,1);

m.dA1dk = zeros(0,0);
m.dA2dk = zeros(0,0);
m.dA3dk = zeros(0,0);
m.dA4dk = zeros(0,0);
m.dA5dk = zeros(0,0);
m.dA6dk = zeros(0,0);
m.dadk  = zeros(0,0);
m.dA1dk_fk_x  = zeros(0,0);
m.dA2dk_fk_xx = zeros(0,0);
m.dA3dk_fk_ux = zeros(0,0);
m.dA4dk_fk_xu = zeros(0,0);
m.dA5dk_fk_uu = zeros(0,0);
m.dA6dk_fk_u  = zeros(0,0);

m.f = @(t,x,u)(zeros(0,1));

m.dfdx = @(t,x,u)(zeros(0,0));
m.dfdu = @(t,x,u)(zeros(0,0));
m.dfdk = @(t,x,u)(zeros(0,0));

m.d2fdx2  = @(t,x,u)(zeros(0,0));
m.d2fdu2  = @(t,x,u)(zeros(0,0));
m.d2fdk2  = @(t,x,u)(zeros(0,0));
m.d2fdudx = @(t,x,u)(zeros(0,0));
m.d2fdxdu = @(t,x,u)(zeros(0,0));
m.d2fdkdx = @(t,x,u)(zeros(0,0));
m.d2fdxdk = @(t,x,u)(zeros(0,0));
m.d2fdkdu = @(t,x,u)(zeros(0,0));
m.d2fdudk = @(t,x,u)(zeros(0,0));

m.S = zeros(0,0);

m.D1 = zeros(0,0);
m.D2 = zeros(0,0);
m.D3 = zeros(0,0);
m.D4 = zeros(0,0);
m.D5 = zeros(0,0);
m.D6 = zeros(0,0);
m.d  = zeros(0,0);

m.dD1dk = zeros(0,0);
m.dD2dk = zeros(0,0);
m.dD3dk = zeros(0,0);
m.dD4dk = zeros(0,0);
m.dD5dk = zeros(0,0);
m.dD6dk = zeros(0,0);
m.dddk  = zeros(0,0);
m.dD1dk_rk_x  = zeros(0,0);
m.dD2dk_rk_xx = zeros(0,0);
m.dD3dk_rk_ux = zeros(0,0);
m.dD4dk_rk_xu = zeros(0,0);
m.dD5dk_rk_uu = zeros(0,0);
m.dD6dk_rk_u  = zeros(0,0);

m.r = @(t,x,u)(zeros(0,1));

m.drdx = @(t,x,u)(zeros(0,0));
m.drdu = @(t,x,u)(zeros(0,0));
m.drdk = @(t,x,u)(zeros(0,0));

m.d2rdx2  = @(t,x,u)(zeros(0,0));
m.d2rdu2  = @(t,x,u)(zeros(0,0));
m.d2rdk2  = @(t,x,u)(zeros(0,0));
m.d2rdudx = @(t,x,u)(zeros(0,0));
m.d2rdxdu = @(t,x,u)(zeros(0,0));
m.d2rdkdx = @(t,x,u)(zeros(0,0));
m.d2rdxdk = @(t,x,u)(zeros(0,0));
m.d2rdkdu = @(t,x,u)(zeros(0,0));
m.d2rdudk = @(t,x,u)(zeros(0,0));

m.v       = @(t,x,u)(zeros(0,1));
m.dvdx    = @(t,x,u)(zeros(0,0));
m.dvdu    = @(t,x,u)(zeros(0,0));
m.d2vdx2  = @(t,x,u)(zeros(0,0));
m.d2vdu2  = @(t,x,u)(zeros(0,0));
m.d2vdudx = @(t,x,u)(zeros(0,0));
m.d2vdxdu = @(t,x,u)(zeros(0,0));

m.y       = @(t,x,u)(zeros(0,1));

m.dydx    = @(t,x,u)(zeros(0,0));
m.dydu    = @(t,x,u)(zeros(0,0));
m.dydk    = @(t,x,u)(zeros(0,0));

m.d2ydx2  = @(t,x,u)(zeros(0,0));
m.d2ydu2  = @(t,x,u)(zeros(0,0));
m.d2ydk2  = @(t,x,u)(zeros(0,0));
m.d2ydudx = @(t,x,u)(zeros(0,0));
m.d2ydxdu = @(t,x,u)(zeros(0,0));
m.d2ydkdx = @(t,x,u)(zeros(0,0));
m.d2ydxdk = @(t,x,u)(zeros(0,0));
m.d2ydkdu = @(t,x,u)(zeros(0,0));
m.d2ydudk = @(t,x,u)(zeros(0,0));

m.Ready  = true;
m.add    = struct;
m.Update = @(k,x0,q)(InitializeModel(name));

m.add.nv  = 0;
m.add.nk  = 0;
m.add.ns  = 0;
m.add.nu  = 0;
m.add.nx  = 0;
m.add.nr  = 0;
m.add.ny  = 0;
m.add.nz  = 0;
