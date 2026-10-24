// Sanity check for the Kerr-Schild inverse metric used in the shader.
// Mirrors invMetric() in index.html. Run: node test_metric.mjs
const eta = [[-1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];

function invMetric(X, a){                 // X=[t,x,y,z], M=1
  const [ , x,y,z] = X, a2=a*a, rho2=x*x+y*y+z*z;
  const r2 = 0.5*((rho2-a2)+Math.sqrt((rho2-a2)**2 + 4*a2*z*z));
  const r  = Math.sqrt(r2);
  const f  = 2*r2*r/(r2*r2 + a2*z*z);
  const l  = [1,(r*x+a*y)/(r2+a2),(r*y-a*x)/(r2+a2),z/r];
  const lu = [-l[0],l[1],l[2],l[3]];
  const g=[];
  for(let i=0;i<4;i++){g[i]=[]; for(let j=0;j<4;j++) g[i][j]=eta[i][j]-f*lu[i]*lu[j];}
  return g;
}
function mul(g,p){ return g.map(row=>row.reduce((s,v,j)=>s+v*p[j],0)); }
function quad(g,p){ return mul(g,p).reduce((s,v,i)=>s+v*p[i],0); }
function assert(c,m){ if(!c){console.error('FAIL:',m); process.exit(1);} }

// 1) Far from the hole g^{mu nu} -> Minkowski diag(-1,1,1,1)
for(const a of [0,0.5,0.9]){
  const g=invMetric([0,1e5,0,1e3],a);   // f=2M/r ~ 2e-5, must approach Minkowski
  assert(Math.abs(g[0][0]+1)<1e-3 && Math.abs(g[1][1]-1)<1e-3, `asymptotic flat a=${a}`);
}

// 2) A null covariant momentum stays null: build p from a direction, solve p_0, check g^{ab}pa pb=0
for(const a of [0,0.6,0.95]){
  const X=[0,15,4,6], g=invMetric(X,a), dir=[0.3,-0.7,0.64];
  const p=[0,...dir];
  const A=g[0][0], B=2*(g[0][1]*p[1]+g[0][2]*p[2]+g[0][3]*p[3]);
  const C=quad(g,[0,...dir]);
  p[0]=(-B-Math.sqrt(B*B-4*A*C))/(2*A);
  assert(Math.abs(quad(g,p))<1e-9, `null condition a=${a} (got ${quad(g,p)})`);
}

// 3) Schwarzschild (a=0) is spherically symmetric: g^{tt} depends only on radius
{ const r=8;
  const g1=invMetric([0,r,0,0],0), g2=invMetric([0,0,r,0],0), g3=invMetric([0,0,0,r],0);
  assert(Math.abs(g1[0][0]-g2[0][0])<1e-9 && Math.abs(g1[0][0]-g3[0][0])<1e-9,'a=0 isotropy');
  // g^{tt} = -(1+2M/r) in Kerr-Schild; check numerically
  assert(Math.abs(g1[0][0] - (-(1+2/r)))<1e-9, 'g^{tt} = -(1+2M/r)');
}
// 4) Matrix-free hot-loop path must equal the mat4 path (guards the refactor).
//    ksPrim -> f, lup;  dx = eta p - f phi lup  must equal g*p ;  H_scalar == quad
function ksPrim(X,a){
  const [ ,x,y,z]=X, a2=a*a, rho2=x*x+y*y+z*z;
  const r2=Math.max(0.5*((rho2-a2)+Math.sqrt((rho2-a2)**2+4*a2*z*z)),1e-6), r=Math.sqrt(r2);
  const f=2*r2*r/(r2*r2+a2*z*z);
  return {f, lup:[-1,(r*x+a*y)/(r2+a2),(r*y-a*x)/(r2+a2),z/r]};
}
for(const a of [0,0.6,0.95]){
  const X=[0,12,-3,5], p=[0.4,0.3,-0.7,0.5], g=invMetric(X,a);
  const {f,lup}=ksPrim(X,a), phi=lup.reduce((s,v,i)=>s+v*p[i],0);
  const dxScalar=[-p[0]-f*phi*lup[0], p[1]-f*phi*lup[1], p[2]-f*phi*lup[2], p[3]-f*phi*lup[3]];
  const dxMat=mul(g,p);
  for(let i=0;i<4;i++) assert(Math.abs(dxScalar[i]-dxMat[i])<1e-9, `dx scalar==mat4 a=${a} comp ${i}`);
  const Hscalar = -p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+p[3]*p[3] - f*phi*phi;
  assert(Math.abs(Hscalar - quad(g,p))<1e-9, `H scalar==quad a=${a}`);
}
// 5) Taub-NUT: verify the analytic inverse metric really inverts g_{mu nu},
//    the null condition holds, and l=0 reduces to Schwarzschild (spherical).
// Leguery gauge (paper eq 3.1): NUT potential N = 2l(1 - cos th), string on theta=pi only.
function tnCov(r,th,l){ const R2=r*r+l*l, f=(r*r-2*r-l*l)/R2, s=Math.sin(th), N=2*l*(1-Math.cos(th));
  const g=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]];
  g[0][0]=-f; g[3][3]=R2*s*s - f*N*N; g[0][3]=g[3][0]=-f*N; g[1][1]=1/f; g[2][2]=R2; return g; }
function tnInv(r,th,l){ const R2=r*r+l*l, f=(r*r-2*r-l*l)/R2, s=Math.sin(th), s2=s*s, N=2*l*(1-Math.cos(th));
  const g=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]];
  g[0][0]=-1/f+N*N/(R2*s2); g[1][1]=f; g[2][2]=1/R2; g[3][3]=1/(R2*s2);
  g[0][3]=g[3][0]=-N/(R2*s2); return g; }
function mm(a,b){ return a.map((row,i)=>b[0].map((_,j)=>row.reduce((s,_,k)=>s+a[i][k]*b[k][j],0))); }

for(const [r,th,l] of [[10,1.1,0.5],[7,0.8,1.0],[20,2.0,0.3]]){
  const I = mm(tnCov(r,th,l), tnInv(r,th,l));                       // must be identity
  for(let i=0;i<4;i++) for(let j=0;j<4;j++)
    assert(Math.abs(I[i][j]-(i===j?1:0))<1e-7, `TN g*g^-1=I at (${r},${th},${l}) [${i}][${j}]=${I[i][j]}`);
  // null momentum stays null under the inverse metric
  const g=tnInv(r,th,l), pr=0.5, pth=-0.3, pph=0.7;
  const A=g[0][0], B=2*g[0][3]*pph, C=g[1][1]*pr*pr+g[2][2]*pth*pth+g[3][3]*pph*pph;
  const pt=(-B-Math.sqrt(B*B-4*A*C))/(2*A);
  const p=[pt,pr,pth,pph];
  let H=0; for(let i=0;i<4;i++)for(let j=0;j<4;j++) H+=g[i][j]*p[i]*p[j];
  assert(Math.abs(H)<1e-9, `TN null condition (${r},${th},${l}) got ${H}`);
}
// l=0 Taub-NUT == Schwarzschild in spherical coords
{ const r=9, th=1.3, g=tnInv(r,th,0), f=1-2/r;
  assert(Math.abs(g[0][0]-(-1/f))<1e-9 && Math.abs(g[1][1]-f)<1e-9
      && Math.abs(g[2][2]-1/(r*r))<1e-9 && Math.abs(g[3][3]-1/(r*r*Math.sin(th)**2))<1e-9
      && Math.abs(g[0][3])<1e-12, 'TN(l=0) == Schwarzschild spherical'); }

console.log('all metric checks passed');
