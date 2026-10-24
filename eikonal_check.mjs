// Eikonal cross-check: does the Taub-NUT scalar QNM track its photon sphere?
//   QNM (from the continued fraction)  vs  (q+1/2)*Omega_c - i (n+1/2)*lambda
// Units match the paper's code: M = 1/2.   Run: node eikonal_check.mjs
const cadd=(a,b)=>[a[0]+b[0],a[1]+b[1]], csub=(a,b)=>[a[0]-b[0],a[1]-b[1]];
const cmul=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]];
const cinv=a=>{const d=a[0]*a[0]+a[1]*a[1];return[a[0]/d,-a[1]/d];};
const cdiv=(a,b)=>cmul(a,cinv(b)), smul=(s,a)=>[s*a[0],s*a[1]], cabs=a=>Math.hypot(a[0],a[1]);
const RR=x=>[x,0], CIj=[0,1];

// --- Taub-NUT scalar continued fraction (mirror of the shader / SpaceTimeDefinition.py) ---
function cfTN(w,l,Cas,depth){
  const l2=l*l,M=0.5,rt=Math.sqrt(M*M+l2),rp=M+rt,rm=M-rt,rp2=rp*rp,rm2=rm*rm,s=(rp2+l2)/(2*rt),dr=rp-rm,dr2=dr*dr;
  const iw=cmul(CIj,w),w2=cmul(w,w),iws=smul(s,iw),omw=csub(RR(1),iw);
  const K1=(-2*rm*rp2-4*l2*l2-2*l2*rp*rm-4*l2*rp2)/dr2,K2=(rm2*rp2+2*rm2*rp+6*l2*l2+2*l2*rm2+4*l2*rm*rp+2*l2*rp2)/dr2;
  const K3=(-2*rm*rm2*rp-4*l2*l2-4*l2*rm2-2*l2*rm*rp)/dr2,K4=(rm2*rm2+l2*l2+2*l2*rm2)/dr2;
  const P0=smul((rp2+l2*l2+2*l2*rp2)/dr2,w2);
  const P1=csub(cadd(cadd(iw,smul(dr,iw)),smul(K1+4*l2,w2)),cadd(RR(1),RR(Cas)));
  const P2=csub(csub(csub(csub(cmul(csub(RR(3),iw),omw),smul(2*dr,cadd(w2,iw))),smul(2,csub(smul(4*l2,w2),RR(Cas)))),smul(dr2,w2)),smul(-K2,w2));
  const P3=cadd(cadd(cadd(cmul(csub(smul(2,iw),RR(3)),omw),smul(dr,cadd(smul(2,w2),iw))),csub(smul(4*l2,w2),RR(Cas))),smul(K3,w2));
  const P4=cadd(cmul(omw,omw),smul(K4,w2));
  const tiny=[1e-30,0];
  let L1p0=RR(1),L1p1=RR(1),L1p2=RR(1),L1p3=RR(1),L2p0=RR(1),L2p1=RR(1),L2p2=RR(1),aP=RR(1),f=[0,0],Cl=[0,0],Dl=[0,0];
  for(let i=0;i<=depth;i++){const n=i,t=x=>csub(RR(n+x),iws);
    const A=cadd(P0,cmul(t(1),cadd(t(0),RR(1))));
    const B=cadd(P1,cmul(t(0),cadd(smul(-4,t(-1)),cadd(cadd(smul(-2,omw),smul(2*dr,iw)),RR(-4)))));
    const Cc=cadd(P2,cmul(t(-1),cadd(smul(6,t(-2)),cadd(csub(smul(6,omw),smul(4*dr,iw)),RR(6)))));
    const Dd=cadd(P3,cmul(t(-2),cadd(smul(-4,t(-3)),cadd(cadd(smul(-6,omw),smul(2*dr,iw)),RR(-4)))));
    const E=cadd(P4,cmul(t(-3),cadd(t(-4),cadd(smul(2,omw),RR(1)))));
    let L10,L11,L12,L13;
    if(i>=3){const q=cdiv(E,L1p3);L10=A;L11=csub(B,cmul(L1p0,q));L12=csub(Cc,cmul(L1p1,q));L13=csub(Dd,cmul(L1p2,q));}else{L10=A;L11=B;L12=Cc;L13=Dd;}
    let al,be,ga;
    if(i>=2){const q=cdiv(L13,L2p2);al=L10;be=csub(L11,cmul(L2p0,q));ga=csub(L12,cmul(L2p1,q));}else{al=L10;be=L11;ga=L12;}
    if(i===0){f=cabs(be)<1e-30?tiny:be;Cl=f;Dl=[0,0];}else{const an=smul(-1,cmul(aP,ga)),bn=be;
      Dl=cadd(bn,cmul(an,Dl));if(cabs(Dl)<1e-30)Dl=tiny;Cl=cadd(bn,cdiv(an,Cl));if(cabs(Cl)<1e-30)Cl=tiny;Dl=cinv(Dl);f=cmul(f,cmul(Cl,Dl));}
    L1p0=L10;L1p1=L11;L1p2=L12;L1p3=L13;L2p0=al;L2p1=be;L2p2=ga;aP=al;}
  return f;}
function newton(w0,l,Cas){let w=w0.slice();for(let k=0;k<80;k++){const f=cfTN(w,l,Cas,110),h=[1e-6,0];
  const fp=cdiv(csub(cfTN(cadd(w,h),l,Cas,110),cfTN(csub(w,h),l,Cas,110)),[2e-6,0]);if(cabs(fp)<1e-30)break;
  const st=cdiv(f,fp);w=csub(w,st);if(cabs(st)<1e-10)break;}return cabs(cfTN(w,l,Cas,110))<1e-2?w:null;}

// --- Equatorial photon sphere of Taub-NUT (M=1/2, theta=pi/2 so NUT potential N=2l) ---
function geq(r,l){const M=0.5,R2=r*r+l*l,f=(r*r-2*M*r-l*l)/R2;
  return {f, gtt:-1/f+4*l*l/R2, gtp:-2*l/R2, gpp:1/R2};}       // inverse-metric comps
function bBranch(r,l,sgn){const {gtt,gtp,gpp}=geq(r,l);
  return (gtp+sgn*Math.sqrt(gtp*gtp-gpp*gtt))/gpp;}             // critical impact param b(r)
function hval(r,l,b){const {gtt,gtp,gpp}=geq(r,l);return -gtt+2*gtp*b-gpp*b*b;}
function photon(l,sgn){
  // r_ph = extremum of b(r) along a branch (db/dr = 0)
  const rh=0.5+Math.sqrt(0.25+l*l); let rlo=rh+1e-3, rhi=40, rph=null;
  const db=r=>(bBranch(r+1e-5,l,sgn)-bBranch(r-1e-5,l,sgn))/2e-5;
  let prev=db(rlo);
  for(let r=rlo+0.002;r<rhi;r+=0.002){const cur=db(r);
    if(prev*cur<0){ // bisection
      let a=r-0.002,c=r;for(let k=0;k<60;k++){const m=(a+c)/2; (db(a)*db(m)<0)?c=m:a=m;} rph=(a+c)/2; break;}
    prev=cur;}
  if(rph===null) return null;
  const b=bBranch(rph,l,sgn), {gtt,gtp,gpp,f}=geq(rph,l);
  const Omega=(-gtp+gpp*b)/(-gtt+gtp*b);                        // coord angular velocity
  const tdot=-gtt+gtp*b;                                        // dt/dlambda per E
  const R=r=>geq(r,l).f*hval(r,l,b);                            // rdot^2/E^2 = f*h
  const d=5e-4, Rpp=(R(rph+d)-2*R(rph)+R(rph-d))/(d*d);
  const lam=Math.sqrt(Math.max(Rpp,0)/(2*tdot*tdot));           // Lyapunov exponent
  return {rph,Omega:Math.abs(Omega),lam};
}

// --- Run: track the fundamental (q=2, C=6, n=0) scalar QNM vs the photon sphere ---
const C=6, q=2;
console.log("  l     QNM omega                 |   Omega_c(pro/retro) lambda(pro/retro) | eikonal (q+.5)Om -i(.5)lam   ratio Im/(.5 lam_avg)");
let seed=[0.967288,-0.193518];                                 // Schwarzschild fundamental
for(const l of [0,0.25,0.5,0.75,1.0,1.5,2.0]){
  const w=newton(seed,l,C); if(w) seed=w;
  const pp=photon(l,+1), pr=photon(l,-1);
  const OmAvg=(pp.Omega+pr.Omega)/2, lamAvg=(pp.lam+pr.lam)/2;
  const eikR=(q+0.5)*OmAvg, eikI=0.5*lamAvg;
  const ratio = w? (-w[1])/(0.5*lamAvg) : NaN;
  const fmt=(x,d=4)=>x.toFixed(d).padStart(7);
  console.log(`${fmt(l,2)}  ${w?`${fmt(w[0])} ${w[1]<0?'-':'+'} ${fmt(Math.abs(w[1]))}i`:'   (lost)     '}  |  ${fmt(pp.Omega)}/${fmt(pr.Omega)}  ${fmt(pp.lam)}/${fmt(pr.lam)} |  ${fmt(eikR)} - ${fmt(eikI)}i        ${fmt(ratio,3)}`);
}
console.log("\nn=0 fundamental, scalar, C=q(q+1)=6.  M=1/2 units.  Eikonal uses the pro/retro-averaged photon orbit.");
