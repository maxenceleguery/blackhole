// Validate the QNM continued-fraction shader math (mirror of qnm.html GLSL) against
// known Leaver QNM values. Run: node test_qnm.mjs
const C=(a,b)=>[a,b];
const cadd=(a,b)=>[a[0]+b[0],a[1]+b[1]], csub=(a,b)=>[a[0]-b[0],a[1]-b[1]];
const cmul=(a,b)=>[a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]];
const cinv=a=>{const d=a[0]*a[0]+a[1]*a[1];return [a[0]/d,-a[1]/d];};
const cdiv=(a,b)=>cmul(a,cinv(b));
const smul=(s,a)=>[s*a[0],s*a[1]];
const cabs=a=>Math.hypot(a[0],a[1]);
const R=x=>[x,0]; const CI=[0,1];
const assert=(c,m)=>{if(!c){console.error('FAIL:',m);process.exit(1);}};

// ---- Schwarzschild 3-term coefficients (eps=-1 scalar), ell multipole ----
function coefsSchw(n,w,ell){
  const eps=-1, miw=cmul([0,-1],w);
  const al=cadd(cadd(cadd(R(n*n), smul(n,cadd(smul(2,miw),R(2)))), smul(2,miw)), R(1));
  const be=smul(-1, cadd(cadd(cadd(cadd(R(2*n*n), smul(n,cadd(smul(8,miw),R(2)))),
                 smul(8,cmul(miw,miw))), smul(4,miw)), R(ell*(ell+1)-eps)));
  const ga=cadd(cadd(cadd(R(n*n), smul(4*n,miw)), smul(4,cmul(miw,miw))), R(-eps-1));
  return [al,be,ga];
}
// Lentz CF for a 3-term recurrence; QNM <=> |CF| -> 0
function logCF3(w,ell,depth){
  const tiny=[1e-30,0];
  let [al,be,ga]=coefsSchw(0,w,ell);
  let f=cabs(be)<1e-30?tiny:be, Cl=f, Dl=[0,0], aPrev=al;
  for(let i=1;i<=depth;i++){
    [al,be,ga]=coefsSchw(i,w,ell);
    const an=smul(-1,cmul(aPrev,ga)), bn=be;
    Dl=cadd(bn,cmul(an,Dl)); if(cabs(Dl)<1e-30) Dl=tiny;
    Cl=cadd(bn,cdiv(an,Cl)); if(cabs(Cl)<1e-30) Cl=tiny;
    Dl=cinv(Dl); f=cmul(f,cmul(Cl,Dl)); aPrev=al;
  }
  return Math.log(cabs(f)+1e-30);
}
// ---- Taub-NUT 5-term, fused 5->3 reduction + Lentz (mirror of logCF_TN) ----
function logCF_TN(w,l,Cas,depth){
  const l2=l*l, M=0.5, rt=Math.sqrt(M*M+l2), rp=M+rt, rm=M-rt;
  const rp2=rp*rp, rm2=rm*rm, s=(rp2+l2)/(2*rt), dr=rp-rm, dr2=dr*dr;
  const iw=cmul(CI,w), w2=cmul(w,w), iws=smul(s,iw), omw=csub(R(1),iw);
  const K1=(-2*rm*rp2-4*l2*l2-2*l2*rp*rm-4*l2*rp2)/dr2;
  const K2=(rm2*rp2+2*rm2*rp+6*l2*l2+2*l2*rm2+4*l2*rm*rp+2*l2*rp2)/dr2;
  const K3=(-2*rm*rm2*rp-4*l2*l2-4*l2*rm2-2*l2*rm*rp)/dr2;
  const K4=(rm2*rm2+l2*l2+2*l2*rm2)/dr2;
  const P0=smul((rp2+l2*l2+2*l2*rp2)/dr2, w2);
  const P1=csub(cadd(cadd(iw,smul(dr,iw)), smul(K1+4*l2,w2)), cadd(R(1),R(Cas)));
  const P2=csub(csub(csub(csub(cmul(csub(R(3),iw),omw), smul(2*dr,cadd(w2,iw))),
             smul(2,csub(smul(4*l2,w2),R(Cas)))), smul(dr2,w2)), smul(-K2,w2));
  const P3=cadd(cadd(cadd(cmul(csub(smul(2,iw),R(3)),omw), smul(dr,cadd(smul(2,w2),iw))),
             csub(smul(4*l2,w2),R(Cas))), smul(K3,w2));
  const P4=cadd(cmul(omw,omw), smul(K4,w2));
  const tiny=[1e-30,0];
  let L1p0=R(1),L1p1=R(1),L1p2=R(1),L1p3=R(1),L2p0=R(1),L2p1=R(1),L2p2=R(1),aPrev=R(1);
  let f=[0,0],Cl=[0,0],Dl=[0,0];
  for(let i=0;i<=depth;i++){
    const n=i;
    const t=x=>csub(R(n+x),iws);
    const A=cadd(P0, cmul(t(1), cadd(t(0),R(1))));
    const B=cadd(P1, cmul(t(0), cadd(smul(-4,t(-1)), cadd(cadd(smul(-2,omw),smul(2*dr,iw)),R(-4)))));
    const Cc=cadd(P2, cmul(t(-1), cadd(smul(6,t(-2)), cadd(csub(smul(6,omw),smul(4*dr,iw)),R(6)))));
    const Dd=cadd(P3, cmul(t(-2), cadd(smul(-4,t(-3)), cadd(cadd(smul(-6,omw),smul(2*dr,iw)),R(-4)))));
    const E=cadd(P4, cmul(t(-3), cadd(t(-4), cadd(smul(2,omw),R(1)))));
    let L10,L11,L12,L13;
    if(i>=3){const q=cdiv(E,L1p3);L10=A;L11=csub(B,cmul(L1p0,q));L12=csub(Cc,cmul(L1p1,q));L13=csub(Dd,cmul(L1p2,q));}
    else{L10=A;L11=B;L12=Cc;L13=Dd;}
    let al,be,ga;
    if(i>=2){const q=cdiv(L13,L2p2);al=L10;be=csub(L11,cmul(L2p0,q));ga=csub(L12,cmul(L2p1,q));}
    else{al=L10;be=L11;ga=L12;}
    if(i===0){f=cabs(be)<1e-30?tiny:be;Cl=f;Dl=[0,0];}
    else{const an=smul(-1,cmul(aPrev,ga)),bn=be;
      Dl=cadd(bn,cmul(an,Dl));if(cabs(Dl)<1e-30)Dl=tiny;
      Cl=cadd(bn,cdiv(an,Cl));if(cabs(Cl)<1e-30)Cl=tiny;
      Dl=cinv(Dl);f=cmul(f,cmul(Cl,Dl));}
    L1p0=L10;L1p1=L11;L1p2=L12;L1p3=L13;L2p0=al;L2p1=be;L2p2=ga;aPrev=al;
  }
  return Math.log(cabs(f)+1e-30);
}

// Known Leaver Schwarzschild ell=2 scalar QNMs (2M=1 units) from the paper
const modes=[[0.967288,-0.193518],[0.927701,-0.591208]];
for(const [wr,wi] of modes){
  const atMode=logCF3([wr,wi],2,120);
  const off  =logCF3([wr+0.3,wi-0.3],2,120);
  console.log(`Schw  omega=${wr}${wi}i : log|CF|=${atMode.toFixed(2)}  (off-mode ${off.toFixed(2)})`);
  assert(atMode < off-3, `Schwarzschild CF should dip sharply at the QNM ${wr}${wi}i`);
}
// Taub-NUT must reproduce the Schwarzschild mode as l -> 0 (C=6=q(q+1), q=2)
{ const atMode=logCF_TN([0.967288,-0.193518],1e-4,6,120);
  const off  =logCF_TN([1.3,-0.5],1e-4,6,120);
  console.log(`TN l->0 omega=0.9673-0.1935i : log|CF|=${atMode.toFixed(2)}  (off ${off.toFixed(2)})`);
  assert(atMode < off-2, 'Taub-NUT(l->0) must recover the Schwarzschild QNM'); }

console.log('all QNM checks passed');
