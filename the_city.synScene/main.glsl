//#define WATER

#ifdef KODELIFE
const float
  neon_towers=.2
, motion_blur=.5
;
#endif

const float
  TAU=2.*PI
, MISS=-1000.
, SKY_LIMIT=50.
;


float freq(float x) {
#ifdef KODELIFE
  return exp(-2.*fract(TIME+x));
#else
  return (texture(syn_Spectrum, x).y);
#endif
}

float g_seed;

// License: Unknown, author: 0b5vr, found: https://www.shadertoy.com/view/ss3SD8
float random(){
  float i = ++g_seed;
  return fract(sin((i)*114.514)*1919.810);
}

// License: Unknown, author: Unknown, found: don't remember
float hash(float co) {
  return fract(sin(co*12.9898) * 13758.5453);
}

// License: Unknown, author: Unknown, found: don't remember
float hash(vec2 co) {
  return fract(sin(dot(co.xy ,vec2(12.9898,58.233))) * 13758.5453);
}

vec3 uniform_lambert(vec3 X, vec3 Y, vec3 Z){
  float
    p=PI*2.*random()
  , cost=sqrt(random())
  , sint=sqrt(1.12-cost*cost)
  ;
  return cos(p)*sint*X+sin(p)*sint*Y+cost*Z;
}

vec3 noisy_ray_dir(vec2 p, vec3 X, vec3 Y, vec3 Z) {
  p += sqrt(2.)/RENDERSIZE.y*(-1.+2.*vec2(random(),random()));
  return normalize(-p.x*X+p.y*Y+2.*Z);
}

float iray_box(vec3 ro, vec3 ird, vec3 boxSize, out vec3 NZ)  {
  vec3
    n = ird*ro
  , k = abs(ird)*boxSize
  , t1 = -n - k
  , t2 = -n + k
  ;
  float
    tN = max( max( t1.x, t1.y ), t1.z )
  , tF = min( min( t2.x, t2.y ), t2.z )
  ;
  if( tN>tF || tF<0.0) return MISS;
  NZ = step(vec3(tN),t1)*-sign(ird);
  return tN;
}


float fbm(vec2 p) {
  const mat2
    PP = mat2(1.2, 1.6, -1.6, 1.2)
  ;

  float h = dot(sin(p), cos(p*1.618).yx);
  p *= PP;
  h += 0.5 * dot(sin(p), cos(p*1.618).yx);

  return h;
}



vec4 doPass0() {
  bool
    isr
  , x0
  , x1
  ;
  float
    j
  , n =1.
  , d0
  , d1
  , d2
  , F
  , FB=.1*TIME
  , FH0=hash(floor(FB)+123.4)
  , FH1=fract(8667.*FH0)
  , FO
  , FD=FH1>.5?1.:-1.
  , FT=fract(FB)
  , SPD=.5*TIME
  , H0
  , H1
  , H2
  , H3
  , H
  , bi
  , si
  , ti
  , tz
  , z
  , MX
  , A
  , xi
  , SIN
  , W
  ;
  FT=FD>0.?FT:1.-FT;
  FO=smoothstep(.6,.4,FT);

  vec2
    p =2.*_uvc
  , q =_uv
  , NN
  , CC
  , S
  ;

  g_seed=fract(hash(p)+float(FRAMECOUNT)/1337.0);
  vec3
    ro=vec3(.5,2,SPD-2.)
  , la=vec3(.5,1.,SPD)
  ;
  vec3
    Z =normalize(la-ro)
  , X =normalize(cross(Z,vec3(0,1,0)))
  , Y =cross(X,Z)
  , col=vec3(0)
  , P
  , G
  , PP
  , PN
  , IPN
  , NX
  , NY
  , NZ
  , R
  , L
  , pcol=texture(syn_FinalPass,q).xyz
  , FL
  , FP=la+vec3(0,0,mix(-4.,30.,FT))
  , XX
  , XY
  , XZ
  ;
  vec4
    M
  ;

  FL=FO*mix(vec3(.1,.1,1)*.5,vec3(1,.1,.5),step(hash(floor(TIME*19.)),.5));

  PP=ro;
  PN=noisy_ray_dir(p,X,Y,Z);
  IPN=1./PN;
  tz=0.;
  for(j=0.;j<120.;++j) {
    NN=floor(PP.xz+.5);
    CC=PP.xz-NN;
    S=(sign(PN.xz)*.5-CC)*IPN.xz;
    MX=min(S.x,S.y)+1e-3;
    ti=(3.-PP.y)*IPN.y;
    si=(IPN.x>0.?(.8-PP.x):(.2-PP.x))*IPN.x;
    bi=(-PP.y)*IPN.y;
    H0=hash(NN);
    H1=fract(6977.*H0);
    H2=fract(7577.*H0);
    H3=fract(8677.*H0);
    H=.5+.5*H0;
    xi=iray_box(PP-vec3(NN.x,H-.02,NN.y),IPN,vec3(.1+.1*H1,H,.2+.1*H),XZ);
    z=1e3;
    if(ti>0.)           { z=ti; NZ=vec3(0,-1,0);}
    if(bi>0.&&bi<z)     { z=bi; NZ=vec3(0, 1,0);}
    x0=ro!=PP||si>0.&&si<z;
    if(x0&&xi>0.&&xi<z) { z=xi; NZ=XZ;  }

    NX = NZ.yzx;
    NY = NZ.zxy;

    if(x0&&MX<z) {
      // Step to next cell
      PP=PP+PN*MX;
      tz+=MX;
      continue;
    }

    tz+=z;

    P=PP+PN*z;
    SIN=sin((20.*TAU)*P.y);
#ifdef WATER
    W=fbm(P.xz*2.);
#endif

    G=P;
    G.xz-=FP.xz-vec2(-.1*FD,0);
    d0=dot(G,G);
    col+=(4e-3/max(d0,5e-4))*FL;

    G=P;
    G.xz-=vec2(.5+.1*FD,.5);
    d1=smoothstep(2.,0.,FD*(FP.z-P.z));
    col+=
       (1e-2*FO*(FD*P.z<FD*FP.z?1.:0.)*d1/max(length(G.xy),1e-3))
      *(.3*smoothstep(2.,0.,FP.z-P.z)+vec3(1,0.1,.25))
      ;

    G=P;
    G.xz-=.5;
    G.z-=floor(G.z+.5);
    d1=min(abs(G.x)-.2,abs(G.z)-.1);

    d2=min(abs(d1),max(abs(G.x),G.z))-.01;

    isr=z==xi&&(SIN)>.0;
    // If we are too far away or lost enough energy
    x0=z==ti||A<1e-1||length(P-ro)>20.;
    // Neon Skyscraper
    x1=H1<neon_towers&&isr&&(P.y<H*2.*freq(H3));
    if(x0||x1) {
      A*=step(tz,100.);
      if(!x0&&x1) {
        col+=(A*(1.-dot(NZ,PN)))*(1.+sin(P.z+P.y+TAU*H2+vec3(2,0,3)));
      }
      // Reset current pos and ray
      PP=ro;
      PN=noisy_ray_dir(p,X,Y,Z);
      IPN=1./PN;
      A=1.;
      ++n;
      continue;
    }

    F=1.+dot(PN,NZ);
    F*=F;
    if(!(z==xi&&(SIN)>0.)) {
      F*=F;
      F*=F;
      F*=F;
    }

    R=reflect(PN,NZ);
    L=uniform_lambert(NX,NY,NZ);

#ifdef WATER
    if(random()<F||(bi==z&&W<.0)) {
#else
    if(random()<F) {
#endif
      PN=R;
      H0=.75;
    } else {
      PN=L;
      H0=d2<0.
        ? .9
        : d1<0.
          ? .05
          : .3
      ;
    }
    A*=H0;

    IPN=1./PN;

    PP=P+1e-2*NZ;
  }

  col/=n;
  col*=.5;
  col=tanh(col);
  col=max(col,0.);
  col=mix(col,pcol*pcol,motion_blur);
  col=sqrt(col);
#ifndef KODELIFE
  M=_loadMedia();
  col=mix(col,M.xyz,M.w*media_transparency*media_multiplier);
#endif
  return vec4(col,1);
}

vec4 renderMain() {
  return doPass0();
}