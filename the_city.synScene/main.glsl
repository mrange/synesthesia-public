#ifdef KODELIFE
const float
  neon_towers =.2
, motion_blur =.4
, height      =2.
;

const vec2
  reflection    =vec2(.0625,.75)
, neon_color    =vec2(2,3)
, path_control  =vec2(1,.33);
;
#endif

const float
  TAU         =2.*PI
, MISS        =-1000.
, SKY_LIMIT   =50.
;

vec3 path(float speed) {
  return vec3(.5+path_control.x*sin(path_control.y*speed),height,speed);
}

vec3 dpath(float speed) {
  return vec3(path_control.x*path_control.y*cos(path_control.y*speed),0,1);
}

vec3 ddpath(float speed) {
  return vec3(-path_control.x*path_control.y*path_control.y*sin(path_control.y*speed),0,0);
}

float freq(float x) {
#ifdef KODELIFE
  return exp(-2.*fract(TIME+x));
#else
  return height_mul*(textureLod(syn_Spectrum, x, 0.0).y-fft_distinct)/(1.0-fft_distinct);
#endif
}

// License: Unknown, author: Unknown, found: don't remember
float hash(float co) {
  return fract(sin(co*12.9898) * 13758.5453);
}

float g_seed;

// License: Unknown, author: 0b5vr, found: https://www.shadertoy.com/view/ss3SD8
float random(){
  float i = ++g_seed;
  return fract(sin((i)*114.514)*1919.810);
}

// License: Unknown, author: Unknown, found: don't remember
float hash(vec2 co) {
  return fract(sin(dot(co.xy ,vec2(12.9898,58.233))) * 13758.5453);
}

vec3 uniform_lambert(vec3 X, vec3 Y, vec3 Z){
  float
    p=TAU*random()
  , cost=sqrt(random())
  , sint=sqrt(1.-cost*cost)
  ;
  return cos(p)*sint*X+sin(p)*sint*Y+cost*Z;
}

vec3 noisy_ray_dir(vec2 p, vec3 X, vec3 Y, vec3 Z) {
  p += 1.41/RENDERSIZE.y*(-1.+2.*vec2(random(),random()));
  return normalize(-p.x*X+p.y*Y+2.*Z);
}

float iray_box_(vec3 ro, vec3 ird, vec3 boxSize, out vec3 NZ)  {
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

float iray_box(vec3 ro, vec3 ird, vec3 boxSize, out vec3 NZ) {
  // Branchless
  vec3
    n = ro * ird
  , k = abs(ird) * boxSize
  , t1 = -n - k
  , t2 = -n + k
  ;

  float
    tN = max( max( t1.x, t1.y ), t1.z )
  , tF = min( min( t2.x, t2.y ), t2.z )
  ;

  float hit = step(tN, tF) * step(0.0, tF);
  NZ = step(vec3(tN), t1) * -sign(ird);
  return mix(MISS, tN, hit);
}

// License: Unknown, author: Claude Brezinski, found: https://mathr.co.uk/blog/2017-09-06_approximating_hyperbolic_tangent.html
vec3 tanh_approx(vec3 x) {
  //  Found this somewhere on the interwebs
  //  return tanh(x);
  vec3 x2 = x*x;
  return clamp(x*(27.0 + x2)/(27.0+9.0*x2), -1.0, 1.0);
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
#ifdef KODELIFE
  , speed=.5*TIME
#endif
  , H0
  , H
  , bi
  , ti
  , zi=0.
  , z
  , MX
  , A=1.
  , xi
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
    ro=path(speed)
  , Z =normalize(dpath(speed)+vec3(0,-.6,0))
  , X =normalize(cross(Z,vec3(0,1,0)+2.*ddpath(speed)))
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
  , FP=vec3(.5,1.,speed+mix(-4.,30.,FT))
  , XX
  , XY
  , XZ
  , NEON
  ;

  vec4
    M
  , HH
  ;

  FL=FO*mix(vec3(.1,.1,1)*.5,vec3(1,.1,.5),step(hash(floor(TIME*19.)),.5));

  PP=ro;
  PN=noisy_ray_dir(p,X,Y,Z);
  IPN=1./PN;
  for(j=0.;j<110.;++j) {
    NN=floor(PP.xz+.5);
    CC=PP.xz-NN;
    S=(sign(PN.xz)*.5-CC)*IPN.xz;
    MX=min(S.x,S.y)+1e-2;
    ti=(3.-PP.y)*IPN.y;
    bi=(-PP.y)*IPN.y;
    H0=hash(NN);
    HH=fract(vec4(5711,6977,7577,8677)*H0);
    H=.5+.5*H0;
    xi=iray_box(PP-vec3(NN.x,H-.02,NN.y),IPN,vec3(.1+.1*HH.x,H,.2+.1*HH.y),XZ);
    z=1e3;
    if(ti>0.)           { z=ti; NZ=vec3(0,-1,0);}
    if(bi>0.&&bi<z)     { z=bi; NZ=vec3(0, 1,0);}
    if(xi>0.&&xi<z)     { z=xi; NZ=XZ;  }

    NX = NZ.yzx;
    NY = NZ.zxy;

    if(MX<z) {
      // Step to next cell
      PP=PP+PN*MX;
      continue;
    }

    if(zi==0.) {
      zi=(PP-ro).x*IPN.x-5e-2;
    }


    P=PP+PN*z;
    NEON=sin((20.*TAU)*P);
    G=P;
    G.xz-=FP.xz-vec2(-.1*FD,0);
    d0=dot(G,G);
    col+=(A*5e-3/max(d0,5e-4))*FL;

    G=P;
    G.xz-=vec2(.5+.1*FD,.5);
    d1=smoothstep(2.,0.,FD*(FP.z-P.z));
    col+=
       (A*1.5e-2*FO*(FD*P.z<FD*FP.z?1.:0.)*d1/max(length(G.xy),1e-3))
      *(.3*smoothstep(2.,0.,FP.z-P.z)+vec3(1,0.1,.25))
      ;

    G=P;
    G.xz-=.5;
    G.z-=floor(G.z+.5);

    G.x=G.x-floor(G.x+.5);
    d1=min(abs(G.x)-.2,abs(G.z)-.1);
    d2=min(abs(d1),max(abs(G.x),G.z))-.01;

    isr=z==xi&&(HH.z>0.5?NEON.x*NEON.z>0.0:NEON.y>0.);
    x0=z==ti||A<1e-1||length(P-ro)>20.;
    x1=HH.w<neon_towers&&isr&&(P.y<H*2.*freq(HH.x+HH.y));
    if(x0||x1) {
      if(x1) {
        col+=(A*(1.-dot(NZ,PN)))*(1.+sin(P.z+P.y+TAU*(HH.z+HH.w)+vec3(neon_color.x,0,neon_color.y)));
      }
      // Reset current pos and ray
      PN=noisy_ray_dir(p,X,Y,Z);
      PP=ro+PN*zi;
      IPN=1./PN;
      A=1.;
      ++n;
      continue;
    }

    F=1.+dot(PN,NZ);
    F*=F;
    F=mix(reflection.x,reflection.y,F);
    R=reflect(PN,NZ);
    L=uniform_lambert(NX,NY,NZ);

    if(isr) {
      PN=R;
      H0=F;
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
  col*=0.5;
  col=tanh_approx(col);
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