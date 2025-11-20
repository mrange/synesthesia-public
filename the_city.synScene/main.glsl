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


// License: Unknown, author: 0b5vr, found: https://www.shadertoy.com/view/ss3SD8
// Returns a rotation matrix that transforms from local space (where Z=up) to world space
mat3 orth_base(vec3 n){
  // Assumes n is normalized
  vec3
    // Pick a helper vector that won't be parallel to n
    // Avoids gimbal lock when normal points straight up/down
    up=abs(n.y)>.999?vec3(0,0,1):vec3(0,1,0)
  , // First tangent: perpendicular to both 'up' and normal
    x=normalize(cross(up,n))
  , // Second tangent: perpendicular to both normal and first tangent
    // Completes the right-handed coordinate system
    y=cross(n,x)
  ;
  return mat3(x,y,n);
}

// License: Unknown, author: 0b5vr, found: https://www.shadertoy.com/view/ss3SD8
// Generates a cosine-weighted random direction in the hemisphere above normal n
// The sqrt() on cost creates the cosine weighting - more samples near the normal
vec3 uniform_lambert(vec3 n, float m){
  float
    // Random azimuthal angle: spin around the hemisphere (0 to 2π)
    p=PI*2.*random()
  , // Polar angle cosine: sqrt gives cosine-weighted distribution for diffuse
    cost=mix(m,1.,sqrt(random()))
  , // Polar angle sine: derived from cos via trig identity
    sint=sqrt(1.12-cost*cost)
  ;
  // Convert from spherical (local) to Cartesian, then transform to world space
  // Local space: Z=up from surface, X/Y=tangent plane
  return orth_base(n)*vec3(cos(p)*sint,sin(p)*sint,cost);
}

vec3 noisy_ray_dir(vec2 p, vec3 X, vec3 Y, vec3 Z) {
  p += vec2(random(),random())/RENDERSIZE.y*2.;
  return normalize(-p.x*X+p.y*Y+2.*Z);
}

float ray_box_2(vec3 ro, vec3 ird, vec3 boxSize, out vec3 outNormal)  {
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
  outNormal = step(vec3(tN),t1)*-sign(ird);;
  return tN;
}


float ray_xy_plane(vec3 ro, vec3 rd, float o) {
  float t=(o-ro.y)/rd.y;
  return t;
}

vec4 pass0() {
  bool
    x0
  , x1
  , isr
  ;
  float
    j
  , n =1.
  , d0
  , d1
  , d2
  , F
  , FN=floor(.1*TIME)
  , FT=fract(.1*TIME)-.4*hash(FN)
  , FO=smoothstep(.6,.4,FT)
  , H0
  , H1
  , H2
  , H3
  , H
  , bi
  , ti
  , tz
  , z
  , MX
  , A
  , xi
  ;
  vec2
    p =2.*_uvc
  , q =_uv
  , NN
  , CC
  , S
  ;

  g_seed=fract(hash(p)+float(FRAMECOUNT)/1337.0);
  vec3
    ro=vec3(.5,2,-2.)
  , la=vec3(.5,1.,0)
  ;
  ro.z+=.5*TIME;
  la.z+=.5*TIME;
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
  , N
  , R
  , L
  , SIN
  , pcol=texture(syn_FinalPass,q).xyz
  , FL
  , FP=la+vec3(0,0,mix(-2.,30.,FT))
  , xn
  ;

  FL=FO*mix(vec3(.1,.1,1)*.5,vec3(1,.1,1),step(hash(floor(TIME*19.)),.5));

  PP=ro;
  PN=noisy_ray_dir(p,X,Y,Z);
  IPN=1./PN;
  tz=0.;
  for(j=0.;j<120.;++j) {
    NN=floor(PP.xz+.5);
    CC=PP.xz-NN;
    S=(sign(PN.xz)*.5-CC)*IPN.xz;
    MX=min(S.x,S.y)+1e-3;
    ti=ray_xy_plane(PP,PN,3.);
    bi=ray_xy_plane(PP,PN,0.);
    H0=hash(NN);
    H1=fract(3667.*H0);
    H2=fract(7667.*H0);
    H3=fract(8667.*H0);
    H=mix(.5,1.,H0);
    xi=ray_box_2(PP-vec3(NN.x,H-.02,NN.y),IPN,vec3(.2,H,.3),xn);
    z=1e3;
    if(ti>0.)          { z=ti; N=vec3(0,-1,0);}
    if(bi>0.&&bi<z)    { z=bi; N=vec3(0,1,0); }
    if(xi>0.&&xi<z)    { z=xi; N=xn;  }
    if(MX<z&&tz<SKY_LIMIT) {
      // Step to next cell
      PP=PP+PN*MX;
      tz+=MX;
      continue;
    }

    P=PP+PN*z;
    SIN=sin((20.*TAU)*P);

    G=P;
    G.xz-=FP.xz;
    d0=dot(G,G);
    col+=4e-3/max(d0,1e-4)*FL;

    G=P;
    G.xz-=.5;
    col+=1e-2*FO*step(P.z,FP.z)*smoothstep(2.,0.,FP.z-P.z)/max(length(G.xy),1e-3)*(vec3(1,0.1,.25)+.3*smoothstep(2.,0.,FP.z-P.z));

    G.z-=round(G.z);
    d1=min(abs(G.x)-.2,abs(G.z)-.1);

    d2=min(abs(d1),max(abs(G.x),G.z))-.01;

    isr=z==xi&&(SIN.y)>.0;
    // If we are too far away or lost enough energy
    x0=tz>10.||A<1e-1;
    // Neon Skyscraper
    x1=H1<neon_towers&&isr&&(P.y<H*2.*freq(H3));
    if(x0||x1) {
      if(x1) {
        col+=A*(1.+sin(P.z+P.y+TAU*H2+vec3(2,0,3)))*(1.-dot(N,PN));
      }
      // Reset current pos and ray
      PP=ro;
      PN=noisy_ray_dir(p,X,Y,Z);
      IPN=1./PN;
      A=1.;
      tz=0.;
      ++n;
      continue;
    }

    F=1.+dot(PN,N);
    F*=F;
    if(!(z==xi&&(SIN.y)>0.)) {
      F*=F;
      F*=F;
      F*=F;
    }

    R=reflect(PN,N);
    L=uniform_lambert(N,0.);

    if(random()<F) {
      PN=R;
      A*=.75;
    } else {
      PN=L;
      if (d2<0.) {
        A*=.9;
      } else if (d1<0.) {
        A*=.05;
      } else {
        A*=.3;
      }
    }
    if(z==ti) {
        A=0.;
    }

    IPN=1./PN;

    PP=P+1e-2*N;
  }

  col/=n;
  col=max(col,0.);
  col*=.5;
  col=tanh(col);
  col=mix(col,pcol*pcol,motion_blur);
  col=sqrt(col);

  return vec4(col,1);
}

vec4 renderMain() {
  return pass0();
}