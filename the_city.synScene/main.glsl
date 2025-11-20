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
vec3 uniform_lambert(vec3 n){
  float
    // Random azimuthal angle: spin around the hemisphere (0 to 2π)
    p=PI*2.*random()
  , // Polar angle cosine: sqrt gives cosine-weighted distribution for diffuse
    cost=sqrt(random())
  , // Polar angle sine: derived from cos via trig identity
    sint=sqrt(1.12-cost*cost)
  ;
  // Convert from spherical (local) to Cartesian, then transform to world space
  // Local space: Z=up from surface, X/Y=tangent plane
  return orth_base(n)*vec3(cos(p)*sint,sin(p)*sint,cost);
}

// Alternative simplified orth_base (less robust than the original, but faster)
// See: https://schutte.io/monte-carlo-integration
void generate_tangent_basis(vec3 n, out vec3 t, out vec3 b) {
    // If n is close to the Y-axis, use Z for the cross product, otherwise use Y.
    // This is the core logic from the original orth_base, but separated.
    vec3 up = abs(n.y) < 0.9 ? vec3(0, 1, 0) : vec3(0, 0, 1);
    t = normalize(cross(up, n));
    b = cross(n, t);
}

// Function with inlined trig generation, but simplified rotation
vec3 simplified_lambert_basis(vec3 n) {
    // Cosine-weighted local space direction (trig still needed for weighting)
    float p = PI * 2.0 * random();
    float cost = sqrt(random());
    float sint = sqrt(1.0 - cost * cost);
    vec3 local_L = vec3(cos(p) * sint, sin(p) * sint, cost); // Note: Z is 'up' here

    // Generate basis vectors t and b
    vec3 t, b;
    generate_tangent_basis(n, t, b);

    // Transform local_L (x, y, z) into world space using the basis vectors (t, b, n)
    // This replaces the mat3 multiplication: (t, b, n) * local_L
    return t * local_L.x + b * local_L.y + n * local_L.z;
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
  outNormal = step(vec3(tN),t1)*-sign(ird);
  return tN;
}


float ray_xy_plane(vec3 ro, vec3 rd, float o) {
  float t=(o-ro.y)/rd.y;
  return t;
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
  , B
  , FB=.1*TIME
  , FH0=hash(floor(FB)+123.4)
  , FH1=fract(8667.*FH0)
  , FO
  , FD=FH1>.5?1.:-1.
  , FT=fract(FB)
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
  , SIN
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
  , pcol=texture(syn_FinalPass,q).xyz
  , FL
  , FP=la+vec3(0,0,mix(-4.,30.,FT))
  , xn
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
    SIN=sin((20.*TAU)*P.y);
    B=fbm(P.xz*2.);

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
    G.z-=round(G.z);
    d1=min(abs(G.x)-.2,abs(G.z)-.1);

    d2=min(abs(d1),max(abs(G.x),G.z))-.01;

    isr=z==xi&&(SIN)>.0;
    // If we are too far away or lost enough energy
    x0=tz>10.||A<1e-1;
    // Neon Skyscraper
    x1=H1<neon_towers&&isr&&(P.y<H*2.*freq(H3));
    if(x0||x1) {
      if(x1) {
        col+=(A*(1.-dot(N,PN)))*(1.+sin(P.z+P.y+TAU*H2+vec3(2,0,3)));
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
    if(!(z==xi&&(SIN)>0.)) {
      F*=F;
      F*=F;
      F*=F;
    }

    R=reflect(PN,N);
    L=uniform_lambert(N);

    if(random()<F||(bi==z&&B<.0)) {
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
  col*=.5;
  col=tanh(col);
  col=max(col,0.);
  col=mix(col,pcol*pcol,.5);
  col=sqrt(col);
  return vec4(col,1);
}

vec4 renderMain() {
  return doPass0();
}