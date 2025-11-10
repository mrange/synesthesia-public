// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
const vec4 hsv2rgb_K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
vec3 hsv2rgb(vec3 c) {
  vec3 p = abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www);
  return c.z * mix(hsv2rgb_K.xxx, clamp(p - hsv2rgb_K.xxx, 0.0, 1.0), c.y);
}
// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
//  Macro version of above to enable compile-time constants
#define HSV2RGB(c)  (c.z * mix(hsv2rgb_K.xxx, clamp(abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www) - hsv2rgb_K.xxx, 0.0, 1.0), c.y))


#define ROT(a)      mat2(cos(a), sin(a), -sin(a), cos(a))
#define SCA(a)      vec2(sin(a), cos(a))

const float
  max_marches_1 = 90.
, tolerance_1   = 1e-3
, max_depth_1   = 100.
, normal_eps_1  = 1e-3
, TAU=2.*PI
, top_plane=9.
;

float length4(vec2 p) {
  p*=p;
  return sqrt(sqrt(dot(p,p)));
}

float length4(vec3 p) {
  p*=p;
  return sqrt(sqrt(dot(p,p)));
}

float segment_y(vec3 p) {
  float
    d0=length4(p)
  , d1=length4(p.xz)
  ;
  return p.y>0.?d0:d1;
}

float fbm(vec2 p) {
  const float
#ifdef KODELIFE
    off=.4
#else
    off=.2
#endif
  ;
  return (texture(t_fbm,2e-3*(p-.5*vec2(-1,1)*TIME)+.5).x-off)*11.;
}

float fbm2(vec2 p) {
  float 
    a=0.
  , aa=1.
  , i
  , d=0.
  ;
  for (i=0.;i<3.;++i) {
    a+=aa*sin(p.x);
    d+=aa;
    aa*=.5;
    p*=mat2(6,8,-8,6)/5.;
    p+=1.234;
  }
  return a/d;
}


float freq(float x, vec2 o) {
#ifdef KODELIFE
  return smoothstep(.0,.9,sin(TAU*x*TIME));
#else  
  float f=texture(syn_Spectrum,o.x+o.y*x).z;
  f*=f;
  f*=f;
  f*=3.5;
  return f;
#endif
}

// License: Unknown, author: Unknown, found: don't remember
float hash(vec2 co) {
  return fract(sin(dot(co.xy ,vec2(12.9898,58.233))) * 13758.5453);
}

vec2 hexagon(vec2 p, vec2 r) {
  p=p.yx;
  const vec3 
    k = 0.5*vec3(-sqrt(3.0), 1, sqrt(4.0/3.0))
  ;
  p = abs(p);
  p -= 2.0*min(dot(k.xy,p),0.0)*k.xy;
  vec2
    p0=p
  , p1=p
  ;
  p0-=vec2(clamp(p0.x, -k.z*r.x, k.z*r.x), r.x);
  p1.x=abs(p1.x);
  p1-=r.y*vec2(sqrt(1./3.),1);
  return vec2(length(p0)*sign(p0.y),length(p1));
}

vec3 hexagon(vec3 p, vec3 r) {
  vec2 
    d   = hexagon(p.xz,r.xy)
  , w0  = vec2(d.x,abs(p.y)-r.z)
  ;
  return vec3(
    min(max(w0.x,w0.y),0.)+length(max(w0,0.))
  , d
  )
  ;
}


// License: Unknown, author: Martijn Steinrucken, found: https://www.youtube.com/watch?v=VmrIDyYiJBA
vec2 hextile(inout vec2 p) {
  // See Art of Code: Hexagonal Tiling Explained!
  // https://www.youtube.com/watch?v=VmrIDyYiJBA
  const vec2 sz       = vec2(1, sqrt(3.));
  const vec2 hsz      = 0.5*sz;

  vec2 p1 = mod(p, sz)-hsz;
  vec2 p2 = mod(p - hsz, sz)-hsz;
  vec2 p3 = dot(p1, p1) < dot(p2, p2) ? p1 : p2;
  vec2 n = ((p3 - p + hsz)/sz);
  p = p3;

  n -= vec2(0.5);
  // Rounding to make hextile 0,0 well behaved
  return round(n*2.0)*0.5;
}

float nearest_hex_wall(vec2 p, vec2 rd) {
  const mat3x2
    PROJ = mat3x2(
        vec2(1 ,  0       )
      , vec2(.5, .8660253 )
      , vec2(.5, -.8660253)
      )
  ;

  vec3
    pp  = p*PROJ
  , prd = rd*PROJ
  , ird = 1./prd
  , dro = (sign(prd)*.5)-pp
  , dt  = dro*ird
  ;

  return min(min(dt.x, dt.y), dt.z);
}

// In
vec2
  g_hsrd
, g_ird
, g_rd
;


// Out
vec2 
  g_g
, g_C
;
vec3
  g_col
;
float df_1(vec3 p) {
  if(p.y>top_plane) {
    return p.y-top_plane+1.;
  }
  vec3
    p0=p
  , p1=p
  , d0
  ;
  vec2
    n
  , c=p.xz
  ;
  ;
  n=hextile(c);
  p0.xz=c;
  float
    h0=hash(n)
  , h1=fract(8667.*h0)
  , h
  , d1=p1.y
  , d
  , f
  , F
  ;
  
  vec3
    col=white_col
  ;

  h=fbm(n);
  F=fbm2(0.23*n);
  f=smoothstep(bouncy_islands.x,bouncy_islands.y,abs(F));
  h=mix(h,freq(h1,F>0.?red_freq:black_freq),f);
  col=mix(col,F>0.?red_col:black_col,f);
  g_col=col;

  // Cool bug
  // g_col=
  d0=hexagon(p0,vec3(.40,.45,h))-vec3(0.05,0,0);
  g_col=g_col;
  g_g=d0.yz;
  d=d0.x;

  float
    cd=1e-3+nearest_hex_wall(c,g_rd);
  ;
  d=min(d,cd);

  d=min(d,d1);
  return d;
}

float ray_march_1(vec3 ro, vec3 rd, float initz) {
  g_hsrd=sign(rd.xz)*.5;
  g_rd=rd.xz;
  g_ird=1./(rd.xz);

  float
    d
  , i
  , z=initz
  ;
  for (i=0.;i<max_marches_1;++i) {
    d=df_1(ro+rd*z);
    if(d<tolerance_1||z>max_depth_1) {
      break;
    }
    z+=d;
  }

  return z;
}

vec3 normal_1(vec3 p) {
  vec2
    e=vec2(normal_eps_1,0)
  ;
  return normalize(vec3(
      df_1(p+e.xyy)-df_1(p-e.xyy)
    , df_1(p+e.yxy)-df_1(p-e.yxy)
    , df_1(p+e.yyx)-df_1(p-e.yyx)
    ));
}

vec3 render1(vec3 ro, vec3 rd) {
  const vec3
    ld=normalize(vec3(-1,.5,1))
  , lc=HSV2RGB(vec3(.58,.0 ,1.))
  , sc=HSV2RGB(vec3(.58,.0,.2))
  , sky=HSV2RGB(vec3(.58,.0,1.))
  ;
  float
    i
  , z
  , Z
  , D
  , dl
  , ds
  , sl
  ;
  vec2
    g
  ;
  vec3
    col=sky
  , bcol
  , p
  , n
  , r
  ;
  D=(top_plane-ro.y)/rd.y;
  z=ray_march_1(ro,rd,D);
  g=g_g;
  bcol=g_col;
  p=ro+rd*z;
  n=normal_1(p);
  r=reflect(rd,n);
  Z=ray_march_1(p+1e-1*n,ld,1e-1);
  dl=max(dot(n,ld),0.);
  dl=sqrt(dl);
  sl=pow(max(dot(r,ld),0.),40.);
  ds=(1.+n.y)*.5;
  if (z<max_depth_1) {
    col=vec3(0);
    bcol=mix(.1*sign(bcol), bcol,mix(.25,1.,smoothstep(0.02,0.04,min(g.y,abs(g.x-.05+mix(0.04,0.03,n.y))))));
    bcol=mix(vec3(1),bcol,smoothstep(.07,.06,g.x));
    bcol*=min(
      dot(n,ld)>0.?1.:0.25
    , Z<max_depth_1?mix(1.,.25,exp(-.05*Z)):1.
    );
    col+=bcol;
  }

  col=mix(sky,col,exp(-.05*max(z-.5*max_depth_1,0.)));

  return col;
}

vec4 fpass0() {
  const vec3
  , Z=normalize(vec3(-1,-1,1))
  , X=normalize(cross(Z,vec3(0,1,0)))
  , Y=cross(X,Z)
  ;
  vec2
    p=2.*_uvc
  ;
  vec3
    ro=vec3(0,20.,TIME)
  , rd =normalize(-p.x*X+p.y*Y+2.*Z)
  , col;
  vec4
    pcol
  ;

  col=render1(ro,rd);
  col=sqrt(col);
  pcol=texture(syn_FinalPass,_uv);
  col=mix(col,pcol.xyz,.3);
  return vec4(col,1.);
}

vec4 renderMain() {
  return fpass0();
}
