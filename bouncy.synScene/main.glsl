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
    off=.25
#endif
  ;
  return (texture(t_fbm,0.002*(p-TIME)+.5).x-off)*12.;
}
float freq(float x) {
#ifdef KODELIFE
  return smoothstep(.0,.9,sin(TAU*x*TIME));
#else  
  return texture(syn_Spectrum,x).z;
#endif
}

// License: Unknown, author: Unknown, found: don't remember
float hash(vec2 co) {
  return fract(sin(dot(co.xy ,vec2(12.9898,58.233))) * 13758.5453);
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
float hexagon(vec2 p, float r) {
  p=p.yx;
  const vec3 k = 0.5*vec3(-sqrt(3.0), 1.0, sqrt(4.0/3.0));
  p = abs(p);
  p -= 2.0*min(dot(k.xy,p),0.0)*k.xy;
  p -= vec2(clamp(p.x, -k.z*r, k.z*r), r);
  return length(p)*sign(p.y);
}

float hexagon(vec3 p, vec2 r) {
  float d = hexagon(p.xz,r.x);
  vec2 w = vec2(d,abs(p.y)-r.y);
  return min(max(w.x,w.y),0.)+length4(max(w,0.));
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

float nearest_hex_wall_(vec2 p, vec2 rd) {
  const vec2
    N1 = vec2(1 , .0)
  , N2 = vec2(.5, sqrt(3.)/2.)
  , N3 = vec2(.5,-sqrt(3.)/2.)
  ;

  vec3
    pp  = vec3(dot(p, N1), dot(p, N2), dot(p, N3))
  , prd = vec3(dot(rd, N1), dot(rd, N2), dot(rd, N3))
  , ird = 1. / prd
  , dro = (sign(prd) * 0.5) - pp
  , dt  = dro * ird
  ;

  return min(min(dt.x, dt.y), dt.z);
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

vec2
  g_hsrd
, g_ird
, g_rd
;

const float
  top_plane=9.
;

float df_1(vec3 p) {
  if(p.y>top_plane) {
    return p.y-top_plane+1.;
  }
  vec3
    p0=p
  , p1=p
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
  , h2=fract(3667.*h0)
  , h=h0>.2
    ?fbm(n)
    :4.*(freq(h1)-.5)
  , d0=hexagon(p0,vec2(.3,h))-.15
  , d1=p1.y
  , d=d0
  ;

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
    ld=normalize(vec3(1,.25,1))
  , lc=HSV2RGB(vec3(.58,.25 ,1.5))
  , sc=HSV2RGB(vec3(.58,.5,.1))
  , sky=HSV2RGB(vec3(.58,.25,.9))
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
  vec3
    col=sky
  , p
  , n
  , r
  ;
  D=(top_plane-ro.y)/rd.y;
  z=ray_march_1(ro,rd,D);
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
    col+=dl*lc*mix(1.0,0.1,exp(-.125*Z*Z));
    col+=ds*sc;
    col+=sl;
  }

  col=mix(sky,col,exp(-.05*max(z-.5*max_depth_1,0.)));

  return col;
}

vec4 fpass0() {
  const vec3
  , Z=normalize(vec3(0,-1,2))
  , X=normalize(cross(Z,vec3(0,1,0)))
  , Y=cross(X,Z)
  ;
  vec2
    p=2.*_uvc
  ;
  vec3
    ro=vec3(0,14.,.5*TIME)
  , rd =normalize(-p.x*X+p.y*Y+2.*Z)
  , col;


  col=render1(ro,rd);
//  col=tanh(col);
  col=sqrt(col);
  return vec4(col,1.);
}

#ifdef KODELIFE

float segment(vec2 p, vec2 a, vec2 b ) {
  vec2 pa = p-a, ba = b-a;
  float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0);
  return length( pa - ba*h );
}

vec4 fpass1() {
  vec2
    p=2.*_uvc
  , mp0=-1.+2.*abs(mouse.xy)
  , mp1=-1.+2.*abs(mouse.zw)
  , c=p
  , n
  , ro
  , rd
  , p1
  ;
  mp0.x*=RENDERSIZE.x/RENDERSIZE.y;
  mp1.x*=RENDERSIZE.x/RENDERSIZE.y;
  p.y*=-1.;
  n=hextile(c);

  ro=mp1;
  rd=normalize(mp0-mp1);

  float
    aa=sqrt(2.)/RENDERSIZE.y
  , dray=segment(p,mp0,mp1)-0.0075
  , dhex=abs(hexagon(c,.5))-aa
  , hd=nearest_hex_wall(ro,rd)
  , dhd=segment(p,ro,ro+rd*hd)-0.0075
  , dh0=hexagon(mp0,.5)
  , dh1=hexagon(mp1,.5)
  ;

  vec3
    col=vec3(0)
  ;

  col=mix(col,vec3(.5),smoothstep(aa,-aa,dhex));
  if(dh0<0.&&dh1<0.) {
    col=mix(col,vec3(0.25,0,1),smoothstep(aa,-aa,dhd));
  }
  col=mix(col,vec3(1,0,.25),smoothstep(aa,-aa,dray));
  col=sqrt(col);
  return vec4(col,1);
}

#endif

vec4 renderMain() {
  return fpass0();
}
