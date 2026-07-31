#define RGB(r,g,b)  pow((vec3(float(r),float(g),float(b))/255.),vec3(2))
#define SC(a)       vec2(sin(a), cos(a))
#define ROT(a)      mat2(cos(a), sin(a), -sin(a), cos(a))

const float
  TAU=2.*PI
, BPM=120.
, TF=120./60.
, BPMF=TF*TAU
, MAX_ITER=77.
, TOL     =1e-4
, EPS     =5e-3
, MAX_DIST=9.
;

const vec3
  GOLD=RGB(247,228,126)
, BLUE=RGB(30,34,120)
;

float beat() {
#ifdef KODELIFE
  return exp(-2.*fract(TIME*TF));
#else
  return smoothstep(lo_beat,1.01,syn_BassHits);
#endif

}

float isotri(vec2 p, vec2 q) {
  p.x = abs(p.x);
  vec2 a = p - q*clamp( dot(p,q)/dot(q,q), 0.0, 1.0 );
  vec2 b = p - q*vec2( clamp( p.x/q.x, 0.0, 1.0 ), 1.0 );
  float s = -sign( q.y );
  vec2 d = min( vec2( dot(a,a), s*(p.x*q.y-p.y*q.x) ),
                vec2( dot(b,b), s*(p.y-q.y)  ));
  return -sqrt(d.x)*sign(d.y);
}

float bend(vec2 p, vec2 sc, float r) {
  p.x = abs(p.x);
  return sc.y*p.x>sc.x*p.y ? dot(p-sc*r,sc) : length(p)-r;
}

// License: MIT, author: Inigo Quilez, found: https://www.iquilezles.org/www/articles/smin/smin.htm
float pmin(float a, float b, float k) {
  float h = clamp(0.5+0.5*(b-a)/k, 0.0, 1.0);
  return mix(b, a, h) - k*h*(1.0-h);
}

// License: CC0, author: Mårten Rånge, found: https://github.com/mrange/glsl-snippets
float pmax(float a, float b, float k) {
  return -pmin(-a, -b, k);
}


float decode(vec4 c) {
  const vec4 
    decoder=vec4(1,1./255.,1./65025.,0)
  ;
  const float 
    halfoff=+.5/(255.)+.5/(65025.)
  ;
  
  return dot(c,decoder)-halfoff;
}

float dlinTex(sampler2D tex, vec2 p) {
  vec2 
    sz=vec2(textureSize(tex, 0))
  , tx=clamp(p*sz-.5,vec2(0), sz-2.)
  , ntx=floor(tx)
  , ftx=fract(tx)
  ;
  ivec2
    itx=ivec2(ntx)
  ;
  
  vec4
    c00 = texelFetch(tex, itx+ivec2(0,0), 0)
  , c01 = texelFetch(tex, itx+ivec2(0,1), 0)
  , c10 = texelFetch(tex, itx+ivec2(1,0), 0)
  , c11 = texelFetch(tex, itx+ivec2(1,1), 0)
  ;

  float
    d00 = decode(c00)
  , d01 = decode(c01)
  , d10 = decode(c10)
  , d11 = decode(c11)
  , d   = mix(
      mix(d00, d01, ftx.y)
    , mix(d10, d11, ftx.y)
    , ftx.x
    )
  ;


  return -1.+2.*d;
}

float L4(vec2 p) {
  return sqrt(length(p*p));
}

float L4(vec3 p) {
  return sqrt(length(p*p));
}

float dssg(sampler2D tex, vec2 p) {
  p.x-=.006;
  p.y-=.005;
  const float D=4.;
  p*=.5;
  p+=.5;
#ifdef KODELIFE  
  p.y=1.-p.y;
#endif  
  return D*dlinTex(tex, p);
}

float dssg1(vec2 p) {
  return dssg(t_dssg1, p);
}

float dssg1(vec3 p, vec3 D) {
  float 
    d=dssg1(p.xy)+D.y
  ;
  vec2 
    w=vec2(d, abs(p.z)-D.x+D.y)
  ;
  return min(max(w.x,w.y),0.)+L4(max(w,0.))-D.z;
}

float dssg2(vec2 p) {
  return dssg(t_dssg2, p);
}

float dssg2(vec3 p, vec3 D) {
  float 
    d=dssg2(p.xy)+D.y
  ;
  vec2 
    w=vec2(d, abs(p.z)-D.x+D.y)
  ;
  return min(max(w.x,w.y),0.)+L4(max(w,0.))-D.z;
}

float dssg0(vec2 p) {
  const float
    a2=51.4
  , r2=.0675
  ;
  const mat2
    rot2=ROT(-radians(a2))
  ;
  const vec2
    sc2=SC(radians(a2))
  ;
  p.x=abs(p.x);
  vec2
    p0=p
  , p1=p
  , p2=p
  ;
  const float off=.54;
  p0.x-=(p0.y>0. && p0.x>.5*off)?off:0.;
  p0.y=-abs(p0.y);
  p0.y+=.77;
  p1.y=abs(p1.y)-.67;
  p2-=vec2(0.84-r2,-0.67+r2);
  p2.y=-p2.y ;
  p2*=rot2;
  float
    d0=isotri(p0,vec2(.45*1.5,.75))
  , d1=p1.y
  , d2=bend(p2, sc2, r2)
  , d=d1
  ;

//  if(p.y>.0) d0=max(d0,d2);  
  d=pmin(d,d0,.05);
  if(!(p.y<-.3&&p.x<.3))
    d=max(d,d2);
  return d;
}

float dssg0(vec3 p, vec3 D) {
  float 
    d=dssg0(p.xy)+D.y
  ;
  vec2 
    w=vec2(d, abs(p.z)-D.x+D.y)
  ;
  return min(max(w.x,w.y),0.)+L4(max(w,0.))-D.z;
}

vec3 flash(vec3 col, vec2 p, float aa) {
  float 
    Z=mix(1.,1.1,beat())
  , BEAT=beat()
  ;
  vec2
    dp=p/Z
  ;
  dp*=ROT(-radians(logo_rot));
  float
    d1=dssg1(dp)*Z
  , d2=dssg2(dp)*Z
  , d3=dssg0(dp)*Z
  , d4=dssg0(dp-.05*vec2(1,-1))*Z
  , aa2=aa*mix(1.,80.,BEAT*BEAT)*Z
  ;
  ;
  dp/=Z;  
  float 
  ;
  col=mix(col, BLUE, smoothstep(aa,-aa,d4*.01));
  col=mix(col, BLUE, smoothstep(aa,-aa,d3));
  col=mix(col, GOLD, smoothstep(aa,-aa,abs(d3)-.005*Z));
  col=mix(col, GOLD, smoothstep(aa,-aa,d2));
  col=mix(col, GOLD, smoothstep(aa2,-aa2,d1));
  col.y-=3.*d2*smoothstep(aa,-aa,d3-.005*Z)/Z;
  col+=(1.+sin(p.y-TIME+10.*d3+vec3(0,1,2)))*sin(d1*300.)*smoothstep(.2*BEAT,.0,abs(d1));
  return col;
}


mat2 g_rot0;
mat2 g_rot1;
vec4 g_hit;

float df(vec3 p) {
//  return length(p)-1.;
  p.zy*=g_rot0;
  p.yx*=g_rot0;
  p.xz*=g_rot1;
  const float Z=1.2;
  p.x*=p.z<0.?1.:-1.;
  float
    d0=dssg0((p-vec3(0,0,0.0))/Z,vec3(.04,.0,.0))*Z-.0
  , d1=dssg1(p,vec3(.05,.01,.015))
  , d2=dssg2(p,vec3(.04,.03,.035))
  , d3=abs(p.z)-.01
  , d =min(d1,d2)
  ;
  d3=max(d0+.01,d3);
  d0=pmax(d0, .02-d,.02);
  d3=min(d0,d3);
  d=min(d,d3);
  g_hit = vec4(d,d0,d1,d2);
  return d;
}

// License: MIT, author: IQ, found: https://www.shadertoy.com/view/4ds3zn
vec3 nf(vec3 p) {
  const vec2 
    e=vec2(EPS,-EPS)
  ;
  return normalize(
    e.xyy*df(p+e.xyy)
  + e.yyx*df(p+e.yyx)
  + e.yxy*df(p+e.yxy)
  + e.xxx*df(p+e.xxx)
  );
}


float raymarch(vec3 ro, vec3 rd, out float oi) {
  float i,d,z=0.;
  
  for(i=0.;i<77.;++i) {
    d=df(z*rd+ro);
    if(d<TOL||z>MAX_DIST)
      break;
    z+=d;
  }
  
  oi=i;
  
  return min(MAX_DIST, z);
}

vec3 threed(vec3 col, vec2 p, float aa) {
  vec3 
    P
  , RD=normalize(vec3(p,2))
  , RO=vec3(0,0,beat()*sin(TIME+dot(p,p))-logo_distance)
  , N
  , R
  , LD0=normalize(vec3(1,1,-3))
  , LD1=normalize(vec3(-2,1,-3))
  , mat=vec3(0)
  , hcol
  ;
  float
    z
  , i
  , d0
  , d1
  , s0
  , s1
  ;
  vec4
    hit
  ;
  
  g_rot0=ROT(radians(172.+logo_rot));
  vec2 S=sin(p+TIME);
  g_rot1=ROT(logo_speed+.5*logo_distort*p.y*p.x+.2*S.x+.3*S.y*S.x+p.y*sin(logo_speed*logo_distort));
  z=raymarch(RO,RD, i);
  hit=g_hit;
  P=z*RD+RO;
  N=nf(P);
  R=reflect(RD,N);
  d0=
    .6*(max(dot(N,LD0),0.))
  + .1*(.5*(1.-N.y))
  ;
  d1=
    .6*(max(dot(N,LD1),0.))
  + .1*(.5*(1.-N.y))
  ;
  d0*=d0;
  d1*=d1;
  d0*=.1;
  d1*=.1;
  s0=pow(max(dot(R,LD0),0.),40.);
  s1=pow(max(dot(R,LD1),0.),40.);

  if(z<MAX_DIST) {
    if (g_hit.x==g_hit.y) {
      mat=(BLUE);
    } else if (g_hit.x==g_hit.z) {
      mat=GOLD;
    } else if (g_hit.x==g_hit.w) {
      mat=vec3(1);
    } else {
      mat=vec3(0);
    }
    
    hcol=
      s0*mat*vec3(1,2,3)
    + s1*mat*vec3(3,2,1)
    + d0*mat*vec3(1,2,3)
    + d1*mat*vec3(3,2,1)
    ;
    hcol*=2.;
    col=hcol;
  }
  return col;
}

vec4 renderMain() {

  vec2
    p=2.*_uvc
  ;
  
  float
    aa=sqrt(2.)/RENDERSIZE.y
  , l2=dot(vec2(p.x*p.y, dot(p,p)),spectrum_shape)
  , n
  , c
  , h=0.
  ;

  vec3
    col=vec3(0)
    ;
  vec4
    spec=texture(syn_Spectrum,mix(.1,.5,.5+.5*sin(-TIME+l2)))
  ;
  
  col=
    2e-2*BLUE/max(dot(p,p),3e-1)
  ;
  
  if(spectrum_level<1.)
    col+=(1.+sin(vec3(2,1,0)+l2-.707*TIME))*smoothstep(.0,.1,mix(spec.x, spec.y, spectrum_picker)-spectrum_level);
  //col*=0.;


  n=clamp(floor((p.y+beat_bar_size*.5)/beat_bar_size+.5),-1.,2.);
  c=p.y-beat_bar_size*n+beat_bar_size*.5;
  
  switch (int(n))
  {
  case -1:
    h=syn_BassHits;
    break;
  case 0:
    h=syn_MidHits;
    break;
  case 1:
    h=syn_MidHighHits;
    break;
  case 2:
    h=syn_HighHits;
    break;
  }
  
  if(beat_bar_size>0.&&beat_bar_size*.4>abs(c))
    col+=(1.+sin(n+vec3(2,1,0)+TIME+p.y+.5*p.x))*smoothstep(.2,1.,h);
  
  if(effect<.5)
    col=flash(col,p,aa);
  else
    col=threed(col,p,aa);
  col=tanh(col);
  col=sqrt(col)-.03;
  return vec4(col,1);
}

