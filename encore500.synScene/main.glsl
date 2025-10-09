#define RGB(a)      (vec3(a*a)/(255.*255.))

const float
  TAU         = 2.*PI
, MaxDistance = 30.
, MaxIter     = 65.
, MinDistance = 1e-3
, NormalEps   = 1e-3
;

const vec2
  nstripe=normalize(vec2(1,1.06))
;

const vec3
  amiga_blue  =RGB(ivec3(0x00,0x34,0x9A))
, amiga_blue2 =RGB(ivec3(0x00,0x49,0xCC))
, amiga_orange=RGB(ivec3(0xFF,0x80,0x18))
, amiga_white =RGB(ivec3(0xE4,0xE4,0xE4))
, amiga_black =RGB(ivec3(0x11,0x12,0x13))
, bb_dim      =vec3(1e3,1e3,10)
, light_dir   =normalize(vec3(-1,1,-1))
;


float length8(vec3 p) {
  p*=p;
  p*=p;
  return pow(dot(p,p),1./8.);
}

void rot(inout vec2 p, float a) {
  float
    c=cos(a)
  , s=sin(a)
  ;
  p=vec2(c*p.x+s*p.y,-s*p.x+c*p.y);
}

float hash(vec3 r)  {
  return fract(sin(dot(r.xy,vec2(1.38984*sin(r.z),1.13233*cos(r.z))))*653758.5453);
}

vec3 hash3(vec3 p) {
  // Claude suggested this. Not sure if it's any good, probably claude stolen from somewhere
  p = vec3(
    dot(p, vec3(127.1, 311.7, 74.7)),
    dot(p, vec3(269.5, 183.3, 246.1)),
    dot(p, vec3(113.5, 271.9, 124.6))
  );
  return fract(sin(p) * 43758.5453123);
}
/*
vec3 hash3(vec3 r) {
  return fract(sin(dot(r.xy,vec2(1.38984*sin(r.z),1.13233*cos(r.z))))*653758.5453);
}
*/


vec3 shash3(vec3 p) {
  return -1.+2.*hash3(p);
}

float box(vec2 p, vec2 b) {
  vec2 d = abs(p)-b;
  return length(max(d,0.)) + min(max(d.x,d.y),0.);
}


float badbox(vec3 p, vec3 d) {
  p=abs(p)-d;
  return max(max(p.x,p.y),p.z);
}

float dbox(vec3 c, vec3 n) {
  vec3 tn=clamp(n,-bb_dim, bb_dim);
  if(tn!=n) return 1e3;
  vec3 s0=shash3(n);
  float
    h0=fract(s0.x+s0.y)
  , h1=mix(-.2,.2,fract(s0.x+s0.z))*TIME
  , h2=mix(-.2,.2,fract(s0.y+s0.z))*TIME
  ;
  if(h0<.666) return 1e3;
  if(abs(6.+dot(n,normalize(vec3(-nstripe.y,nstripe.x,1))))>1.) return 1e3;
  vec3
    off=.4*shash3(n)
  , bc=c-off
  ;
  rot(bc.xy,TAU*h1);
  rot(bc.xz,TAU*h2);

//  return length8(bc)-.05;
  return badbox(bc,vec3(.05));
}

// in
vec3 g_hsrd ;
vec3 g_ird  ;

// out
vec3 g_N    ;
float df(vec3 p) {
  float dbb=badbox(p,bb_dim+1.);
  if(dbb>1.) return dbb;
  const float EPS=1e-3;
  vec3
    n=round(p)
  , c=p-n
  ;

  g_N=n;
  float
    d=dbox(c,n);
  ;


  vec3
    dro=g_hsrd-c
  , dt=dro*g_ird
  ;

  float
    nd=EPS+min(min(dt.x,dt.y),dt.z);
  ;
  d=min(d,nd);

  return d;
}


float raymarch(vec3 ro, vec3 rd, float tinit, out float iter) {
  float
    z=tinit
  , d
  , i
  , N=1e3
  , Z=z
  ;
  vec3
    p=ro+z*rd
  ;
  for(i=0.;i<MaxIter;++i) {
    d=df(p);
    if(d<MinDistance||z>MaxDistance) {
      break;
    }
    if (d<N) {
      N=d;
      Z=z;
    }
    z+=d;
    p+=d*rd;
  }

  if(i==MaxIter) {
    z=MaxDistance;
  }
  iter=i;
  return z;
}

vec3 nf(vec3 p) {
  vec2 e=vec2(NormalEps,0);
  return normalize(vec3(
    df(p+e.xyy)-df(p-e.xyy)
  , df(p+e.yxy)-df(p-e.yxy)
  , df(p+e.yyx)-df(p-e.yyx)
  ));
}

float dstripe(vec2 p) {
  float
    ds=dot(nstripe,p)
  , d0=abs(ds)-1.
  , d1=.125-abs(ds-.56)
  , d2=.03-abs(ds+.8)
  , d=max(max(d0,d1),d2)
  ;
  return d;
}

vec3 stripe(vec3 col, vec2 p, float aa, out float d) {
  float
    Z2=.13-.015*p.y
  , ds=dstripe((p-vec2(.59,0))/Z2)*Z2
  ;
  col=mix(col,amiga_white,smoothstep(aa,-aa,ds));
  d=ds;
  return col;
}

float dencore(vec2 p) {
  const vec2
    bd=vec2(.465,.4/3.)
  , nn=normalize(vec2(1,1.775))
  ;
  vec2
    c=p-vec2(.5,0)
  , c0=c
  , c1=c
  ;

  float
    nx=clamp(round(c.x),-3.,2.)
  , ny=clamp(round(c.y*3.),-1.,1.)
  ;
  c0-=vec2(nx,ny/3.);
  c1-=vec2(nx,0);

  float
    d0=box(c0, bd)
  , d1=box(abs(c1)-vec2(1./3.,0), bd.yx)
  , d3=max(d0,1./6.-abs(c1.y))
  , d
  , dc=min(d3,max(d1,c1.x))
  , de=d0
  , dn=min(d1,max(abs(c1.x)-.125,abs(dot(c1,nn))-.125))
  , do_=min(d3,d1)
  , dr=max(d0,min(d1,-1./6.-c1.y))
  ;


  if(nx==-3.||nx==2.) d=de;
  if(nx==-2.) d=dn;
  if(nx==-1.) d=dc;
  if(nx== 0.) d=do_;
  if(nx== 1.) d=dr;

  return d;
}

float d500(vec2 p) {
  vec2
    c=p
  ;

  float
    nx=clamp(round(c.x),-1.,1.)
  ;
  c.x-=nx;
  float
    db=box(c,vec2(.465))
  , ds=min(.03-abs(abs(c.y)-1./6.),sign(c.y)*c.x+.205)
  , d=1e3
  , d0=max(db,-db-.27)
  , d5=max(db,ds)
  ;
  d=d0;
  if(nx==-1.) d=d5;
  return d;
}

float hash(float co) {
  return fract(sin(co*12.9898) * 13758.5453);
}

vec2 hash2(vec2 p) {
  p = vec2(dot (p, vec2 (127.1, 311.7)), dot (p, vec2 (269.5, 183.3)));
  return fract(sin(p)*43758.5453123);
}

vec3 logo(vec3 col, vec2 p, float aa) {
  const float
    Z0=.266
  , Z1=.355
  ;
  const vec2
    gsz=vec2(10,1)/20.
  ;
  vec2
    n=round(p/gsz)
  , h0=hash2(n-round(TIME*10.)*.1234)
  ;
  float
    h1=fract(h0.x+h0.y)
  , h2=hash(round(TIME*.5))
  ;
  if(h1>.9&&h2>.9)
    p+=gsz*vec2(.5,1)*(-1.+2.*h0);
  float
    d0=dencore((p-vec2(-.663,-.197))/Z0)*Z0
  , d1=d500((p-vec2(-.663,-.505))/Z1)*Z1
  , d=min(d0,d1)
  ;
  col=mix(col,amiga_orange,smoothstep(aa,-aa,d));
  return col;
}


vec3 pass0() {
  vec2
    r=RENDERSIZE
  , p=2.*_uvc
  ;
  float
    i
  , z
  , aa=sqrt(2.)/r.y
  , sd
  ;
  const vec3
  , Z=vec3(0,0,1)
  , X=normalize(cross(vec3(0,1,0),Z))
  , Y=cross(X,Z)
  ;
  vec3
    ro=vec3(nstripe*TIME*.3,-12)
  , rd=normalize(p.x*X+p.y*Y+2.*Z)
  , bkg=amiga_blue*(1.-.1*dot(nstripe,p))
  , col=bkg
  ;

  col=stripe(col,p,aa, sd);
  g_hsrd=sign(rd)*.5;
  g_ird=1./rd;


  z=raymarch(ro,rd,0.,i);
  if(z<MaxDistance&&(sd>.0||z<6.)) {
    vec3
      N3=g_N
    , p3=ro+z*rd
    , n3=nf(p3)
    , r3=reflect(rd,n3)
    ;
    float
      h0=hash(N3)
    , spe=pow(max(dot(r3,light_dir),0.),10.)
    ;
    if(h0<.25) {
      col=amiga_orange;
    } else if (h0<.5) {
      col=amiga_blue2;
    } else if (h0<.75) {
      col=amiga_white;
    } else {
      col=amiga_black;
    }
    col*=mix(1.,.5,i/MaxIter);
    col+=spe;
    col=mix(bkg, col, exp(-.002*z*z));
  }
  // Surprisingly nice!
  // col=vec3(i/MaxIter);
  return col;
}

vec3 gb(sampler2D pp, vec2 n, vec2 dir) {
  vec2 
    q=_uv
  , p=2.*_uvc
  ;
  vec3 
    col=texture(pp,q).xyz
  ;

  float
    s=max(2.+1.*dot(n,p),.001)
  , s2=2.*s*s
  , w
  , ws=1.
  ;

  for(float i=1.;i<9.;++i) {
    w=exp(-(i*i)/s2);
    vec2 off=dir*i;

    col+=w*(texture(pp,q-off).xyz+texture(pp,q+off).xyz);
    ws+=2.*w;
  }
  col/=ws;
  return col;
}

vec3 pass1() {
  vec3 col=vec3(0);
  vec2
    r=RENDERSIZE
  , q=_uv
  , p=2.*_uvc
  , n=normalize(vec2(1.,1.06))
  , dir=n/r
  ;
  float
    s=max(2.+1.*dot(n,p),.001)
  , s2=2.*s*s
  , w
  , ws=1.
  ;
  col+=texture(passA,q).xyz;
  for(float i=1.;i<9.;++i) {
    w=exp(-(i*i)/s2);
    vec2 off=dir*i;

    col+=w*(texture(passA,q-off).xyz+texture(passA,q+off).xyz);
    ws+=2.*w;
  }
  col/=ws;
  return col;
}

vec3 pass2() {
  vec3 col=vec3(0);
  vec2
    r=RENDERSIZE
  , q=_uv
  , p=2.*_uvc
  , n=normalize(vec2(1.,1.06))
  , dir=vec2(n.y,-n.x)/r
  ;
  float
    s=max(2.+1.*dot(n,p),.001)
  , s2=2.*s*s
  , w
  , ws=1.
  ;
  col+=texture(passB,q).xyz;
  for(float i=1.;i<9.;++i) {
    w=exp(-(i*i)/s2);
    vec2 off=dir*i;

    col+=w*(texture(passB,q-off).xyz+texture(passB,q+off).xyz);
    ws+=2.*w;
  }
  col/=ws;
  return col;
}

vec3 pass3() {
  vec2
    r=RENDERSIZE
  , q=_uv
  , p=2.*_uvc
  , n=normalize(vec2(1.,1.06))
  , dir=n/r
  , off=(3.+1.*dot(n,p))*dir
  ;

  float
    aa=sqrt(2.)/r.y
  ;

  vec4 pcol=texture(syn_FinalPass, q);
  vec3 col=vec3(0);
  col+=vec3(
    texture(passC, q-off).x
  , texture(passC, q).y
  , texture(passC, q+off).z
  );
  col=logo(col,p,aa);
  col=sqrt(col);
  col=mix(col,pcol.xyz,.5);
  return col;
}

vec4 renderMain() {
  vec3 col=vec3(0);
  switch(PASSINDEX)
  {
  case 0:
    col=pass0();
    break;
  case 1:
    col=pass1();
    break;
  case 2:
    col=pass2();
    break;
  default:
    col=pass3();
  };
  return vec4(col,1);
}
