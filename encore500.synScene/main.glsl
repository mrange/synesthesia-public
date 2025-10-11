#define RGB(a)      (vec3(a*a)/(255.*255.))
#define OKRGB(a)    LINEARTOOKLAB((vec3(a*a)/(255.*255.)))
#define ROT(a)      mat2(cos(a), sin(a), -sin(a), cos(a))

#ifdef KODELIFE
const float
    beat_speed      =60.
,   blurriness      =.2
,   color_distortion=.2
,   glitch_freq     =.9
,   glitch_level    =.9
,   glitch_size     =vec2(10,1)/20.
,   motion_blur     =.5
,   show_beat       =0.
,   volume_control  =0.
;
#endif

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
  amiga_bkg   =RGB(ivec3(0x00,0x34,0x9A))
, amiga_blue  =RGB(ivec3(0x00,0x49,0xCC))
, amiga_orange=RGB(ivec3(0xFF,0x80,0x18))
, amiga_white =RGB(ivec3(0xE4,0xE4,0xE4))
, amiga_black =RGB(ivec3(0x11,0x12,0x13))
, real_black  =RGB(ivec3(0x01,0x01,0x01))
, real_orange =RGB(ivec3(0xFF,0x80,0x18))
, real_white  =RGB(ivec3(0xFD,0xFD,0xFD))
, bb_dim      =vec3(1e3,1e3,10)
, light_dir   =normalize(vec3(-1,1,-1))
;

const mat3
  OKLAB_M1=mat3(
    0.4122214708, 0.5363325363, 0.0514459929
  , 0.2119034982, 0.6806995451, 0.1073969566
  , 0.0883024619, 0.2817188376, 0.6299787005
  )
, OKLAB_M2=mat3(
    0.2104542553,  0.7936177850, -0.0040720468
  , 1.9779984951, -2.4285922050,  0.4505937099
  , 0.0259040371,  0.7827717662, -0.8086757660
  )
, OKLAB_N1 = mat3(
    1,  0.3963377774,  0.2158037573
  , 1, -0.1055613458, -0.0638541728
  , 1, -0.0894841775, -1.2914855480
  )
, OKLAB_N2 = mat3(
     4.0767416621, -3.3077115913,  0.2309699292
  , -1.2684380046,  2.6097574011, -0.3413193965
  , -0.0041960863, -0.7034186147,  1.7076147010
  )
;

#define LINEARTOOKLAB(c) (OKLAB_M2*pow(OKLAB_M1*(c), vec3(1./3.)))
vec3 linearToOklab(vec3 c) {
  return OKLAB_M2*pow(OKLAB_M1*(c), vec3(1./3.));
}

#define OKLABTOLINEAR(c) (OKLAB_N2*pow(OKLAB_N1*(c),vec3(3)))
vec3 oklabToLinear(vec3 c) {
  return OKLAB_N2*pow(OKLAB_N1*(c),vec3(3));
}

vec3 cool(float t) {
  const vec3
    ok_red   = OKRGB(ivec3(0xF0,0x00,0x0A))
  , ok_orange= OKRGB(ivec3(0xFB,0x82,0x17))
  , ok_white = OKRGB(ivec3(0xFF,0xFF,0xFF))
  , ok_blue  = OKRGB(ivec3(0x00,0x40,0xFA))
  ;

  vec3
    ok_color
  , ok_from
  , ok_to
  ;

  float
    sub
  ;

  if(t < 1./3.) {
    ok_from=ok_red    ;
    ok_to  =ok_orange ;
    sub    =0.        ;
  } else if (t < 2./3.){
    ok_from=ok_orange ;
    ok_to  =ok_white  ;
    sub    =1./3.     ;
  } else {
    ok_from=ok_white  ;
    ok_to  =ok_blue   ;
    sub    =2./3.     ;
  }

  ok_color = mix(ok_from, ok_to, smoothstep(0., .33, t-sub));

  return oklabToLinear(ok_color);
}

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
  if(any(notEqual(tn, n))) return 1e3;
  vec3 s0=shash3(n);
  float
    h0=fract(s0.x+s0.y)
  ;
  if(h0<.6667) return 1e3;
  float
    h1=mix(-.2,.2,fract(s0.x+s0.z))*TIME
  , h2=mix(-.2,.2,fract(s0.y+s0.z))*TIME
  ;
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
  ;
  vec3
    p=ro+z*rd
  ;
  for(i=0.;i<MaxIter;++i) {
    d=df(p);
    if(d<MinDistance||z>MaxDistance) {
      break;
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

vec3 stripe(vec3 col, vec3 fg, vec2 p, float aa, out float d) {
  float
    Z2=.13-.015*p.y
  , ds=dstripe((p-vec2(.59,0))/Z2)*Z2
  ;
  col=mix(col,fg,smoothstep(aa,-aa,ds));
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
float d2(vec2 p) {
  float
    db=box(p,vec2(.465))
  , ds=min(.03-abs(abs(p.y)-1./6.),-sign(p.y)*p.x+.205)
  , d =max(db,ds)
  ;
  return d;
}

float dencore500(vec2 p) {
  const float
    Z0=.266
  , Z1=.355
  ;
  float
    d0=dencore((p-vec2(0., .134))/Z0)*Z0
  , d1=d500((p-vec2(0.,-.174))/Z1)*Z1
  ;

  return min(d0, d1);
}

float dheart(vec2 p, out vec2 cellid) {
  const mat2
    R=ROT(radians(-45.))
  ;
  p*=R;
  vec2
    n=clamp(round(p),vec2(-1), vec2(1))
  , c=p-n
  ;
  float
    d0=box(c,vec2(.465))
  , d1=-max(.5+p.x,.5-p.y)
  , d=max(d0,d1)
  ;
  cellid=n;
  return d;
}

float hash(float co) {
  return fract(sin(co*12.9898) * 13758.5453);
}

vec2 hash2(vec2 p) {
  p = vec2(dot (p, vec2 (127.1, 311.7)), dot (p, vec2 (269.5, 183.3)));
  return fract(sin(p)*43758.5453123);
}


//#define DEBUG

vec3 logo(vec3 col, vec2 p, float aa) {
#ifdef DEBUG
  vec2 ap=abs(p);

#endif
  p.y -= 0.283;
  const float
    Z0 = 0.53*1.11
  , ZE = Z0
  , ZH =.2
  ;
  const vec2
    off=vec2(.89,-.425)
  , bsz=5.2*vec2(1,85./588.)
  ;
  vec3
    bcol=real_orange
    ;
  vec2
    n=round(p/glitch_size)
  , h0=hash2(n-round(TIME*10.)*.1234)
  , cellid
  , ph
  ;
  float
    h1=fract(h0.x+h0.y)
  , h2=hash(round(TIME*.5))
  , bz=mix(9.,1.,volume_control)
  , Z2=.10*Z0*bz
  ;
  if(h1>glitch_level&&h2>glitch_freq)
    p+=glitch_size*vec2(.5,1)*(-1.+2.*h0);
  ph=(p-vec2(-off.x, off.y))/ZH;
  ph*=mix(vec2(1.),vec2(.9,.95), show_beat*smoothstep(.5,.9,sin(beat_speed*(TAU/60.)*TIME+ph.y*.5)));
  float
    de=dencore500((p-vec2(-off.x,-off.y+.024))/ZE)*ZE
  , dh=dheart(ph,cellid)*ZH
  , d2=d2((p-vec2(-off.x,-off.y-.27))/Z2)*Z2
  , db=box(p-vec2(-off.x,-off.y-.27),Z2*bsz)
  ;

  col=mix(
    col
  , mix((cellid.x==1.||cellid.y==-1.)?amiga_white:amiga_orange, cool(.43+ph.y*.2), volume_control)
  , smoothstep(aa,-aa,dh)
  );

  col=mix(
    col
  , mix(amiga_orange, real_white, volume_control)
  , smoothstep(aa,-aa,de)
  );

  bcol=mix(
    bcol
  , real_black
  , smoothstep(aa,-aa,d2)
  );

  col=mix(
    col
  , mix(col, bcol, volume_control)
  , smoothstep(aa,-aa,db)
  );
#ifdef DEBUG
  col=mix(
    col
  , real_white
  , smoothstep(aa,-aa,min(ap.x,ap.y))
  );
#endif

  return col;
}

vec3 gb(sampler2D pp, vec2 dir) {
  vec2
    q=_uv
  , p=2.*_uvc
  ;
  vec3
    col=texture(pp,q).xyz
  ;

  float
    s=max(2.+1.*dot(nstripe,p),.001)
  , s2=blurriness*10.*s*s
  , w
  , ws=1.
  ;

  for(float i=1.;i<19.;++i) {
    w=exp(-(i*i)/s2);
    vec2 off=dir*i;

    col+=w*(texture(pp,clamp(q-off, 0., 1.)).xyz+texture(pp,clamp(q+off, 0., 1.)).xyz);
    ws+=2.*w;
  }
  col/=ws;
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
  , box_1
  , box_2
  , bkg_1=amiga_bkg*(1.-.1*dot(nstripe,p))
  , bkg_2=real_black
  , grd_2=cool(.25+.5*dot(nstripe,p))
  , str_1=amiga_white
  , str_2=grd_2
  , str  =mix(str_1, str_2, volume_control)
  , bkg  =mix(bkg_1, bkg_2, volume_control)
  , col  =bkg
  ;

  col=stripe(col, str, p , aa, sd);

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
      box_1=amiga_orange;
    } else if (h0<.5) {
      box_1=amiga_blue;
    } else if (h0<.75) {
      box_1=amiga_white;
    } else {
      box_1=amiga_black;
    }
    box_2=grd_2;
    col=mix(box_1, box_2, volume_control);
    col*=mix(1.,.5,i/MaxIter);
    col+=spe;
    col=mix(bkg, col, exp(-.002*z*z));
  }
  // Surprisingly nice!
  // col=vec3(i/MaxIter);
  return col;
}

vec3 pass1() {
  return gb(passA,vec2(1,0)/RENDERSIZE);
}

vec3 pass2() {
  return gb(passB,vec2(0,1)/RENDERSIZE);
}

vec3 pass3() {
  vec2
    r=RENDERSIZE
  , q=_uv
  , p=2.*_uvc
  , dir=nstripe/r
  , off=(3.+1.*dot(nstripe,p))*dir*color_distortion*5.
  ;

  float
    aa=sqrt(2.)/r.y
  ;

  vec4 pcol=texture(syn_FinalPass, q);
  vec3 col=vec3(0);
  col+=vec3(
    texture(passC, clamp(q-off,0.,1.)).x
  , texture(passC, clamp(q,0.,1.)).y
  , texture(passC, clamp(q+off,0.,1.)).z
  );
  col=logo(col,p,aa);
  col=sqrt(col);
  col=mix(col,pcol.xyz,motion_blur);
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
    break;
  };
  return vec4(col,1);
}
