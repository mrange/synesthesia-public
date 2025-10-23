#define ROT(a)      mat2(cos(a), sin(a), -sin(a), cos(a))
#define SC(a)       vec2(sin(a),cos(a))


// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
const vec4 hsv2rgb_K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
vec3 hsv2rgb(vec3 c) {
  vec3 p = abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www);
  return c.z * mix(hsv2rgb_K.xxx, clamp(p - hsv2rgb_K.xxx, 0.0, 1.0), c.y);
}
// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
//  Macro version of above to enable compile-time constants
#define HSV2RGB(c)  (c.z * mix(hsv2rgb_K.xxx, clamp(abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www) - hsv2rgb_K.xxx, 0.0, 1.0), c.y))

#ifdef KODELIFE
const float
  flickerness     =.6
, motion_blur     =.4
, rotation_speed  = .2
, sea_level       =-.66
, tilt_control    =-.4
;

const vec2
  camera        =vec2(3.,.55)
, media_control =vec2(10.,.5)
;

const vec3
, bars_col        = HSV2RGB(vec3(0.85, 0.90, .5))
, horiz_col       = HSV2RGB(vec3(0.06, 0.90, .5))
, sun_col         = HSV2RGB(vec3(0.03, 0.90, .5))
, top_box_col     = HSV2RGB(vec3(0.60, 0.90, .5))

, border_col      = HSV2RGB(vec3(0.64,.98 ,.5))
, air_col         = HSV2RGB(vec3(0.0 ,.98 ,.5))
, on_col          = HSV2RGB(vec3(0.1 ,.0  ,.5))

/*
, bars_col        = vec3(0.50, 0.05, 0.45)
, horiz_col       = vec3(0.50, 0.21, 0.05)
, sun_col         = vec3(0.50, 0.13, 0.05)
, top_box_col     = vec3(0.05, 0.23, 0.50)

, border_col      = vec3(0.01, 0.09, 0.50)
, air_col         = vec3(0.50, 0.01, 0.01)
, on_col          = vec3(0.50, 0.50, 0.50)
*/
;
#endif

const float
  max_distance_1  = 14.
, max_iteration_1 = 70.
, norm_eps_1      = 1e-3
, tolerance_1     = 1e-4
, pi              = acos(-1.)
, tau             = 2.*pi
;

const vec3
  sun_dir         = normalize(vec3(1,.2,1))
;

const vec4
  planet_dim      = vec4(50.*sun_dir-vec3(0,3,-2), 10.)
;

float freq(float x) {
#ifdef KODELIFE
  return .5*sin(2.*6.28*x)+exp(-fract(TIME))-.5;
#else
  x=.5+.45*sin(x*tau);
  return texture(syn_Spectrum ,x).y-.5;
#endif
}

// IQ's ray sphere intersect: https://iquilezles.org/articles/intersectors
vec2 raysphere(vec3 ro, vec3 rd, vec4 sph) {
  vec3
    oc = ro-sph.xyz
  ;
  float
    b = dot(oc, rd)
  , c = dot(oc, oc)-sph.w*sph.w
  , h = b*b-c
  ;
  if (h < 0.) return vec2(-1.0);
  h = sqrt(h);
  return vec2(-b-h,-b+h);
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
float arc(vec2 p, vec2 sc, float ra, float rb ) {
  // sc is the sin/cos of the arc's aperture
  p.x = abs(p.x);
  return ((sc.y*p.x>sc.x*p.y) ? length(p-sc*ra) :
                                abs(length(p)-ra)) - rb;
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

float length4(vec2 p) {
  p*=p;
  return pow(dot(p,p),.25);
}

float length8(vec2 p) {
  p*=p;
  p*=p;
  return pow(dot(p,p),.125);
}

float length4(vec3 p) {
  p*=p;
  return pow(dot(p,p),.25);
}

float length8(vec3 p) {
  p*=p;
  p*=p;
  return pow(dot(p,p),.125);
}


// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
float box(vec2 p, vec2 b) {
  vec2 d = abs(p)-b;
  return length(max(d,0.0)) + min(max(d.x,d.y),0.0);
}

float box(vec2 p, vec2 b, float r) {
  vec2 d = abs(p)-b+r;
  return length(max(d,0.0)) + min(max(d.x,d.y),0.0)-r;
}

float box(vec3 p, vec2 b, vec2 r) {
  float d = box(p.xy,b,r.x);
  vec2 w = vec2(d, abs(p.z) - r.y);
  return min(max(w.x,w.y),0.) + length4(max(w,0.));
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
float box(vec2 p, vec2 b, vec4 r) {
  r.xy = (p.x>0.0)?r.xy : r.zw;
  r.x  = (p.y>0.0)?r.x  : r.y;
  vec2 q = abs(p)-b+r.x;
  return min(max(q.x,q.y),0.0) + length(max(q,0.0)) - r.x;
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
float horseshoe(vec2 p, vec2 c, float r, vec2 w) {
  p.x = abs(p.x);
  float l = length(p);
  p = mat2(-c.x, c.y, c.y, c.x)*p;
  p = vec2((p.y>0.0 || p.x>0.0)?p.x:l*sign(-c.x),
           (p.x>0.0)?p.y:l );
  p = vec2(p.x,abs(p.y-r))-w;
  return length(max(p,0.0)) + min(0.0,max(p.x,p.y));
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
float segment(vec2 p, vec2 a, vec2 b) {
  vec2 pa = p-a, ba = b-a;
  float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
  return length( pa - ba*h );
}

float segmenty(vec2 p, vec2 d) {
  p.y=abs(p.y)-d.x;
  float
    d0=length(p)
  , d1=abs(p.x);

  return (p.y>0.?d0:d1)-d.y;
}

float segmentx(vec2 p, vec2 d) {
  p.x=abs(p.x)-d.x;
  float
    d0=length(p)
  , d1=abs(p.y);

  return (p.x>0.?d0:d1)-d.y;
}

float box(vec3 p, vec3 b, vec3 bb) {
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float torus(vec3 p, vec2 t) {
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length4(q)-t.y;
}

// License: Unknown, author: Unknown, found: don't remember
float hash(float co) {
  return fract(sin(co*12.9898) * 13758.5453);
}

vec3 onair(vec2 p) {
  const float
    RB=0.19
  ;
  const float
    AN=10.635
  ;
  const vec2
    SCO=SC(radians(165.))
  , SCA=SC(radians(-15.))
  , SCN=SC(radians(-AN))
  ;
  const mat2
    ROTN=ROT(radians(AN))
  ;
  float
    d
  , da
  , db
  , di
  , dn
  , do_
  , dr
  , sn
  , don
  , dair
  ;
  vec2
    pa=p
  , pb=p
  , pi=p
  , pn=p
  , po=p
  , pr=p
  ;
  pa-=vec2(0.26,0.19);
  pa.y=-pa.y;
  pi-=vec2(.52,0);
  pn-=vec2(-0.22,-0.0);
  sn=pn.x;
  pn.x=abs(pn.x);
  pn-=vec2(0.11,0.0);
  if(sn>0.) {
    pn.y-=-.2;
  } else {
    pn.y=.2-pn.y;
  }
  pn*=ROTN;
  po-=vec2(-.675,0);
  pr-=vec2(.75,0.11);
  da=min(
      segmentx(pa-vec2(0,0.26), vec2(.1,.0))
    , horseshoe(pa,SCA,.04,vec2(0.45,0))
    );
  db=max(abs(box(pb,vec2(1.,0.4),RB)),.075-max(pb.y, abs(pb.x)));
  di=segmenty(pi,vec2(0.23,0));
  dn=horseshoe(pn,SCN,0.03,vec2(.43,0));
  do_=arc(po,SCO,.23,0.);
  dr=
    min(min(
        abs(box(pr,vec2(0.1,0.12),vec4(0.1,0.1,.04,0)))
      , segmenty(pr-vec2(-0.1,-.14),vec2(0.2,0)))
    , segment(pr,vec2(0,-0.12), vec2(0.11,-.34))
    );

  don=min(do_,dn);
  dair=min(min(da,di),dr);

  return vec3(db, min(do_,dn), min(min(da,di),dr));
}

vec3 onair_(vec3 p) {
  vec4 d = vec4(onair(p.xy),abs(p.z));
  return vec3(
    min(max(d.x,d.w),0.) + length(max(d.xw,0.))
  , min(max(d.y,d.w),0.) + length(max(d.yw,0.))
  , min(max(d.z,d.w),0.) + length(max(d.zw,0.))
  );
}

vec3 onair(vec3 p) {
  float
    z2 = p.z*p.z
  ;

  vec3
    d2d  = onair(p.xy)
  , d2d2 = d2d*d2d
  ;

  return sqrt(vec3(
    d2d.x <= 0.0 ? z2 : (d2d2.x + z2)
  , d2d.y <= 0.0 ? z2 : (d2d2.y + z2)
  , d2d.z <= 0.0 ? z2 : (d2d2.z + z2)
  ));
}

vec3 g_G;
vec4 g_H;
float df1(vec3 p) {
  vec3
    d0
  , p0=p
  , p1=p
  , p3=p
  ;
  vec2
    S=sin(p.xz*vec2(23.0,13.0))
  ;
  p0.z=-abs(p0.z);
  p0.z+=.07;
  d0=onair(p0);
  float
    d1=box(p1,vec2(1.1,.5),vec2(.29,.01))-.01
  , d2=p.y-sea_level+1e-3/(1.+0.1*dot(p.xz,p.xz))*S.x*S.y
  , d3=torus(p3.yzx,1.41*vec2(1.,.02))
  , d=min(min(d0.x,d0.y),d0.z)-.02
  ;
  d1=min(d1,d3);
  d1=min(d1,d2);
  g_G=min(g_G,d0);
  d=min(d,d1);
  g_H=vec4(d0,d1);
  return d;
}

vec3 nf1(vec3 p) {
  const vec2
    e=vec2(norm_eps_1,0)
  ;
  return normalize(vec3(
    df1(p+e.xyy)-df1(p-e.xyy)
  , df1(p+e.yxy)-df1(p-e.yxy)
  , df1(p+e.yyx)-df1(p-e.yyx)
  ));
}

float raymarch1(vec3 ro, vec3 rd, float initz) {
  float
    i
  , d
  , z=initz
  , D=1e3
  , Z=initz
  ;

  for(i=0.;i<max_iteration_1;++i) {
    vec3
      p=ro+rd*z
    ;
    d=df1(p);
    if(d<D) {
      D=d;
      Z=z;
    }
    if(d<tolerance_1||z>max_distance_1) {
      break;
    }
    z+=d;
  }

  if(i==max_iteration_1) {
//    z=Z;
  }
  return z;
}

vec4 tcol(vec2 p) {
  vec2
    sz=vec2(textureSize(syn_Media,0))
  , S
  ;
  p.x*=sz.y/sz.x;
  S=step(abs(p),vec2(.5));
  p+=.5;

#ifdef KODELIFE
  p.y=1.-p.y;
#endif
  vec4
    tcol=texture(syn_Media, p)
  ;
  tcol.xyz*=tcol.xyz;
  tcol.w*=S.x*S.y;
  return tcol;
}

vec3 render0(vec3 ro, vec3 rd) {
  vec3 col = vec3(0.0);

  float
    tp   = (9.-ro.y)/(rd.y)
  ;

  if (tp>0.) {
    vec3 pos  = ro + tp*rd;
    vec2 pp = pos.xz;
    float db = box(pp, vec2(5.0, 9.0))-3.0;

    col += rd.y*rd.y*smoothstep(0.25, 0.0, db)*top_box_col;
    col += 2e-1*exp(-0.5*max(db, 0.0))*top_box_col;
    col += 5e-2*max(-db, 0.0)*sqrt(top_box_col);
    col *= 4.;
  }

  vec2
    P=vec2(pi+atan(rd.x,rd.z),rd.y)/pi
  , S=raysphere(ro,rd,planet_dim)
  ;

  vec4
    tcol=tcol(media_control.x*P-vec2(media_control.x*.25,media_control.y))
  ;
  col += 8e-7/pow(1.01-dot(rd,sun_dir),4.)*mix(0.,1.,exp(-.3*(S.y-S.x)))*sun_col;
  col += 4e-3/abs(rd.y)*horiz_col;
  const float ZZ=.01,CC=2.;
  float
    N=round(P.x/ZZ)
  , D=1e3
  , j
  ;
  P.x-=N*ZZ;
  for(j=-CC;j<=CC;++j) {
    D=min(D, box(P+vec2(-ZZ*j,0),vec2(ZZ*.1,(freq((N+j)*ZZ))*.04))-ZZ*.2);
  }

  col+=1./max(D,1e-3)*(bars_col*2E-3-.02*D);
  col = mix(col,tcol.xyz*tcol.xyz,tcol.w);

  return col;
}

vec3 render1(vec3 ro, vec3 rd) {
  bool
    hit
  ;
  vec3
    col
  , n
  , p
  , r
  , rr
  , rcol
  , G
  , S=vec3(0)
  , tcol=vec3(0)
  , mcol=vec3(1)
  , col_on
  , col_air
  , col_border
  ;

  float
    d
  , f
  , z
  , Z
  , L
  , g
  , i
  , h0=hash(floor(9.*TIME))
  , h1=fract(6047.*h0)
  , h2=fract(7907.*h0)
  , H0=hash(floor(.5*TIME))
  , H1=fract(6047.*H0)
  , H2=fract(7907.*H0)
  ;
  vec4
    H
  ;
  col_on=(h0>H0-flickerness?2.:0.)*on_col;
  col_air=(h1>H1-flickerness?2.:0.)*air_col;
  col_border=(h2>H2-flickerness?2.:0.)*border_col;
  for(i=0.;i<3.;++i) {
    col=vec3(0.);
    Z=(sea_level-ro.y)/rd.y;
    g_G=vec3(1e3);
    z=raymarch1(ro,rd,0.01);
    G=g_G;
    H=g_H;
    p=ro+rd*z;
    n=nf1(p);

    if (z<=max_distance_1) {
      hit = true;
    } else if(Z>=max_distance_1){
      hit=true;
      n=vec3(0,1,0);
      // TODO: Unsure why *Z doesn't work
      p=ro+rd*max_distance_1;
    } else {
      hit = false;
    }

    r=reflect(rd,n);
    rr=refract(rd,n,.2);
    f=1.+dot(rd,n);
    f*=f;
    g=.05/(1.01-(dot(rr,rd)));
    g*=g;
    if(hit) {
      d=min(min(min(H.x,H.y),H.z),H.w);
      if(d==H.x) {
        S=col_border;
      } else if(d==H.y) {
        S=col_on;
      } else if(d==H.z) {
        S=col_air;
      } else {
        f=mix(.3,1.,f);
      }
      col+=clamp(S*g,0.,3.);
    }
    G*=G;

    col += 5e-4/G.x*col_border;
    col += 5e-4/G.y*col_on;
    col += 5e-4/G.z*col_air;
    tcol+=col*mcol;
    if(!hit) {
      tcol+=mcol*render0(ro,rd);
      break;
    }
    mcol*=f;
    rd=r;
    ro=p+tolerance_1*1e1*(n+rd*1e1);
  }

  return tcol;
}

vec4 renderMain() {
  vec2
    r =RENDERSIZE.xy
  , q =_uv
  , p =2.*_uvc
  ;
  vec3
    ro  =vec3(0,camera.y,-camera.x)
  , up  = normalize(vec3(-tilt_control,1,0))
  ;
#ifdef KODELIFE
  mat2
//    R=ROT(radians(3.+180.))
    R=ROT(rotation_speed*TIME)
  ;
#else
  mat2
    R=ROT(rotation_speed)
  ;
#endif
  ro.xz *= R;
  up.xz *= R;
  vec3
    la  =vec3(0,0.0,0)
  , Z   =normalize(la-ro)
  , X   =normalize(cross(Z,up))
  , Y   =cross(X,Z)
  , rd  =normalize(-p.x*X+p.y*Y+2.*Z)
  , col =vec3(0)
  ;

  col=render1(ro,rd);
  col -=vec3(2,3,1)*3e-3*(.25+length(p));
  col=tanh(col);
  col=clamp(col,0.,1.);
  vec4
    pcol=texture(syn_FinalPass,q)
  ;
  col=mix(col,pcol.xyz*pcol.xyz,motion_blur);
  col=sqrt(col);
  return vec4(col,1.);
}
