// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
const vec4 hsv2rgb_K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);

// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
//  Macro version of above to enable compile-time constants
#define HSV2RGB(c)  (c.z * mix(hsv2rgb_K.xxx, clamp(abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www) - hsv2rgb_K.xxx, 0.0, 1.0), c.y))

// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
vec3 hsv2rgb(vec3 c) {
  vec3 p = abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www);
  return c.z * mix(hsv2rgb_K.xxx, clamp(p - hsv2rgb_K.xxx, 0.0, 1.0), c.y);
}

// License: Unknown, author: Unknown, found: don't remember
float hash(vec2 co) {
  return fract(sin(dot(co.xy ,vec2(12.9898,58.233))) * 13758.5453);
}

// License: Unknown, author: Claude Brezinski, found: https://mathr.co.uk/blog/2017-09-06_approximating_hyperbolic_tangent.html
vec3 tanh_approx(vec3 x) {
  vec3 x2 = x*x;
  return clamp(x*(27.0 + x2)/(27.0+9.0*x2), -1.0, 1.0);
}


// License: MIT, author: Inigo Quilez, found: https://www.iquilezles.org/www/articles/spherefunctions/spherefunctions.htm
float ray_sphere(vec3 ro, vec3 rd, vec4 sph) {
  vec3 oc = ro - sph.xyz;
  float 
    b = dot(oc, rd)
  , c = dot(oc, oc)- sph.w*sph.w
  , h = b*b-c
  ;
  h = sqrt(h);
  return -b-h;
}

float doctahedron(vec3 p, float s) {
  p = abs(p);
  return (p.x+p.y+p.z-s)*0.57735027;
}

float dfm(vec3 p) {
   float
    d=p.y+.6
  , a=1.
  ;

  vec2
    D=vec2(0)
  , P=p.xz*.23
  ;

  vec4
    o
  ;

  for(int j=0;j<7;++j) {
    o=cos(P.xxyy+vec4(11,0,11,0));
    p=o.yxx*o.zwz;
    D+=p.xy;
    d-=a*(1.+p.z)/(1.+3.*dot(D,D));
    P*=mat2(1.2,1.6,-1.6,1.2);
    a*=.55;
  }
  return d;
}

vec3 path(float z) {
  return vec3(2.*vec2(3.1,.7)*cos(z*.5*vec2(.11,.07)),z)+vec3(25,3,0);
}

vec3 dpath(float z) {
  float dt=.05;
  return (path(z+dt)-path(z-dt))/(2.*dt);
}

vec3 ddpath(float z) {
  float dt=.05;
  return (dpath(z+dt)-dpath(z-dt))/(2.*dt);
}

float dfo(vec3 p, out vec4 oo) {
  const float
    ZZ=10.;

  vec2
    n=floor(p.xz/ZZ+.5)
  ;
  p.xz-=n*ZZ;

  float
    h0=hash(n)
  , h1=fract(8677.*h0)
  , h2=fract(9677.*h0)
#ifdef KODELIFE
  , s=sqrt(h1)
#else
//  , s=smoothstep(.7,1.,textureLod(syn_Spectrum,.9*h1+.05,0).z)
  , s=syn_BassLevel
#endif
  , h=.3*ZZ*h0*h0+0.1
  , d0=doctahedron(p,h)
  , d
  ;
  oo=vec4(1e3,0,0,0);
  if(h2<.66) return 1e3;
  d=d0;
  oo=vec4(d,h0,h,s);
  return d;
}

float df(vec3 p, out vec4 oo) {
  p.y=abs(p.y);

  float
    d0=dfm(p)
  , d1=dfo(p,oo)
  , d
  ;
  d=d0;
  d=min(d,d1);
  return d;
}

const float
  TAU=2.*PI
, OFF=.7
;

const vec3
    BY=HSV2RGB(vec3(.05+OFF,.7,.8))
  , BM=HSV2RGB(vec3(.95+OFF,.6,.3))
  , BS=HSV2RGB(vec3(.55+OFF,.3,2.))
  , BO=HSV2RGB(vec3(.82+OFF,.6,2.))
  ;

float surface(float x) {
  x/=100.;
  float 
    a=1.
  , h=0.
  ;
  
  for(int i=0;i<5;++i) {
    h+=a*sin(x);
    x*=2.03;
    x+=123.4;
    a*=.55;
  }
  
  return abs(h);
}

vec4 renderMain() {
  float
      d=1.
    , z=0.
    , T=TIME*4.
    ;
  vec2
      p2=_uvc*2.
  ;
  vec3
      O=vec3(0)
    , p
    , P=path(T)
    , ZZ=normalize(dpath(T)+vec3(0,-0.1,0))
    , XX=normalize(cross(ZZ,vec3(0,1,0)+ddpath(T)))
    , YY=cross(XX,ZZ)
    , R=normalize(-p2.x*XX+p2.y*YY+2.*ZZ)
    , Y=(1.+R.x)*BY
    , S=(1.-R.y*R.y)*BS*Y
    ;
  vec4
      oo
    ;

  for(int i=0;i<50&&d>1e-5&&z<2e2;++i) {
    p=z*R+P;
    d=df(p,oo);
    if(p.y>0.) {
      O+=BM+d*Y;
    } else {
      O+=S;
      oo.x*=1e2;
    }

    O+=(0.5+2.*smoothstep(0.7,1.,oo.w))*smoothstep(oo.z*.78,oo.z*.8,abs(p.y))*1e-0/max(oo.x+oo.x*oo.x*oo.x*oo.x*9.,1e-2)*BO;

    z+=d*.7;
  }

  O*=9E-3;
  O=tanh_approx(O);
  
  if(R.y>0.0) {
    vec3
      N
    , H
    ;
    vec4
      S=vec4(P+vec3(-700,300,1000),400.)
    ;
    float 
      si=ray_sphere(P,R,S)
    ;
    
    H=tanh_approx(hsv2rgb(vec3(OFF-.4*R.y,.5+1.*R.y,3./(1.+800.*R.y*R.y*R.y))));
    if(si>0.) {
      p=P+R*si;
      N=normalize(p-S.xyz);
      H+=
          max(dot(N,normalize(vec3(1,-0.5,4))),0.)
        * smoothstep(0.0,0.2,R.y)
        * smoothstep(1.0,.89,1.+dot(R,N))
        * surface(p.y)
        ;
    }
    O*=H;
  }

  O-=.04*vec3(1,2,0)*(length(-1.+2.*_uv)+.2);
  O=max(O,0.);

#ifdef KODELIFE 
  O=sqrt(O);
#else
  vec4 M=_loadMedia();
  O=sqrt(O);
  O=mix(O,M.xyz,(p2.y+.5)*M.w);
#endif  
  return vec4(O,1);
}