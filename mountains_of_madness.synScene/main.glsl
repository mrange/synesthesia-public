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
  , s=smoothstep(.7,1.,textureLod(syn_Spectrum,.9*h1+.05,0).z)
  , h=.2*ZZ*h0
  , d0=doctahedron(p,h)
  , d1=length(p.xz)-.05*s+.04
  , d=ZZ*.25
  ;
  if(h2<.2) return d;
  d=d0;
//  d=min(d,d1);
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
;

vec4 renderMain() {
  float
      d=1.
    , z=0.
    ;

  const float OFF=0.75;
  const vec3
      BY=HSV2RGB(vec3(.05+OFF,.7,.8))
    , BM=HSV2RGB(vec3(.95+OFF,.6,.3))
    , BS=HSV2RGB(vec3(.55+OFF,.3,2.))
    , BO=HSV2RGB(vec3(.82+OFF,.6,2))
    , BB=HSV2RGB(vec3(.22+OFF,.5,1))
    ;
  vec3
      O=vec3(0)
    , p
    , R=normalize(vec3(_uvc-vec2(0,.5),1))
    , P=vec3(25.,4,TIME)
    , Y=(1.+R.x)*BY
    , S=(1.-R.y*R.y)*BS*Y
    ;
  vec4 
      oo
    ;

  for(int i=0;i<50&&d>1e-5&&z<9e3;++i) {
    p=z*R+P;
    d=df(p,oo);
    if(p.y>0.) {
      O+=BM+d*Y;
    } else {
      O+=S;
      oo.x*=1e3;
    }
    O+=step(oo.z*.8,abs(p.y))*1e-1/max(oo.x,1e-2)*BO;
    O+=step(abs(p.y),oo.w*oo.z*.8)*1e-1/max(oo.x,1e-2)*BB;

    z+=d*.7;
  }

  O*=9E-3;
  O=tanh_approx(O);
  O=sqrt(O);
  return vec4(O,1);
}