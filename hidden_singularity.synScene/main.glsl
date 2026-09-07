// License: Unknown, author: Claude Brezinski, found: https://mathr.co.uk/blog/2017-09-06_approximating_hyperbolic_tangent.html
vec3 tanh_approx(vec3 x) {
  //  Found this somewhere on the interwebs
  //  return tanh(x);
  vec3 x2 = x*x;
  return clamp(x*(27.0 + x2)/(27.0+9.0*x2), -1.0, 1.0);
}

vec4 renderMain() {
  vec2
    r=RENDERSIZE
  , P=2.*_uvc
  ;
  
  float
    d
  , i=0.
  , k
  , K
  , T=syn_BassTime
  , z=0.
  , f=7.+4.*sin(T*PI/16.)
  ;
  
  vec4 
    p
  , o=vec4(0)
  , O
  , L=vec4(syn_BassTime, syn_MidTime, syn_MidHighTime, syn_HighTime)
  , S=sin(L*PI/32.)
  ;
  
  vec3
    R=normalize(tan(PI*(L.yzw)/32.))
  , D=normalize(vec3(_xy-.5*r,r.y))
  ;
  
  for(
      i=0.
    ; i<77.&&z<49.
    ; ++i
  )
    p=z*vec4(D,dot(D,vec3(1,shape_dist)))-vec4(0,0,9,0)
  , p.xyz=R*dot(R,p.xyz)+cross(R,p.xyz)
  , p*=k=7./dot(p,p)
  , p+=S
  , p*=K=f/dot(p,p)
  , k*=K
  , O=1.+sin(log2(k)+T+vec4(0,1,2,0))
  , z+=d=.7*abs((length(fract(p.xy+.5)-.5)-.03*p.w*p.w-.01)/k)+1e-3
  , o+=O.w/d*O+20.*K*syn_BassHits*vec4(1,2,3,0)
  ;
  
  o/=2e4;
  o+=flash_bang*vec4(1,2,3,0)*(1.+length(P));
  return vec4(tanh_approx(syn_FadeInOut*o.xyz),1);
}
