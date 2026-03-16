#define ROT(a)      mat2(cos(a), sin(a), -sin(a), cos(a))

const float
  TAU=2.*PI
;

float fft(float x) {
#ifdef KODELIFE
  return .5+.5*sin(x*.2);
#else
  return textureLod(syn_Spectrum,.05+.85*(.5+.5*sin(x*.2)),0).y;
#endif
}

#ifdef KODELIFE
mat3 rotationFromAxisAngle(vec3 axis, float angle) {
  float
    c = cos(angle)
  , s = sin(angle)
  , t = 1. - c
  ;
  vec3
    a = normalize(axis)
  ;
  return mat3(
    t*a.x*a.x + c,      t*a.x*a.y - s*a.z,  t*a.x*a.z + s*a.y,
    t*a.x*a.y + s*a.z,  t*a.y*a.y + c,      t*a.y*a.z - s*a.x,
    t*a.x*a.z - s*a.y,  t*a.y*a.z + s*a.x,  t*a.z*a.z + c
  );
}
#endif


mat3 angle() {
#ifdef KODELIFE
  return rotationFromAxisAngle(vec3(1,0,0),TIME);
#else
  return mat3(u_rot_x, u_rot_y, u_rot_z);
#endif
}

// License: Unknown, author: Unknown, found: don't remember
float hash(float co) {
  return fract(sin(co*12.9898) * 13758.5453);
}

// License: Unknown, author: Unknown, found: don't remember
float hash(vec3 r)  {
  return fract(sin(dot(r.xy,vec2(1.38984*sin(r.z),1.13233*cos(r.z))))*653758.5453);
}
vec3 sphere_pos() {
#ifdef KODELIFE
  float
    B=fract(TIME*.5)-.5
  ;
  B*=B;
  B=(1.-4.*B)*3.-3.;
  return vec3(0,B,0);
#else
  return u_sphere_pos;
#endif
}

float hash(vec2 co) {
  return fract(sin(dot(co.xy ,vec2(12.9898,58.233))) * 13758.5453);
}

vec3 point_on_sphere(vec2 r) {
  r=vec2(PI*2.*r.x, 2.*r.y-1.);
  return vec3(sqrt(1. - r.y * r.y) * vec2(cos(r.x), sin(r.x)), r.y);
}

vec3 uniform_lambert(vec2 r, vec3 n) {
  return normalize(n*(1.001) + point_on_sphere(r)); // 1.001 required to avoid NaN
}

vec3 noisy_ray_dir(vec2 r, vec2 uv, vec3 cam_right, vec3 cam_up, vec3 cam_fwd) {
  uv += (-1. + 2.*r) / RENDERSIZE.y;
  return normalize(-uv.x*cam_right + uv.y*cam_up + 2.*cam_fwd);
}

// License: Unknown, author: Unknown, found: don't remember
vec2 hash2(float co) {
  return fract(sin(co*vec2(12.9898,78.233))*43758.5453);
}

vec4 amiga(mat3 R, vec3 p) {
  p*=R;
  vec2 pp=.5+.5*vec2(p.z,.25*p.y/p.x);
  vec4 tcol=textureLod(passAmiga, pp, 0.);
  return tcol;
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/articles/intersectors/
float ray_unitsphere(vec3 ro, vec3 rd) {
  float
    b=dot(ro, rd)
  , c=dot(ro, ro)-1.
  , h=b*b-c
  ;
  if(h<.0) return -1.;
  return -b-sqrt(h);
}

float ray_isphere4(vec3 ro, vec3 rd) {
  float
    b = dot(ro, rd)
  , c = dot(ro, ro)-16.
  , h = b*b-c
  , t
  ;
  if (h<0.) return -1.;

  return -b+sqrt(h);
}

float L4(vec3 p) {
  return sqrt(length(p*p));
}

vec4 pass_main() {
  vec3
    ro        = vec3(3,1,3.)
  ;
  ro.xz*=ROT(rotation_speed_0);
  ro.xy*=ROT(rotation_speed_1);
  vec3
    la        = vec3(0,-2,0.)
  , cam_fwd   = normalize(la - ro)
  , cam_right = normalize(cross(cam_fwd, vec3(0,1,0)))
  , cam_up    = cross(cam_right, cam_fwd)
  ;

  vec2 p=2.*_uvc,xy=_xy;

  float
    samples = 0.
  , fresnel
  , fade_out=.005
  , t
  , throughput
  , t_sphere
  , t_isphere
  , seed
  , h0
  , h1
  , d
  , f
  ;

  vec2
    r
  ;

  mat3 R=angle();

  seed = fract(hash(p) + TIME/1337.);

  vec3
    sphere_center = sphere_pos()
  , col           = vec3(0)
  , pos
  , prev_pos
  , prev_normal
  , iprev_normal
  , normal
  , reflect_dir
  , diffuse_dir
  ;

  vec4
    acol
  ;


  bool
    missed
  , hit_amiga
  , hit_grid
  , hit_fft
  ;

  prev_pos    = ro;
  prev_normal = noisy_ray_dir(r, p, cam_right, cam_up, cam_fwd);
  throughput  = 1.;

  for(int i=0; i<50; ++i) {
    ++seed;
    r=hash2(seed);
    iprev_normal=1./prev_normal;

    t_sphere  = ray_unitsphere(prev_pos - sphere_center, prev_normal);
    t_isphere = ray_isphere4(prev_pos, prev_normal);

    t = 1e3;
    if(t_sphere>0. && t_sphere<t)   { t=t_sphere      ; normal = (prev_pos - sphere_center + prev_normal*t); }
    if(t_isphere>0. && t_isphere<t) { t=t_isphere     ; normal = -.25*(prev_pos + prev_normal*t_isphere); }

    pos = prev_pos + prev_normal*t;

    vec3
      n = floor(pos+.5)
    , c = abs(pos-n)
    ;
    h0=hash(n+.123);
    h1=fract(8667.*h0);

    fresnel = 1. + dot(prev_normal, normal);
    f=fft(h0)*main_bass_gain*mix(min_bass_gain,max_bass_gain,syn_BassHits);
    f*=f;
    f+=.2;
    f=clamp(f,0.,1.);
    d=(L4(c)-.4*f);

    missed      = t==1e3 || throughput<1e-1;
    hit_amiga   = t==t_sphere ? (acol=amiga(R, pos-sphere_center), true) : false;
    f=smoothstep(low_bass_edge,high_bass_edge,f);
    
    hit_amiga = hit_amiga && r.y*r.y>fresnel;
    hit_grid    = t==t_isphere && (c.x<.01||c.y<.01||c.z<.01) && r.x>.5;
    hit_fft     = t==t_isphere && (abs(d)<.01 || f*5e-4/max(d*mix(d,abs(d),step(.8,f)),1e-4)>r.x) ;
    if(i==0 && missed) {
      break;
    }

    if(missed || hit_amiga || hit_grid ||hit_fft) {
      throughput/=(1.+fade_out*t*t);

      vec3 tc=vec3(0);
      if (hit_grid) {
        tc += .5*(1.+sin(3.+.5*pos.y+vec3(2,1,0)));
      }

      if (hit_fft) {
        tc += (1.+f+sin(pos.y+TIME+0.*h0+vec3(0,1,2)))-d;
      }

      if (hit_amiga) {
        tc += acol.xyz*sqrt(max(0.,-dot(prev_normal,normal)));
      }

      col += throughput*tc;

      prev_pos    = ro;
      prev_normal = noisy_ray_dir(r, p, cam_right, cam_up, cam_fwd);
      throughput  = 1.;
      ++samples;
      continue;
    }

    fresnel *= fresnel;
    fresnel *= fresnel;
    fresnel *= fresnel;

    reflect_dir = reflect(prev_normal, normal);
    diffuse_dir = uniform_lambert(r, reflect_dir);
    diffuse_dir = normalize(mix(diffuse_dir, reflect_dir, .9));

    if(
         r.x < fresnel
      || t==t_sphere
      ) {
      prev_normal = reflect_dir;
      throughput *= .9;
    } else {
      prev_normal = diffuse_dir;
      throughput *= .3;
    }

    throughput/=(1.+fade_out*t*t);

    prev_pos = pos + 1e-3*normal;
  }

  col /= max(samples, 1.);
  return vec4(col, 1.);
}

vec4 pass_amiga() {
  if(FRAMECOUNT==1) {
    vec2
      q=_uv
    , p=-1.+2.*q
    , S=vec2(acos(p.x),atan(p.y*4.))
    ;

    float
      f=sin(8.*S.y)*sin(8.*S.x)
    ;

    vec3
      amiga=mix(vec3(1.), vec3(1,.01,.01), step(f,0.))
    ;

    return vec4(amiga,step(f,.0));
  } else {
    return texelFetch(passAmiga, ivec2(_xy), 0);
  }
}


float dot2(vec3 p) {
  return dot(p,p);
}

vec3 denoise(ivec2 xy) {
  const int
    MAX=3
  ;
  const float
    DIV=1./float(MAX*MAX)
  ;
  vec3
    center  = texelFetch(passMain, xy, 0).xyz
  , s
  , sum     = center
  ;
  float
    rangeW
  , spatialW
  , w
  , weight = 1.0
  ;
  for(int dy = -MAX; dy <= MAX; ++dy)
  for(int dx = -MAX; dx <= MAX; ++dx) {
    if(dx == 0 && dy == 0) continue;

    s = texelFetch(passMain, xy + ivec2(dx, dy), 0).xyz;
    w=exp(-float(dx*dx + dy*dy) * DIV-dot2(s - center) *1E1);
    sum += s * w;
    weight += w;
  }

  return sum / weight;
}


vec4 pass_denoise() {
  ivec2
    xy=ivec2(_xy)
  ;
  vec3
    col=denoise(xy)
  , pcol=texelFetch(passDenoise, xy,0).xyz
  ;
  col=mix(col,pcol,.1);
  return vec4(col,1);
}

// License: Unknown, author: Matt Taylor (https://github.com/64), found: https://64.github.io/tonemapping/
vec3 aces_approx(vec3 v) {
  const float
    a = 2.51
  , b = 0.03
  , c = 2.43
  , d = 0.59
  , e = 0.14
  ;
  v = max(v, 0.);
  v *= .6;
  return clamp((v*(a*v+b))/(v*(c*v+d)+e), 0., 1.);
}

// From: https://www.shadertoy.com/view/XtlSD7
vec2 crt_distort(vec2 q) {
  q = _uv*2. - 1;
  vec2 
    o = crt_effect*q.yx/vec2(6,4)
  ;
  q = q + q*o*o;
  return q*.5 + .5;
}

// From: https://www.shadertoy.com/view/XtlSD7
float vig(vec2 q) {    
  float 
    v = q.x*q.y*(1.0 - q.x)*(1.0 - q.y)
  ;
  v = clamp(pow(16.*v, .3), 0., 1.);
  return v;
}


vec4 pass_post() {
  ivec2
    xy=ivec2(_xy)
  ;
  
  vec2
    q=crt_distort(_uv)
  , s=step(abs(q-.5),vec2(.5))
  ;
  
  vec3
    col=textureLod(passDenoise, q, 0.).xyz
  ;
  col *= mix(1.,vig(q),crt_effect);
  col -=1e-2*vec3(2,3,1);
  col=max(col,0.);
  col *= 1.5;
  col=aces_approx(col);
  col=sqrt(col)-.05;

  col*=mix(1.,1.5+.5*sin(_xy.y*TAU/6.), crt_effect);

#ifndef KODELIFE
  vec4 mcol=_loadMedia();
  col=mix(col,mcol.xyz,media_opaque*mix(dot(mcol.xyz,vec3(0.299, 0.587, 0.114)), mcol.w, mix_mode));
#endif


  return vec4(col,1);
}

vec4 renderMain() {
  switch(PASSINDEX) {
  case 0:
    return pass_amiga();
  case 1:
    return pass_main();
  case 2:
    return pass_denoise();
  default:
    return pass_post();
  }
}