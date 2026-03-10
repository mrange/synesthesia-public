#define ROT(a)      mat2(cos(a), sin(a), -sin(a), cos(a))

const float
  TAU=2.*PI
;

vec2 angle() {
#ifdef KODELIFE
  return vec2(TIME,0.324*TIME);
#else
  return u_angle;
#endif
}

// License: Unknown, author: Unknown, found: don't remember
float hash(float co) {
  return fract(sin(co*12.9898) * 13758.5453);
}

vec3 sphere_pos() {
#ifdef KODELIFE
  float
    B=fract(TIME*.5)-.5
  ;
  B*=B;
  B=(1.-4.*B)*4.;
  return vec3(0,B,0);
#else
  return u_sphere_pos;
#endif
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

float ray_plane(vec3 ro, vec3 rd, vec4 p) {
  return -(dot(ro,p.xyz)+p.w)/dot(rd,p.xyz);
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

vec3 amiga(mat2 R0, mat2 R1, vec3 p) {
  p.xy*=R1;
  p.zx*=R0;
  vec3 s=sign(p);
  vec2 pp=vec2(.5+.5*(p.z),.5+.5*(p.y/p.x)*.25);
  return textureLod(passAmiga, pp,0.).xyz;
}

float L4(vec3 p) {
  return sqrt(length(p*p));
}
float L4(vec2 p) {
  return sqrt(length(p*p));
}

float segment(vec2 p) {
  float
    d0=L4(p)
  , d1=abs(p.x)
  ;
  return p.y > 0.?d0:d1;
}

vec4 pass_main() {
  const float off=4.;

  const vec3
    ro        = vec3(-6,4,0.)
  , la        = vec3(0,3,0.)
  , cam_fwd   = normalize(la - ro)
  , cam_right = normalize(cross(cam_fwd, vec3(0,1,0)))
  , cam_up    = cross(cam_right, cam_fwd)
  ;

  const vec4
    floor_0=vec4(normalize(vec3(0,1,1)),1)
  , floor_1=vec4(normalize(vec3(0,1,-1)),1)
  ;

  vec2 p=2.*_uvc,xy=_xy;

  float
    samples = 0.
  , fresnel
  , fade_out=.005
  , t
  , throughput
  , t_wall_0
  , t_wall_1
  , t_floor_0
  , t_floor_1
  , t_roof
  , t_sphere
  , seed
  , fd
  , fn
  , fh
  , iz
  ;

  vec2
    r
  , a=angle()
  , wpos
  ;

  mat2 R0=ROT(a.x);
  mat2 R1=ROT(a.y);

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
  , acol
  ;


  bool
    missed
  , hit_amiga
  , hit_fft
  ;

  prev_pos    = ro;
  prev_normal = noisy_ray_dir(r, p, cam_right, cam_up, cam_fwd);
  throughput  = 1.;

  for(int i=0; i<100; ++i) {
    ++seed;
    r=hash2(seed);
    iprev_normal=1./prev_normal;

    t_floor_0 = ray_plane(prev_pos, prev_normal, floor_0);
    t_floor_1 = ray_plane(prev_pos, prev_normal, floor_1);
    t_wall_0  = ( off - prev_pos.z)*iprev_normal.z;
    t_wall_1  = (-off - prev_pos.z)*iprev_normal.z;
    t_roof    = (6. - prev_pos.y)*iprev_normal.y;
    t_sphere  = ray_unitsphere(prev_pos - sphere_center, prev_normal);

    t = 1e3;
    if(t_floor_0>0.&& t_floor_0<t){ t=t_floor_0 ; normal=floor_0.xyz; }
    if(t_floor_1>0.&& t_floor_1<t){ t=t_floor_1 ; normal=floor_1.xyz; }
    if(t_wall_0>0. && t_wall_0<t) { t=t_wall_0  ; normal=vec3(0,0,-1); }
    if(t_wall_1>0. && t_wall_1<t) { t=t_wall_1  ; normal=vec3(0,0,1); }
    if(t_roof>0.&& t_roof<t)      { t=t_roof    ; normal=vec3(0,-1,0); }
    if(t_sphere>0. && t_sphere<t) { t=t_sphere  ; normal=normalize(prev_pos+prev_normal*t_sphere-sphere_center);}

    pos = prev_pos + prev_normal*t;

    wpos = pos.xy+vec2(TIME,0);
    fn=floor(wpos.x+.5);
    fh=.5+.5*sin(fn*.2);
    wpos.x -= fn;
    
    fd=segment(wpos-vec2(0,2.3+2.*textureLod(syn_Spectrum,.05+.85*fh,0).y))-.4;

    missed      = t==1e3 || throughput<1e-1;
    hit_amiga   = t==t_sphere ? (acol=amiga(R0, R1, pos-sphere_center), true) : false;
    hit_fft     = (t==t_wall_0||t==t_wall_1)&&fd<0.;
    hit_amiga   = hit_amiga && acol.z<0.1 && r.x>.33;

    if(i==0 && missed) {
      break;
    }

    if(missed || hit_amiga || hit_fft) {
      throughput/=(1.+fade_out*t*t);
      if (hit_amiga) {
        col += throughput*acol*sqrt(max(0.,-dot(prev_normal,normal)));
      }
      if (hit_fft) {
        col += throughput*(fd+1.+sin(2.+3.*fh+vec3(2,1,0)));
      }

      prev_pos    = ro;
      prev_normal = noisy_ray_dir(r, p, cam_right, cam_up, cam_fwd);
      throughput  = 1.;
      ++samples;
      continue;
    }

    fresnel = 1. + dot(prev_normal, normal);
    fresnel *= fresnel;
    fresnel *= fresnel;
    fresnel *= fresnel;

    reflect_dir = reflect(prev_normal, normal);
    diffuse_dir = uniform_lambert(r, reflect_dir);
    diffuse_dir = normalize(mix(diffuse_dir, reflect_dir, (t==t_floor_0||t==t_floor_1)&&abs(sin(.5*pos.x+TIME)-pos.y+.5)<.5?.8:.0));
    if(
        r.x < fresnel
      || t==t_sphere
      || t==t_wall_0
      || t==t_wall_1
    //  || (t==t_floor_0||t==t_floor_1)&&abs(sin(.5*pos.x+TIME)-pos.y+.5)<.5
      ) {
      prev_normal = reflect_dir;
      throughput *= .7;
    } else {
      prev_normal = diffuse_dir;
      throughput *= .5;
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
  
    return vec4(amiga,1);
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
  col=mix(col,pcol,.3);
  return vec4(col,1);
}

vec4 pass_post() {
  ivec2
    xy=ivec2(_xy)
  ;
  vec3 
    col=texelFetch(passDenoise, xy,0).xyz
  ;
  col=max(col,0.);
  col *= 1.5;
  col=tanh(col);
  col=sqrt(col)-.05;
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