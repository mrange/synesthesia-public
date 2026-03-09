#define ROT(a)      mat2(cos(a), sin(a), -sin(a), cos(a))

const float
  TAU=2.*PI
;

vec3 sphere_pos() {
  return u_sphere_pos;
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

vec4 pMain() {
  const float off=4.;

  const vec3
    ro        = vec3(-5,3,0.)
  , la        = vec3(0,2,0.)
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
  , t_sphere
  , seed
  ;

  vec2
    r
  ;

  mat2 R0=ROT(angle0);
  mat2 R1=ROT(angle1);

  seed = fract(hash(p) + TIME/1337.);


  vec3
    sphere_center = sphere_pos()
  , col           = vec3(0)
  , pos
  , prev_pos
  , prev_normal
  , normal
  , reflect_dir
  , diffuse_dir
  , acol
  ;


  bool
    missed
  , hit_amiga
  ;

  prev_pos    = ro;
  prev_normal = noisy_ray_dir(r, p, cam_right, cam_up, cam_fwd);
  throughput  = 1.;

  for(int i=0; i<100; ++i) {
    ++seed;
    r=hash2(seed);

    t_floor_0 = ray_plane(prev_pos, prev_normal, floor_0);
    t_floor_1 = ray_plane(prev_pos, prev_normal, floor_1);
    t_wall_0  = ( off - prev_pos.z) / prev_normal.z;
    t_wall_1  = (-off - prev_pos.z) / prev_normal.z;
    t_sphere  = ray_unitsphere(prev_pos - sphere_center, prev_normal);

    t = 1e3;
    if(t_floor_0>0.&& t_floor_0<t){ t=t_floor_0;  normal=floor_0.xyz; }
    if(t_floor_1>0.&& t_floor_1<t){ t=t_floor_1;  normal=floor_1.xyz; }
    if(t_wall_0>0. && t_wall_0<t) { t=t_wall_0; normal=vec3(0,0,-1); }
    if(t_wall_1>0. && t_wall_1<t) { t=t_wall_1; normal=vec3(0,0,1); }
    if(t_sphere>0. && t_sphere<t) { t=t_sphere; normal=normalize(prev_pos+prev_normal*t_sphere-sphere_center);}

    pos = prev_pos + prev_normal*t;

    missed      = t==1e3 || throughput<1e-1;
    hit_amiga   = t==t_sphere ? (acol=amiga(R0, R1, pos-sphere_center), true) : false;
    //hit_amiga   = hit_amiga && acol.z<0.1 && r.x>.5;

    if(i==0 && missed) {
      break;
    }

    if(missed || hit_amiga) {
      throughput/=(1.+fade_out*t*t);
      if (hit_amiga) {
        col += throughput*acol*sqrt(max(0.,-dot(prev_normal,normal)));
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

    reflect_dir = reflect(prev_normal, normal);
    diffuse_dir = uniform_lambert(r, normal);

    if(
        r.x < fresnel
      || t==t_sphere
      ) {
      prev_normal = reflect_dir;
      throughput *= .8;
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

vec4 pAmiga() {
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


vec4 pDenoise() {
  ivec2
    xy=ivec2(_xy)
  ;
  vec3 
    col=denoise(xy)
  , pcol=texelFetch(passDenoise, xy,0).xyz
  ;
  col=mix(col,pcol,0.2);
  return vec4(col,1);
}

vec4 pPost() {
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
    return pAmiga();
  case 1:
    return pMain();
  case 2:
    return pDenoise();
  default:
    return pPost();
  }
}