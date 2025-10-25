// This file is released under CC0 1.0 Universal (Public Domain Dedication).
// To the extent possible under law, Mårten Rånge has waived all copyright
// and related or neighboring rights to this work.
// See <https://creativecommons.org/publicdomain/zero/1.0/> for details.

//#define IMAGE_SAMPLE

vec4 image(vec2 p, vec2 tz) {
  float
    fade = smoothstep(15., 1., abs(p.x))*smoothstep(9., 1., p.y)
  ;
  p *= tz;
  p.x+=.1*TIME-.125*
  floor(p.y);
  p.x = fract(p.x);
  if (fade < 0.1 || p.y < 00.0)
    return vec4(0.0);
  p.y = fract(p.y);
#ifdef KODELIFE
  p.y=1.-p.y;
#endif
  vec4 col = texture(syn_Media, p);
  col.xyz*=col.xyz;
  return vec4(col.xyz, col.w*fade);
}

void rot(inout vec2 p, float a) {
  float
    c=cos(a)
  , s=sin(a)
  ;
  p = vec2(c*p.x+s*p.y,c*p.y-s*p.x);
}


vec4 alphaBlend(vec4 back, vec4 front) {
  // Based on: https://en.wikipedia.org/wiki/Alpha_compositing
  float w = front.w + back.w*(1.0-front.w);
  vec3 xyz = (front.xyz*front.w + back.xyz*back.w*(1.-front.w))/w;
  return w > 0.0 ? vec4(xyz, w) : vec4(0.0);
}

float freq(float x) {
#ifdef KODELIFE
  return 0.;
#else
  return texture(syn_Spectrum,.5+.5*sin(x)).y;
#endif
}

vec4 renderMain() {
  const float
    maxIterF = 80.
  ;

#ifdef KODELIFE
  float
    tm = .5*TIME
#else
  float
    tm = dot(speed_control,vec2(TIME, syn_BassTime))
#endif
  , i
  , l
  ;
  vec2
    q   = _uv
  , p   = 2.*_uvc
  , op  = p
  , tsz = vec2(textureSize(syn_Media,0))
  , tz  = vec2(tsz.y/tsz.x, 1)
  , c = vec2(-0.76, 0.15)
  , z
  ;

  rot(p,.3*tm);
  p = vec2(0.5, -0.05) + p*0.75*pow(0.9, 20.0*(0.5+0.5*cos(0.3*sqrt(2.0)*tm)));

  vec3
    ss = mix(vec3(0.2, 0.2, 0.5), vec3(0.2,-0.2,1.0), 2.2 + 1.25*sin(tm/2.0))
  , col = vec3(0)
  ;

  z=p;
  rot(c, 0.2*sin(tm*sqrt(3.0)/6.0));

  for(i=0.; i<maxIterF; ++i) {
    float
      re2  = z.x*z.x
    , im2  = z.y*z.y
    , reim = z.x*z.y
    ;
    if ((re2+im2)>20.) break;

    z = vec2(re2 - im2, 2.*reim) + c;

    float shade = smoothstep(maxIterF, 0., i);
#ifdef IMAGE_SAMPLE
    vec4 img = image(ss.xy + ss.z*z,tz);
    img.w*=tanh(shade);
    col=alphaBlend(img, col);
#endif
    shade /= re2+im2+0.01;
    col+=shade*3e-3*(1.+sin(vec3(4.*z, 0.2*i)+vec3(0,color_control)))/abs(z.y+.1*(freq(2.*z.x)-.5));
//    bg+=shade*3e-3*(1.+sin(vec3(4.*z, 0.2*i)+vec3(0,1,2)))/abs(z.y+sin(TIME+10.*z.x)-.2);
//    bg+=shade*3e-3*(1.+sin(vec3(5.*z, 0.1*i)+vec3(2,1,0)))/abs(z.y*z.y-.2);
//    bg+=shade*3e-3*(1.+sin(vec3(5.*z, 0.1*i)+vec3(2,1,0)))/abs(z.x*z.x-.3);
//  bg+=shade*3e-3*(1.+sin(vec3(6.*z, 0.1*i)+vec3(1,2,0)))/abs(z.x*z.y-.3);
  }
  l = i - log2(log2(dot(z,z)));

  col /= 9.;
  if (!isnan(l)) {
    col += (0.5 + 0.5*cos(l*.5 + 3.+ vec3(0,2.*color_control)))/(1.+l*l);
  }
  col = tanh(col);
  col = mix(col, 1.-col, invert_color);

  col = max(col,0.);
  vec4 pcol = texture(syn_FinalPass, q);
  col.xyz = mix(col.xyz, pcol.xyz, mix(motion_blur_intensity.x, motion_blur_intensity.y, smoothstep(motion_blur_limits.x, motion_blur_limits.y, length(p))));

#ifdef KODELIFE
#else
  vec4 mcol=_loadMedia(media_warp.x*normalize(p)*freq(media_warp.y*length(p)));
  col=mix(col,mcol.xyz,mcol.w*media_opacity*media_multiplier);
#endif
  
  return vec4(col.xyz, 1);
}

