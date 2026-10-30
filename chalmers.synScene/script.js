function mix(a,b,c) {
  return (1-c)*a+c*b;
}

function smoothstep(edge0, edge1, x) {
  var t=(x-edge0)/(edge1-edge0);
  t=t<0?0:(t > 1 ? 1 : t);
  return t*t*(3-2*t);
}

function vdot(a,b) {
  var s=0;
  for(var i=0;i<a.length;++i) {
    s+=a[i]*b[i];
  }
  return s;
}

function vdot2(a) {
  var s=0;
  for(var i=0;i<a.length;++i) {
    s+=a[i]*a[i];
  }
  return s;
}

function vl2(a) {
  return Math.sqrt(vdot2(a));
}

function vzero(a) {
  const c = new Float32Array(a.length);
  return c;
}

function vscale(s,a) {
  const c = vzero(a);
  for(var i=0;i<c.length;++i) {
    c[i]=s*a[i];
  }
  return c;
}

function vdiv(a,b) {
  const c = vzero(a);
  for(var i=0;i<c.length;++i) {
    c[i]=a[i]/b[i];
  }
  return c;
}

function vsin(a) {
  const c = vzero(a);
  for(var i=0;i<c.length;++i) {
    c[i]=Math.sin(a[i]);
  }
  return c;
}

function vsmoothstep(a,b,x) {
  const c = vzero(a);
  for(var i=0;i<c.length;++i) {
    c[i]=smoothstep(a[i],b[i],x[i]);
  }
  return c;
}

function vnormalize(a) {
  return vscale(1/vl2(a),a);
}

function v3(a,b,c) {
  const v = new Float32Array(3);
  a=a===undefined?0:a;
  b=b===undefined?a:b;
  c=c===undefined?b:c;
  v[0]=a;
  v[1]=b;
  v[2]=c;
  return v;
}

var vperiods = v3(11,13,17);

var t =0;
function update(dt) {
  const tt=speed*mix(2.*syn_BassLevel*dt,dt, speed_mix);
  t+=tt;
  const vrotor  = vnormalize(vsin(vdiv(v3(t+100),vperiods)));
  setUniform("effect_time", t);
  setUniform("u_rotor", vrotor[0], vrotor[1], vrotor[2]);
}
