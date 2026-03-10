const N0 = [0, 1/Math.SQRT2,  1/Math.SQRT2];
const N1 = [0, 1/Math.SQRT2, -1/Math.SQRT2];
const GRAVITY = [0, -9.8, 0];
const RESTITUTION = .98;
const RADIUS = 1.0;
const RANDOM_STRENGTH = 1.2;
const SPIN_DECAY = .5; // higher = faster decay; half-life ~= ln(2)/SPIN_DECAY

var pos = [0, 5, 0];
var vel = [0, 0, 3];
var angle0 = 0;
var angle1 = 0;
var spin0 = 0; // angular velocity for angle0
var spin1 = 0; // angular velocity for angle1

function dot(a, b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

function reflect(v, n) {
  const d = dot(v, n);
  return [v[0] - 2*d*n[0], v[1] - 2*d*n[1], v[2] - 2*d*n[2]];
}

function randomTangentKick(normal) {
  const up = Math.abs(normal[1]) < 0.9 ? [0,1,0] : [1,0,0];
  const d0 = dot(up, normal);
  const t0 = [up[0]-d0*normal[0], up[1]-d0*normal[1], up[2]-d0*normal[2]];
  const l0 = Math.sqrt(dot(t0,t0));
  const T0 = [t0[0]/l0, t0[1]/l0, t0[2]/l0];
  const T1 = [
    normal[1]*T0[2] - normal[2]*T0[1],
    normal[2]*T0[0] - normal[0]*T0[2],
    normal[0]*T0[1] - normal[1]*T0[0],
  ];
  const r0 = (Math.random()-0.5)*2;
  const r1 = (Math.random()-0.5)*2;
  const r2 = (Math.random()-0.5)*2;
  const r3 = (Math.random()-0.5)*2;
  return [
    0,
    RANDOM_STRENGTH * (r0*T0[1] + r1*T1[1]),
    RANDOM_STRENGTH * (r2*T0[2] + r3*T1[2]),
  ];
}

function collidePlane(normal, w) {
  const dist = dot(pos, normal) + w;
  if (dist < RADIUS) {
    const penetration = RADIUS - dist;
    pos[0] += normal[0] * penetration;
    pos[1] += normal[1] * penetration;
    pos[2] += normal[2] * penetration;
    if (dot(vel, normal) < 0) {
      vel = reflect(vel, normal);
      vel[0] *= RESTITUTION;
      vel[1] *= RESTITUTION;
      vel[2] *= RESTITUTION;
      const kick = randomTangentKick(normal);
      vel[0] += kick[0];
      vel[1] += kick[1];
      vel[2] += kick[2];

      spin0 += 4*kick[1];
      spin1 += 4*kick[2];
    }
  }
}

function update(dt) {
  vel[0] += GRAVITY[0] * dt;
  vel[1] += GRAVITY[1] * dt;
  vel[2] += GRAVITY[2] * dt;
  pos[0] += vel[0] * dt;
  pos[1] += vel[1] * dt;
  pos[2] += vel[2] * dt;
  collidePlane(N0, 1);
  collidePlane(N1, 1);

  // Exponential decay of spin
  const decay = Math.exp(-SPIN_DECAY * dt);
  spin0 *= decay;
  spin1 *= decay;

  // Integrate angles (wrap to avoid float drift over time)
  angle0 = (angle0 + spin0 * dt) % (2 * Math.PI);
  angle1 = (angle1 + spin1 * dt) % (2 * Math.PI);

  setUniform('u_sphere_pos', pos[0], pos[1], pos[2]);
  setUniform('u_angle', angle0, angle1);
}