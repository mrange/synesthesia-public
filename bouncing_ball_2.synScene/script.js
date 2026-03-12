var pos = [0, 2, 0];
var vel = [1, 0, 0.5];
var angVel = [0, 0, 0];

// Rotation stored as axis-angle: [ax, ay, az, angle]
var rotAxis = [1, 0, 0];
var rotAngle = 0.0;

const INNER_R = 1.0;
const OUTER_R = 4.0;
const BOUNCE_R = OUTER_R - INNER_R;
const GRAVITY = 9.8;
const DAMPING = 1.0;
const ANG_DAMPING = 0.995;

function len(v) {
  return Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}
function cross(a, b) {
  return [
    a[1]*b[2] - a[2]*b[1],
    a[2]*b[0] - a[0]*b[2],
    a[0]*b[1] - a[1]*b[0]
  ];
}
// Compose two axis-angle rotations into one using Rodrigues
function composeRotations(ax1, an1, ax2, an2) {
  // Convert both to quaternions, multiply, convert back
  var s1 = Math.sin(an1/2), c1 = Math.cos(an1/2);
  var s2 = Math.sin(an2/2), c2 = Math.cos(an2/2);
  var q1 = [ax1[0]*s1, ax1[1]*s1, ax1[2]*s1, c1];
  var q2 = [ax2[0]*s2, ax2[1]*s2, ax2[2]*s2, c2];
  // Quaternion multiply q1 * q2
  var qx = q1[3]*q2[0] + q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1];
  var qy = q1[3]*q2[1] - q1[0]*q2[2] + q1[1]*q2[3] + q1[2]*q2[0];
  var qz = q1[3]*q2[2] + q1[0]*q2[1] - q1[1]*q2[0] + q1[2]*q2[3];
  var qw = q1[3]*q2[3] - q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2];
  // Back to axis-angle
  var sinHalf = Math.sqrt(qx*qx + qy*qy + qz*qz);
  if (sinHalf < 0.0001) return { axis: [1,0,0], angle: 0 };
  return {
    axis: [qx/sinHalf, qy/sinHalf, qz/sinHalf],
    angle: 2 * Math.atan2(sinHalf, qw)
  };
}

function update(dt) {
  vel[1] -= GRAVITY * dt;

  pos[0] += vel[0] * dt;
  pos[1] += vel[1] * dt;
  pos[2] += vel[2] * dt;

  // Accumulate rotation: compose current rotation with this frame's delta
  var w = len(angVel);
  if (w > 0.0001) {
    var deltaAxis = [angVel[0]/w, angVel[1]/w, angVel[2]/w];
    var deltaAngle = w * dt;
    var r = composeRotations(rotAxis, rotAngle, deltaAxis, deltaAngle);
    rotAxis  = r.axis;
    rotAngle = r.angle;
  }

  var dist = Math.sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
  if (dist > BOUNCE_R) {
    var nx = pos[0] / dist;
    var ny = pos[1] / dist;
    var nz = pos[2] / dist;

    pos[0] = nx * BOUNCE_R;
    pos[1] = ny * BOUNCE_R;
    pos[2] = nz * BOUNCE_R;

    var vDotN = vel[0]*nx + vel[1]*ny + vel[2]*nz;
    vel[0] = (vel[0] - 2*vDotN*nx) * DAMPING;
    vel[1] = (vel[1] - 2*vDotN*ny) * DAMPING;
    vel[2] = (vel[2] - 2*vDotN*nz) * DAMPING;

    // Spin from tangential velocity at bounce
    var vDotN2 = vel[0]*nx + vel[1]*ny + vel[2]*nz;
    var vTan = [
      vel[0] - vDotN2*nx,
      vel[1] - vDotN2*ny,
      vel[2] - vDotN2*nz
    ];
    var spin = cross([nx, ny, nz], vTan);
    angVel[0] = spin[0] / INNER_R;
    angVel[1] = spin[1] / INNER_R;
    angVel[2] = spin[2] / INNER_R;
  }

  angVel[0] *= ANG_DAMPING;
  angVel[1] *= ANG_DAMPING;
  angVel[2] *= ANG_DAMPING;

  // Pass axis (xyz) and angle as 4 values split across two uniforms
  setUniform('u_sphere_pos', pos[0], pos[1], pos[2]);
  setUniform('u_angle', rotAxis[0], rotAxis[1], rotAxis[2], rotAngle);
}