open System
open System.Numerics

let inline fract1 (v : float32) = v - floor v
let inline fract3 (v : Vector3) = Vector3 (fract1 v.X, fract1 v.Y, fract1 v.Z)

// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
let hsv2rgb (h : float32) (s : float32) (v : float32) : Vector3 =
  let zero  = Vector3.Zero
  let one   = Vector3.One
  let three = 3.F*one
  let xyz   = Vector3(3.F, 2.F, 1.F)/3.F
  v*Vector3.Lerp(
      one
    , Vector3.Clamp(
        Vector3.Abs(6.F*(fract3 (Vector3(h)+xyz))-three)-one
      , zero
      , one
      )
    , s
    )

let OKLAB_M1 =
  Matrix4x4(
    0.4122214708f, 0.2119034982f, 0.0883024619f, 0.f,
    0.5363325363f, 0.6806995451f, 0.2817188376f, 0.f,
    0.0514459929f, 0.1073969566f, 0.6299787005f, 0.f,
    0.f          , 0.f          , 0.f          , 1.f)
    |> Matrix4x4.Transpose

let OKLAB_M2 =
  Matrix4x4(
    +0.2104542553f, +1.9779984951f, +0.0259040371f, +0.f,
    +0.7936177850f, -2.4285922050f, +0.7827717662f, +0.f,
    -0.0040720468f, +0.4505937099f, -0.8086757660f, +0.f,
    +0.f          , +0.f          , +0.f          , +1.f)
    |> Matrix4x4.Transpose

let linearToOklab (c : Vector3) : Vector3 =
  let cubeRoot (x : float32) = float32 (Math.Pow(float x, 1.0/3.0))
  let x = Vector3.Transform(c,OKLAB_M1)
  let y = Vector3(cubeRoot x.X, cubeRoot x.Y, cubeRoot x.Z)
  Vector3.Transform(y,OKLAB_M2)

let inline HSV2RGB (h, s, v) =
  let rgb = hsv2rgb (float32 h) (float32 s) (float32 v)
  printfn "HSV2RGB(%.2f,%.2f,%.2f): vec3(%.2f, %.2f, %.2f)" h s v rgb.X rgb.Y rgb.Z

let inline HSV2OKLAB (h, s, v) =
  let rgb = hsv2rgb (float32 h/360.F) (float32 s/100.F) (float32 v/100.F)
  // Square as a fake reverse sRGB
  let ok  = linearToOklab (rgb*rgb)
  printfn "HSV2OKLAB(%.2f,%.2f,%.2f): vec3(%.2f, %.2f, %.2f)" h s v ok.X ok.Y ok.Z


HSV2RGB (0.58 , 0.7 , 1.0)
HSV2RGB (0.63 , 0.9 , 1.0)


HSV2OKLAB(222.,57.,33.)
HSV2OKLAB(218.,45.,45.)
HSV2OKLAB(205.,30.,54.)
HSV2OKLAB(021.,15.,62.)
HSV2OKLAB(038.,30.,68.)
HSV2OKLAB(028.,42.,74.)
HSV2OKLAB(015.,55.,80.)
HSV2OKLAB(035.,50.,50.)
HSV2OKLAB(190.,50.,80.)
