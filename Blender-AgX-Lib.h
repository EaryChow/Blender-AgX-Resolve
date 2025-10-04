/* libDTColorMath
    --------------------------------------

      This was originally from Jed smith, but I added math functions for my DCTLS
*/

/* ##########################################################################
    Constants
    ---------------------------
*/

__CONSTANT__ float pi = 3.14159265358979323846f;
__CONSTANT__ float pi2 = 1.57079632679489661923f;
__CONSTANT__ float pi4 = 0.785398163397448309616f;

/* ##########################################################################
    Custom Structs
    ---------------------------
*/

// Custom 3x3 matrix struct float3x3
typedef struct {
  float3 x, y, z;
} float3x3;

// Struct to hold whitepoint normalization
typedef struct {
  float rec2020, p3d65, p3d60, p3dci, rec709, dcdm;
} wpnorm;


// Struct for chromaticity coordinates of color spaces
typedef struct {
  float2 red; float2 green; float2 blue; float2 white;
} Chromaticities;

/* ##########################################################################
    Color Conversion Matrices
    ---------------------------
*/

#define identity_mtx make_float3x3(make_float3(1.0f, 0.0f, 0.0f), make_float3(0.0f, 1.0f, 0.0f), make_float3(0.0f, 0.0f, 1.0f))

// Color Spaces Coordinates
//AP0, AP1 and P3D60 was changed to D65. Method is, first arrive at RGB D60 to XYZ D65 matrix by "Bradford matrix * RGB to XYZ D60 matrix", and then reverse engineer the D65 version of the primaries coordinates from the matrix.
//This appraoch basically baked the Bradford adaptation matrix into the primaries coordinates relative to D65.
//P3DCI also adapted to D65 coordinates using the same approach.
__CONSTANT__ Chromaticities AP0 =
{ {0.7348552433737107f, 0.2642253252455353f}, {-0.006170912478622476f, 1.0113149590212862f}, {0.01596755925504054f, -0.0642355031285507f}, {0.3127f, 0.329f} };
//__CONSTANT__ Chromaticities AP0 =
//{ {0.7347f, 0.2653f}, {0.0f, 1.0f}, {0.0001f, -0.077f}, {0.32168f, 0.33767f} };
__CONSTANT__ Chromaticities AP1 =
{ {0.7131958876620503f, 0.2926889144633257f}, {0.15950855654177495f, 0.8387885161509631f}, {0.12867299528535026f, 0.04389557116052796f}, {0.3127f, 0.329f} };
__CONSTANT__ Chromaticities REC709_PRI =
{ {0.64f, 0.33f}, {0.3f, 0.6f}, {0.15f, 0.06f}, {0.3127f, 0.329f} };
__CONSTANT__ Chromaticities P3D60_PRI =
{ {0.6793813705302119f, 0.32015311540444413f}, {0.2597770990951463f, 0.6939530572700782f}, {0.14869419443105955f, 0.05860092634705156f}, {0.3127, 0.329f} };
__CONSTANT__ Chromaticities P3D65_PRI =
{ {0.68f, 0.32f}, {0.265f, 0.69f}, {0.15f, 0.06f}, {0.3127f, 0.329f} };
__CONSTANT__ Chromaticities P3DCI_PRI =
{ {0.680701401314748f, 0.3188954328985409f}, {0.28120700362868956f, 0.6741679571843063f}, {0.14883188869747185f, 0.05766741627278289f}, {0.3127f, 0.329f} };
__CONSTANT__ Chromaticities ARRI_ALEXA_WG_PRI =
{ {0.684f, 0.313f}, {0.221f, 0.848f}, {0.0861f, -0.102f}, {0.3127f, 0.329f} };
__CONSTANT__ Chromaticities REC2020_PRI =
{ {0.708f, 0.292f}, {0.17f, 0.797f}, {0.131f, 0.046f}, {0.3127f, 0.329f} };
__CONSTANT__ Chromaticities ARRI_ALEXA_WG4_PRI =
{ {0.7347f,  0.2653f}, {0.1424f, 0.8576f}, {0.0991f,-0.0308f}, {0.3127f, 0.329f} };
__CONSTANT__ Chromaticities SGAMUT3cine_PRI =
{ {0.766f,  0.275f}, {0.225f, 0.8f}, {0.089f, -0.087f}, {0.3127f, 0.329f} };
__CONSTANT__ Chromaticities SGAMUT3_PRI =
{ {0.73f,  0.28f}, {0.14f, 0.855f}, {0.1f, -0.05f}, {0.3127f, 0.329f} };
__CONSTANT__ Chromaticities EGAMUT_PRI =
{ {0.8f,  0.3177f}, {0.18f, 0.9f}, {0.065f, -0.0805f}, {0.3127f, 0.329f} };
__CONSTANT__ Chromaticities DWGINT_PRI =
{ {0.8f,  0.3130f}, {0.1682f, 0.9877f}, {0.079f, -0.1155f}, {0.3127f, 0.329f} };
__CONSTANT__ Chromaticities SONYA6000_PRI =
{ {0.709212f,  0.27432f}, {0.219737f, 0.961661f}, {0.023252f, -0.305702f}, {0.3127f, 0.329f} };
__CONSTANT__ Chromaticities REDWG_PRI =
{ {0.780308f, 0.30425299f}, {0.12159501f, 1.49399403f}, {0.095612f, -0.084589f}, {0.3127f, 0.329f} };
__CONSTANT__ Chromaticities BLACKMAGICWG_PRI =
{ {0.7177215143781699f, 0.3171180860551045f}, {0.2280409960285721f, 0.8615690087709122f}, {0.10058410063338172f, -0.08204520062464217f}, {0.3127f, 0.329f} };
__CONSTANT__ Chromaticities CANONCINEMA_PRI =
{ {0.7400000016535166f, 0.27000000186020623f}, {0.17000000222850284f, 1.1400000018352379f}, {0.07999999816658841f, -0.09999999770823552f}, {0.3127f, 0.329f} };
__CONSTANT__ Chromaticities AGX_LOG_BLENDER_PRI =
{ {0.9838378198847221f, 0.239412646694087f}, {0.09704998110120326f, 0.9900566739456024f}, {0.10813075783923522f, -0.0323324412351137f}, {0.3127f, 0.329f} };

/* Matrix for conversion from CIE 1931 XYZ tristumulus to CIE 2006 LMS or "Truelight LMS", described in:
  "Chromaticity Coordinates for Graphic Arts Based on CIE 2006 LMS with Even Spacing of Munsell Colours" by Richard Kirk
  https://doi.org/10.2352/issn.2169-2629.2019.27.38
*/
#define matrix_xyz_to_truelightlms make_float3x3(make_float3(0.257085f, 0.859943f, -0.031061f), make_float3(-0.394427, 1.175800f, 0.106423f), make_float3(0.064856f, -0.07625f, 0.559067f))

// Whitepoint scaling factors for Truelight LMS
// eg: (1, 1, 1) in RGB (D65 whitepoint) -> XYZ -> TLMS /= catd65 *= catd55 -> XYZ -> Yxy == D65 white
#define catd50 make_float3(1.08221f, 0.883260f, 0.447579f)
#define catd55 make_float3(1.07730f, 0.896467f, 0.500927f)
#define catd60 make_float3(1.07344f, 0.907523f, 0.549381f)
#define catd65 make_float3(1.07046f, 0.916817f, 0.594251f)
#define catd75 make_float3(1.06600f, 0.931715f, 0.670839f)
#define catd93 make_float3(1.06098f, 0.950462f, 0.776150f)

/* ##########################################################################
    Functions
    ---------------------------------
*/
__DEVICE__ Chromaticities make_chromaticities( float2 A, float2 B, float2 C, float2 D) {
  Chromaticities E;
  E.red = A; E.green = B; E.blue = C; E.white = D;
  return E;
}

__DEVICE__ float _radians(float d) {return d * (pi / 180.0f);}
__DEVICE__ float _degrees(float r) {return r * (180.0f / pi);}

__DEVICE__ float3 sqrtf3(float3 a) {
  // For each component of float3 a, compute the square-root
  return make_float3(_sqrtf(a.x), _sqrtf(a.y), _sqrtf(a.z));
}

__DEVICE__ float3 clampf3(float3 a, float mn, float mx) {
  // Clamp each component of float3 a to be between float mn and float mx
  return make_float3(
    _fminf(_fmaxf(a.x, mn), mx),
    _fminf(_fmaxf(a.y, mn), mx),
    _fminf(_fmaxf(a.z, mn), mx)
  );
}

__DEVICE__ float3 maxf3(float b, float3 a) {
  // For each component of float3 a, return max of component and float b
  return make_float3(_fmaxf(a.x, b), _fmaxf(a.y, b), _fmaxf(a.z, b));
}

__DEVICE__ float3 minf3(float b, float3 a) {
  // For each component of float3 a, return min of component and float b
  return make_float3(_fminf(a.x, b), _fminf(a.y, b), _fminf(a.z, b));
}

__DEVICE__ float _sign(float x) {
  // Return the sign of float x
  if (x > 0.0f) return 1.0f;
  if (x < 0.0f) return -1.0f;
  return 0.0f;
}

__DEVICE__ float3 powf3(float3 a, float b) {
  // Raise each component of float3 a to power b
  return make_float3(_powf(a.x, b), _powf(a.y, b), _powf(a.z, b));
}

__DEVICE__ float spowf(float a, float b) {
  // Compute "safe" power of float a, reflected over the origin

  a=_sign(a)*_powf(_fabs(a), b);
  return a;
}

__DEVICE__ float3 spowf3(float3 a, float b) {
  // Compute "safe" power of float3 a, reflected over the origin
  return make_float3(
    _sign(a.x)*_powf(_fabs(a.x), b),
    _sign(a.y)*_powf(_fabs(a.y), b),
    _sign(a.z)*_powf(_fabs(a.z), b)
  );
}

__DEVICE__ float _mixf(float a, float b, float f) {
  // Linear interpolation between float a and float b by factor f. Extrapolates.
  return a * (1.0f - f) + b * f;
}

__DEVICE__ float3 _mixf3(float3 a, float3 b, float f) {
  // Linear interpolation between float3 a and float3 b by factor f. Extrapolates.
  return make_float3(_mixf(a.x, b.x, f), _mixf(a.y, b.y, f), _mixf(a.z, b.z, f));
}

__DEVICE__ float3 _log2f3(float3 RGB) {
  return make_float3(_log2f(RGB.x), _log2f(RGB.y), _log2f(RGB.z));
}

__DEVICE__ float _smoothstepf(float e0, float e1, float x) {
  // return smoothstep of float x between e0 and e1
  x = _clampf((x - e0) / (e1 - e0), 0.0f, 1.0f);
  return x * x * (3.0f - 2.0f * x);
}

__DEVICE__ float3 _smoothstepf3(float e0, float e1, float3 x) {
  // return smoothstep of float3 x between e0 and e1
  return make_float3(_smoothstepf(e0, e1, x.x), _smoothstepf(e0, e1, x.y), _smoothstepf(e0, e1, x.z));
}

__DEVICE__ float chroma(float3 rgb, int norm) {
  // Calculate and return classical chroma. If norm, normalize by mx
  float mx = _fmaxf(rgb.x, _fmaxf(rgb.y, rgb.z));
  float mn = _fminf(rgb.x, _fminf(rgb.y, rgb.z));
  float ch = mx - mn;
  if (norm == 1) ch = mx == 0.0f ? 0.0f : ch / mx;
  return ch;
}

__DEVICE__ float hue(float3 rgb) {
  // Calculate and return hue in degrees between 0 and 6
  float mx = _fmaxf(rgb.x, _fmaxf(rgb.y, rgb.z));
  float mn = _fminf(rgb.x, _fminf(rgb.y, rgb.z));
  float ch = mx - mn;
  float h = 0.0;
  if (ch == 0.0f) h = 0.0f;
  else if (mx == rgb.x) h = _fmod((rgb.y - rgb.z) / ch + 6.0f, 6.0f);
  else if (mx == rgb.y) h = (rgb.z - rgb.x) / ch + 2.0f;
  else if (mx == rgb.z) h = (rgb.x - rgb.y) / ch + 4.0f;
  return h;
}

// Helper function to create a float3x3
__DEVICE__ float3x3 make_float3x3(float3 a, float3 b, float3 c) {
  float3x3 d;
  d.x = a, d.y = b, d.z = c;
  return d;
}

// Multiply float3 vector a and 3x3 matrix m
__DEVICE__ float3 mult_f3_f33(float3 a, float3x3 m) {
  return make_float3(
    m.x.x * a.x + m.x.y * a.y + m.x.z * a.z,
    m.y.x * a.x + m.y.y * a.y + m.y.z * a.z,
    m.z.x * a.x + m.z.y * a.y + m.z.z * a.z
  );
}

// Calculate inverse of 3x3 matrix: https://stackoverflow.com/questions/983999/simple-3x3-matrix-inverse-code-c
__DEVICE__ float3x3 inv_f33(float3x3 m) {
  float d = m.x.x * (m.y.y * m.z.z - m.z.y * m.y.z) -
            m.x.y * (m.y.x * m.z.z - m.y.z * m.z.x) +
            m.x.z * (m.y.x * m.z.y - m.y.y * m.z.x);
  float id = 1.0f / d;
  float3x3 c = identity_mtx;
  c.x.x = id * (m.y.y * m.z.z - m.z.y * m.y.z);
  c.x.y = id * (m.x.z * m.z.y - m.x.y * m.z.z);
  c.x.z = id * (m.x.y * m.y.z - m.x.z * m.y.y);
  c.y.x = id * (m.y.z * m.z.x - m.y.x * m.z.z);
  c.y.y = id * (m.x.x * m.z.z - m.x.z * m.z.x);
  c.y.z = id * (m.y.x * m.x.z - m.x.x * m.y.z);
  c.z.x = id * (m.y.x * m.z.y - m.z.x * m.y.y);
  c.z.y = id * (m.z.x * m.x.y - m.x.x * m.z.y);
  c.z.z = id * (m.x.x * m.y.y - m.y.x * m.x.y);
  return c;
}

__DEVICE__ float3x3 transpose_f33( float3x3 A) {
  float3x3 B = A;
  A.x=make_float3(B.x.x,B.y.x,B.z.x);
  A.y=make_float3(B.x.y,B.y.y,B.z.y);
  A.z=make_float3(B.x.z,B.y.z,B.z.z);

  return A;
}

__DEVICE__ float3x3 mult_f33_f33( float3x3 A, float3x3 B) {
  A = transpose_f33(A);
  float3x3 C = B;
  B.x= mult_f3_f33(A.x,C);
  B.y= mult_f3_f33(A.y,C);
  B.z= mult_f3_f33(A.z,C);
  B = transpose_f33(B);

  return B;
}

__DEVICE__ float3x3 RGBtoXYZ( Chromaticities N) {
  float3x3 M = make_float3x3(
    make_float3(N.red.x/N.red.y, N.green.x / N.green.y, N.blue.x / N.blue.y),
    make_float3(1.0, 1.0, 1.0),
    make_float3(
      (1-N.red.x-N.red.y) / N.red.y, (1-N.green.x-N.green.y) / N.green.y, (1-N.blue.x-N.blue.y)/N.blue.y
    )
  );
  float3 wh = make_float3(
    N.white.x / N.white.y, 1.0, (1-N.white.x-N.white.y) / N.white.y
  );
  wh = mult_f3_f33(wh, inv_f33(M));
  M = make_float3x3(
    make_float3(M.x.x*wh.x , M.x.y*wh.y , M.x.z*wh.z),
    make_float3(M.y.x*wh.x, M.y.y*wh.y, M.y.z*wh.z),
    make_float3(M.z.x*wh.x,M.z.y*wh.y,M.z.z*wh.z)
  );
  return M;
}

__DEVICE__ float3x3 XYZtoRGB( Chromaticities N) {
  float3x3 M = inv_f33(RGBtoXYZ(N));
  return M;
}

__DEVICE__ Chromaticities Insetcalc(Chromaticities N,float cpr,float cpg,float cpb){
  float3 scale = make_float3(
    1.0 / _powf(1.0 - cpr, 2.0),
    1.0 / _powf(1.0 - cpg, 2.0),
    1.0 / _powf(1.0 - cpb, 2.0)
  );

  Chromaticities adj = make_chromaticities(
    make_float2(
      (N.red.x - N.white.x) * scale.x + N.white.x,
      (N.red.y-N.white.y) * scale.x + N.white.y
    ),
    make_float2(
      (N.green.x - N.white.x) * scale.y + N.white.x,
      (N.green.y - N.white.y) * scale.y + N.white.y
    ),
    make_float2(
      (N.blue.x - N.white.x) * scale.z + N.white.x,
      (N.blue.y - N.white.y) * scale.z + N.white.y
    ),
    make_float2(N.white.x, N.white.y)
  );

  return adj;
}

__DEVICE__ float3x3 Insetcalcmatrix(Chromaticities N,float cpr,float cpg,float cpb){
  float3 scale = make_float3(1/_powf(1-cpr,2),1/_powf(1-cpg,2),1/_powf(1-cpb,2));

  Chromaticities adj = make_chromaticities(make_float2((N.red.x-N.white.x)*scale.x+N.white.x,(N.red.y-N.white.y)*scale.x+N.white.y),
      make_float2((N.green.x-N.white.x)*scale.y+N.white.x,(N.green.y-N.white.y)*scale.y+N.white.y),
      make_float2((N.blue.x-N.white.x)*scale.z+N.white.x,(N.blue.y-N.white.y)*scale.z+N.white.y),make_float2(N.white.x,N.white.y));

  float3x3 In2XYZ = RGBtoXYZ(N);
  float3x3 XYZ2Adj = XYZtoRGB(adj);

  float3x3 RGBtoAdj = mult_f33_f33(In2XYZ,XYZ2Adj);

  return RGBtoAdj;
}

__DEVICE__ Chromaticities Primaries2Moment(Chromaticities N){
  float2 momr = make_float2((N.red.x-N.white.x)/N.red.y,(N.red.y-N.white.y)/N.red.y);
  float2 momg = make_float2((N.green.x-N.white.x)/N.green.y,(N.green.y-N.white.y)/N.green.y);
  float2 momb = make_float2((N.blue.x-N.white.x)/N.blue.y,(N.blue.y-N.white.y)/N.blue.y);

  Chromaticities M = make_chromaticities(momr,momg,momb,N.white);

  return M;
}

__DEVICE__ Chromaticities CenterPrimaries(Chromaticities N){
  N.red.x = N.red.x-N.white.x;
  N.red.y = N.red.y-N.white.y;
  N.green.x = N.green.x-N.white.x;
  N.green.y = N.green.y-N.white.y;
  N.blue.x = N.blue.x-N.white.x;
  N.blue.y = N.blue.y-N.white.y;

  return N;
}

__DEVICE__ Chromaticities DeCenterPrimaries(Chromaticities N){
  N.red.x = N.red.x+N.white.x;
  N.red.y = N.red.y+N.white.y;
  N.green.x = N.green.x+N.white.x;
  N.green.y = N.green.y+N.white.y;
  N.blue.x = N.blue.x+N.white.x;
  N.blue.y = N.blue.y+N.white.y;

  return N;
}

__DEVICE__ Chromaticities ScalePrim(Chromaticities N,float rs,float gs,float bs){
  N = CenterPrimaries(N);

  N.red = make_float2(N.red.x*rs,N.red.y*rs);
  N.green = make_float2(N.green.x*gs,N.green.y*gs);
  N.blue = make_float2(N.blue.x*bs,N.blue.y*bs);

  N = DeCenterPrimaries(N);

  return N;
}

__DEVICE__ float2 cartesian_to_polar2(float2 a) {
  float2 b = a;
  b.y = _atan2f(a.y,a.x);

  return make_float2(_sqrtf(a.x*a.x+ a.y*a.y),b.y);
}

__DEVICE__ float2 polar_to_cartesian2(float2 a) {

  return make_float2(a.x * _cosf(a.y), a.x * _sinf(a.y));
}

__DEVICE__ Chromaticities RotatePrimary(Chromaticities N,float rrot,float grot,float brot){
  //rotatation parameter excepted in degrees, but internally transformed to radians

  N = CenterPrimaries(N);
  N.red = cartesian_to_polar2(N.red);
  N.green = cartesian_to_polar2(N.green);
  N.blue = cartesian_to_polar2(N.blue);

  rrot = _radians(rrot);
  grot = _radians(grot);
  brot = _radians(brot);

  N.red.y = N.red.y + rrot;
  N.green.y = N.green.y + grot;
  N.blue.y = N.blue.y + brot;

  N.red = polar_to_cartesian2(N.red);
  N.green = polar_to_cartesian2(N.green);
  N.blue = polar_to_cartesian2(N.blue);

  N = DeCenterPrimaries(N);

  return N;
}

__DEVICE__ float2 Line_equation(float2 a, float2 b)
{
    float dx = b.x - a.x;
    if (fabsf(dx) < 1e-6f) {
        /* vertical line → slope = +∞, intercept = x */
        const float kInf = 1e30f * 1e30f;
        return make_float2(kInf, a.x);
    }
    float m = (b.y - a.y) / dx;
    float c = a.y - m * a.x;
    return make_float2(m, c);
}

__DEVICE__ Chromaticities PrimariesLines(Chromaticities N){
    Chromaticities M = N;

    N.red = Line_equation(M.red,M.white);
    N.green = Line_equation(M.green,M.white);
    N.blue = Line_equation(M.blue,M.white);

    return N;
}

__DEVICE__ Chromaticities Polygon(Chromaticities N){
  Chromaticities M = N;

  N.red = Line_equation(M.red,M.green);
  N.green = Line_equation(M.red,M.blue);
  N.blue = Line_equation(M.blue,M.green);

  return N;
}

__DEVICE__ float2 intersection(float2 l1, float2 l2) {
  float m1 = l1.x, c1 = l1.y;
  float m2 = l2.x, c2 = l2.y;

  bool inf1 = isfinite(m1) == false;
  bool inf2 = isfinite(m2) == false;

  float x, y;
  if (inf1 && inf2) {
    // two verticals (parallel) → just default (0,0)
    x = y = 0.0f;
  }
  else if (inf1) {
    // first is vertical: x = c1
    x = c1;
    y = m2 * x + c2;
  }
  else if (inf2) {
    // second is vertical: x = c2
    x = c2;
    y = m1 * x + c1;
  }
  else {
    x = (c2 - c1) / (m1 - m2);
    y = m1 * x + c1;
  }
  return make_float2(x, y);
}

__DEVICE__ Chromaticities InsetPrimaries(Chromaticities N, float cpr, float cpg, float cpb, float ored, float og, float ob, float achromatic_rotate = 0, float achromatic_outset = 0) {
  Chromaticities M = N;

  Chromaticities original_N = N;
  Chromaticities scaled_N = ScalePrim(N, 4.0, 4.0, 4.0);

  N = RotatePrimary(scaled_N, ored, og, ob);

  M = Polygon(M);

  float2 redline = ored > 0 ? M.red : M.green;
  float2 greenline = og > 0 ? M.blue : M.red;
  float2 blueline = ob > 0 ? M.green : M.blue;

  // compute the line eqns from each rotated‐primary → whitepoint
  float2 wp = original_N.white;
  float2 lr = Line_equation( make_float2(N.red.x,   N.red.y),   wp );
  float2 lg = Line_equation( make_float2(N.green.x, N.green.y), wp );
  float2 lb = Line_equation( make_float2(N.blue.x,  N.blue.y),  wp );

  // intersect each with the chosen triangle edge
  float2 Pr = intersection(lr, redline);
  float2 Pg = intersection(lg, greenline);
  float2 Pb = intersection(lb, blueline);

  // assign back into N
  N.red   = Pr;
  N.green = Pg;
  N.blue  = Pb;

  cpr = 1 - cpr;
  cpg = 1 - cpg;
  cpb = 1 - cpb;

  N = ScalePrim(N, cpr, cpg, cpb);

  // --- Start of Achromatic Tinting ---
  float2 original_white = Polygon(original_N).white;
  const float arbitrary_scale = 4.0f;

  // Scale & rotate the achromatic point
  float2 scaled_achromatic = make_float2(
      original_white.x,
      original_white.y * arbitrary_scale
  );
  float achromatic_rotate_radians = _radians(achromatic_rotate);
  float dx = scaled_achromatic.x - original_white.x;
  float dy = scaled_achromatic.y - original_white.y;
  float2 rotated_achromatic = make_float2(
      original_white.x + dx * cos(achromatic_rotate_radians) - dy * sin(achromatic_rotate_radians),
      original_white.y + dx * sin(achromatic_rotate_radians) + dy * cos(achromatic_rotate_radians)
  );

  // Build the infinite achromatic ray and triangle edges
  float2 la = Line_equation(rotated_achromatic, original_white);
  float2 e1 = Line_equation(original_N.red,   original_N.green);
  float2 e2 = Line_equation(original_N.green, original_N.blue);
  float2 e3 = Line_equation(original_N.blue,  original_N.red);

  // Find their intersections
  float2 i1 = intersection(la, e1);
  float2 i2 = intersection(la, e2);
  float2 i3 = intersection(la, e3);

  // Pick the first intersect that lies on both the achromatic ray segment (t∈[0,1])
  // and the edge segment (u∈[0,1]), else fallback to white
  float2 hull_achromatic = original_white;
  float2 A = rotated_achromatic, B = original_white;
  auto on_segment = [&](float2 P, float2 C, float2 D){
    float t = (fabs(B.x-A.x)>fabs(B.y-A.y))
      ? (P.x-A.x)/(B.x-A.x)
      : (P.y-A.y)/(B.y-A.y);
    float u = (fabs(D.x-C.x)>fabs(D.y-C.y))
      ? (P.x-C.x)/(D.x-C.x)
      : (P.y-C.y)/(D.y-C.y);
    return (t >= 0.0f && t <= 1.0f && u >= 0.0f && u <= 1.0f);
  };
  if      (on_segment(i1, original_N.red,   original_N.green)) hull_achromatic = i1;
  else if (on_segment(i2, original_N.green, original_N.blue )) hull_achromatic = i2;
  else if (on_segment(i3, original_N.blue,  original_N.red  )) hull_achromatic = i3;

  // Move whitepoint towards hull_achromatic by achromatic_outset
  float2 interp = make_float2(
      (original_white.x - hull_achromatic.x) * (1.0f - achromatic_outset),
      (original_white.y - hull_achromatic.y) * (1.0f - achromatic_outset)
  );
  N.white.x = hull_achromatic.x + interp.x;
  N.white.y = hull_achromatic.y + interp.y;

  return N;
}


__DEVICE__ float3x3 RGBtoRGB(Chromaticities N,Chromaticities M){
  float3x3 In2XYZ = RGBtoXYZ(N);
  float3x3 XYZ2Out = XYZtoRGB(M);

  float3x3 rgbtorgb = mult_f33_f33(In2XYZ,XYZ2Out);

  return rgbtorgb;
}

__DEVICE__ float3 XYZ_2_xyY( float3 XYZ) {
  float3 xyY;
  //float divisor = (_fabs(XYZ.x) + _fabs(XYZ.y) + _fabs(XYZ.z));
  float divisor = ((XYZ.x) + XYZ.y + XYZ.z);
  //if (divisor == 0.0f) divisor = 1e-10f;
  xyY.x =divisor == 0.0f? 0.0f:(XYZ.x) / divisor;
  xyY.y = divisor == 0.0f? 0.0f:(XYZ.y)/ divisor;
  //xyY.z = _fabs(XYZ.y);
  xyY.z = XYZ.y;
  return xyY;
}

__DEVICE__ float3 xyY_2_XYZ( float3 xyY) {
  float3 XYZ;
  //XYZ.x = xyY.x * xyY.z / _fmaxf( xyY.y, 1e-10f);
  XYZ.x = xyY.y==0.0f? 0.0f: xyY.x * xyY.z / xyY.y;
  XYZ.y = xyY.z;
  //XYZ.z = (1.0f - xyY.x - xyY.y) * xyY.z / _fmaxf( xyY.y, 1e-10f);
  XYZ.z = xyY.y==0.0f? 0.0f: (1.0f - xyY.x - xyY.y) * xyY.z /( xyY.y);
  return XYZ;
}

// Blender-AgX's Compensate Low side
__DEVICE__ float3 compensate_low_side(float3 rgb, int use_heuristic_lerp, Chromaticities working_chrom) {
    // Hardcoded Rec.2020 luminance coefficients (2015 CMFs)
    const float3 luminance_coeffs = make_float3(0.2589235355689848f, 0.6104985346066525f, 0.13057792982436284f);
    Chromaticities rec2020 = REC2020_PRI;

    // Convert RGB to Rec.2020 for luminance calculation
    float3x3 working_to_rec2020 = RGBtoRGB(working_chrom, rec2020);
    float3 rgb_rec2020 = mult_f3_f33(rgb, working_to_rec2020);

    // Calculate original luminance Y
    float Y = rgb_rec2020.x * luminance_coeffs.x +
              rgb_rec2020.y * luminance_coeffs.y +
              rgb_rec2020.z * luminance_coeffs.z;

    // Calculate inverse RGB in working space
    float max_rgb = fmaxf(rgb.x, fmaxf(rgb.y, rgb.z));
    float3 inverse_rgb = make_float3(max_rgb - rgb.x, max_rgb - rgb.y, max_rgb - rgb.z);

    // Calculate max of the inverse
    float max_inv_rgb = fmaxf(inverse_rgb.x, fmaxf(inverse_rgb.y, inverse_rgb.z));

    // Convert inverse RGB to Rec.2020 for Y calculation
    float3 inverse_rec2020 = mult_f3_f33(inverse_rgb, working_to_rec2020);
    float Y_inverse = inverse_rec2020.x * luminance_coeffs.x +
                      inverse_rec2020.y * luminance_coeffs.y +
                      inverse_rec2020.z * luminance_coeffs.z;

    // Calculate compensation values
    float y_compensate = (max_inv_rgb - Y_inverse + Y);
    if (use_heuristic_lerp) {
        float Y_clipped = _clampf(powf(Y, 0.08f), 0.0f, 1.0f);
        y_compensate = y_compensate + Y_clipped * (Y - y_compensate);
    }

    // Offset to avoid negatives
    float min_rgb = fminf(rgb.x, fminf(rgb.y, rgb.z));
    float offset = fmaxf(-min_rgb, 0.0f);
    float3 rgb_offset = make_float3(rgb.x + offset, rgb.y + offset, rgb.z + offset);

    // Calculate max of the offseted RGB
    float max_offset = fmaxf(rgb_offset.x, fmaxf(rgb_offset.y, rgb_offset.z));

    // Calculate new luminance after offset
    float3 offset_rec2020 = mult_f3_f33(rgb_offset, working_to_rec2020);
    float Y_new = offset_rec2020.x * luminance_coeffs.x +
                  offset_rec2020.y * luminance_coeffs.y +
                  offset_rec2020.z * luminance_coeffs.z;

    // Calculate the inverted RGB offset
    float3 inverse_offset = make_float3(max_offset - rgb_offset.x,
                                       max_offset - rgb_offset.y,
                                       max_offset - rgb_offset.z);

    // Calculate max of the inverse
    float max_inv_offset = fmaxf(inverse_offset.x, fmaxf(inverse_offset.y, inverse_offset.z));

    float3 inverse_offset_rec2020 = mult_f3_f33(inverse_offset, working_to_rec2020);
    float Y_inverse_offset = inverse_offset_rec2020.x * luminance_coeffs.x +
                             inverse_offset_rec2020.y * luminance_coeffs.y +
                             inverse_offset_rec2020.z * luminance_coeffs.z;

    float Y_new_compensate = (max_inv_offset - Y_inverse_offset + Y_new);
    if (use_heuristic_lerp) {
        float Y_new_clipped = _clampf(powf(Y_new, 0.08f), 0.0f, 1.0f);
        Y_new_compensate = Y_new_compensate + Y_new_clipped * (Y_new - Y_new_compensate);
    }

    // Adjust luminance ratio
    float ratio = (Y_new_compensate > y_compensate) ? (y_compensate / Y_new_compensate) : 1.0f;
    return make_float3(rgb_offset.x * ratio, rgb_offset.y * ratio, rgb_offset.z * ratio);
}

// Blender-AgX's Circular Chromaticity Angle Mixing
__DEVICE__ float3 RGB_to_HSV(float3 rgb) {
    float cmax = fmaxf(rgb.x, fmaxf(rgb.y, rgb.z));
    float cmin = fminf(rgb.x, fminf(rgb.y, rgb.z));
    float delta = cmax - cmin;

    float h = 0.0f;
    if (delta != 0.0f) {
        if (cmax == rgb.x) h = fmodf((rgb.y - rgb.z)/delta + 6.0f, 6.0f);
        else if (cmax == rgb.y) h = (rgb.z - rgb.x)/delta + 2.0f;
        else h = (rgb.x - rgb.y)/delta + 4.0f;
    }
    h /= 6.0f; // Normalize to 0-1

    float s = cmax == 0.0f ? 0.0f : delta/cmax;
    return make_float3(h, s, cmax);
}

__DEVICE__ float3 HSV_to_RGB(float3 hsv) {
    float h = hsv.x * 6.0f;
    float c = hsv.z * hsv.y;
    float x = c * (1.0f - fabsf(fmodf(h, 2.0f) - 1.0f));

    float3 rgb;
    if (h < 1.0f) rgb = make_float3(c, x, 0.0f);
    else if (h < 2.0f) rgb = make_float3(x, c, 0.0f);
    else if (h < 3.0f) rgb = make_float3(0.0f, c, x);
    else if (h < 4.0f) rgb = make_float3(0.0f, x, c);
    else if (h < 5.0f) rgb = make_float3(x, 0.0f, c);
    else rgb = make_float3(c, 0.0f, x);

    float m = hsv.z - c;
    return make_float3(rgb.x + m, rgb.y + m, rgb.z + m);
}

__DEVICE__ float lerp_chromaticity_angle(float h1, float h2, float t) {
    float delta = h2 - h1;
    if (delta > 0.5f) delta -= 1.0f;
    else if (delta < -0.5f) delta += 1.0f;
    float lerped = h1 + t * delta;
    return lerped - floorf(lerped);
}

/* ##########################################################################
    Transfer Functions
    ---------------------------------
*/

__DEVICE__ float3 encode_inverse_EOTF(float3 rgb, float EOTF) {
  rgb.x = _powf(rgb.x, 1.0f / EOTF);
  rgb.y = _powf(rgb.y, 1.0f / EOTF);
  rgb.z = _powf(rgb.z, 1.0f / EOTF);

  return rgb;
}

__DEVICE__ float3 log2lin(float3 rgb, int tf, float generic_log2_min_expo = -10, float generic_log2_max_expo = 6.5) {
  if (tf == 0) return rgb;
  else if (tf == 1) { // ACEScct
    rgb.x = rgb.x > 0.155251141552511f ? _powf(2.0f, rgb.x * 17.52f - 9.72f) : (rgb.x - 0.0729055341958355f) / 10.5402377416545f;
    rgb.y = rgb.y > 0.155251141552511f ? _powf(2.0f, rgb.y * 17.52f - 9.72f) : (rgb.y - 0.0729055341958355f) / 10.5402377416545f;
    rgb.z = rgb.z > 0.155251141552511f ? _powf(2.0f, rgb.z * 17.52f - 9.72f) : (rgb.z - 0.0729055341958355f) / 10.5402377416545f;
  } else if (tf == 2) { // Arri V3 LogC EI 800
    rgb.x = rgb.x > 0.149658f ? (_powf(10.0f, (rgb.x - 0.385537f) / 0.24719f) - 0.052272f) / 5.555556f : (rgb.x - 0.092809f) / 5.367655f;
    rgb.y = rgb.y > 0.149658f ? (_powf(10.0f, (rgb.y - 0.385537f) / 0.24719f) - 0.052272f) / 5.555556f : (rgb.y - 0.092809f) / 5.367655f;
    rgb.z = rgb.z > 0.149658f ? (_powf(10.0f, (rgb.z - 0.385537f) / 0.24719f) - 0.052272f) / 5.555556f : (rgb.z - 0.092809f) / 5.367655f;
  } else if (tf == 3) { // Red Log3G10
    rgb.x = rgb.x > 0.0f ? (_powf(10.0f, rgb.x / 0.224282f) - 1.0f) / 155.975327f - 0.01f : (rgb.x / 15.1927f) - 0.01f;
    rgb.y = rgb.y > 0.0f ? (_powf(10.0f, rgb.y / 0.224282f) - 1.0f) / 155.975327f - 0.01f : (rgb.y / 15.1927f) - 0.01f;
    rgb.z = rgb.z > 0.0f ? (_powf(10.0f, rgb.z / 0.224282f) - 1.0f) / 155.975327f - 0.01f : (rgb.z / 15.1927f) - 0.01f;
  } else if (tf == 4) { // Sony SLog3
    rgb.x = rgb.x >= 171.2102946929f / 1023.0f ? _powf(10.0f, ((rgb.x * 1023.0f - 420.0f) / 261.5f)) * (0.18f + 0.01f) - 0.01f : (rgb.x * 1023.0f - 95.0f) * 0.01125000f / (171.2102946929f - 95.0f);
    rgb.y = rgb.y >= 171.2102946929f / 1023.0f ? _powf(10.0f, ((rgb.y * 1023.0f - 420.0f) / 261.5f)) * (0.18f + 0.01f) - 0.01f : (rgb.y * 1023.0f - 95.0f) * 0.01125000f / (171.2102946929f - 95.0f);
    rgb.z = rgb.z >= 171.2102946929f / 1023.0f ? _powf(10.0f, ((rgb.z * 1023.0f - 420.0f) / 261.5f)) * (0.18f + 0.01f) - 0.01f : (rgb.z * 1023.0f - 95.0f) * 0.01125000f / (171.2102946929f - 95.0f);
  } else if (tf == 5) { // Filmlight T-Log
    rgb.x = rgb.x < 0.075f ? (rgb.x - 0.075f) / 16.18437649f : _expf((rgb.x - 0.55201266f) / 0.09232903f) - 0.00570482f;
    rgb.y = rgb.y < 0.075f ? (rgb.y - 0.075f) / 16.18437649f : _expf((rgb.y - 0.55201266f) / 0.09232903f) - 0.00570482f;
    rgb.z = rgb.z < 0.075f ? (rgb.z - 0.075f) / 16.18437649f : _expf((rgb.z - 0.55201266f) / 0.09232903f) - 0.00570482f;
  } else if (tf == 6) { // DaVinci Intermediate
    rgb.x = rgb.x <= 0.02740668f ? rgb.x / 10.44426855f : _powf(2.0f, (rgb.x / 0.07329248f) - 7.0f) - 0.0075f;
    rgb.y = rgb.y <= 0.02740668f ? rgb.y / 10.44426855f : _powf(2.0f, (rgb.y / 0.07329248f) - 7.0f) - 0.0075f;
    rgb.z = rgb.z <= 0.02740668f ? rgb.z / 10.44426855f : _powf(2.0f, (rgb.z / 0.07329248f) - 7.0f) - 0.0075f;
  } else if (tf == 7) { // Blackmagic Film Gen5
    rgb.x = rgb.x < 0.13388378f ? (rgb.x - 0.09246575342465753f) / 8.283605932402494f : _expf((rgb.x - 0.5300133392291939f) / 0.08692876065491224f) - 0.005494072432257808f;
    rgb.y = rgb.y < 0.13388378f ? (rgb.y - 0.09246575342465753f) / 8.283605932402494f : _expf((rgb.y - 0.5300133392291939f) / 0.08692876065491224f) - 0.005494072432257808f;
    rgb.z = rgb.z < 0.13388378f ? (rgb.z - 0.09246575342465753f) / 8.283605932402494f : _expf((rgb.z - 0.5300133392291939f) / 0.08692876065491224f) - 0.005494072432257808f;
  } else if (tf == 8) { // CanonLog3
    rgb.x = rgb.x < 0.097465473f ?-( _powf( 10.0f, ( 0.12783901f - rgb.x ) / 0.36726845f ) - 1.0f ) / 14.98325f : 0.097465473f <= rgb.x && rgb.x <= 0.15277891f?(rgb.x - 0.12512219f) / 1.9754798f:( _powf( 10.0f, ( rgb.x - 0.12240537f ) / 0.36726845f ) - 1.0f ) / 14.98325f;
    rgb.y = rgb.y < 0.097465473f ?-( _powf( 10.0f, ( 0.12783901f - rgb.y ) / 0.36726845f ) - 1.0f ) / 14.98325f : 0.097465473f <= rgb.y && rgb.y <= 0.15277891f?(rgb.y - 0.12512219f) / 1.9754798f:( _powf( 10.0f, ( rgb.y - 0.12240537f ) / 0.36726845f ) - 1.0f ) / 14.98325f;
    rgb.z = rgb.z < 0.097465473f ?-( _powf( 10.0f, ( 0.12783901f - rgb.z ) / 0.36726845f ) - 1.0f ) / 14.98325f : 0.097465473f <= rgb.z && rgb.z <= 0.15277891f?(rgb.z - 0.12512219f) / 1.9754798f:( _powf( 10.0f, ( rgb.z - 0.12240537f ) / 0.36726845f ) - 1.0f ) / 14.98325f;
    rgb.x = rgb.x*0.9f;
    rgb.y = rgb.y*0.9f;
    rgb.z = rgb.z*0.9f;
  } else if (tf == 9){  // Arri LogC 4
    const float a = (_powf(2.0f, 18.0f) - 16.0f) / 117.45f;
    const float b = (1023.0f - 95.0f) / 1023.0f;
    const float c = 95.0f / 1023.f;
    const float s = (7.f * _logf(2.0f) * _powf(2.0f, 7.0f - 14.0f * c / b)) / (a * b);
    const float t = (_powf(2.0f, 14.0f * ((-1.0f * c) / b) + 6.0f) - 64.0f) / a;

    rgb.x = rgb.x < 0.0f ? rgb.x * s + t : (_powf(2.0f, (14.0f * (rgb.x - c) / b + 6.0f)) - 64.0f) / a;
    rgb.y = rgb.y < 0.0f ? rgb.y * s + t : (_powf(2.0f, (14.0f * (rgb.y - c) / b + 6.0f)) - 64.0f) / a;
    rgb.z = rgb.z < 0.0f ? rgb.z * s + t : (_powf(2.0f, (14.0f * (rgb.z - c) / b + 6.0f)) - 64.0f) / a;
  } else if (tf == 10){ //PureLog2 -10 stop under 0.18 and 6.5 over
    float mx = 6.5;
    float mn = -10;

    rgb.x = 0.18*_powf(2,(rgb.x*(mx-mn)+mn));
    rgb.y = 0.18*_powf(2,(rgb.y*(mx-mn)+mn));
    rgb.z = 0.18*_powf(2,(rgb.z*(mx-mn)+mn));
  } else if (tf == 11){  // Flog2
    const float a = 5.555556f;
    const float b = 0.064829f;
    const float c = 0.245281f;
    const float d = 0.384316f;
    const float e = 8.799461f;
    const float f = 0.092864f;
    //const float cut1 = 0.000889f; // Should be equal to ((cut2 - f) / e)
    const float cut2 = 0.100686685370811f; // should be equal to (e * cut1 + f)

    rgb.x = rgb.x>=cut2?((_exp10f((rgb.x - d) / c) - b) / a):((rgb.x - f) / e);
    rgb.y = rgb.y>=cut2?((_exp10f((rgb.y - d) / c) - b) / a):((rgb.y - f) / e);
    rgb.z = rgb.z>=cut2?((_exp10f((rgb.z - d) / c) - b) / a):((rgb.z - f) / e);
  } else if (tf == 12){  // BMFilm
    const float p1 = 0.06031746;
    const float p2 = 0.00712130;
    const float p3 = 0.20123443;
    const float p4 = 0.72309772;
    const float p5 = 0.42929868;
    const float p6 = 0.02476426;
    const float p7 = 0.76154283;


    rgb.x = rgb.x<=p1?rgb.x*p3-p2:(spowf(10,(rgb.x-p4)/p5)-p6)/p7;
    rgb.y = rgb.y<=p1?rgb.y*p3-p2:(spowf(10,(rgb.y-p4)/p5)-p6)/p7;
    rgb.z = rgb.z<=p1?rgb.z*p3-p2:(spowf(10,(rgb.z-p4)/p5)-p6)/p7;
  } else if (tf == 13){ //PureLog2 -10 stop under 0.18 and 15 over (AgX Log Blender)
    float mx = 15;
    float mn = -10;

    rgb.x = 0.18*_powf(2,(rgb.x*(mx-mn)+mn));
    rgb.y = 0.18*_powf(2,(rgb.y*(mx-mn)+mn));
    rgb.z = 0.18*_powf(2,(rgb.z*(mx-mn)+mn));
  } else if (tf == 14){ // User controlled PureLog2
    float mx = generic_log2_max_expo;
    float mn = generic_log2_min_expo;

    rgb.x = 0.18*_powf(2,(rgb.x*(mx-mn)+mn));
    rgb.y = 0.18*_powf(2,(rgb.y*(mx-mn)+mn));
    rgb.z = 0.18*_powf(2,(rgb.z*(mx-mn)+mn));
  }
  return rgb;
}

__DEVICE__ float3 lin2log(float3 rgb, int tf, float generic_log2_min_expo = -10, float generic_log2_max_expo = 6.5) {
  float log_floor = log2lin(make_float3(0.0f), tf, generic_log2_min_expo, generic_log2_max_expo).x;
  rgb = maxf3(log_floor, rgb);
  if (tf == 0) return rgb;
  else if (tf == 1) { // ACEScct
    rgb.x = rgb.x > 0.0078125f ? (_log2f(rgb.x) + 9.72f) / 17.52f : 10.5402377416545f * rgb.x + 0.0729055341958355f;
    rgb.y = rgb.y > 0.0078125f ? (_log2f(rgb.y) + 9.72f) / 17.52f : 10.5402377416545f * rgb.y + 0.0729055341958355f;
    rgb.z = rgb.z > 0.0078125f ? (_log2f(rgb.z) + 9.72f) / 17.52f : 10.5402377416545f * rgb.z + 0.0729055341958355f;
  } else if (tf == 2) { // Arri V3 LogC EI 800
    rgb.x = rgb.x > 0.010591f ? 0.24719f * _log10f(5.555556f * rgb.x + 0.052272f) + 0.385537f : 5.367655f * rgb.x + 0.092809f;
    rgb.y = rgb.y > 0.010591f ? 0.24719f * _log10f(5.555556f * rgb.y + 0.052272f) + 0.385537f : 5.367655f * rgb.y + 0.092809f;
    rgb.z = rgb.z > 0.010591f ? 0.24719f * _log10f(5.555556f * rgb.z + 0.052272f) + 0.385537f : 5.367655f * rgb.z + 0.092809f;
  } else if (tf == 3) { // Red Log3G10
    rgb.x = rgb.x > -0.01f ? 0.224282f * _log10f(((rgb.x + 0.01f) * 155.975327f) + 1.0f) : (rgb.x + 0.01f) * 15.1927f;
    rgb.y = rgb.y > -0.01f ? 0.224282f * _log10f(((rgb.y + 0.01f) * 155.975327f) + 1.0f) : (rgb.y + 0.01f) * 15.1927f;
    rgb.z = rgb.z > -0.01f ? 0.224282f * _log10f(((rgb.z + 0.01f) * 155.975327f) + 1.0f) : (rgb.z + 0.01f) * 15.1927f;
  } else if (tf == 4) { // Sony SLog3
    rgb.x = rgb.x >= 0.01125f ? (420.0f + _log10f((rgb.x + 0.01f) / (0.18f + 0.01f)) * 261.5f) / 1023.0f : (rgb.x * (171.2102946929f - 95.0f) / 0.01125000f + 95.0f) / 1023.0f;
    rgb.y = rgb.y >= 0.01125f ? (420.0f + _log10f((rgb.y + 0.01f) / (0.18f + 0.01f)) * 261.5f) / 1023.0f : (rgb.y * (171.2102946929f - 95.0f) / 0.01125000f + 95.0f) / 1023.0f;
    rgb.z = rgb.z >= 0.01125f ? (420.0f + _log10f((rgb.z + 0.01f) / (0.18f + 0.01f)) * 261.5f) / 1023.0f : (rgb.z * (171.2102946929f - 95.0f) / 0.01125000f + 95.0f) / 1023.0f;
  } else if (tf == 5) { // Filmlight T-Log
    rgb.x = rgb.x < 0.0f ? 16.18437649f * rgb.x + 0.075f : _logf(rgb.x + 0.00570482f) * 0.09232903f + 0.55201266f;
    rgb.y = rgb.y < 0.0f ? 16.18437649f * rgb.y + 0.075f : _logf(rgb.y + 0.00570482f) * 0.09232903f + 0.55201266f;
    rgb.z = rgb.z < 0.0f ? 16.18437649f * rgb.z + 0.075f : _logf(rgb.z + 0.00570482f) * 0.09232903f + 0.55201266f;
  } else if (tf == 6) { // DaVinci Intermediate
    rgb.x = rgb.x <= 0.00262409f ? rgb.x * 10.44426855f : (_log2f(rgb.x + 0.0075f) + 7.0f) * 0.07329248f;
    rgb.y = rgb.y <= 0.00262409f ? rgb.y * 10.44426855f : (_log2f(rgb.y + 0.0075f) + 7.0f) * 0.07329248f;
    rgb.z = rgb.z <= 0.00262409f ? rgb.z * 10.44426855f : (_log2f(rgb.z + 0.0075f) + 7.0f) * 0.07329248f;
  } else if (tf == 7) { // Blackmagic Film Gen5
    rgb.x = rgb.x <= 0.004999993693740552f ? rgb.x * 8.283611088773256f + 0.09246580021201303f : 0.08692875224330131f * _logf(rgb.x + 0.0054940711907293955f) + 0.5300133837514731f;
    rgb.y = rgb.y <= 0.004999993693740552f ? rgb.y * 8.283611088773256f + 0.09246580021201303f : 0.08692876065491224f * _logf(rgb.y + 0.0054940711907293955f) + 0.5300133837514731f;
    rgb.z = rgb.z <= 0.004999993693740552f ? rgb.z * 8.283611088773256f + 0.09246580021201303f : 0.08692876065491224f * _logf(rgb.z + 0.0054940711907293955f) + 0.5300133837514731f;
  } else if (tf == 8) { // CanonLog3
    rgb.x = rgb.x/0.9f ;
    rgb.y = rgb.y/0.9f ;
    rgb.z = rgb.z/0.9f ;
    rgb.x = rgb.x < -0.014f ? -0.36726845f * _log10f( 1.0f - 14.98325f * rgb.x ) + 0.12783901f: -0.014f <= rgb.x && rgb.x <= 0.014f?1.9754798f * rgb.x + 0.12512219f:0.36726845f * _log10f( 14.98325f * rgb.x + 1.0f ) + 0.12240537f;
    rgb.y = rgb.y < -0.014f ? -0.36726845f * _log10f( 1.0f - 14.98325f * rgb.y ) + 0.12783901f: -0.014f <= rgb.y && rgb.y <= 0.014f?1.9754798f * rgb.y + 0.12512219f:0.36726845f * _log10f( 14.98325f * rgb.y + 1.0f ) + 0.12240537f;
    rgb.z = rgb.z < -0.014f ? -0.36726845f * _log10f( 1.0f - 14.98325f * rgb.z ) + 0.12783901f: -0.014f <= rgb.z && rgb.z <= 0.014f?1.9754798f * rgb.z + 0.12512219f:0.36726845f * _log10f( 14.98325f * rgb.z + 1.0f ) + 0.12240537f;
  } else if (tf == 9){  // Arri LogC 4
    const float a = (_powf(2.0f, 18.0f) - 16.0f) / 117.45f;
    const float b = (1023.0f - 95.0f) / 1023.0f;
    const float c = 95.0f / 1023.f;
    const float s = (7.f * _logf(2.0f) * _powf(2.0f, 7.0f - 14.0f * c / b)) / (a * b);
    const float t = (_powf(2.0f, 14.0f * ((-1.0f * c) / b) + 6.0f) - 64.0f) / a;

    rgb.x = rgb.x >= t ? ((_log2f(a * rgb.x + 64.f) - 6.f) / 14.f) * b + c : (rgb.x - t) / s;
    rgb.y = rgb.y >= t ? ((_log2f(a * rgb.y + 64.f) - 6.f) / 14.f) * b + c : (rgb.y - t) / s;
    rgb.z = rgb.z >= t ? ((_log2f(a * rgb.z + 64.f) - 6.f) / 14.f) * b + c : (rgb.z - t) / s;
  } else if (tf == 10) { // Was CanonLog2, now Normalized Log2
    // total_exposure = maximum_ev - minimum_ev

    //   in_od = numpy.asarray(in_od)
    //   in_od[in_od <= 0.0] = numpy.finfo(float).eps

    //   output_log = numpy.clip(
    //       numpy.log2(in_od / in_middle_grey),
    //       minimum_ev,
    //       maximum_ev
    //   )
    // return as_numeric((output_log - minimum_ev) / total_exposure)

    rgb = _log2f3(rgb / 0.18f);
    rgb = clampf3(rgb, -10.0f, 6.5f);

    rgb = (rgb + 10.0f) / 16.5f;
  } else if (tf == 11){  // Flog2
    const float a = 5.555556f;
    const float b = 0.064829f;
    const float c = 0.245281f;
    const float d = 0.384316f;
    const float e = 8.799461f;
    const float f = 0.092864f;
    const float cut1 = 0.000889f; // Should be equal to ((cut2 - f) / e)
    //const float cut2 = 0.100686685370811f; // should be equal to (e * cut1 + f)

    rgb.x = rgb.x>=cut1?(c * _log10f(a * rgb.x + b) + d):(e * rgb.x + f);
    rgb.y = rgb.y>=cut1?(c * _log10f(a * rgb.y + b) + d):(e * rgb.y + f);
    rgb.z = rgb.z>=cut1?(c * _log10f(a * rgb.z + b) + d):(e * rgb.z + f);
  } else if (tf == 14) { // User Controlled Generic Log2 Curve
    rgb = _log2f3(rgb / 0.18f);
    rgb = clampf3(rgb, generic_log2_min_expo, generic_log2_max_expo);
    rgb = (rgb + _fabs(generic_log2_min_expo)) / (fabs(generic_log2_min_expo)+fabs(generic_log2_max_expo));
  }
  // else if (tf == 12) {
  //     // BMFilm Placeholder
  // }
  return rgb;
}

//Based on Jed Smith Sigmoid
__DEVICE__ float tonescale(float in, float sp, float tp, float Pslope, float px, float py,float s0=1.0,float t0=0.0)
{
  //calculate Shoulder
  //float s0= 1.0;
  //float t0= 0.0;
  float ss =spowf(((spowf((Pslope*((s0-px)/(1-py))),sp)-1)*(spowf(Pslope*(s0-px),-sp))),-1/sp);
  float ms = Pslope*(in-px)/ss;
  float fs = ms/spowf(1+(spowf(ms,sp)),1/sp);

  //calculate Toe
  float ts =spowf(((spowf((Pslope*((px-t0)/(py))),tp)-1)*(spowf(Pslope*(px-t0),-tp))),-1/tp);
  float mr = (Pslope*(in-px))/-ts;
  float ft = mr/spowf(1+(spowf(mr,tp)),1/tp);

  in = in>=px? ss*fs+py:-ts*ft+py;

  return in;
}

// Pure LintoLog
__DEVICE__ float LintoLog(float in,float mn, float mx)
{
  float offs = _powf(2,mn);
  in = _logf(in/0.18+offs)/_logf(2);
  //in = _fmaxf(in, mn);
  in = (in-mn)/(mx-mn);

  return in;
}

  __DEVICE__ float LogtoLin(float in,float mn, float mx)
{
  float offs = _powf(2,mn);
  in = 0.18*(_powf(2,(in*(mx-mn)+mn))-offs);

  return in;
}

__DEVICE__ float ShoulderSigmoid(float in,float sp,float Pslength,float slx){
  float px =  0.5;
  // float py =  0.5;
  // float tlx = 0.0;
  // float tly = 0.0;
  float sly = 1.0;

  float stx = Pslength/_sqrtf(2) + px;
  float sty = stx;

  float ss = spowf(((spowf((slx-stx)/(sly-sty),sp)-1)*spowf(slx-stx,-sp)),-1/sp);

  float sx = (in-stx)/ss;
  float fs = sx/spowf(1+spowf(sx,sp),1/sp);

  float fss = ss*fs + sty;

  float Out = in > stx ? fss : in;

  return Out;
}

__DEVICE__ float3 eotf_hlg(float3 rgb, int inverse, float sdr_peak_nits=1000.0) {
  // Aply the HLG Forward or Inverse EOTF. Implements the full ambient surround illumination model
  // ITU-R Rec BT.2100-2 https://www.itu.int/rec/R-REC-BT.2100
  // ITU-R Rep BT.2390-8: https://www.itu.int/pub/R-REP-BT.2390
  // Perceptual Quantiser (PQ) to Hybrid Log-Gamma (HLG) Transcoding: https://www.bbc.co.uk/rd/sites/50335ff370b5c262af000004/assets/592eea8006d63e5e5200f90d/BBC_HDRTV_PQ_HLG_Transcode_v2.pdf

  const float HLG_Lw = 1000.0f;
  // const float HLG_Lb = 0.0f;
  const float HLG_Ls = 5.0f;
  const float h_a = 0.17883277f;
  const float h_b = 1.0f - 4.0f * 0.17883277f;
  const float h_c = 0.5f - h_a * _logf(4.0f * h_a);
  const float h_g = 1.2f * _powf(1.111f, _log2f(HLG_Lw / 1000.0f)) * _powf(0.98f, _log2f(_fmaxf(1e-6f, HLG_Ls) / 5.0f));

  if (inverse == 1) {
    rgb /= (HLG_Lw / sdr_peak_nits);
    float Yd = 0.2627f * rgb.x + 0.6780f * rgb.y + 0.0593f * rgb.z;
    // HLG Inverse OOTF
    rgb = rgb * _powf(Yd, (1.0f - h_g) / h_g);
    // HLG OETF
    rgb.x = rgb.x <= 1.0f / 12.0f ? _sqrtf(3.0f * rgb.x) : h_a * _logf(12.0f * rgb.x - h_b) + h_c;
    rgb.y = rgb.y <= 1.0f / 12.0f ? _sqrtf(3.0f * rgb.y) : h_a * _logf(12.0f * rgb.y - h_b) + h_c;
    rgb.z = rgb.z <= 1.0f / 12.0f ? _sqrtf(3.0f * rgb.z) : h_a * _logf(12.0f * rgb.z - h_b) + h_c;
  } else {
    // HLG Inverse OETF
    rgb.x = rgb.x <= 0.5f ? rgb.x * rgb.x / 3.0f : (_expf((rgb.x - h_c) / h_a) + h_b) / 12.0f;
    rgb.y = rgb.y <= 0.5f ? rgb.y * rgb.y / 3.0f : (_expf((rgb.y - h_c) / h_a) + h_b) / 12.0f;
    rgb.z = rgb.z <= 0.5f ? rgb.z * rgb.z / 3.0f : (_expf((rgb.z - h_c) / h_a) + h_b) / 12.0f;
    // HLG OOTF
    float Ys = 0.2627f * rgb.x + 0.6780f * rgb.y + 0.0593f * rgb.z;
    rgb = rgb * _powf(Ys, h_g - 1.0f);
    rgb *= (HLG_Lw / sdr_peak_nits);
  }
  return rgb;
}

__DEVICE__ float3 eotf_pq(float3 rgb, int inverse, int jz, float sdr_peak_nits = 10000.0) {
  // Apply the ST-2084 PQ Forward or Inverse EOTF
  // Normalized such that input display linear light code value 1.0 equals 10,000 nits
  // ITU-R Rec BT.2100-2 https://www.itu.int/rec/R-REC-BT.2100
  // ITU-R Rep BT.2390-9 https://www.itu.int/pub/R-REP-BT.2390

  float Lp = 10000.0f;
  const float m1 = 2610.0f / 16384.0f;
  float m2 = 2523.0f / 32.0f;
  const float c1 = 107.0f / 128.0f;
  const float c2 = 2413.0f / 128.0f;
  const float c3 = 2392.0f / 128.0f;

  // Custom values for JzAzBz colorspace
  if (jz == 1) {
    m2 *= 1.7f;
    Lp = 10000.0f;
  }

  if (inverse == 1) {
    rgb /= (Lp / sdr_peak_nits);
    rgb = spowf3(rgb, m1);
    // Prevent shitting of the bed when there are negatives, for JzAzBz conversion
    rgb.x = _sign(rgb.x) * _powf((c1 + c2 * _fabs(rgb.x)) / (1.0f + c3 * _fabs(rgb.x)), m2);
    rgb.y = _sign(rgb.y) * _powf((c1 + c2 * _fabs(rgb.y)) / (1.0f + c3 * _fabs(rgb.y)), m2);
    rgb.z = _sign(rgb.z) * _powf((c1 + c2 * _fabs(rgb.z)) / (1.0f + c3 * _fabs(rgb.z)), m2);
  } else {
    rgb = spowf3(rgb, 1.0f / m2);
    rgb.x = _sign(rgb.x) * _powf((_fabs(rgb.x) - c1) / (c2 - c3 * _fabs(rgb.x)), 1.0f / m1) * Lp;
    rgb.y = _sign(rgb.y) * _powf((_fabs(rgb.y) - c1) / (c2 - c3 * _fabs(rgb.y)), 1.0f / m1) * Lp;
    rgb.z = _sign(rgb.z) * _powf((_fabs(rgb.z) - c1) / (c2 - c3 * _fabs(rgb.z)), 1.0f / m1) * Lp;
    rgb *= (Lp / sdr_peak_nits);
  }
  return rgb;
}

__DEVICE__ float3 sRGB_piecewise_transfer_function(float3 rgb, bool decode) {
  if (decode) {
    rgb.x = rgb.x <= 0.04045f ? rgb.x / 12.92f : _powf((rgb.x + 0.055f) / 1.055f, 2.4f);
    rgb.y = rgb.y <= 0.04045f ? rgb.y / 12.92f : _powf((rgb.y + 0.055f) / 1.055f, 2.4f);
    rgb.z = rgb.z <= 0.04045f ? rgb.z / 12.92f : _powf((rgb.z + 0.055f) / 1.055f, 2.4f);
  } else {
    rgb.x = rgb.x <= 0.0031308f ? 12.92f * rgb.x : 1.055f * _powf(rgb.x, 1.0f / 2.4f) - 0.055f;
    rgb.y = rgb.y <= 0.0031308f ? 12.92f * rgb.y : 1.055f * _powf(rgb.y, 1.0f / 2.4f) - 0.055f;
    rgb.z = rgb.z <= 0.0031308f ? 12.92f * rgb.z : 1.055f * _powf(rgb.z, 1.0f / 2.4f) - 0.055f;
  }
  return rgb;
}
