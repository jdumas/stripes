uniform vec3 eye;
uniform vec3 light;
varying vec3 position;
varying vec3 normal;
varying vec3 color;
varying vec4 texcoord;
varying vec3 barycoords;
varying vec2 param;
varying float index;
uniform sampler2D stripeTex;
uniform float stripeAmount;

#define M_PI 3.1415926535897932384626433832795

float diffuse( vec3 N, vec3 L )
{
   return max( 0., dot( N, L ));
}

float specular( vec3 N, vec3 L, vec3 E )
{
   const float shininess = 8.;
   vec3 R = 2.*dot(L,N)*N - L;
   return pow( max( 0., dot( R, E )), shininess );
}

float fresnel( vec3 N, vec3 E )
{
   const float sharpness = 10.;
   float NE = max( 0., dot( N, E ));
   return pow( sqrt( 1. - NE*NE ), sharpness );
}

void main()
{
   vec3 N = normalize( normal );
   vec3 L = normalize( light - position );
   vec3 E = normalize( eye - position );
   vec3 R = 2.*dot(L,N)*N - L;
   vec3 one = vec3( 1., 1., 1. );
   const float ambient = 0.3;

   gl_FragColor.rgb = (ambient+diffuse(N,L))*color + .5*specular(N,L,E)*one + .5*fresnel(N,E)*one;
   gl_FragColor.a = .1;
}

