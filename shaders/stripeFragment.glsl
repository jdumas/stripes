uniform vec3 eye;
uniform vec3 light;
varying vec3 position;
varying vec3 normal;
varying vec3 color;
varying vec4 texcoord;
varying vec3 barycoords;
uniform sampler2D stripeTex;

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

   const float f = 1.; // controls line frequency
   const float s = 20.; // controls line sharpness
   const float w = .93; // controls line width

   float theta = texcoord.x;
   float phi = texcoord.y;
   float nu = texcoord.z;
   float nv = texcoord.w;
   
   // compute 1Log_n
   vec3 b = barycoords;
   float lLog = 0.;
        if( b.z <= b.x && b.z <= b.y ) lLog = (M_PI/3.) * ( 1. + (b.y-b.x)/(1.-3.*b.z));
   else if( b.x <= b.y && b.x <= b.z ) lLog = (M_PI/3.) * ( 3. + (b.z-b.y)/(1.-3.*b.x));
   else                                lLog = (M_PI/3.) * ( 5. + (b.x-b.z)/(1.-3.*b.y));

   // adjust texture coordinates
   theta += lLog * nu;
   phi   += lLog * nv;

   vec2 t = vec2( theta/(2.*M_PI), phi/(2.*M_PI) );
   vec3 c = texture2D( stripeTex, t ).rgb;
   gl_FragColor.rgb = diffuse(N,L)*c + .5*specular(N,L,E)*one + .5*fresnel(N,E)*one;

   //float u = 1./(1.+exp(s*(cos(f*( theta ))-w)));
   //float v = 1./(1.+exp(s*(cos(f*( phi   ))-w)));
   //
   //u = 1.-u;
   //v = 1.-v;
   //vec3 c = (1.-u)*(1.-v)*vec3( 1.,1.,1.) + u*vec3(0.,0.,1.) + v*vec3(1.,0.,0.); // red and blue principal lines; don't multiply by u,v

   //gl_FragColor.rgb = diffuse(N,L)*c + .5*specular(N,L,E)*one*u*v + .5*fresnel(N,E)*one;
   //gl_FragColor.rgb = diffuse(N,L)*color + .5*specular(N,L,E)*one*u*v + .5*fresnel(N,E)*one;
   gl_FragColor.a = 1.;
}


// u = 1.-u; v = 1.-v; vec3 c = (1.-u)*(1.-v)*vec3( 1.,1.,1.) + u*vec3(0.,0.,1.) + v*vec3(1.,0.,0.); // red and blue principal lines; don't multiply by u,v
