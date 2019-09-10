varying vec3 position;
varying vec3 normal;
varying vec3 color;
varying vec4 texcoord;
varying vec3 barycoords;
varying vec2 param;
varying float index;

void main()
{	
   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;

   position = gl_Vertex.xyz;
   normal = gl_Normal.xyz;
   color = gl_Color.xyz;
   texcoord = gl_MultiTexCoord0.xyzw;
   barycoords = gl_MultiTexCoord1.xyz;
   index = gl_MultiTexCoord2.x;
   param = gl_MultiTexCoord3.xy;
}
