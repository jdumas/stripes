varying vec3 position;
varying vec3 normal;
varying vec3 color;
varying vec4 texcoord;
varying vec3 barycoords;

void main()
{	
   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;

   position = gl_Vertex.xyz;
   normal = gl_Normal.xyz;
   color = gl_Color.xyz;
   texcoord = gl_MultiTexCoord0.xyzw;
   barycoords = gl_MultiTexCoord1.xyz;
}
