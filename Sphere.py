import bpy
import math
from mathutils import Vector

''' Input '''
n = 12    # Number of rings (vertical, theta)
m = 24    # Number of segments (horizontal, phi)
r = 2.0   # Radius
c = 0.7   # Ratio 0 < c <= 1.0  (c = 1.0 ==> Perfect sphere)
loc = Vector((2,3,4))  # Position of the shpere origin

# Create a new mesh
ob_name = "Sphere"
if ob_name in bpy.data.objects:
    bpy.data.objects.remove(bpy.data.objects[ob_name])

mesh_name = ob_name + "_mesh"
if mesh_name in bpy.data.meshes:
    bpy.data.meshes.remove(bpy.data.meshes[mesh_name])
mesh = bpy.data.meshes.new(mesh_name)

# Create a new object with the mesh
ob = bpy.data.objects.new(ob_name, mesh)

''' ############ Vertices ############## '''
verts = [] 
''' Top'''
verts.append(Vector((0,0,c*r)))

''' In case cutting equally in z direction
# The ratio (c) was not implemented in this version.
for i in range(n-1):

    z = 1.0 - (i+1)*2/n
    theta = math.acos(z)

    for j in range(m):
        phi = j * 2 * math.pi / m

        x = r * math.sin(theta) * math.cos(phi)
        y = r * math.sin(theta) * math.sin(phi)
    
        verts.append(Vector((x,y,z)))
'''
for i in range(n-1):

    theta = (i+1) * math.pi / n
    z = c*r * math.cos(theta)

    for j in range(m):
        phi = j * 2 * math.pi / m

        x = r * math.sin(theta) * math.cos(phi)
        y = r * math.sin(theta) * math.sin(phi)
    
        verts.append(Vector((x,y,z)))

''' Bottom '''
verts.append(Vector((0,0,-c*r)))

''' Translate to the location '''
verts = [loc + v for v in verts]

''' ############ Edges ############## '''
edges = []  # Do not set edges

''' ############ Faces ############## '''
faces = []
idx = [0]*m
preset = 1
for i in range(n):
    idx_pre = idx
    if i < n-1:
        idx = [preset + i for i in range(m)]
    else:
        idx = [preset]*m

    for j in range(m):

        k = j + 1
        if k == m:
            k = 0

        faces.append((idx_pre[j], idx[j], idx[k], idx_pre[k]))

    preset += m

# Add it to the mesh
mesh.from_pydata(verts, edges, faces)

# Link the object to the first collection
bpy.data.collections[0].objects.link(ob)