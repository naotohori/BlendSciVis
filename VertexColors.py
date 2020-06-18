import bpy

# Delete the cube
bpy.data.objects['Cube'].select_set(True)
bpy.ops.object.delete()

# Create a new mesh
ob_name = "triangle"
mesh_name = ob_name + "_mesh"

objects = bpy.data.objects
meshes = bpy.data.meshes
collection = bpy.data.collections[0]

if mesh_name in meshes:
    meshes.remove(meshes[mesh_name])
mesh = meshes.new(mesh_name)

if ob_name in bpy.data.objects:
    objects.remove(objects[ob_name])
ob = objects.new(ob_name, mesh)
collection.objects.link(ob)


# Define some geometry
verts = [ (0,0,0), (0,2,0), (0,1,2) , (0,3,2) ]
edges = [ (0,1), (1,2), (2,0), (1,3), (3, 2) ]
faces = [ (0,1,2), (1,3,2) ]

mesh.from_pydata(verts, edges, faces)

# Check indices
for poly in mesh.polygons:
    for i in range(len(poly.vertices)):
        print('Polygon index: ', list(mesh.polygons).index(poly),
        mesh.vertices.data.vertices[poly.vertices[i]].co,
        'Mesh index: ', poly.vertices[i],
        'Loop index: ', poly.loop_indices[i]) 


# Colors
color_layer = mesh.vertex_colors.new(name='vert_colors')
print(color_layer)

# ONE COLOR PER VERTEX
vertex_colors = [[1,0,0,1], [0,1,0,1], [0,1,0,1], [1,0,0,1]]
for poly in mesh.polygons:
    for i in range(len(poly.vertices)):
        color_layer.data[poly.loop_indices[i]].color = vertex_colors[poly.vertices[i]]
        print( vertex_colors[poly.vertices[i]])

# Check if the color is assigned
for poly in mesh.polygons:
    for i in range(len(poly.vertices)):
        print(list(color_layer.data[poly.loop_indices[i]].color))
        
print('Number of entries in the color layer')
print('len(color_layer.data) = %i' % len(color_layer.data))
print('Total number of vertices')
print('len(mesh.vertices) = %i' % len(mesh.vertices))


# ONE COLOR PER POLYGON
#vertex_colors_per_loop_index = [[1,0,0,1],[1,0,0,1],[1,0,0,1], [0,1,0,1], [0,1,0,1], [0,1,0,1]]
#
#for poly in mesh.polygons:
#    for i in range(len(poly.vertices)):
#        color_layer.data[poly.loop_indices[i]].color = vertex_colors_per_loop_index[poly.loop_indices[i]]
        


# Make material
triangle_material_name = "triangle_mat"
mat = bpy.data.materials.new(triangle_material_name)

mat.use_nodes = True
nodes = mat.node_tree.nodes
links = mat.node_tree.links

# Clear default nodes
nodes.clear()

# Create an output for the shader
node_output = nodes.new(type='ShaderNodeOutputMaterial')
node_output.location = 400, 300

shader = nodes.new(type='ShaderNodeBsdfPrincipled')
shader.location = 0, 300 # Location in the node window
#shader.inputs['Base Color'].default_value = (1,0,0,1)
links.new(shader.outputs['BSDF'], node_output.inputs['Surface'])

input_vertex_color = nodes.new(type='ShaderNodeVertexColor')
input_vertex_color.location = -400, 300 # Location in the node window
input_vertex_color.layer_name = 'vert_colors'
links.new(input_vertex_color.outputs['Color'], shader.inputs['Base Color'])

mesh.materials.append( mat )