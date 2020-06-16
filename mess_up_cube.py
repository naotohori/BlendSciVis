import bpy

#bpy.ops.object.editmode_toggle()
## To see what heppens if the mode is toggled.
## When the script starts, the mode is "OBJECT". If toggled to "EDIT" mode, the code below won't work.

bpy.data.meshes[0].vertices[0].co.x += 2

bpy.ops.wm.save_as_mainfile(filepath='a_cube_messed_up.blend')
