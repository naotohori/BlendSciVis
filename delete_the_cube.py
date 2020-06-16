import bpy

bpy.data.objects['Cube'].select_set(True)
bpy.ops.object.delete()

#bpy.ops.wm.save_as_mainfile(filepath='A_new_world_without_the_cube.blend')
# This returns an error
# "ERROR (bke.bpath): /home/sources/buildbot-x86_64-slave/linux_glibc217_x86_64_cmake/blender.git/source/blender/blenkernel/intern/bpath.c:233 BKE_bpath_relative_convert: basedir='', this is a bug"
# although it can successfully generates a file and finishes.

bpy.ops.wm.save_mainfile(filepath='A_new_world_without_the_cube.blend')
# This one does not return the error.
