import os
import tempfile
import webbrowser
import threading
import pickle
tmp = []
GraphVizExt = ["bmp",
"canon",
"dot",
"gv",
"xdot",
"xdot1.2",
"xdot1.4",
"cgimage",
"cmap",
"eps",
"exr",
"fig",
"gd",
"gd2",
"gif",
"gtk",
"ico",
"imap",
"cmapx",
"imap_np",
"cmapx_np",
"ismap",
"jp2",
"jpg",
"jpeg",
"jpe",
"json",
"json0",
"dot_json",
"xdot_json",
"pct",
"pict",
"pdf",
"pic",
"plain",
"plain-ext",
"png",
"pov",
"ps",
"ps2",
"psd",
"sgi",
"svg",
"svgz",
"tga",
"tif",
"tiff",
"tk",
"vml",
"vmlz",
"vrml",
"wbmp",
"webp",
"xlib",
"x11"]
def save_graphviz(graphviz_str, filename, extension=None, directory=None, verbose=False):
    (fd,file_dot) = tempfile.mkstemp(suffix=".dot",
                                     dir=directory)
    file_out = filename    
    l = file_out.split(".")
    if extension == None:
        extension = l[len(l)-1]
    if extension not in GraphVizExt:
        raise ValueError("Not a valid Graphviz output format")
    f = open(fd, "w")
    f.write(graphviz_str)
    f.close()
    compil = 'dot -T{} {} -o {}'.format(extension,
                                        file_dot,
                                        file_out)
    if verbose:
        print("Graphviz Compilation:", compil)
    os.system(compil)
    
def view_graphviz(graphviz_str,save_to_file=None, extension="svg", directory=None, verbose=False):
    if (save_to_file==None):
        (fd,file_svg) = tempfile.mkstemp(suffix="."+extension,
                                         dir=directory)
        open(fd).close()
    else:
        file_svg = save_to_file
    save_graphviz(graphviz_str, file_svg,
                  extension=extension,
                  directory=directory,
                  verbose=verbose)
    webbrowser.open(file_svg)
    
def delete_file(s):
     os.system('rm %s'% s)

def save(obj,file_name):
    f = open(file_name,"wb")
    pickle.dump(obj,f)
    f.close()
def load(file_name):
    f = open(file_name,"rb")
    obj = pickle.load(f)
    f.close()
    return obj

def gcd(a,b):
    """Compute the greatest common divisor of a and b"""
    while b > 0:
        a, b = b, a % b
    return a
    
def lcm(a, b):
    """Compute the lowest common multiple of a and b"""
    return a * b / gcd(a, b)
