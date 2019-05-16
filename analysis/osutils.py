import os

def ensure_folder_exists(filename):
    if os.path.isdir(filename):
        pass
    elif os.path.isfile(filename):
        ensure_folder_exists(os.path.dirname(filename))
    else:
        head, tail = os.path.split(filename)
        if head and not os.path.isdir(head):
            ensure_folder_exists(head)
        #print "_mkdir %s" % repr(newdir)
        if tail:
            os.mkdir(filename)