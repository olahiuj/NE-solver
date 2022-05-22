import os


def nash(in_path, out_path):
    # load file
    # get NE
    # write file
    pass

if __name__ == '__main__':
    for f in os.listdir('input'):
        if f.endswith('.nfg'):
            nash('input/'+f, 'output/'+f.replace('nfg','ne'))