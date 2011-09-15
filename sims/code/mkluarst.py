#! /usr/bin/env python
import optparse

def main():
    parser = optparse.OptionParser()
    parser.add_option("-l", "--lua", dest="lua", help="Name of Lua program")
    options, args = parser.parse_args()

    options, args = parser.parse_args()
    luaFile = options.lua

    if luaFile == None:
        print "Must provide a Lua file to process!"
        exit(1)

    # contruct name
    outName = luaFile[:-3] + "rst"
    fh = open(outName, "w")
    fh.writelines(".. literalinclude:: %s\n" % luaFile[luaFile.find("/")+1:])
    fh.writelines("  :language: lua\n")

    fh.close()

if __name__ == "__main__" : main()

