#! /usr/bin/env python
import glob

def main():
    luaList = glob.glob("s*/*.lua")
    for luaFile in luaList:

        outName = luaFile[:luaFile.find(".lua")] + ".rst"
        print "Working on %s (outname is %s) ..." % (luaFile, outName)
        fh = open(outName, "w")
        fh.writelines(".. literalinclude:: %s\n" % luaFile[luaFile.find("/")+1:])
        fh.writelines("  :language: %s\n" % "lua")
        fh.close() 

if __name__ == "__main__" : main()
