#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2019-2026  CEA, EDF
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#
# Author : Anthony Geay (EDF R&D)

def WriteInTmpDir(func):
    def decoratedFunc(*args,**kwargs):
        import tempfile,os,platform
        ret = None
        currentDir = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            ret = func(*args,**kwargs)
            # On Windows, if one deletes the temporary directory while busy
            # with it, it raises an error. So we change back to the original
            if platform.system() == "Windows" :
                os.chdir(currentDir)
        return ret
    return decoratedFunc

def RetriveBaseLine(imgFile):
    import os,sys
    try:
        baselineIndex = sys.argv.index('-B')+1
        baselinePath = sys.argv[baselineIndex]
    except:
        print("Could not get baseline directory. Test failed.")
        exit(1)
    baseline_file = os.path.join(baselinePath, imgFile)
    return os.path.abspath(baseline_file)

