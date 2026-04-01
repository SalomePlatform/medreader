// Copyright (C) 2022-2026  CEA, EDF
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#include "vtkMEDReader.h"

#include <functional>

class vtkMutableDirectedGraph;
class vtkCollection;

using SILModifierFunc = std::function< void(vtkMutableDirectedGraph *, vtkCollection *) >;

namespace
{
    struct CallbackWrapper
    {
        SILModifierFunc func;
    };

    class CallbackForSILManager
    {
    public:
       CallbackForSILManager() = default;
       void setReader( vtkMEDReader *mr ) { _mr = mr; }
       void setCallbackInfo( SILModifierFunc func, vtkCollection *collection )
       {
        _cbw.func = func;
        if( _mr )
            _mr->setSILCallback( &_cbw, collection );
       }
       ~CallbackForSILManager()
       {
        if( _mr )
            _mr->setSILCallback( nullptr, nullptr);
       }
    private:
        vtkMEDReader *_mr = nullptr;
        CallbackWrapper _cbw;
    };
}
