// Copyright (c) 2015 Davide Rambaldi.
// All rights reserved.
//
// This file is part of ofxMol.
//
// ofxMol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ofxMol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ofxMol.  If not, see <http://www.gnu.org/licenses/>.
//
//
// Additional permission under GNU GPL version 3 section 7
//
// If you modify this Library, or any covered work, by linking or
// combining it with ofxMol (or a modified version of that library), the
// licensors of this Library grant you additional permission to convey
// the resulting work. Corresponding Source for a non-source form of
// such a combination shall include the source code for the parts of CGAL
// used as well as that of the covered work.
//
//
//
// Author(s)     :  Davide Rambaldi

#pragma once

#include "ofMain.h"
#include <ESBTL/default.h>
#include <ESBTL/coarse_creators.h>
#include <ESBTL/coarse_classifier.h>

namespace OfxMol
{
    class Coarse_Atom
    {
    public:
        Coarse_Atom() : _atom(ESBTL::Default_system_with_coarse_grain::Residue::Coarse_atom()) {}
        Coarse_Atom(ESBTL::Default_system_with_coarse_grain::Residue::Coarse_atom eatom);
        ~Coarse_Atom() {}
        
        const ofVec3f position() const { return ofVec3f(_atom.x(),_atom.y(),_atom.z()); }
        
        ofFloatColor getColor() const
        {
            return color;
        }
        
        float radius() const
        {
            return _radius;
        }
        
        bool is_backbone()
        {
            return _is_backbone;
        }
        
        //! Return a sphere primitive
        ofSpherePrimitive sphere(int resolution = 24);
        
    protected:
        ESBTL::Default_system_with_coarse_grain::Residue::Coarse_atom _atom;
        ofFloatColor color;
        ESBTL::Generic_classifier<ESBTL::Radius_of_coarse_atom<double,ESBTL::Default_system_with_coarse_grain::Residue::Coarse_atom> > radius_classifier;
        float _radius;
        bool _is_backbone;
    };
}