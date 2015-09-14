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

#include "ofxMol/Coarse_Atom.h"

namespace OfxMol
{
    typedef ESBTL::Color_of_atom<ESBTL::Default_system_with_coarse_grain::Residue::Coarse_atom> OfxMol_Coarse_Atom_color;
    
    Coarse_Atom::Coarse_Atom(ESBTL::Default_system_with_coarse_grain::Residue::Coarse_atom eatom)
    {
        _atom = eatom;
        _radius = radius_classifier.get_properties(eatom).value();
        OfxMol_Coarse_Atom_color color_of = OfxMol_Coarse_Atom_color();
        std::vector<string> rgb = ofSplitString(color_of(eatom), ",");
        color.set( ofToFloat(rgb[0]), ofToFloat(rgb[1]), ofToFloat(rgb[2]));
    }
    
    ofSpherePrimitive Coarse_Atom::sphere(int resolution)
    {
        ofSpherePrimitive sphere;
        sphere.setResolution(resolution);
        sphere.setRadius(radius());
        sphere.setPosition(position());
        return sphere;
    }
}