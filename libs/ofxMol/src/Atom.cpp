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

#include "ofxMol/Atom.h"

namespace OfxMol
{
    typedef ESBTL::Color_of_atom<ESBTL::Default_system_with_coarse_grain::Residue::Atom> OfxMol_Atom_color;
    
    Atom::Atom(ESBTL::Default_system_with_coarse_grain::Atom eatom)
    {
        _atom = eatom;
        OfxMol_Atom_color color_of = OfxMol_Atom_color();
        std::vector<string> rgb = ofSplitString(color_of(eatom), ",");
        _color.set( ofToFloat(rgb[0]), ofToFloat(rgb[1]), ofToFloat(rgb[2]));
        _name = ESBTL::get_atom_name(eatom);
         _is_backbone = ESBTL::is_backbone(eatom);
        _radius = radius_classifier.get_properties(eatom).value();
    }
    
    std::string Atom::log()
    {
        std::string msg;
        msg += "Name: [" + name() + "]";
        msg += " residue name: [" + residue_name() + "]";
        msg += " element [" + element() + "]";
        msg += " position [" + ofToString(position()) + "]";
        msg += " backbone: [" + ofToString(is_backbone()? "yes" : "no") + "]";
        msg += " radius: [" + ofToString(radius()) + "]";
        msg += " color: [" + ofToString(getColor()) + "]";
        msg += "\n";
        return msg;
    }
    
    ofSpherePrimitive Atom::sphere(int resolution)
    {
        return makeSpherePrimitive(resolution, radius(), position());
    }
    
    ofSpherePrimitive Atom::sphere(float radius, int resolution)
    {
        return makeSpherePrimitive(resolution, radius, position());
    }
    
    ofSpherePrimitive Atom::makeSpherePrimitive(int resolution, float radius, ofVec3f position)
    {
        ofSpherePrimitive sphere;
        sphere.setResolution(resolution);
        sphere.setRadius(radius);
        sphere.setPosition(position);
        return sphere;
    }
}






