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
    class Atom
    {
    public:
        Atom() : _atom(ESBTL::Default_system_with_coarse_grain::Atom()),
            _color(ofColor()), _name(""), _is_backbone(false), _radius(0.0f) {}
        
        Atom( ESBTL::Default_system_with_coarse_grain::Atom eatom);
        ~Atom(){}
        
        //! Logger
        std::string log();
        //! Return a sphere primitive
        ofSpherePrimitive sphere(int resolution = 24);
        //! Return a sphere primitive with given radius
        ofSpherePrimitive sphere(float radius, int resolution = 24);
        
        const ofVec3f position() const { return ofVec3f(_atom.x(),_atom.y(),_atom.z()); }
        const double occupancy() const { return _atom.occupancy(); }
        const std::string residue_name() const { return _atom.residue_name(); }
        void setColor(ofFloatColor new_color)
        {
            _color = new_color;
        }
        
        ofFloatColor getColor() const
        {
            return _color;
        }
        
        std::string name() const
        {
            return _name;
        }
        
        std::string element() const
        {
            return _atom.element();
        }
        
        bool is_backbone() const
        {
            return _is_backbone;
        }
        
        /*
         * Function filling default radius of atoms. Current implementation uses radii
         * from Tsai J, Taylor R, Chothia C, Gerstein M. J Mol Biol. 1999 Jul 2;290(1):253-66.
         * see table: <ESBTL/properties/Tsai_jmb_99_radii.h>
         * Default is: 1.8
         */
        double radius() const
        {
            return _radius;
        }
        
    protected:
        ESBTL::Default_system_with_coarse_grain::Atom _atom;
        ESBTL::Generic_classifier<ESBTL::Radius_of_atom<double,ESBTL::Default_system_with_coarse_grain::Residue::Atom> > radius_classifier;
        ofFloatColor _color;
        std::string _name;
        bool _is_backbone;
        double _radius;
        ofSpherePrimitive makeSpherePrimitive(int resolution, float radius, ofVec3f position);
    };
}






